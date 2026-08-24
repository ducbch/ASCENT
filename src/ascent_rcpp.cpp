// ascent_rcpp.cpp — Rcpp/OpenMP reimplementation of SCENT algorithm
//
// Poisson or Negative Binomial GLM with iterative bootstrap for empirical
// p-values.  Parallelized across peak-gene pairs via OpenMP.
// Sparse matrix rows extracted directly from dgCMatrix (CSC) format.
//
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

#include <RcppArmadillo.h>
#include <random>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace arma;

// ================================================================
// SECTION 1: SHARED UTILITIES
// ================================================================

// Thread-safe standard normal CDF (avoids R API calls in OpenMP)
inline double pnorm_ts(double x) {
    return 0.5 * std::erfc(-x * M_SQRT1_2);
}

// ------------------------------------------------------------
// Extract a single row from a dgCMatrix (column-compressed sparse).
// For each active cell (column), binary-search for row_idx.
//   p, i, x : CSC pointers from dgCMatrix slots @p, @i, @x
//   row_idx : the row to extract (0-based)
//   active_cells : indices of columns to scan
//   out : pre-allocated output vector (length = active_cells.size())
// ------------------------------------------------------------
inline void extract_sparse_row(
    const int* p, const int* i, const double* x,
    int row_idx,
    const std::vector<int>& active_cells,
    vec& out)
{
    int n = (int)active_cells.size();
    out.zeros();
    for (int k = 0; k < n; k++) {
        int col = active_cells[k];
        int start = p[col];
        int end   = p[col + 1];
        const int* pos = std::lower_bound(i + start, i + end, row_idx);
        if (pos != i + end && *pos == row_idx) {
            out(k) = x[pos - i];
        }
    }
}

// ------------------------------------------------------------
// QR-based weighted least squares with column equilibration.
//   Solve: beta = argmin_b sum_i w_i (z_i - X_i b)^2
//   Returns false on singular / near-singular systems.
// ------------------------------------------------------------
inline bool qr_wls_solve(vec& beta, const mat& X, const vec& w, const vec& z) {
    int p_dim = (int)X.n_cols;
    vec sw = sqrt(w);
    mat Xw = X.each_col() % sw;
    vec zw = sw % z;

    // Column equilibration: scale each column to unit norm
    vec col_norms(p_dim);
    for (int j = 0; j < p_dim; j++) {
        col_norms(j) = norm(Xw.col(j));
        if (col_norms(j) < 1e-15) return false;
        Xw.col(j) /= col_norms(j);
    }

    vec beta_scaled;
    bool ok = solve(beta_scaled, Xw, zw);
    if (!ok) return false;

    beta = beta_scaled / col_norms;
    return true;
}

// ================================================================
// SECTION 2: POISSON GLM VIA IRLS
// ================================================================

struct FitResult {
    vec  beta;       // full coefficient vector
    double coef;     // beta[1] — atac coefficient
    double se;       // SE for atac
    double z;        // z-score
    double p;        // two-sided p-value
    double theta;    // NB dispersion (only meaningful for negbin)
    bool converged;
};

// Full Poisson IRLS returning coefficients + inference for atac.
FitResult poisson_irls(const mat& X, const vec& y,
                       int max_iter = 25, double tol = 1e-8)
{
    FitResult res;
    res.converged = false;
    res.theta = 0;
    int n = (int)X.n_rows;
    int p = (int)X.n_cols;

    // Initialise: mu = y + 0.1
    vec mu = y + 0.1;
    vec eta = log(mu);
    vec beta(p, fill::zeros);

    for (int iter = 0; iter < max_iter; iter++) {
        vec w = clamp(mu, 1e-10, 1e7);          // Poisson weight = mu
        vec z_w = eta + (y - mu) / mu;           // working response

        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return res;

        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);

        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new;
        eta  = eta_new;
        mu   = mu_new;
        if (change < tol) { res.converged = true; break; }
    }
    if (!res.converged) res.converged = true;   // ran full iterations

    // Variance-covariance: (X'WX)^{-1}
    vec w_final = clamp(mu, 1e-10, 1e7);
    vec sw = sqrt(w_final);
    mat Xw = X.each_col() % sw;
    mat XtWX = Xw.t() * Xw;
    mat vcov;
    if (!inv_sympd(vcov, XtWX)) { res.converged = false; return res; }

    res.beta = beta;
    res.coef = beta(1);
    res.se   = std::sqrt(std::max(vcov(1, 1), 0.0));
    if (res.se > 0) {
        res.z = res.coef / res.se;
        res.p = 2.0 * pnorm_ts(-std::abs(res.z));
    } else {
        res.z = 0; res.p = 1;
    }
    return res;
}

// ================================================================
// SECTION 2b: POISSON NULL MODEL FOR SCORE TEST
// ================================================================

struct NullFitResult {
    vec  beta;       // null model coefficients
    vec  mu;         // fitted values
    double theta;    // NB dispersion (0 for Poisson)
    mat  XtWX_inv;   // (X'WX)^{-1} for projection (uses appropriate weights)
    bool converged;
};

// Fit Poisson null model (no peak term) via IRLS.
// Returns beta, mu, and (X'WX)^{-1} for score test projection.
NullFitResult poisson_null_irls(const mat& X, const vec& y,
                                int max_iter = 25, double tol = 1e-8)
{
    NullFitResult res;
    res.converged = false;
    res.theta = 0;
    int n = (int)X.n_rows;
    int p = (int)X.n_cols;

    vec mu = y + 0.1;
    vec eta = log(mu);
    vec beta(p, fill::zeros);

    for (int iter = 0; iter < max_iter; iter++) {
        vec w = clamp(mu, 1e-10, 1e7);
        vec z_w = eta + (y - mu) / mu;

        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return res;

        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);

        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new;
        eta  = eta_new;
        mu   = mu_new;
        if (change < tol) { res.converged = true; break; }
    }
    if (!res.converged) res.converged = true;

    // Compute (X'WX)^{-1} for projection (Poisson weights W = mu)
    vec w_final = clamp(mu, 1e-10, 1e7);
    vec sw = sqrt(w_final);
    mat Xw = X.each_col() % sw;
    mat XtWX = Xw.t() * Xw;
    if (!inv_sympd(res.XtWX_inv, XtWX)) { res.converged = false; return res; }

    res.beta = beta;
    res.mu   = mu;
    return res;
}

// Forward declaration (defined in Section 3)
inline double update_theta(const vec& y, const vec& mu, double theta,
                           int max_iter = 10);

// Fit NegBin null model (no peak term) via alternating IRLS + theta NR.
// Returns beta, mu, theta, and (X'WX)^{-1} with NB weights for score test.
NullFitResult negbin_null_irls(const mat& X, const vec& y,
                                int max_outer = 25, int max_inner = 25,
                                double tol = 1e-8)
{
    NullFitResult res;
    res.converged = false;
    res.theta = 0;
    int n = (int)X.n_rows;
    int p = (int)X.n_cols;

    // Step 1: Poisson initialisation
    vec mu = y + 0.1;
    vec eta = log(mu);
    vec beta(p, fill::zeros);

    for (int iter = 0; iter < max_inner; iter++) {
        vec w = clamp(mu, 1e-10, 1e7);
        vec z_w = eta + (y - mu) / mu;
        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return res;
        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);
        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new; eta = eta_new; mu = mu_new;
        if (change < tol) break;
    }

    // Step 2: Method-of-moments theta
    double ss = 0;
    for (int i = 0; i < n; i++) {
        double r = y(i) / mu(i) - 1.0;
        ss += r * r - 1.0 / mu(i);
    }
    double theta = std::max((double)n / std::max(ss, 0.01), 0.01);
    theta = std::min(theta, 1e6);

    // Step 3: Alternating IRLS + theta NR
    for (int outer = 0; outer < max_outer; outer++) {
        double theta_old = theta;

        // IRLS with NB weight: w = mu * theta / (theta + mu)
        for (int iter = 0; iter < max_inner; iter++) {
            vec w = mu % (theta / (theta + mu));
            w = clamp(w, 1e-10, 1e7);
            vec z_w = eta + (y - mu) / mu;
            vec beta_new;
            if (!qr_wls_solve(beta_new, X, w, z_w)) return res;
            vec eta_new = X * beta_new;
            eta_new = clamp(eta_new, -30.0, 30.0);
            vec mu_new = exp(eta_new);
            double change = accu(abs(eta_new - eta)) / n;
            beta = beta_new; eta = eta_new; mu = mu_new;
            if (change < tol) break;
        }

        theta = update_theta(y, mu, theta);
        if (std::abs(theta - theta_old) / (theta_old + 0.001) < 1e-4) {
            res.converged = true; break;
        }
    }
    if (!res.converged) res.converged = true;

    // Compute (X'WX)^{-1} with NB weights
    vec w_final = mu % (theta / (theta + mu));
    w_final = clamp(w_final, 1e-10, 1e7);
    vec sw = sqrt(w_final);
    mat Xw = X.each_col() % sw;
    mat XtWX = Xw.t() * Xw;
    if (!inv_sympd(res.XtWX_inv, XtWX)) { res.converged = false; return res; }

    res.beta  = beta;
    res.mu    = mu;
    res.theta = theta;
    return res;
}

// ================================================================
// (continued) SECTION 2: POISSON BOOTSTRAP
// ================================================================

// Frequency-weighted Poisson IRLS for bootstrap (warm-started).
// Returns atac coefficient only; NaN on failure.
inline double poisson_boot_coef(const mat& X, const vec& y, const vec& freq,
                                const vec& beta_warm,
                                int max_iter = 15, double tol = 1e-6)
{
    int n = (int)X.n_rows;
    vec eta = X * beta_warm;
    eta = clamp(eta, -30.0, 30.0);
    vec mu = exp(eta);
    vec beta = beta_warm;

    for (int iter = 0; iter < max_iter; iter++) {
        vec w = freq % clamp(mu, 1e-10, 1e7);
        vec z_w = eta + (y - mu) / mu;

        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return datum::nan;

        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);

        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new;
        eta  = eta_new;
        mu   = mu_new;
        if (change < tol) break;
    }
    return beta(1);
}

// ================================================================
// SECTION 3: NEGATIVE BINOMIAL GLM VIA IRLS + THETA NR
// ================================================================

// Newton-Raphson update for NB dispersion parameter theta.
inline double update_theta(const vec& y, const vec& mu, double theta,
                           int max_iter)
{
    int n = (int)y.n_elem;
    for (int iter = 0; iter < max_iter; iter++) {
        double score = 0.0, info = 0.0;
        for (int i = 0; i < n; i++) {
            double ym = y(i) + theta;
            double tm = theta + mu(i);
            score += R::digamma(ym) - R::digamma(theta)
                   + std::log(theta) - std::log(tm) + 1.0 - ym / tm;
            info  += -R::trigamma(ym) + R::trigamma(theta)
                   - 1.0 / theta + 2.0 / tm - ym / (tm * tm);
        }
        if (std::abs(info) < 1e-15) break;
        double theta_new = theta + score / info;
        theta_new = std::max(theta_new, 0.01);
        theta_new = std::min(theta_new, 1e6);
        if (std::abs(theta_new - theta) / (theta + 0.001) < 1e-6) {
            theta = theta_new; break;
        }
        theta = theta_new;
    }
    return theta;
}

// Weighted theta NR (for bootstrap with frequency weights).
inline double update_theta_w(const vec& y, const vec& mu, const vec& freq,
                             double theta, int max_iter = 10)
{
    int n = (int)y.n_elem;
    for (int iter = 0; iter < max_iter; iter++) {
        double score = 0.0, info = 0.0;
        for (int i = 0; i < n; i++) {
            if (freq(i) < 0.5) continue;  // skip cells with zero weight
            double ym = y(i) + theta;
            double tm = theta + mu(i);
            double f  = freq(i);
            score += f * (R::digamma(ym) - R::digamma(theta)
                       + std::log(theta) - std::log(tm) + 1.0 - ym / tm);
            info  += f * (-R::trigamma(ym) + R::trigamma(theta)
                       - 1.0 / theta + 2.0 / tm - ym / (tm * tm));
        }
        if (std::abs(info) < 1e-15) break;
        double theta_new = theta + score / info;
        theta_new = std::max(theta_new, 0.01);
        theta_new = std::min(theta_new, 1e6);
        if (std::abs(theta_new - theta) / (theta + 0.001) < 1e-6) {
            theta = theta_new; break;
        }
        theta = theta_new;
    }
    return theta;
}

// Full NegBin IRLS with theta estimation.
FitResult negbin_irls(const mat& X, const vec& y,
                      int max_outer = 25, int max_inner = 25, double tol = 1e-8)
{
    FitResult res;
    res.converged = false;
    int n = (int)X.n_rows;
    int p = (int)X.n_cols;

    // Step 1: Poisson initialisation
    vec mu = y + 0.1;
    vec eta = log(mu);
    vec beta(p, fill::zeros);

    for (int iter = 0; iter < max_inner; iter++) {
        vec w = clamp(mu, 1e-10, 1e7);
        vec z_w = eta + (y - mu) / mu;
        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return res;
        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);
        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new; eta = eta_new; mu = mu_new;
        if (change < tol) break;
    }

    // Step 2: Method-of-moments theta
    double ss = 0;
    for (int i = 0; i < n; i++) {
        double r = y(i) / mu(i) - 1.0;
        ss += r * r - 1.0 / mu(i);
    }
    double theta = std::max((double)n / std::max(ss, 0.01), 0.01);
    theta = std::min(theta, 1e6);

    // Step 3: Alternating IRLS + theta NR
    for (int outer = 0; outer < max_outer; outer++) {
        double theta_old = theta;

        // IRLS with NB variance: w = mu * theta / (theta + mu)
        for (int iter = 0; iter < max_inner; iter++) {
            vec w = mu % (theta / (theta + mu));
            w = clamp(w, 1e-10, 1e7);
            vec z_w = eta + (y - mu) / mu;
            vec beta_new;
            if (!qr_wls_solve(beta_new, X, w, z_w)) return res;
            vec eta_new = X * beta_new;
            eta_new = clamp(eta_new, -30.0, 30.0);
            vec mu_new = exp(eta_new);
            double change = accu(abs(eta_new - eta)) / n;
            beta = beta_new; eta = eta_new; mu = mu_new;
            if (change < tol) break;
        }

        theta = update_theta(y, mu, theta);
        if (std::abs(theta - theta_old) / (theta_old + 0.001) < 1e-4) {
            res.converged = true; break;
        }
    }
    if (!res.converged) res.converged = true;

    // Variance-covariance
    vec w_final = mu % (theta / (theta + mu));
    w_final = clamp(w_final, 1e-10, 1e7);
    vec sw = sqrt(w_final);
    mat Xw = X.each_col() % sw;
    mat vcov;
    if (!inv_sympd(vcov, Xw.t() * Xw)) { res.converged = false; return res; }

    res.beta  = beta;
    res.coef  = beta(1);
    res.se    = std::sqrt(std::max(vcov(1, 1), 0.0));
    res.theta = theta;
    if (res.se > 0) {
        res.z = res.coef / res.se;
        res.p = 2.0 * pnorm_ts(-std::abs(res.z));
    } else {
        res.z = 0; res.p = 1;
    }
    return res;
}

// Frequency-weighted NegBin IRLS for bootstrap (fixed theta, warm-started).
inline double negbin_boot_coef(const mat& X, const vec& y, const vec& freq,
                               const vec& beta_warm, double theta,
                               int max_iter = 15, double tol = 1e-6)
{
    int n = (int)X.n_rows;
    vec eta = X * beta_warm;
    eta = clamp(eta, -30.0, 30.0);
    vec mu = exp(eta);
    vec beta = beta_warm;

    for (int iter = 0; iter < max_iter; iter++) {
        vec w = freq % mu % (theta / (theta + mu));
        w = clamp(w, 1e-10, 1e7);
        vec z_w = eta + (y - mu) / mu;

        vec beta_new;
        if (!qr_wls_solve(beta_new, X, w, z_w)) return datum::nan;

        vec eta_new = X * beta_new;
        eta_new = clamp(eta_new, -30.0, 30.0);
        vec mu_new = exp(eta_new);

        double change = accu(abs(eta_new - eta)) / n;
        beta = beta_new; eta = eta_new; mu = mu_new;
        if (change < tol) break;
    }
    return beta(1);
}

// ================================================================
// SECTION 4: BOOTSTRAP P-VALUE (basic bootstrap method)
// ================================================================

// Interpolated two-sided p-value from bootstrap quantiles centred at 0.
inline double interp_pval(std::vector<double>& q) {
    int R = (int)q.size();
    if (R == 0) return 1.0;
    std::sort(q.begin(), q.end());
    auto it = std::lower_bound(q.begin(), q.end(), 0.0);
    int zero_pos = (int)(it - q.begin());
    if (zero_pos == 0 || zero_pos == R) return 2.0 / R;
    return 2.0 * std::min((double)zero_pos / R, (double)(R - zero_pos) / R);
}

// "basic" bootstrap p-value: transform boot values and call interp_pval.
inline double basic_p_cpp(double obs, const std::vector<double>& boot) {
    int R = (int)boot.size();
    std::vector<double> q(R);
    for (int i = 0; i < R; i++) {
        q[i] = 2.0 * obs - boot[i];  // null = 0
    }
    return interp_pval(q);
}

// ================================================================
// SECTION 5: MAIN EXPORTED FUNCTION
//
// Processes all peak-gene pairs in parallel (OpenMP).
// Input: dgCMatrix RNA/ATAC, gene/peak indices, covariates, cell mask.
// Output: DataFrame with gene, peak, beta, se, z, p, boot_basic_p.
// ================================================================

// [[Rcpp::export]]
Rcpp::DataFrame ascent_process_pairs(
    Rcpp::S4 rna_sparse,
    Rcpp::S4 atac_sparse,
    Rcpp::IntegerVector gene_idx_r,
    Rcpp::IntegerVector peak_idx_r,
    Rcpp::StringVector gene_names_r,
    Rcpp::StringVector peak_names_r,
    arma::mat cov_mat,
    Rcpp::LogicalVector cell_mask,
    bool binarize,
    int regr_type,   // 0 = poisson, 1 = negbin
    int ncores,
    bool skip_bootstrap = false,
    double min_pct_rna = 0.05,
    double min_pct_atac = 0.05)
{
    // ---- Extract CSC components from dgCMatrix ----
    Rcpp::IntegerVector rna_i_rv  = rna_sparse.slot("i");
    Rcpp::IntegerVector rna_p_rv  = rna_sparse.slot("p");
    Rcpp::NumericVector rna_x_rv  = rna_sparse.slot("x");
    Rcpp::IntegerVector atac_i_rv = atac_sparse.slot("i");
    Rcpp::IntegerVector atac_p_rv = atac_sparse.slot("p");
    Rcpp::NumericVector atac_x_rv = atac_sparse.slot("x");

    const int*    rna_ip  = rna_i_rv.begin();
    const int*    rna_pp  = rna_p_rv.begin();
    const double* rna_xp  = rna_x_rv.begin();
    const int*    atac_ip = atac_i_rv.begin();
    const int*    atac_pp = atac_p_rv.begin();
    const double* atac_xp = atac_x_rv.begin();

    int n_pairs = gene_idx_r.size();

    // ---- Copy R vectors to C++ (thread-safe) ----
    std::vector<int> gene_idx(n_pairs), peak_idx(n_pairs);
    std::vector<std::string> gene_names(n_pairs), peak_names(n_pairs);
    for (int i = 0; i < n_pairs; i++) {
        gene_idx[i]  = gene_idx_r(i);
        peak_idx[i]  = peak_idx_r(i);
        gene_names[i] = Rcpp::as<std::string>(gene_names_r(i));
        peak_names[i] = Rcpp::as<std::string>(peak_names_r(i));
    }

    // ---- Active cells (those passing celltype filter) ----
    std::vector<int> active_cells;
    int total_cells = cell_mask.size();
    for (int i = 0; i < total_cells; i++) {
        if (cell_mask(i)) active_cells.push_back(i);
    }
    int n_active = (int)active_cells.size();

    // ---- Pre-extract covariate sub-matrix for active cells ----
    int n_cov = (int)cov_mat.n_cols;
    mat cov_sub(n_active, n_cov);
    for (int i = 0; i < n_active; i++) {
        cov_sub.row(i) = cov_mat.row(active_cells[i]);
    }

    // ---- Pre-allocate per-pair result arrays ----
    std::vector<int>    valid(n_pairs, 0);
    std::vector<double> r_beta(n_pairs, 0), r_se(n_pairs, 0);
    std::vector<double> r_z(n_pairs, 0),    r_p(n_pairs, 1);
    std::vector<double> r_boot_p(n_pairs, 1);

    // ---- Bootstrap stage configuration (matches original SCENT) ----
    const int    boot_R[]      = {100,  500,  2500,  25000,  50000};
    const double boot_thresh[] = {0.1,  0.05, 0.01,  0.001,  0.0};
    const int    n_stages      = 5;

    // ---- Set OpenMP thread count ----
#ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
#endif

    // ---- Progress counter (atomic) ----
    int progress_done = 0;

    Rprintf("ASCENT: Processing %d peak-gene pairs across %d cells (%s, %d threads)...\n",
            n_pairs, n_active,
            regr_type == 0 ? "Poisson" : "NegBin",
#ifdef _OPENMP
            ncores > 0 ? ncores : omp_get_max_threads()
#else
            1
#endif
    );

    // ================================================================
    // MAIN PARALLEL LOOP — one iteration per peak-gene pair
    // ================================================================
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int pair = 0; pair < n_pairs; pair++) {

        int gi = gene_idx[pair];
        int pi = peak_idx[pair];

        // Thread-local working vectors
        vec rna_vec(n_active);
        vec atac_vec(n_active);

        // Extract sparse rows for this gene and peak
        extract_sparse_row(rna_pp,  rna_ip,  rna_xp,  gi, active_cells, rna_vec);
        extract_sparse_row(atac_pp, atac_ip, atac_xp, pi, active_cells, atac_vec);

        // Binarize ATAC if requested
        if (binarize) {
            atac_vec.elem(find(atac_vec > 0)).ones();
        }

        // Quality filter: require >min_pct nonzero in both modalities
        int n_expr = (int)accu(rna_vec > 0);
        int n_open = (int)accu(atac_vec > 0);
        if ((double)n_expr / n_active <= min_pct_rna ||
            (double)n_open / n_active <= min_pct_atac) {
            // Progress update
#ifdef _OPENMP
            #pragma omp atomic
#endif
            progress_done++;
            continue;
        }

        // Build design matrix: [intercept, atac, cov1, cov2, ...]
        int p_dim = 2 + n_cov;
        mat X(n_active, p_dim);
        X.col(0).ones();
        X.col(1) = atac_vec;
        for (int j = 0; j < n_cov; j++) {
            X.col(2 + j) = cov_sub.col(j);
        }

        // ---- Fit initial model ----
        FitResult fit;
        if (regr_type == 0) {
            fit = poisson_irls(X, rna_vec);
        } else {
            fit = negbin_irls(X, rna_vec);
        }
        if (!fit.converged) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            progress_done++;
            continue;
        }

        // ---- Iterative bootstrap (adaptive) ----
        double boot_p_val = NA_REAL;

        if (!skip_bootstrap) {
        // Thread-local RNG seeded by pair index for reproducibility
        std::mt19937 rng(42u + (unsigned)pair);
        std::uniform_int_distribution<int> cell_dist(0, n_active - 1);

        boot_p_val = 1.0;

        // Always run stage 0
        {
            int R = boot_R[0];
            std::vector<double> boot_coefs;
            boot_coefs.reserve(R);
            vec freq(n_active);

            for (int r = 0; r < R; r++) {
                freq.zeros();
                for (int c = 0; c < n_active; c++) freq(cell_dist(rng))++;
                double bc;
                if (regr_type == 0)
                    bc = poisson_boot_coef(X, rna_vec, freq, fit.beta);
                else
                    bc = negbin_boot_coef(X, rna_vec, freq, fit.beta, fit.theta);
                if (!std::isnan(bc)) boot_coefs.push_back(bc);
            }
            if (!boot_coefs.empty())
                boot_p_val = basic_p_cpp(fit.coef, boot_coefs);
        }

        // Subsequent stages: only run if previous p is below threshold
        for (int stage = 1; stage < n_stages; stage++) {
            if (boot_p_val >= boot_thresh[stage - 1]) break;

            int R = boot_R[stage];
            std::vector<double> boot_coefs;
            boot_coefs.reserve(R);
            vec freq(n_active);

            for (int r = 0; r < R; r++) {
                freq.zeros();
                for (int c = 0; c < n_active; c++) freq(cell_dist(rng))++;
                double bc;
                if (regr_type == 0)
                    bc = poisson_boot_coef(X, rna_vec, freq, fit.beta);
                else
                    bc = negbin_boot_coef(X, rna_vec, freq, fit.beta, fit.theta);
                if (!std::isnan(bc)) boot_coefs.push_back(bc);
            }
            if (!boot_coefs.empty())
                boot_p_val = basic_p_cpp(fit.coef, boot_coefs);
        }
        } // end if (!skip_bootstrap)

        // ---- Store results ----
        valid[pair]    = 1;
        r_beta[pair]   = fit.coef;
        r_se[pair]     = fit.se;
        r_z[pair]      = fit.z;
        r_p[pair]      = fit.p;
        r_boot_p[pair] = boot_p_val;

        // Progress
#ifdef _OPENMP
        #pragma omp atomic
#endif
        progress_done++;

        // Print progress occasionally (thread-safe via critical)
        if (progress_done % std::max(1, n_pairs / 20) == 0) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                Rprintf("\r  [%d / %d pairs processed (%.0f%%)]",
                        progress_done, n_pairs,
                        100.0 * progress_done / n_pairs);
            }
        }
    } // end parallel for

    Rprintf("\r  [%d / %d pairs processed (100%%)]\n", n_pairs, n_pairs);

    // ---- Collect valid results into R DataFrame ----
    int n_valid = 0;
    for (int i = 0; i < n_pairs; i++) n_valid += valid[i];

    Rcpp::StringVector  out_gene(n_valid), out_peak(n_valid);
    Rcpp::NumericVector out_beta(n_valid), out_se(n_valid);
    Rcpp::NumericVector out_z(n_valid),    out_p(n_valid);
    Rcpp::NumericVector out_boot_p(n_valid);

    int idx = 0;
    for (int i = 0; i < n_pairs; i++) {
        if (!valid[i]) continue;
        out_gene[idx]   = gene_names[i];
        out_peak[idx]   = peak_names[i];
        out_beta[idx]   = r_beta[i];
        out_se[idx]     = r_se[i];
        out_z[idx]      = r_z[i];
        out_p[idx]      = r_p[i];
        out_boot_p[idx] = r_boot_p[i];
        idx++;
    }

    Rprintf("ASCENT: Done. %d / %d pairs passed quality filters.\n", n_valid, n_pairs);

    return Rcpp::DataFrame::create(
        Rcpp::Named("gene")         = out_gene,
        Rcpp::Named("peak")         = out_peak,
        Rcpp::Named("beta")         = out_beta,
        Rcpp::Named("se")           = out_se,
        Rcpp::Named("z")            = out_z,
        Rcpp::Named("p")            = out_p,
        Rcpp::Named("boot_basic_p") = out_boot_p,
        Rcpp::Named("stringsAsFactors") = false
    );
}


// ================================================================
// Joint IRLS refinement of the peak coefficient.
//
// The one-step score beta (U / I_info) undershoots the MLE for large
// effects. This refines it with a few joint IRLS iterations of the full
// model [X_null | peak], warm-started at the one-step estimate, holding
// theta fixed (NegBin). Each iteration solves the (p_null + 1) system via
// a Schur complement so only the peak coefficient is profiled out.
// Mirrors fasthurdle's refine_beta_joint_ztnb, adapted to a single-part
// log-link GLM. 3 iterations suffice for <1% error; intended to be called
// only on significant pairs (|z| > 2).
// ================================================================
static double refine_beta_glm(double beta_init,
                              const vec& x,        // peak (n)
                              const vec& y,        // gene counts (n)
                              const mat& X_null,   // n x p_null
                              vec& beta_null,      // in/out, length p_null
                              int regr_type, double theta,
                              int max_iter = 3)
{
    const int p = (int)X_null.n_cols;
    const double th = (regr_type == 1) ? std::max(theta, 1e-2) : 0.0;
    double beta_test = beta_init;

    for (int iter = 0; iter < max_iter; iter++) {
        vec eta = X_null * beta_null + beta_test * x;
        eta = clamp(eta, -30.0, 30.0);   // guard against exp() overflow
        vec mu = exp(eta);
        vec W  = (regr_type == 0) ? mu : (mu % (th / (th + mu)));
        // Working response (log link, no offset): z = eta + (y - mu)/mu
        vec z  = eta + (y - mu) / mu;

        // Schur-complement solve for beta_test, profiling out beta_null.
        vec Wx   = W % x;
        mat WX   = X_null.each_col() % W;
        mat XtWX = X_null.t() * WX;          // p x p
        vec XtWx = X_null.t() * Wx;          // p
        vec XtWz = X_null.t() * (W % z);     // p

        mat XtWX_inv;
        if (!inv_sympd(XtWX_inv, XtWX)) {
            if (!inv(XtWX_inv, XtWX)) break;
        }
        double xWz = dot(Wx, z);
        double xWx = dot(Wx, x);
        double xt_inv_xz = dot(XtWx, XtWX_inv * XtWz);
        double xt_inv_xx = dot(XtWx, XtWX_inv * XtWx);
        double denom = xWx - xt_inv_xx;
        if (denom <= 0 || !std::isfinite(denom)) break;

        double beta_new = (xWz - xt_inv_xz) / denom;
        beta_null = XtWX_inv * (XtWz - XtWx * beta_new);

        if (std::abs(beta_new - beta_test) < 1e-6 * (std::abs(beta_test) + 1e-8)) {
            beta_test = beta_new;
            break;
        }
        beta_test = beta_new;
    }
    (void)p;
    return std::isfinite(beta_test) ? beta_test : beta_init;
}

// ================================================================
// SECTION 6: SCORE TEST (model-based score test + refined coefficient)
//
// Parallelised across genes (not pairs). For each gene, fits a null
// GLM (without peak term) via IRLS, then computes the model-based
// score statistic for every linked peak.
// ================================================================

// [[Rcpp::export]]
Rcpp::DataFrame ascent_score_pairs(
    Rcpp::S4 rna_sparse,
    Rcpp::S4 atac_sparse,
    Rcpp::IntegerVector gene_idx_r,
    Rcpp::IntegerVector peak_idx_r,
    Rcpp::StringVector gene_names_r,
    Rcpp::StringVector peak_names_r,
    arma::mat cov_mat,
    Rcpp::LogicalVector cell_mask,
    bool binarize,
    int regr_type,   // 0 = poisson (negbin: future)
    int ncores,
    double min_pct_rna = 0.05,
    double min_pct_atac = 0.05,
    bool skip_bootstrap = false)  // FALSE (default) = bootstrap every pair (like the Wald test)
{
    // ---- Extract CSC components from dgCMatrix ----
    Rcpp::IntegerVector rna_i_rv  = rna_sparse.slot("i");
    Rcpp::IntegerVector rna_p_rv  = rna_sparse.slot("p");
    Rcpp::NumericVector rna_x_rv  = rna_sparse.slot("x");
    Rcpp::IntegerVector atac_i_rv = atac_sparse.slot("i");
    Rcpp::IntegerVector atac_p_rv = atac_sparse.slot("p");
    Rcpp::NumericVector atac_x_rv = atac_sparse.slot("x");

    const int*    rna_ip  = rna_i_rv.begin();
    const int*    rna_pp  = rna_p_rv.begin();
    const double* rna_xp  = rna_x_rv.begin();
    const int*    atac_ip = atac_i_rv.begin();
    const int*    atac_pp = atac_p_rv.begin();
    const double* atac_xp = atac_x_rv.begin();

    int n_pairs = gene_idx_r.size();

    // ---- Copy R vectors to C++ (thread-safe) ----
    std::vector<int> gene_idx(n_pairs), peak_idx(n_pairs);
    std::vector<std::string> gene_names(n_pairs), peak_names(n_pairs);
    for (int i = 0; i < n_pairs; i++) {
        gene_idx[i]  = gene_idx_r(i);
        peak_idx[i]  = peak_idx_r(i);
        gene_names[i] = Rcpp::as<std::string>(gene_names_r(i));
        peak_names[i] = Rcpp::as<std::string>(peak_names_r(i));
    }

    // ---- Active cells (those passing celltype filter) ----
    std::vector<int> active_cells;
    int total_cells = cell_mask.size();
    for (int i = 0; i < total_cells; i++) {
        if (cell_mask(i)) active_cells.push_back(i);
    }
    int n_active = (int)active_cells.size();

    // ---- Pre-extract covariate sub-matrix for active cells ----
    int n_cov = (int)cov_mat.n_cols;
    mat cov_sub(n_active, n_cov);
    for (int i = 0; i < n_active; i++) {
        cov_sub.row(i) = cov_mat.row(active_cells[i]);
    }

    // ---- Build null design matrix: [intercept, covariates] ----
    int p_null = 1 + n_cov;
    mat X_null(n_active, p_null);
    X_null.col(0).ones();
    for (int j = 0; j < n_cov; j++) {
        X_null.col(1 + j) = cov_sub.col(j);
    }

    // ---- Group pairs by gene ----
    std::map<int, std::vector<int>> gene_to_pairs;
    for (int i = 0; i < n_pairs; i++) {
        gene_to_pairs[gene_idx[i]].push_back(i);
    }

    std::vector<int> unique_gene_rows;
    std::vector<std::vector<int>> gene_pair_lists;
    for (auto& kv : gene_to_pairs) {
        unique_gene_rows.push_back(kv.first);
        gene_pair_lists.push_back(std::move(kv.second));
    }
    int n_genes = (int)unique_gene_rows.size();

    // ---- Pre-allocate per-pair result arrays ----
    std::vector<int>    valid(n_pairs, 0);
    std::vector<double> r_beta(n_pairs, 0);        // refined coefficient estimate
    std::vector<double> r_se(n_pairs, 0), r_z(n_pairs, 0);           // model-based
    std::vector<double> r_score_U(n_pairs, 0);
    std::vector<double> r_score_V(n_pairs, 0);
    std::vector<double> r_score_T(n_pairs, 0), r_score_p(n_pairs, 1);         // model-based
    std::vector<double> r_boot_p(n_pairs, NA_REAL);  // bootstrap p (NA if not bootstrapped)

    // Adaptive bootstrap schedule (same as the Wald path)
    const int    boot_R[]      = {100,  500,  2500,  25000,  50000};
    const double boot_thresh[] = {0.1,  0.05, 0.01,  0.001,  0.0};
    const int    n_stages      = 5;

    // ---- Set OpenMP thread count ----
#ifdef _OPENMP
    if (ncores > 0) omp_set_num_threads(ncores);
#endif

    int progress_done = 0;

    Rprintf("ASCENT [score]: %d pairs (%d genes) across %d cells (%s, %d threads)...\n",
            n_pairs, n_genes, n_active,
            regr_type == 0 ? "Poisson" : "NegBin",
#ifdef _OPENMP
            ncores > 0 ? ncores : omp_get_max_threads()
#else
            1
#endif
    );

    // ================================================================
    // MAIN PARALLEL LOOP — one iteration per gene
    // ================================================================
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int gi = 0; gi < n_genes; gi++) {

        int gene_row = unique_gene_rows[gi];
        const std::vector<int>& pair_list = gene_pair_lists[gi];
        int n_peaks_for_gene = (int)pair_list.size();

        // Extract RNA row for this gene
        vec rna_vec(n_active);
        extract_sparse_row(rna_pp, rna_ip, rna_xp, gene_row,
                           active_cells, rna_vec);

        // Gene sparsity filter: require >min_pct_rna nonzero
        int n_expr = (int)accu(rna_vec > 0);
        if ((double)n_expr / n_active <= min_pct_rna) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            progress_done += n_peaks_for_gene;
            continue;
        }

        // Fit null model via IRLS (X_null is read-only, thread-safe)
        NullFitResult null_fit;
        if (regr_type == 0) {
            null_fit = poisson_null_irls(X_null, rna_vec);
        } else {
            null_fit = negbin_null_irls(X_null, rna_vec);
        }
        if (!null_fit.converged) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            progress_done += n_peaks_for_gene;
            continue;
        }

        // Compute weights and weighted residuals
        // Poisson: W = mu, r* = y - mu
        // NegBin:  W = mu*theta/(theta+mu), r* = (y-mu)*theta/(theta+mu)
        vec W, e_star;
        if (regr_type == 0) {
            W      = null_fit.mu;
            e_star = rna_vec - null_fit.mu;
        } else {
            double th = null_fit.theta;
            W      = null_fit.mu % (th / (th + null_fit.mu));
            e_star = (rna_vec - null_fit.mu) % (th / (th + null_fit.mu));
        }

        // For each peak linked to this gene
        for (int pair_idx : pair_list) {
            int pi = peak_idx[pair_idx];

            vec atac_vec(n_active);
            extract_sparse_row(atac_pp, atac_ip, atac_xp, pi,
                               active_cells, atac_vec);

            if (binarize) {
                atac_vec.elem(find(atac_vec > 0)).ones();
            }

            // Peak sparsity filter
            int n_open = (int)accu(atac_vec > 0);
            if ((double)n_open / n_active <= min_pct_atac) {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                progress_done++;
                continue;
            }

            // Compute b_tilde: WLS residual of atac on X_null
            //   Weights: Poisson W=mu, NegBin W=mu*theta/(theta+mu)
            //   XtWX_inv already uses the correct weights from null fit
            vec Wz = W % atac_vec;
            vec v = X_null.t() * Wz;
            vec gamma_proj = null_fit.XtWX_inv * v;
            vec b_tilde = atac_vec - X_null * gamma_proj;

            // Score statistic using weighted residuals
            double U = dot(b_tilde, e_star);

            // Model-based Fisher information: I = sum(b_tilde^2 * W)
            double I_info = dot(b_tilde % b_tilde, W);  // sum(b_tilde^2 * mu)
            if (I_info <= 0) {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                progress_done++;
                continue;
            }

            // One-step Newton-Raphson beta (warm start only): beta = U / I
            double beta_approx = U / I_info;

            // Model-based: se, z, score stat, p
            double se_model    = 1.0 / std::sqrt(I_info);
            double z_model     = U / std::sqrt(I_info);
            double score_T_mod = U * U / I_info;
            double score_p_mod = std::erfc(std::sqrt(score_T_mod / 2.0));

            // Coefficient estimate + bootstrap p-value.
            //  score_beta is ALWAYS the refined estimate (the one-step value is
            //  only the warm start; refinement early-exits in ~1 iteration when
            //  the one-step is already at the MLE, so refining every pair is
            //  cheap). With bootstrap on (default) the same adaptive bootstrap
            //  as the Wald test is then run over every pair, warm-started at the
            //  refined coefficients, giving boot_p. With bootstrap off, boot_p
            //  is NA.
            double beta_final  = beta_approx;
            double boot_p_val  = NA_REAL;

            // Refine beta (all pairs) -> reported estimate + bootstrap warm start.
            vec beta_null_local = null_fit.beta;
            {
                double br = refine_beta_glm(beta_approx, atac_vec, rna_vec,
                                            X_null, beta_null_local,
                                            regr_type, null_fit.theta);
                if (std::isfinite(br)) beta_final = br;
            }

            if (!skip_bootstrap) {
                // Full design [intercept, atac, covariates] + warm-start coefs.
                int p_dim = 2 + n_cov;
                mat Xf(n_active, p_dim);
                Xf.col(0).ones();
                Xf.col(1) = atac_vec;
                for (int j = 0; j < n_cov; j++) Xf.col(2 + j) = cov_sub.col(j);
                vec warm(p_dim);
                warm(0) = beta_null_local(0);                       // intercept
                warm(1) = beta_final;                               // atac (peak)
                for (int j = 0; j < n_cov; j++) warm(2 + j) = beta_null_local(1 + j);
                double th = null_fit.theta;

                // Adaptive bootstrap (per-pair reproducible RNG).
                std::mt19937 rng(42u + (unsigned)pair_idx);
                std::uniform_int_distribution<int> cell_dist(0, n_active - 1);
                boot_p_val = 1.0;
                for (int stage = 0; stage < n_stages; stage++) {
                    if (stage > 0 && boot_p_val >= boot_thresh[stage - 1]) break;
                    int R = boot_R[stage];
                    std::vector<double> bcs; bcs.reserve(R);
                    vec freq(n_active);
                    for (int r = 0; r < R; r++) {
                        freq.zeros();
                        for (int c = 0; c < n_active; c++) freq(cell_dist(rng))++;
                        double bc = (regr_type == 0)
                            ? poisson_boot_coef(Xf, rna_vec, freq, warm)
                            : negbin_boot_coef(Xf, rna_vec, freq, warm, th);
                        if (!std::isnan(bc)) bcs.push_back(bc);
                    }
                    if (!bcs.empty()) boot_p_val = basic_p_cpp(beta_final, bcs);
                }
            }

            valid[pair_idx]         = 1;
            r_beta[pair_idx]        = beta_final;
            r_boot_p[pair_idx]      = boot_p_val;
            r_se[pair_idx]          = se_model;
            r_z[pair_idx]           = z_model;
            r_score_U[pair_idx]     = U;
            r_score_V[pair_idx]     = I_info;
            r_score_T[pair_idx]     = score_T_mod;
            r_score_p[pair_idx]     = score_p_mod;

#ifdef _OPENMP
            #pragma omp atomic
#endif
            progress_done++;

            // Progress reporting
            if (progress_done % std::max(1, n_pairs / 20) == 0) {
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    Rprintf("\r  [%d / %d pairs processed (%.0f%%)]",
                            progress_done, n_pairs,
                            100.0 * progress_done / n_pairs);
                }
            }
        } // end peak loop
    } // end gene loop

    Rprintf("\r  [%d / %d pairs processed (100%%)]\n", n_pairs, n_pairs);

    // ---- Collect valid results into R DataFrame ----
    int n_valid = 0;
    for (int i = 0; i < n_pairs; i++) n_valid += valid[i];

    Rcpp::StringVector  out_gene(n_valid), out_peak(n_valid);
    Rcpp::NumericVector out_beta(n_valid);
    Rcpp::NumericVector out_se(n_valid), out_z(n_valid);
    Rcpp::NumericVector out_U(n_valid), out_V(n_valid);
    Rcpp::NumericVector out_T(n_valid), out_p(n_valid);
    Rcpp::NumericVector out_boot_p(n_valid);

    int idx = 0;
    for (int i = 0; i < n_pairs; i++) {
        if (!valid[i]) continue;
        out_gene[idx]     = gene_names[i];
        out_peak[idx]     = peak_names[i];
        out_beta[idx]     = r_beta[i];
        out_se[idx]     = r_se[i];
        out_z[idx]      = r_z[i];
        out_U[idx]      = r_score_U[i];
        out_V[idx]      = r_score_V[i];
        out_T[idx]      = r_score_T[i];
        out_p[idx]      = r_score_p[i];
        out_boot_p[idx] = r_boot_p[i];
        idx++;
    }

    Rprintf("ASCENT [score]: Done. %d / %d pairs passed quality filters.\n",
            n_valid, n_pairs);

    return Rcpp::DataFrame::create(
        Rcpp::Named("gene")                = out_gene,
        Rcpp::Named("peak")                = out_peak,
        Rcpp::Named("score_beta")          = out_beta,      // refined / MLE estimate
        Rcpp::Named("score_se")            = out_se,
        Rcpp::Named("score_z")             = out_z,
        Rcpp::Named("score_p")             = out_p,
        Rcpp::Named("boot_p")              = out_boot_p,    // NA unless bootstrapped
        Rcpp::Named("score_U")             = out_U,
        Rcpp::Named("score_V")             = out_V,
        Rcpp::Named("score_stat")          = out_T,
        Rcpp::Named("stringsAsFactors")    = false
    );
}
