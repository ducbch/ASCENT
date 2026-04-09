# ASCENT: Technical Deep Dive

ASCENT (Accelerated Single-Cell ENhancer Target) is a high-performance reimplementation of the [SCENT](https://github.com/immunogenomics/SCENT) algorithm in Rcpp/C++ with OpenMP parallelism. This document details the computational bottlenecks in the original R implementation and the optimizations introduced in ASCENT.

---

## 1. Overview of the SCENT Algorithm

SCENT tests for statistical association between chromatin accessibility (ATAC-seq peaks) and gene expression (RNA) across single cells. For each peak-gene pair:

1. **Fit a GLM** (Poisson or Negative Binomial) regressing gene expression on peak accessibility, controlling for covariates (e.g., total UMI counts).
2. **Bootstrap** the coefficient estimate with an adaptive resampling schedule:
   - 100 replicates; if p < 0.1, escalate to 500
   - 500 replicates; if p < 0.05, escalate to 2,500
   - 2,500 replicates; if p < 0.01, escalate to 25,000
   - 25,000 replicates; if p < 0.001, escalate to 50,000
3. **Compute an empirical p-value** from the bootstrap distribution using the basic bootstrap method.

The computational cost is dominated by the bootstrap: a single pair may require up to 50,000 GLM fits. Genome-wide analyses involve tens of thousands of pairs.

---

## 2. Bottlenecks in the Original R Implementation

### 2.1 Sequential pair processing

The original SCENT iterates over peak-gene pairs in a `for` loop. Each pair is fully processed (initial fit + all bootstrap stages) before moving to the next. There is no parallelism across pairs.

```r
for (n in 1:nrow(peak.info)) {
  # ... extract data, fit GLM, run boot::boot ...
}
```

### 2.2 Fork-based bootstrap parallelism

SCENT uses `boot::boot(parallel = "multicore", ncpus = N)` to parallelize bootstrap replicates within a single pair. This calls `parallel::mclapply`, which `fork()`s the R process N times.

**Problems:**

- **Memory duplication.** Each forked child gets a copy-on-write clone of the entire R session (Seurat object, sparse matrices, metadata). As children modify working memory during GLM fitting, the OS must copy modified pages. With 8 cores, observed peak RSS is ~8x the base process size (~3.3 GB base becomes ~26 GB).
- **Fork overhead per pair.** Forking is repeated for every pair at every bootstrap stage. For 1,543 pairs with an average of ~2 bootstrap stages each, that is ~3,000 fork operations.
- **Diminishing returns.** The bootstrap replicates within a single pair are fast individually (one GLM fit each). The fork/join overhead can exceed the computation time, especially at the 100-replicate stage.

### 2.3 Data copying per bootstrap replicate

`boot::boot` with `stype = "i"` generates an index vector and calls `statistic(data[idx, ])` for each replicate. This subsets the data frame every time, creating a new copy of the design matrix, response vector, and all columns.

For a pair requiring 50,000 replicates on 2,286 cells with 4 predictors, that is 50,000 data frame subset operations, each allocating and populating a new object.

### 2.4 Cold-started GLM fits

Each bootstrap replicate calls `glm()` (or `fastglm::fastglm()`) from scratch. The IRLS algorithm starts from a default initialization (mu = y + 0.1), typically requiring 15-25 iterations to converge. There is no mechanism to reuse the converged coefficients from the initial fit as a starting point.

### 2.5 R interpreter overhead in the IRLS loop

R's `glm()` function wraps LAPACK/BLAS linear algebra in several hundred lines of R code:

- Formula parsing and `model.matrix()` construction
- Family object method dispatch (`$linkfun`, `$variance`, `$mu.eta`, `$dev.resids`)
- R-level convergence checking
- `summary.glm()` for variance-covariance extraction

The underlying QR decomposition is compiled (LAPACK), but everything around it is interpreted R. This overhead is negligible for a single fit but becomes significant when multiplied by up to 50,000 replicates per pair across thousands of pairs.

### 2.6 Sparse matrix access via R's S4 dispatch

Extracting a row from a `dgCMatrix` (CSC sparse matrix) in R, e.g., `object@atac[peak_name, ]`, triggers:

- S4 method dispatch for `[`
- Internal conversion logic for the CSC-to-row extraction
- Allocation of a new dense vector
- Name assignment for the result

This is repeated for both the RNA and ATAC matrix for every pair.

### 2.7 Negative Binomial theta estimation

For `regr = "negbin"`, SCENT uses `MASS::glm.nb()`, which internally alternates between IRLS and theta estimation. Each bootstrap replicate re-estimates theta from scratch, even though the dispersion parameter is relatively stable across bootstrap samples of the same pair.

---

## 3. ASCENT Optimizations

### 3.1 Pair-level OpenMP parallelism

ASCENT parallelizes across **pairs**, not within bootstrap replicates:

```cpp
#pragma omp parallel for schedule(dynamic)
for (int pair = 0; pair < n_pairs; pair++) {
    // ... extract data, fit IRLS, run all bootstrap stages ...
}
```

**Advantages:**

- **No fork overhead.** OpenMP threads are created once at loop entry. Thread creation cost is amortized across all pairs.
- **Dynamic scheduling.** Pairs with significant p-values require more bootstrap stages (up to 50,000 replicates). `schedule(dynamic)` ensures threads that finish fast pairs pick up new work immediately, avoiding load imbalance.
- **Shared memory.** All threads read from the same CSC arrays and covariate matrix. Only per-thread working vectors (design matrix, frequency vector, IRLS temporaries) are thread-local.

**Memory impact:** With 8 threads, ASCENT uses ~3.3 GB total. The R implementations with 8 forked processes use ~26 GB.

### 3.2 Frequency-weighted bootstrap

Instead of resampling rows of the data frame, ASCENT generates a **multinomial frequency vector**:

```cpp
vec freq(n_active);
freq.zeros();
for (int c = 0; c < n_active; c++) freq(cell_dist(rng))++;
```

This frequency vector serves as case weights in the IRLS. The design matrix `X` and response `y` are never copied — only the weight vector changes per replicate.

**Impact:** Eliminates 50,000 data frame copies per pair. The frequency vector is `n_cells` doubles (~18 KB for 2,286 cells), reused in-place.

### 3.3 Warm-started IRLS

Each bootstrap replicate initializes from the initial fit's converged coefficients:

```cpp
double poisson_boot_coef(const mat& X, const vec& y, const vec& freq,
                         const vec& beta_warm, ...) {
    vec eta = X * beta_warm;  // start from converged solution
    ...
}
```

Since bootstrap samples are perturbations of the original data, the converged coefficients are a good starting point. Typical convergence: **2-5 iterations** instead of 15-25 from a cold start.

**Impact:** ~5x fewer IRLS iterations per bootstrap replicate.

### 3.4 Direct CSC sparse row extraction

ASCENT extracts rows from `dgCMatrix` by directly operating on the CSC arrays (`@i`, `@p`, `@x`), using binary search within each column's index range:

```cpp
inline void extract_sparse_row(const int* p, const int* i, const double* x,
                                int row_idx, const std::vector<int>& active_cells,
                                vec& out) {
    for (int k = 0; k < n; k++) {
        int col = active_cells[k];
        const int* pos = std::lower_bound(i + p[col], i + p[col+1], row_idx);
        if (pos != i + p[col+1] && *pos == row_idx)
            out(k) = x[pos - i];
    }
}
```

**Advantages:**

- No S4 method dispatch
- No memory allocation (output vector is pre-allocated and reused)
- Only scans active cells (those passing the cell-type filter), not all columns
- Binary search is O(log nnz) per column vs. linear scan

### 3.5 Compiled IRLS loop

The entire IRLS iteration — weight computation, working response, QR-based weighted least squares, eta update, convergence check — is a tight C++ loop with no R API calls:

```cpp
for (int iter = 0; iter < max_iter; iter++) {
    vec w = clamp(mu, 1e-10, 1e7);
    vec z_w = eta + (y - mu) / mu;
    if (!qr_wls_solve(beta_new, X, w, z_w)) return res;
    vec eta_new = X * beta_new;
    eta_new = clamp(eta_new, -30.0, 30.0);
    vec mu_new = exp(eta_new);
    double change = accu(abs(eta_new - eta)) / n;
    beta = beta_new; eta = eta_new; mu = mu_new;
    if (change < tol) break;
}
```

The QR decomposition itself uses the same LAPACK routines as R (via Armadillo's `solve()`), but removing the R-level wrapper overhead around each iteration accumulates significant savings over millions of total IRLS iterations.

### 3.6 Fixed-theta Negative Binomial bootstrap

For Negative Binomial regression, ASCENT estimates theta once from the initial fit (alternating IRLS + Newton-Raphson), then **fixes theta** for all bootstrap replicates:

```cpp
double negbin_boot_coef(..., double theta, ...) {
    // theta is fixed — only beta is re-estimated
    vec w = freq % mu % (theta / (theta + mu));
    ...
}
```

**Why this is statistically valid:**

- **Theta is a nuisance parameter.** The bootstrap is testing whether beta (the peak-gene association) is significantly different from zero. Theta controls the overall overdispersion of the gene's count distribution across all cells — it is a property of the gene, not of the peak-gene relationship. Resampling cells doesn't change the gene's underlying dispersion.
- **Standard practice.** This is the same approach used by `DESeq2` and `edgeR` — estimate dispersion once on the full data, then hold it fixed for all downstream inference. The full-data estimate is more stable and accurate than per-resample estimates.
- **Re-estimation adds noise, not accuracy.** Theta estimated from a bootstrap resample has high variance, especially on sparse single-cell data. This noisy theta propagates into noisier beta estimates, making the bootstrap distribution artificially wider without improving the test's validity.

**Why SCENT re-estimates theta:** The original SCENT uses `MASS::glm.nb()`, which jointly estimates both beta and theta in a single function call. This was the path of least resistance — fixing theta would require extracting it from the initial fit, storing it, and passing it into a separate `glm(..., family = negative.binomial(theta))` call for each bootstrap replicate. More code and bookkeeping for the same statistical conclusion. The computational penalty of re-estimation was not apparent at the scale SCENT was originally designed for.

**Impact:** For a pair requiring 50,000 replicates, SCENT's `glm.nb()` runs ~50-100 IRLS iterations per replicate (alternating beta and theta). ASCENT's fixed-theta bootstrap runs ~2-3 iterations per replicate (warm-started, beta only). This is the single largest source of speedup for NegBin regression.

### 3.7 Thread-safe RNG with reproducible seeds

Each pair gets a deterministic seed (`42 + pair_index`), ensuring reproducibility regardless of thread count:

```cpp
std::mt19937 rng(42u + (unsigned)pair);
```

This avoids R's single-threaded RNG (which cannot be called from OpenMP threads) and ensures results are identical whether run with 1 or 8 threads.

### 3.8 Pre-subsetting to relevant features

Before entering C++, the R wrapper subsets the sparse matrices to only the genes and peaks present in the pair list:

```r
rna_sub  <- object@rna[genes_needed, , drop = FALSE]
atac_sub <- object@atac[peaks_needed, , drop = FALSE]
```

This reduces the size of the CSC arrays passed to C++, improving cache locality and reducing binary search ranges. The original SCENT accesses the full matrix for every pair.

---

## 4. Benchmark Results

**Test configuration:** 1,800 peak-gene pairs stratified by distance to TSS (0-5kb, 5-20kb, 20-50kb) and sparsity (Borderline 5-10%, Sparse 10-25%, Dense 25-100%), 2,286 CD14 monocyte cells, 8 cores.

### Runtime — Poisson (1,800 pairs)

| Method  | Elapsed (sec) | Speedup vs glm |
|---------|--------------|----------------|
| glm     | 12,918       | 1.0x (baseline)|
| fastglm | 7,200        | 1.8x           |
| rcpp    | 597          | **21.6x**      |

### Runtime — Negative Binomial (32 pairs)

| Method  | Elapsed (sec) | Speedup vs glm |
|---------|--------------|----------------|
| glm     | 6,505        | 1.0x (baseline)|
| fastglm | 341          | 19.1x          |
| rcpp    | 82           | **79.3x**      |

The much larger NegBin speedup is primarily due to ASCENT fixing theta (dispersion) during bootstrap rather than re-estimating it on every replicate (see Section 3.6).

### Peak Memory (MaxRSS)

ASCENT uses OpenMP shared-memory threads; memory stays constant regardless of core count. The R methods use `fork()`-based parallelism (`mclapply`), so memory scales linearly with the number of cores — each forked child gets a copy-on-write clone of the entire R session.

| Method  | Memory (GB) | Note |
|---------|-------------|------|
| glm     | 25.8        | ~3.3 GB base x 8 forked processes |
| fastglm | 25.5        | ~3.3 GB base x 8 forked processes |
| rcpp    | 3.3         | constant regardless of core count |

### Numerical Accuracy

| Comparison       | Max absolute beta diff | Max absolute SE diff | Correlation |
|------------------|-----------------|---------------|-------------|
| rcpp vs glm      | 2.24e-08        | 1.54e-05      | 1.000000    |
| fastglm vs glm   | 3.36e-14        | 7.66e-15      | 1.000000    |

Bootstrap p-values differ due to independent RNG streams but show high Spearman correlation across all strata, with >95% concordance on significance calls at alpha = 0.05.

---

## 5. Summary of Improvements

| Bottleneck | SCENT (R) | ASCENT (Rcpp) |
|-----------|-----------|---------------|
| Parallelism | Fork-based within bootstrap | OpenMP across pairs |
| Memory model | N forked processes (copy-on-write) | Shared address space (threads) |
| Bootstrap data | Full data frame copy per replicate | Frequency vector as case weights |
| IRLS initialization | Cold start every replicate | Warm start from initial fit |
| IRLS loop | R interpreter + S4 dispatch | Compiled C++ loop |
| Sparse access | R's `[` with S4 dispatch | Direct CSC binary search |
| NB theta | Re-estimated every replicate | Fixed from initial fit |
| RNG | R's single-threaded RNG via fork | Per-pair `std::mt19937` |
| Matrix scope | Full genome-wide matrices | Pre-subsetted to relevant features |
