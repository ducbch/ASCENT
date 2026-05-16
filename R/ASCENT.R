#' @import methods
#' @import Matrix
#' @importFrom Rcpp sourceCpp
#' @importFrom data.table fread
#' @importFrom stringr str_split_fixed
#' @importFrom Hmisc cut2
#' @useDynLib ASCENT, .registration = TRUE
NULL


# ================================================================
# BOOTSTRAP P-VALUE HELPERS (used by fastglm fallback path)
# ================================================================

#' Interpolated two-sided p-value from bootstrap distribution
#'
#' Computes an empirical p-value by finding where zero falls in the sorted
#' bootstrap distribution and returning twice the smaller tail proportion.
#'
#' @param q Numeric vector. Bootstrap distribution of the test quantity
#'   (typically \code{2 * obs - boot - null}).
#'
#' @return Numeric scalar. Two-sided empirical p-value.
#' @noRd
interp_pval <- function(q) {
  R <- length(q)
  tstar <- sort(q)
  zero <- findInterval(0, tstar)
  if (zero == 0 || zero == R) return(2 / R)
  2 * min(zero / R, (R - zero) / R)
}

#' Basic bootstrap p-value
#'
#' Computes the basic bootstrap p-value following the method of Davison &
#' Hinkley (1997). Transforms the bootstrap distribution using
#' \code{2 * obs - boot - null} and passes it to \code{\link{interp_pval}}.
#'
#' @param obs Numeric scalar. Observed statistic (e.g. ATAC coefficient).
#' @param boot Numeric vector. Bootstrap replicates of the statistic.
#' @param null Numeric scalar. Null hypothesis value (default 0).
#'
#' @return Numeric scalar. Two-sided basic bootstrap p-value.
#' @noRd
basic_p <- function(obs, boot, null = 0) {
  boot <- boot[!is.na(boot)]
  if (length(boot) == 0) return(NA_real_)
  interp_pval(2 * obs - boot - null)
}


# ================================================================
# FASTGLM HELPERS
# ================================================================

#' Newton-Raphson update for negative binomial dispersion (theta)
#'
#' R implementation of the C++ \code{update_theta} function. Given current
#' fitted values \code{mu} from an NB GLM, updates the dispersion parameter
#' \code{theta} using score/information Newton-Raphson steps.
#'
#' @param y Numeric vector. Response (gene expression counts).
#' @param mu Numeric vector. Current fitted values from NB GLM.
#' @param theta Numeric scalar. Current theta estimate.
#' @param max_iter Integer. Maximum Newton-Raphson iterations (default 10).
#'
#' @return Numeric scalar. Updated theta, clamped to [0.01, 1e6].
#' @noRd
.update_theta_r <- function(y, mu, theta, max_iter = 10) {
  for (i in seq_len(max_iter)) {
    ym <- y + theta
    tm <- theta + mu
    score <- sum(digamma(ym) - digamma(theta) + log(theta) - log(tm) + 1 - ym / tm)
    info  <- sum(-trigamma(ym) + trigamma(theta) - 1/theta + 2/tm - ym / tm^2)
    if (abs(info) < 1e-15) break
    theta_new <- theta + score / info
    theta_new <- max(min(theta_new, 1e6), 0.01)
    if (abs(theta_new - theta) / (theta + 0.001) < 1e-6) {
      theta <- theta_new; break
    }
    theta <- theta_new
  }
  theta
}

#' Fit negative binomial GLM via alternating fastglm IRLS and theta NR
#'
#' Alternates between fitting beta via \code{\link[fastglm]{fastglm}} with a
#' fixed \code{MASS::negative.binomial(theta)} family, and updating theta via
#' \code{\link{.update_theta_r}}. Initialises from a Poisson fit.
#'
#' @param X Numeric matrix. Design matrix (with intercept column).
#' @param y Numeric vector. Response (gene expression counts).
#' @param max_outer Integer. Maximum outer iterations (default 25).
#'
#' @return A fastglm fit object with an additional \code{$theta} element.
#' @noRd
.negbin_fastglm <- function(X, y, max_outer = 25) {
  # Poisson initialisation for mu and theta
  fit <- fastglm::fastglm(x = X, y = y, family = poisson())
  mu <- fit$fitted.values
  ss <- sum((y / pmax(mu, 1e-10) - 1)^2 - 1 / pmax(mu, 1e-10))
  theta <- max(length(y) / max(ss, 0.01), 0.01)
  theta <- min(theta, 1e6)

  for (outer in seq_len(max_outer)) {
    theta_old <- theta
    fit <- fastglm::fastglm(x = X, y = y,
                             family = MASS::negative.binomial(theta))
    mu <- fit$fitted.values
    theta <- .update_theta_r(y, mu, theta)
    if (abs(theta - theta_old) / (theta_old + 0.001) < 1e-4) break
  }
  fit$theta <- theta
  fit
}

#' Bootstrap statistic: Poisson association via fastglm
#'
#' Returns the ATAC coefficient and its variance from a Poisson GLM fit
#' with \code{\link[fastglm]{fastglm}}. Designed for use with
#' \code{\link[boot]{boot}}.
#'
#' @param data Data frame with columns \code{exprs}, \code{atac}, and covariates.
#' @param idx Integer vector. Bootstrap indices (default: all rows).
#' @param x_names Character vector. Predictor column names (atac + covariates).
#'
#' @return Numeric vector of length 2: \code{c(coef_atac, var_atac)}.
#' @noRd
.assoc_poisson_fast <- function(data, idx = seq_len(nrow(data)), x_names) {
  d <- data[idx, , drop = FALSE]
  X <- cbind(Intercept = 1, as.matrix(d[, x_names, drop = FALSE]))
  y <- d$exprs
  gg <- fastglm::fastglm(x = X, y = y, family = poisson())
  c(gg$coefficients[2], gg$se[2]^2)
}

#' Bootstrap statistic: negative binomial association via fastglm
#'
#' Returns the ATAC coefficient and its variance from a NegBin GLM fit
#' with \code{\link[fastglm]{fastglm}} using a fixed dispersion parameter
#' \code{theta}. Designed for use with \code{\link[boot]{boot}}.
#'
#' @param data Data frame with columns \code{exprs}, \code{atac}, and covariates.
#' @param idx Integer vector. Bootstrap indices (default: all rows).
#' @param x_names Character vector. Predictor column names (atac + covariates).
#' @param theta Numeric scalar. Fixed NB dispersion parameter from initial fit.
#'
#' @return Numeric vector of length 2: \code{c(coef_atac, var_atac)}.
#' @noRd
.assoc_negbin_fast <- function(data, idx = seq_len(nrow(data)), x_names, theta) {
  d <- data[idx, , drop = FALSE]
  X <- cbind(Intercept = 1, as.matrix(d[, x_names, drop = FALSE]))
  y <- d$exprs
  gg <- fastglm::fastglm(x = X, y = y,
                          family = MASS::negative.binomial(theta))
  c(gg$coefficients[2], gg$se[2]^2)
}


# ================================================================
# FASTGLM FALLBACK ALGORITHM (pure R, mirrors original SCENT logic)
# ================================================================

#' ASCENT algorithm — fastglm fallback path
#'
#' Pure-R implementation using \pkg{fastglm} for GLM fitting and
#' \code{\link[boot]{boot}} with multicore parallelism for bootstrap.
#' Mirrors the original SCENT logic but replaces \code{glm()} with
#' \code{fastglm::fastglm()} for ~2x speedup. For NegBin, uses
#' \code{\link{.negbin_fastglm}} for the initial fit and fixes theta
#' during bootstrap.
#'
#' @param object ASCENT object.
#' @param celltype Character. Cell type to analyse.
#' @param ncores Integer. Number of bootstrap cores.
#' @param regr Character. \code{"poisson"} or \code{"negbin"}.
#' @param bin Logical. Binarise ATAC counts?
#'
#' @return Data frame of results (gene, peak, beta, se, z, p, boot_basic_p).
#' @noRd
.ASCENT_algorithm_fastglm <- function(object, celltype, ncores, regr, bin, bootstrap = TRUE) {
  if (!requireNamespace("fastglm", quietly = TRUE))
    stop("Install fastglm: install.packages('fastglm')")
  if (regr == "negbin" && !requireNamespace("MASS", quietly = TRUE))
    stop("Install MASS for negbin: install.packages('MASS')")
  if (bootstrap && !requireNamespace("boot", quietly = TRUE))
    stop("Install boot: install.packages('boot')")

  res <- data.frame()
  n_pairs <- nrow(object@peak.info)

  for (n in seq_len(n_pairs)) {
    gene      <- object@peak.info[n, 1]
    this_peak <- object@peak.info[n, 2]

    # Extract ATAC for this peak
    atac_target <- data.frame(
      cell = colnames(object@atac),
      atac = as.numeric(object@atac[this_peak, ])
    )
    if (bin && any(atac_target$atac > 0)) {
      atac_target$atac[atac_target$atac > 0] <- 1
    }

    # Extract RNA for this gene
    mrna_target <- object@rna[gene, ]
    df <- data.frame(cell = names(mrna_target), exprs = as.numeric(mrna_target))
    df <- merge(df, atac_target, by = "cell")
    df <- merge(df, object@meta.data, by = "cell")

    # Filter to celltype
    df2 <- df[df[[object@celltypes]] == celltype, ]
    if (nrow(df2) == 0) next

    nonzero_m <- mean(df2$exprs > 0)
    nonzero_a <- mean(df2$atac > 0)
    if (nonzero_m <= 0.05 || nonzero_a <= 0.05) next

    # Predictor names for design matrix
    pred_var <- c("atac", object@covariates)

    # Build design matrix for initial fit
    X_init <- cbind(Intercept = 1, as.matrix(df2[, pred_var, drop = FALSE]))
    y_init <- df2$exprs

    # ---- Initial fit ----
    base_fit <- tryCatch({
      if (regr == "poisson") {
        fastglm::fastglm(x = X_init, y = y_init, family = poisson())
      } else {
        .negbin_fastglm(X_init, y_init)
      }
    }, error = function(e) NULL)

    if (is.null(base_fit)) next
    coef_atac <- base_fit$coefficients[2]
    se_atac   <- base_fit$se[2]
    if (is.na(coef_atac) || is.na(se_atac) || se_atac <= 0) next

    if (regr == "poisson") {
      assoc_fn   <- .assoc_poisson_fast
      boot_extra <- list(x_names = pred_var)
    } else {
      assoc_fn   <- .assoc_negbin_fast
      boot_extra <- list(x_names = pred_var, theta = base_fit$theta)
    }

    z_atac <- coef_atac / se_atac
    p_atac <- 2 * pnorm(-abs(z_atac))
    coefs  <- c(coef_atac, se_atac, z_atac, p_atac)

    # ---- Iterative bootstrap (same schedule as original SCENT) ----
    if (bootstrap) {
      boot_args <- c(list(data = df2, statistic = assoc_fn, stype = "i",
                          parallel = "multicore", ncpus = ncores),
                     boot_extra)

      bs <- do.call(boot::boot, c(boot_args, list(R = 100)))
      p0 <- basic_p(bs$t0[1], bs$t[, 1])

      if (p0 < 0.1) {
        bs <- do.call(boot::boot, c(boot_args, list(R = 500)))
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.05) {
        bs <- do.call(boot::boot, c(boot_args, list(R = 2500)))
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.01) {
        bs <- do.call(boot::boot, c(boot_args, list(R = 25000)))
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.001) {
        bs <- do.call(boot::boot, c(boot_args, list(R = 50000)))
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
    } else {
      p0 <- NA_real_
    }

    out <- data.frame(gene = gene, peak = this_peak,
                      beta = coefs[1], se = coefs[2],
                      z = coefs[3], p = coefs[4],
                      boot_basic_p = p0)
    res <- rbind(res, out)

    if (n %% max(1, n_pairs %/% 20) == 0) {
      message(sprintf("  [%d / %d pairs (%.0f%%)]", n, n_pairs, 100 * n / n_pairs))
    }
  }

  res
}


# ================================================================
# NATIVE GLM FALLBACK (identical to original SCENT logic)
# ================================================================

#' Bootstrap statistic: Poisson association via base R glm
#'
#' Identical to the original SCENT \code{assoc_poisson} function.
#' Returns the ATAC coefficient and its variance from a Poisson GLM fit
#' with \code{\link[stats]{glm}}. Designed for use with
#' \code{\link[boot]{boot}}.
#'
#' @param data Data frame with columns \code{exprs}, \code{atac}, and covariates.
#' @param idx Integer vector. Bootstrap indices (default: all rows).
#' @param formula Formula. GLM formula (e.g. \code{exprs ~ atac + log_nUMI}).
#'
#' @return Numeric vector of length 2: \code{c(coef_atac, var_atac)}.
#' @noRd
.assoc_poisson_glm <- function(data, idx = seq_len(nrow(data)), formula) {
  gg <- glm(formula, family = "poisson", data = data[idx, , drop = FALSE])
  c(coef(gg)["atac"], diag(vcov(gg))["atac"])
}

#' Bootstrap statistic: negative binomial association via MASS::glm.nb
#'
#' Identical to the original SCENT \code{assoc_negbin} function.
#' Returns the ATAC coefficient and its variance from a NegBin GLM fit
#' with \code{\link[MASS]{glm.nb}}. Unlike the fastglm path, theta is
#' re-estimated jointly with beta on every bootstrap replicate.
#' Designed for use with \code{\link[boot]{boot}}.
#'
#' @param data Data frame with columns \code{exprs}, \code{atac}, and covariates.
#' @param idx Integer vector. Bootstrap indices (default: all rows).
#' @param formula Formula. GLM formula (e.g. \code{exprs ~ atac + log_nUMI}).
#'
#' @return Numeric vector of length 2: \code{c(coef_atac, var_atac)}.
#' @noRd
.assoc_negbin_glm <- function(data, idx = seq_len(nrow(data)), formula) {
  tryCatch({
    gg <- MASS::glm.nb(formula, data = data[idx, , drop = FALSE])
    c(coef(gg)["atac"], diag(vcov(gg))["atac"])
  }, error = function(e) c(NA_real_, NA_real_))
}

#' ASCENT algorithm — native glm fallback path
#'
#' Pure-R implementation identical to the original SCENT algorithm.
#' Uses \code{\link[stats]{glm}} for Poisson and \code{\link[MASS]{glm.nb}}
#' for negative binomial. Bootstrap via \code{\link[boot]{boot}} with
#' multicore parallelism. Slowest option — useful as a ground-truth reference.
#'
#' @param object ASCENT object.
#' @param celltype Character. Cell type to analyse.
#' @param ncores Integer. Number of bootstrap cores.
#' @param regr Character. \code{"poisson"} or \code{"negbin"}.
#' @param bin Logical. Binarise ATAC counts?
#'
#' @return Data frame of results (gene, peak, beta, se, z, p, boot_basic_p).
#' @noRd
.ASCENT_algorithm_glm <- function(object, celltype, ncores, regr, bin, bootstrap = TRUE) {
  if (regr == "negbin" && !requireNamespace("MASS", quietly = TRUE))
    stop("Install MASS for negbin: install.packages('MASS')")
  if (bootstrap && !requireNamespace("boot", quietly = TRUE))
    stop("Install boot: install.packages('boot')")

  res <- data.frame()
  n_pairs <- nrow(object@peak.info)

  for (n in seq_len(n_pairs)) {
    gene      <- object@peak.info[n, 1]
    this_peak <- object@peak.info[n, 2]

    atac_target <- data.frame(
      cell = colnames(object@atac),
      atac = as.numeric(object@atac[this_peak, ])
    )
    if (bin && any(atac_target$atac > 0)) {
      atac_target$atac[atac_target$atac > 0] <- 1
    }

    mrna_target <- object@rna[gene, ]
    df <- data.frame(cell = names(mrna_target), exprs = as.numeric(mrna_target))
    df <- merge(df, atac_target, by = "cell")
    df <- merge(df, object@meta.data, by = "cell")

    df2 <- df[df[[object@celltypes]] == celltype, ]
    if (nrow(df2) == 0) next

    nonzero_m <- mean(df2$exprs > 0)
    nonzero_a <- mean(df2$atac > 0)
    if (nonzero_m <= 0.05 || nonzero_a <= 0.05) next

    # Build formula: exprs ~ atac + cov1 + cov2 + ...
    pred_var <- c("atac", object@covariates)
    formula  <- as.formula(paste("exprs", paste(pred_var, collapse = "+"), sep = "~"))

    # ---- Initial fit ----
    base_fit <- tryCatch({
      if (regr == "poisson") {
        glm(formula, family = "poisson", data = df2)
      } else {
        MASS::glm.nb(formula, data = df2)
      }
    }, error = function(e) NULL)

    if (is.null(base_fit)) next
    coefs_tbl <- summary(base_fit)$coefficients
    if (!"atac" %in% rownames(coefs_tbl)) next
    coefs <- coefs_tbl["atac", ]  # Estimate, Std. Error, z value, Pr(>|z|)

    # ---- Iterative bootstrap (same schedule as original SCENT) ----
    if (bootstrap) {
      assoc_fn <- if (regr == "poisson") .assoc_poisson_glm else .assoc_negbin_glm

      bs <- boot::boot(df2, assoc_fn, R = 100, formula = formula,
                       stype = "i", parallel = "multicore", ncpus = ncores)
      p0 <- basic_p(bs$t0[1], bs$t[, 1])

      if (p0 < 0.1) {
        bs <- boot::boot(df2, assoc_fn, R = 500, formula = formula,
                         stype = "i", parallel = "multicore", ncpus = ncores)
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.05) {
        bs <- boot::boot(df2, assoc_fn, R = 2500, formula = formula,
                         stype = "i", parallel = "multicore", ncpus = ncores)
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.01) {
        bs <- boot::boot(df2, assoc_fn, R = 25000, formula = formula,
                         stype = "i", parallel = "multicore", ncpus = ncores)
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
      if (p0 < 0.001) {
        bs <- boot::boot(df2, assoc_fn, R = 50000, formula = formula,
                         stype = "i", parallel = "multicore", ncpus = ncores)
        p0 <- basic_p(bs$t0[1], bs$t[, 1])
      }
    } else {
      p0 <- NA_real_
    }

    out <- data.frame(gene = gene, peak = this_peak,
                      beta = coefs[1], se = coefs[2],
                      z = coefs[3], p = coefs[4],
                      boot_basic_p = p0)
    res <- rbind(res, out)

    if (n %% max(1, n_pairs %/% 20) == 0) {
      message(sprintf("  [%d / %d pairs (%.0f%%)]", n, n_pairs, 100 * n / n_pairs))
    }
  }

  res
}


# ================================================================
# SCORE TEST ALGORITHM (Poisson, HC0 sandwich)
# ================================================================

#' Score test with HC0 sandwich variance (R path)
#'
#' For each gene, fits a null GLM (without peak term), then computes the
#' robust score statistic for each linked peak using the HC0 sandwich
#' estimator. Much faster than Wald + bootstrap because: (1) the null model
#' is amortised across all peaks per gene, and (2) no bootstrap is needed.
#'
#' @param object ASCENT object.
#' @param celltype Character. Cell type to analyse.
#' @param regr Character. \code{"poisson"} or \code{"negbin"}.
#' @param bin Logical. Binarise ATAC counts?
#' @param method Character. \code{"glm"} or \code{"fastglm"} for null model.
#'
#' @return Data frame with columns: gene, peak, beta, se, z, score_U,
#'   score_V, score_stat, score_p.
#' @noRd
.ASCENT_algorithm_score_r <- function(object, celltype, regr, bin,
                                      method = "glm") {
  if (method == "fastglm" && !requireNamespace("fastglm", quietly = TRUE))
    stop("Install fastglm: install.packages('fastglm')")
  if (regr == "negbin" && !requireNamespace("MASS", quietly = TRUE))
    stop("Install MASS for negbin: install.packages('MASS')")

  # Identify cells of the target celltype
  meta <- object@meta.data
  if ("cell" %in% colnames(meta)) rownames(meta) <- meta$cell
  cell_ids <- colnames(object@rna)
  meta_ordered <- meta[cell_ids, , drop = FALSE]
  ct_col <- meta_ordered[[object@celltypes]]
  is_target <- !is.na(ct_col) & ct_col == celltype
  target_cells <- cell_ids[is_target]
  n_cells <- length(target_cells)
  if (n_cells < 10) return(data.frame())

  # Covariate data for target cells
  cov_target <- meta_ordered[target_cells, object@covariates, drop = FALSE]

  # Group pairs by gene for null model amortisation
  pairs <- object@peak.info
  gene_groups <- split(seq_len(nrow(pairs)), pairs[[1]])

  res_list <- vector("list", nrow(pairs))
  res_idx <- 0
  n_genes_done <- 0
  n_genes_total <- length(gene_groups)

  for (gene in names(gene_groups)) {
    pair_indices <- gene_groups[[gene]]

    # Extract gene expression for target cells
    y <- as.numeric(object@rna[gene, target_cells])

    # Gene sparsity filter
    if (mean(y > 0) <= 0.05) { n_genes_done <- n_genes_done + 1; next }

    # Fit null model (no peak term)
    theta <- 0
    if (method == "fastglm") {
      X_null <- cbind(Intercept = 1, as.matrix(cov_target))
      null_fit <- tryCatch({
        if (regr == "poisson") {
          fastglm::fastglm(x = X_null, y = y, family = poisson())
        } else {
          .negbin_fastglm(X_null, y)
        }
      }, error = function(e) NULL)
      if (is.null(null_fit)) { n_genes_done <- n_genes_done + 1; next }
      mu <- null_fit$fitted.values
      if (regr == "negbin") theta <- null_fit$theta
    } else {
      df_null <- cov_target
      df_null$exprs <- y
      null_formula <- as.formula(
        paste("exprs ~", paste(object@covariates, collapse = " + "))
      )
      null_fit <- tryCatch({
        if (regr == "poisson") {
          glm(null_formula, family = "poisson", data = df_null)
        } else {
          MASS::glm.nb(null_formula, data = df_null)
        }
      }, error = function(e) NULL)
      if (is.null(null_fit)) { n_genes_done <- n_genes_done + 1; next }
      mu <- fitted(null_fit)
      X_null <- model.matrix(null_fit)
      if (regr == "negbin") theta <- null_fit$theta
    }

    # Weights and weighted residuals
    # Poisson: W = mu, e_star = y - mu
    # NegBin:  W = mu*theta/(theta+mu), e_star = (y-mu)*theta/(theta+mu)
    if (regr == "poisson") {
      W      <- mu
      e_star <- y - mu
    } else {
      W      <- mu * theta / (theta + mu)
      e_star <- (y - mu) * theta / (theta + mu)
    }

    # Precompute (X'WX)^{-1} for WLS projection
    XtWX_inv <- tryCatch(
      solve(crossprod(X_null, X_null * W)),
      error = function(e) NULL
    )
    if (is.null(XtWX_inv)) { n_genes_done <- n_genes_done + 1; next }

    # Score test each peak linked to this gene
    for (pi in pair_indices) {
      peak <- pairs[[2]][pi]
      z <- as.numeric(object@atac[peak, target_cells])
      if (bin && any(z > 0)) z[z > 0] <- 1
      if (mean(z > 0) <= 0.05) next

      # b_tilde: residual of WLS regression of atac on X_null
      # Weights: Poisson W=mu, NegBin W=mu*theta/(theta+mu)
      b_tilde <- as.numeric(
        z - X_null %*% (XtWX_inv %*% crossprod(X_null, W * z))
      )

      U <- sum(b_tilde * e_star)
      V <- sum(b_tilde^2 * e_star^2)
      if (V <= 0) next

      # One-step Newton-Raphson beta from null: beta = U / I, se = sqrt(V) / I
      I_info <- sum(b_tilde^2 * W)
      if (I_info <= 0) next
      beta_approx <- U / I_info
      se_approx   <- sqrt(V) / I_info
      z_approx    <- U / sqrt(V)

      score_T <- U^2 / V
      score_p <- pchisq(score_T, df = 1, lower.tail = FALSE)

      res_idx <- res_idx + 1
      res_list[[res_idx]] <- data.frame(
        gene = gene, peak = peak,
        beta = beta_approx, se = se_approx, z = z_approx,
        score_U = U, score_V = V,
        score_stat = score_T, score_p = score_p,
        stringsAsFactors = FALSE
      )
    }

    n_genes_done <- n_genes_done + 1
    if (n_genes_done %% max(1, n_genes_total %/% 20) == 0) {
      message(sprintf("  [%d / %d genes processed (%.0f%%)]",
                      n_genes_done, n_genes_total,
                      100 * n_genes_done / n_genes_total))
    }
  }

  if (res_idx == 0) return(data.frame())
  do.call(rbind, res_list[seq_len(res_idx)])
}


# ================================================================
# S4 CLASS DEFINITION
# ================================================================

#' Validate dimensions and feature names of an ASCENT object
#'
#' Checks that the RNA and ATAC matrices have the same number of cells
#' (columns) and that all genes and peaks referenced in \code{peak.info}
#' exist in the corresponding matrix rownames. Also warns if the ATAC
#' matrix has fewer rows than the RNA matrix (suggesting transposed input).
#'
#' This function is used as the validity method for the ASCENT S4 class
#' and is called automatically on object creation.
#'
#' @param object An ASCENT S4 object.
#'
#' @return \code{TRUE} if valid, otherwise a character vector of error messages.
#' @export
check_dimensions <- function(object) {
  errors <- character()

  num_cells_rna  <- ncol(object@rna)
  num_cells_atac <- ncol(object@atac)
  num_genes <- nrow(object@rna)
  num_peaks <- nrow(object@atac)

  if (num_cells_rna != num_cells_atac) {
    msg <- paste("Number of cells in RNA matrix (", num_cells_rna,
                 ") != ATAC matrix (", num_cells_atac, ").",
                 "Both should be (features x cells).")
    errors <- c(errors, msg)
  }

  if (num_peaks < num_genes) {
    warning("Fewer peaks than genes — verify matrices are (features x cells).")
  }

  if (nrow(object@peak.info) > 0) {
    if (!all(object@peak.info[[1]] %in% rownames(object@rna))) {
      errors <- c(errors, "Some genes in peak.info are missing from RNA matrix rownames.")
    }
    if (!all(object@peak.info[[2]] %in% rownames(object@atac))) {
      errors <- c(errors, "Some peaks in peak.info are missing from ATAC matrix rownames.")
    }
  }

  if (length(errors) == 0) TRUE else errors
}


#' Create an ASCENT object
#'
#' Constructs an S4 object holding all data required by
#' \code{\link{ASCENT_algorithm}}: paired RNA and ATAC sparse count matrices,
#' cell metadata with covariates, and a table of peak-gene pairs to test.
#'
#' @slot rna \code{dgCMatrix}. Gene-by-cell RNA count matrix.
#' @slot atac \code{dgCMatrix}. Peak-by-cell ATAC count matrix. Must have the
#'   same number of columns (cells) as \code{rna}.
#' @slot meta.data \code{data.frame}. Cell metadata. Must contain a
#'   \code{cell} column matching the column names of \code{rna}/\code{atac},
#'   the covariate columns specified in \code{covariates}, and the cell-type
#'   column specified in \code{celltypes}.
#' @slot peak.info \code{data.frame}. Two columns: \code{gene} (matching
#'   \code{rownames(rna)}) and \code{peak} (matching \code{rownames(atac)}).
#' @slot peak.info.list \code{list}. Batched version of \code{peak.info},
#'   populated by \code{\link{CreatePeakToGeneList}} for external
#'   parallelisation.
#' @slot covariates \code{character}. Column names in \code{meta.data} to
#'   include as covariates in the GLM (e.g. \code{c("log_nUMI", "log_nATAC")}).
#' @slot celltypes \code{character}. Column name in \code{meta.data}
#'   containing cell-type labels.
#' @slot ASCENT.result \code{data.frame}. Output from
#'   \code{\link{ASCENT_algorithm}} (empty until the algorithm is run).
#'
#' @return An ASCENT S4 object.
#' @export
CreateASCENTObj <- setClass(
  Class = "ASCENT",
  slots = c(
    rna            = "dgCMatrix",
    atac           = "dgCMatrix",
    meta.data      = "data.frame",
    peak.info      = "data.frame",
    peak.info.list = "list",
    covariates     = "character",
    celltypes      = "character",
    ASCENT.result  = "data.frame"
  ),
  validity = check_dimensions
)


# ================================================================
# MAIN ALGORITHM
# ================================================================

#' Run the ASCENT algorithm
#'
#' Tests for statistical association between chromatin accessibility (ATAC-seq
#' peaks) and gene expression (RNA) across single cells using a GLM framework.
#'
#' Two testing strategies are available via the \code{test} argument:
#' \describe{
#'   \item{\code{"wald"}}{(default) Fits a full Poisson or NegBin GLM per
#'     peak-gene pair including the peak as a predictor, then computes an
#'     empirical p-value via the basic bootstrap method with an adaptive
#'     resampling schedule (100 -> 500 -> 2,500 -> 25,000 -> 50,000
#'     replicates).}
#'   \item{\code{"score"}}{Fits a null model per gene (without the peak term),
#'     then uses the robust score test with HC0 sandwich variance for each
#'     linked peak. Much faster than Wald + bootstrap because the null model
#'     is amortised across all peaks per gene and no bootstrap is needed.
#'     Supports both Poisson and negative binomial.}
#' }
#'
#' Pairs where either the gene or peak has <= 5\% nonzero values in the
#' selected cell type are excluded from the output.
#'
#' @param object An ASCENT S4 object created by \code{\link{CreateASCENTObj}}.
#' @param celltype Character. Cell type to analyse. Must match a value in the
#'   \code{celltypes} column of \code{meta.data}.
#' @param ncores Integer. Number of threads. For \code{method = "rcpp"}, this
#'   controls OpenMP threads parallelising across pairs (Wald) or genes
#'   (score). For \code{"fastglm"} and \code{"glm"} with \code{test = "wald"},
#'   this controls bootstrap cores.
#' @param regr Character. Regression family: \code{"poisson"} (default) or
#'   \code{"negbin"}.
#' @param bin Logical. If \code{TRUE} (default), binarise ATAC counts
#'   (nonzero values set to 1).
#' @param method Character. Computational engine:
#'   \describe{
#'     \item{\code{"rcpp"}}{(default) Full C++/OpenMP engine. Fastest.}
#'     \item{\code{"fastglm"}}{Pure-R fallback using \pkg{fastglm}.}
#'     \item{\code{"glm"}}{Pure-R fallback using base \code{glm()}.}
#'   }
#' @param test Character. Testing strategy: \code{"wald"} (default) for
#'   Wald test + adaptive bootstrap, or \code{"score"} for score test +
#'   HC0 sandwich variance (no bootstrap).
#' @param bootstrap Logical. If \code{TRUE} (default), run adaptive bootstrap
#'   to compute \code{boot_basic_p}. If \code{FALSE}, skip bootstrap entirely
#'   (much faster) and return only the asymptotic Wald p-value. Ignored when
#'   \code{test = "score"}.
#'
#' @return The input ASCENT object with the \code{@@ASCENT.result} slot
#'   populated as a \code{data.frame}. For \code{test = "wald"}: columns
#'   \code{gene}, \code{peak}, \code{beta}, \code{se}, \code{z}, \code{p},
#'   \code{boot_basic_p} (NA when \code{bootstrap = FALSE}).
#'   For \code{test = "score"}: columns \code{gene},
#'   \code{peak}, \code{beta}, \code{se}, \code{z}, \code{score_U},
#'   \code{score_V}, \code{score_stat}, \code{score_p}.
#'
#' @examples
#' \dontrun{
#' # Wald test with bootstrap (default)
#' result <- ASCENT_algorithm(
#'   obj,
#'   celltype = "CD14 Mono",
#'   ncores   = 8L,
#'   regr     = "poisson",
#'   bin      = TRUE,
#'   method   = "rcpp",
#'   test     = "wald"
#' )
#'
#' # Score test (faster, no bootstrap)
#' result <- ASCENT_algorithm(
#'   obj,
#'   celltype = "CD14 Mono",
#'   ncores   = 8L,
#'   regr     = "poisson",
#'   bin      = TRUE,
#'   method   = "rcpp",
#'   test     = "score"
#' )
#' }
#'
#' @export
ASCENT_algorithm <- function(object, celltype, ncores = 1L,
                             regr = "poisson", bin = TRUE,
                             method = "rcpp", test = "wald",
                             bootstrap = TRUE) {

  # ---- Validate inputs ----
  stopifnot(inherits(object, "ASCENT"))
  stopifnot(nrow(object@peak.info) > 0)
  regr   <- match.arg(regr, c("poisson", "negbin"))
  method <- match.arg(method, c("rcpp", "fastglm", "glm"))
  test   <- match.arg(test, c("wald", "score"))

  # ---- Score test dispatch ----
  if (test == "score") {
    if (method %in% c("fastglm", "glm")) {
      message(sprintf(
        "ASCENT [%s/score]: %d pairs | celltype='%s' | %s",
        method, nrow(object@peak.info), celltype, regr
      ))
      res <- .ASCENT_algorithm_score_r(object, celltype, regr, bin, method)
      object@ASCENT.result <- res
      return(object)
    }
    # rcpp score test path — falls through to rcpp section below
  }

  # ---- Wald test: R fallback paths ----
  if (test == "wald" && method %in% c("fastglm", "glm")) {
    message(sprintf(
      "ASCENT [%s]: %d pairs | celltype='%s' | %s | bootstrap=%s | %d cores",
      method, nrow(object@peak.info), celltype, regr, bootstrap, ncores
    ))
    res <- if (method == "fastglm") {
      .ASCENT_algorithm_fastglm(object, celltype, ncores, regr, bin, bootstrap)
    } else {
      .ASCENT_algorithm_glm(object, celltype, ncores, regr, bin, bootstrap)
    }
    object@ASCENT.result <- res
    return(object)
  }

  # ---- rcpp path (Wald or score) ----
  regr_int <- if (regr == "poisson") 0L else 1L

  # ---- Subset to relevant features for speed ----
  genes_needed <- unique(object@peak.info[[1]])
  peaks_needed <- unique(object@peak.info[[2]])
  rna_sub  <- object@rna[genes_needed, , drop = FALSE]
  atac_sub <- object@atac[peaks_needed, , drop = FALSE]

  # Ensure dgCMatrix
  if (!is(rna_sub, "dgCMatrix"))  rna_sub  <- as(rna_sub, "dgCMatrix")
  if (!is(atac_sub, "dgCMatrix")) atac_sub <- as(atac_sub, "dgCMatrix")

  # ---- Prepare cell metadata ----
  cell_ids <- colnames(rna_sub)

  # Match meta.data rows to matrix column order
  if ("cell" %in% colnames(object@meta.data)) {
    rownames(object@meta.data) <- object@meta.data$cell
  }
  meta_ordered <- object@meta.data[cell_ids, , drop = FALSE]

  # Cell-type mask
  ct_col    <- meta_ordered[[object@celltypes]]
  cell_mask <- ct_col == celltype

  if (sum(cell_mask) < 10) {
    warning("Fewer than 10 cells for celltype '", celltype, "'. Returning empty result.")
    object@ASCENT.result <- data.frame()
    return(object)
  }

  # ---- Build numeric covariate matrix ----
  cov_mat <- as.matrix(meta_ordered[, object@covariates, drop = FALSE])
  storage.mode(cov_mat) <- "double"

  # ---- Map gene/peak names to 0-based row indices ----
  gene_name_to_idx <- setNames(seq_len(nrow(rna_sub)) - 1L,  rownames(rna_sub))
  peak_name_to_idx <- setNames(seq_len(nrow(atac_sub)) - 1L, rownames(atac_sub))

  pairs    <- object@peak.info
  gene_idx <- as.integer(gene_name_to_idx[pairs[[1]]])
  peak_idx <- as.integer(peak_name_to_idx[pairs[[2]]])

  if (test == "score") {
    message(sprintf(
      "ASCENT [rcpp/score]: %d pairs | %d cells (celltype='%s') | %s | %d threads",
      nrow(pairs), sum(cell_mask), celltype, regr, ncores
    ))

    res <- ascent_score_pairs(
      rna_sparse   = rna_sub,
      atac_sparse  = atac_sub,
      gene_idx_r   = gene_idx,
      peak_idx_r   = peak_idx,
      gene_names_r = as.character(pairs[[1]]),
      peak_names_r = as.character(pairs[[2]]),
      cov_mat      = cov_mat,
      cell_mask    = cell_mask,
      binarize     = bin,
      regr_type    = regr_int,
      ncores       = as.integer(ncores)
    )
  } else {
    message(sprintf(
      "ASCENT [rcpp]: %d pairs | %d cells (celltype='%s') | %s | %d threads",
      nrow(pairs), sum(cell_mask), celltype, regr, ncores
    ))

    res <- ascent_process_pairs(
      rna_sparse   = rna_sub,
      atac_sparse  = atac_sub,
      gene_idx_r   = gene_idx,
      peak_idx_r   = peak_idx,
      gene_names_r = as.character(pairs[[1]]),
      peak_names_r = as.character(pairs[[2]]),
      cov_mat      = cov_mat,
      cell_mask    = cell_mask,
      binarize     = bin,
      regr_type    = regr_int,
      ncores       = as.integer(ncores),
      skip_bootstrap = !bootstrap
    )
  }

  object@ASCENT.result <- res
  return(object)
}


# ================================================================
# PEAK-TO-GENE LIST GENERATION (kept from original SCENT)
# ================================================================

#' Generate cis peak-gene pair lists via bedtools intersection
#'
#' Identifies all ATAC peaks that fall within a genomic window around each
#' gene body (defined by \code{genebed}), then splits the resulting pairs
#' into \code{nbatch} batches for external parallelisation (e.g. submitting
#' separate SLURM jobs per batch).
#'
#' Requires \code{bedtools} to be available on the system PATH.
#'
#' @param object An ASCENT S4 object.
#' @param genebed Character. Path to a BED file defining gene windows
#'   (e.g. gene body +/- 500 kb). Column 4 must contain gene names matching
#'   \code{rownames(object@@rna)}.
#' @param nbatch Integer. Number of batches to split pairs into for external
#'   parallelisation.
#' @param tmpfile Character. Path for a temporary BED file of ATAC peaks
#'   (default \code{"./temporary_atac_peak.bed"}). Deleted after intersection.
#' @param intersectedfile Character. Path for the gzipped bedtools intersect
#'   output (default \code{"./temporary_atac_peak_intersected.bed.gz"}).
#'
#' @return The input ASCENT object with \code{@@peak.info.list} populated as
#'   a named list of data frames, each with \code{gene} and \code{peak}
#'   columns.
#'
#' @examples
#' \dontrun{
#' obj <- CreatePeakToGeneList(
#'   obj,
#'   genebed = "/path/to/GeneBody_500kb_margin.bed",
#'   nbatch  = 10
#' )
#' # Access batch 1
#' obj@@peak.info.list[["1"]]
#' }
#'
#' @export
CreatePeakToGeneList <- function(object,
                                 genebed = "/path/to/GeneBody_500kb_margin.bed",
                                 nbatch,
                                 tmpfile = "./temporary_atac_peak.bed",
                                 intersectedfile = "./temporary_atac_peak_intersected.bed.gz") {

  peaknames   <- rownames(object@atac)
  peaknames_r <- gsub("[_:]", "-", peaknames)

  peak_bed <- data.frame(
    chr   = str_split_fixed(peaknames_r, "-", 3)[, 1],
    start = str_split_fixed(peaknames_r, "-", 3)[, 2],
    end   = str_split_fixed(peaknames_r, "-", 3)[, 3],
    peak  = peaknames
  )

  write.table(peak_bed, tmpfile, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t")

  system(paste("bedtools intersect -a", genebed, "-b", tmpfile,
               "-wa -wb -loj | gzip -c >", intersectedfile))
  unlink(tmpfile)

  d <- data.frame(fread(intersectedfile, sep = "\t"))
  d <- d[d$V5 != ".", ]

  cis.g2p <- d[, c("V4", "V8")]
  colnames(cis.g2p) <- c("gene", "peak")
  cis.g2p <- cis.g2p[cis.g2p$gene %in% rownames(object@rna), ]

  cis.g2p$index       <- seq_len(nrow(cis.g2p))
  cis.g2p$batch_index <- cut2(cis.g2p$index, g = nbatch, levels.mean = TRUE)
  cis.g2p_list <- split(cis.g2p, f = cis.g2p$batch_index)
  cis.g2p_list <- lapply(cis.g2p_list, function(x) x[, c("gene", "peak")])
  names(cis.g2p_list) <- seq_along(cis.g2p_list)

  object@peak.info.list <- cis.g2p_list
  return(object)
}
