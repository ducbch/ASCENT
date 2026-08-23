# ASCENT: Methods and Implementation

ASCENT (Accelerated Single-Cell ENhancer Target) is a high-performance reimplementation of the [SCENT](https://github.com/immunogenomics/SCENT) algorithm in Rcpp/C++ with OpenMP parallelism. This document describes the key design decisions, the score test methodology, and the performance characteristics.

---

## 1. The Problem: Model-Based P-values Are Unreliable

SCENT tests for statistical association between chromatin accessibility (ATAC-seq peaks) and gene expression (RNA) across single cells, using a generalized linear model (GLM; Poisson or Negative Binomial) for each peak-gene pair.

A standard Wald test computes a p-value from the model-based standard error (SE), which is derived from the Fisher information matrix (a measure of how much information the data carry about the coefficient). However, single-cell count data frequently violates the variance the model assumes:

- **Poisson** assumes the variance equals the mean, but real single-cell RNA data is overdispersed -- its variance is much larger than its mean -- because of biological heterogeneity, dropout, and technical noise.
- **Negative Binomial** allows for overdispersion through a dispersion parameter (theta), but the variance it assumes may still not match the true data-generating process.

When the variance is underestimated, model-based standard errors are too small, z-statistics are inflated, and **p-values are anti-conservative** (they call too many false positives). This is a well-known problem in count regression for genomics.

---

## 2. SCENT's Solution: Bootstrap

To obtain reliable p-values without trusting the model-based SE, SCENT uses an **adaptive bootstrap**:

1. **Fit the GLM** to get the observed coefficient beta.
2. **Resample cells with replacement** and re-fit the GLM on each bootstrap sample to build an empirical distribution of beta.
3. **Compute an empirical p-value** from the bootstrap distribution using the basic bootstrap method.

The adaptive schedule escalates the number of replicates for promising pairs (100 -> 500 -> 2,500 -> 25,000 -> 50,000), concentrating computation on the pairs most likely to be significant.

This approach is statistically robust -- the bootstrap distribution captures the true sampling variability regardless of model misspecification. But it is computationally expensive: a single pair may require up to 50,000 GLM fits. Genome-wide analyses involve hundreds of thousands of pairs.

---

## 3. ASCENT Wald Test: Accelerating the Bootstrap

ASCENT preserves SCENT's bootstrap approach but replaces the computational engine. The key design changes:

### OpenMP parallelism across pairs

SCENT parallelizes *within* each pair -- it uses `boot::boot(parallel = "multicore")` to distribute bootstrap replicates across forked R processes. Each fork duplicates the entire R session (copy-on-write), so memory scales linearly with core count.

ASCENT instead parallelizes *across* pairs using OpenMP threads with dynamic scheduling. All threads share the same memory space -- the sparse matrices, covariate matrix, and metadata are read-only shared data. Only small per-thread working vectors (design matrix, weights, and IRLS temporaries -- IRLS, iteratively reweighted least squares, is the standard algorithm for fitting a GLM) are thread-local. Memory stays constant regardless of core count.

### Fixed-theta Negative Binomial bootstrap

For NegBin regression, SCENT uses `MASS::glm.nb()` which jointly re-estimates both beta and theta (dispersion) on every bootstrap replicate. This is computationally expensive (~50-100 IRLS iterations per replicate due to alternating beta-theta optimization) and statistically unnecessary -- theta is a nuisance parameter that characterizes the gene's overdispersion, not the peak-gene association being tested.

ASCENT estimates theta once from the initial fit, then holds it fixed for all bootstrap replicates. This is the same approach used by DESeq2 and edgeR. Combined with warm-starting IRLS from the initial fit's converged coefficients, bootstrap replicates converge in 2-3 iterations instead of 50-100. This is the largest single source of speedup for NegBin regression (82x vs 22x for Poisson).

### bootstrap = FALSE option

When `bootstrap = FALSE`, ASCENT skips the adaptive bootstrap entirely and returns only the model-based Wald p-value. This provides ~16x speedup over the default bootstrap and is the recommended setting for NegBin regression (see Section 5).

---

## 4. ASCENT Score Test: Analytic Association Testing

### 4.1 Motivation: testing at the null without a per-pair full fit

Recall the core problem (Section 1): model-based SEs are unreliable when the assumed variance function doesn't match real single-cell data. SCENT's solution is bootstrap -- resample the data thousands of times to empirically measure the true variability of beta, bypassing the model-based SE entirely, at a cost of up to 50,000 GLM fits per pair.

The score test uses a different approach: it evaluates the association analytically at the null model, so it never fits a full GLM per peak and gives a p-value in closed form (details below). With a correctly specified negative binomial model this analytic p-value is close to calibrated on its own; the adaptive bootstrap remains available as an option (Section 4.6) for a resampling-based p-value.

### 4.2 Score test output columns

The score test outputs an analytic model-based p-value and a refined coefficient estimate per pair, plus an optional bootstrap p-value:

| Column | Formula | Description |
|--------|---------|-------------|
| `score_beta` | $\hat{\beta}_{\text{MLE}}$ | Coefficient estimate for the peak effect, refined to the full-model maximum likelihood estimate (MLE); see Section 4.3 |
| `score_se` | $1/\sqrt{I}$ | Model-based standard error |
| `score_z` | $U/\sqrt{I}$ | z-statistic |
| `score_p` | $\Pr(\chi^2_1 > U^2/I)$ | Model-based analytic p-value |
| `boot_p` | — | Bootstrap p-value over every pair when `bootstrap = TRUE`; NA otherwise. See Section 4.6 |
| `score_U` | $U = \sum_i \tilde{b}_i\, e_i$ | Score statistic |
| `score_V` | $I = \sum_i \tilde{b}_i^{2}\, W_i$ | Fisher information |
| `score_stat` | $U^{2}/I$ | Chi-squared(1) test statistic |

Here $\tilde{b}$ is the peak accessibility with the covariates projected out (its part not explained by the covariates); $e$ is the gene's working residual under the null model (roughly, observed minus fitted); and $W$ is the GLM weight ($\mu$ for Poisson, $\mu\theta/(\theta+\mu)$ for negative binomial, where $\mu$ is the fitted mean and $\theta$ the dispersion). Intuitively, $U$ measures how strongly the peak co-varies with the part of expression the null model leaves unexplained.

### 4.3 How beta is computed: score vs Wald

The **Wald beta** is the full maximum likelihood estimate. The GLM is fitted with IRLS, which iterates 10-25 times until convergence, fully optimizing the log-likelihood.

The **score beta** is a one-step approximation. Starting from the null model (no peak effect), ASCENT computes a single Newton-Raphson update -- equivalent to one IRLS iteration with the peak variable added to the design matrix. This is much cheaper: one matrix solve instead of 10-25 iterative solves.

**Why this works:** Near the null, the log-likelihood is approximately quadratic, so one Newton step lands close to the true MLE. For the vast majority of peak-gene pairs (small-to-moderate effects), the score beta and Wald beta are nearly identical.

**Where they diverge:** For large effect sizes, the log-likelihood surface is non-quadratic far from the null. A single Newton step from the null uses the curvature *at the null*, which overestimates the curvature at the true optimum. The result is that the score beta **undershoots** the Wald beta -- it is smaller in magnitude.

**How `score_beta` avoids the undershoot:** ASCENT refines the estimate for **every** pair. It runs up to 3 joint IRLS iterations of the full model (covariates plus the peak term), warm-started at the one-step estimate, holding theta fixed for NegBin. Each iteration re-estimates the peak coefficient while re-adjusting the covariate coefficients around it; this is done efficiently by a Schur complement -- a standard block-matrix shortcut that solves the small update over the covariates plus the peak while reusing the already-fitted null design, rather than refitting the whole model from scratch. The refinement early-exits once converged, so a pair whose one-step value is already at the MLE (any small effect) costs ~1 iteration, while a large effect takes the full three. This is the same technique used by the `fasthurdle` score test (`refine_beta_joint_ztnb`), adapted to a single-part log-link GLM.

Empirically, these iterations recover the full-model MLE to within ~1%. On real PBMC CD14 Monocyte data (1,000 peak-gene pairs), the reported `score_beta` matches the full-GLM Wald MLE to a mean absolute error of 0.000 (Poisson) and 0.001 (NegBin, where ASCENT holds theta fixed while `glm.nb` re-estimates it jointly), with an overall correlation of 0.994.

### 4.4 Reusing the null model across peaks

When multiple peaks are tested against the same gene, the null model (gene ~ covariates) is **fitted once and reused** for all peaks linked to that gene. The C++ implementation groups pairs by gene and parallelizes at the gene level with OpenMP. In a typical genome-wide analysis, this reduces the number of null model fits from number of peak-gene pairs to number of genes.

### 4.5 Why the analytic score test is fast

With `bootstrap = FALSE`, the score test avoids the Wald test's dominant cost -- the bootstrap (up to 50,000 IRLS fits x 2-5 iterations per pair). It replaces it with:

- 1 null model fit (~10-25 IRLS iterations), reused across all peaks for the same gene
- 1 score computation per peak (matrix operations, no full-model fit)
- a short beta refinement per peak (up to 3 joint IRLS iterations, warm-started at the one-step estimate and early-exiting once converged, so small effects cost ~1 iteration)

This keeps the total IRLS iterations at tens of thousands (score) versus tens of millions (Wald + bootstrap) -- a reduction of three to four orders of magnitude.

With `bootstrap = TRUE` (the default) the score test additionally bootstraps every pair (Section 4.6), which makes its runtime comparable to the Wald test; the speed advantage above applies to the analytic-only (`bootstrap = FALSE`) mode.

### 4.6 Optional bootstrap for the score test

The `bootstrap` parameter applies to both tests. With `bootstrap = TRUE` (the default), the score test runs the **same adaptive bootstrap as the Wald test** (100 -> 500 -> ... -> 50,000 replicates) over **every** pair, producing a robust `boot_p` alongside the analytic `score_p`. The bootstrap is warm-started at the pair's refined coefficients (Section 4.3), so its per-replicate refits converge in a few IRLS iterations. `score_beta` remains the refined estimate.

With `bootstrap = FALSE`, the score test returns only the fast analytic result (`score_p` plus the refined `score_beta`), with `boot_p = NA`.

An important caveat on cost: enabling the bootstrap makes the score test's runtime **comparable to the Wald+bootstrap test**, not faster. The bootstrap cost is not spread evenly across pairs -- it concentrates in the significant ones (the pairs whose small p-values drive the adaptive schedule up to 25,000-50,000 replicates), while the null pairs terminate at the first stage at negligible cost. Bootstrapping every pair therefore costs about the same as the Wald test, and restricting it to only the significant pairs yields little additional saving. The score test's speed advantage thus comes from being **bootstrap-free** (`bootstrap = FALSE`): a correctly specified negative binomial model yields an analytic `score_p` that is close to calibrated (very slightly inflated relative to NegBin Wald), so the bootstrap is optional rather than required. The bootstrap is implemented in the rcpp engine only.

---

## 5. Speed Reference

The "Score (no boot)" column is the analytic-only score test (`bootstrap = FALSE`). With the default `bootstrap = TRUE`, the score test also bootstraps every pair, so its runtime is comparable to the NegBin Wald+Boot column (Section 4.6).

| Scenario | Score (no boot) | NegBin Wald (no boot) | NegBin Wald+Boot |
|----------|-----------------|----------------------|------------------|
| 2k cells, 1.8k pairs | ~9 sec | ~8 sec | 13 min |
| 35k cells, 200k pairs | ~4 hours | ~20 hours | ~15 days |
