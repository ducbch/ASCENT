# ASCENT: Technical Deep Dive

ASCENT (Accelerated Single-Cell ENhancer Target) is a high-performance reimplementation of the [SCENT](https://github.com/immunogenomics/SCENT) algorithm in Rcpp/C++ with OpenMP parallelism. This document describes the key design decisions, the score test methodology, and the empirical calibration analysis.

---

## 1. The Problem: Model-Based P-values Are Unreliable

SCENT tests for statistical association between chromatin accessibility (ATAC-seq peaks) and gene expression (RNA) across single cells, using a GLM (Poisson or Negative Binomial) for each peak-gene pair.

A standard Wald test computes a p-value from the model-based standard error (derived from the Fisher information matrix). However, single-cell count data frequently violates the assumed variance function:

- **Poisson** assumes Var(y) = mu, but real scRNA-seq data is overdispersed (Var(y) >> mu) due to biological heterogeneity, dropout, and technical noise.
- **Negative Binomial** accounts for overdispersion via a dispersion parameter theta, but the assumed variance Var(y) = mu + mu^2/theta may still not match the true data-generating process.

When the variance is underestimated, model-based SEs are too small, z-statistics are inflated, and **p-values are anti-conservative** (too many false positives). This is a well-known problem in count regression for genomics.

---

## 2. SCENT's Solution: Bootstrap

To obtain reliable p-values without trusting the model-based SE, SCENT uses an **adaptive bootstrap**:

1. **Fit the GLM** to get the observed coefficient beta.
2. **Resample cells with replacement** and re-fit the GLM on each bootstrap sample to build an empirical distribution of beta.
3. **Compute an empirical p-value** from the bootstrap distribution using the basic bootstrap method.

The adaptive schedule escalates the number of replicates for promising pairs (100 -> 500 -> 2,500 -> 25,000 -> 50,000), concentrating computation where it matters.

This approach is statistically robust -- the bootstrap distribution captures the true sampling variability regardless of model misspecification. But it is computationally expensive: a single pair may require up to 50,000 GLM fits. Genome-wide analyses involve hundreds of thousands of pairs.

---

## 3. ASCENT Wald Test: Accelerating the Bootstrap

ASCENT preserves SCENT's bootstrap approach but replaces the computational engine. The key design changes:

### OpenMP parallelism across pairs

SCENT parallelizes *within* each pair -- it uses `boot::boot(parallel = "multicore")` to distribute bootstrap replicates across forked R processes. Each fork duplicates the entire R session (copy-on-write), so memory scales linearly with core count.

ASCENT instead parallelizes *across* pairs using OpenMP threads with dynamic scheduling. All threads share the same memory space -- the sparse matrices, covariate matrix, and metadata are read-only shared data. Only small per-thread working vectors (design matrix, weights, IRLS temporaries) are thread-local. Memory stays constant regardless of core count.

### Fixed-theta Negative Binomial bootstrap

For NegBin regression, SCENT uses `MASS::glm.nb()` which jointly re-estimates both beta and theta (dispersion) on every bootstrap replicate. This is computationally expensive (~50-100 IRLS iterations per replicate due to alternating beta-theta optimization) and statistically unnecessary -- theta is a nuisance parameter that characterizes the gene's overdispersion, not the peak-gene association being tested.

ASCENT estimates theta once from the initial fit, then holds it fixed for all bootstrap replicates. This is the same approach used by DESeq2 and edgeR. Combined with warm-starting IRLS from the initial fit's converged coefficients, bootstrap replicates converge in 2-3 iterations instead of 50-100. This is the largest single source of speedup for NegBin regression (82x vs 22x for Poisson).

### bootstrap = FALSE option

When `bootstrap = FALSE`, ASCENT skips the adaptive bootstrap entirely and returns only the model-based Wald p-value. This provides ~16x speedup over the default bootstrap and is the recommended setting for NegBin regression (see Section 6).

---

## 4. ASCENT Score Test: An Alternative to Bootstrap

### 4.1 Motivation: sandwich SE as an analytical alternative to bootstrap

Recall the core problem (Section 1): model-based SEs are unreliable because the assumed variance function doesn't match real single-cell data. SCENT's solution is bootstrap -- resample the data thousands of times to empirically measure the true variability of beta, bypassing the model-based SE entirely.

The **HC0 sandwich standard error** attempts to solve the same problem analytically. Instead of assuming a variance function, it uses the observed squared residuals to estimate the true variability of each observation. Where the model-based SE asks "how variable would beta be *if the model were correct*?", the sandwich SE asks "how variable is beta *given what the residuals actually look like*?" -- the same question the bootstrap answers, but without resampling.

This is a well-established approach in econometrics (White, 1980) and increasingly used in genomics.

### 4.2 Score test output columns

The score test outputs multiple p-value variants for each pair:

| Column | Formula | Description |
|--------|---------|-------------|
| `score_p` | U^2 / I_info | Model-based (no sandwich correction) |
| `score_p_HC0` | U^2 / V_HC0 | HC0 sandwich (raw) |
| `score_p_HC1` | U^2 / (V_HC0 * n/(n-p)) | HC1 sandwich (degrees-of-freedom corrected) |
| `score_p_HC2` | U^2 / sum(be^2 / (1-h_ii)) | HC2 sandwich (leverage-adjusted) |
| `score_p_HC3` | U^2 / sum(be^2 / (1-h_ii)^2) | HC3 sandwich (aggressive leverage correction) |

Where U is the score statistic, I_info = sum(b_tilde^2 * W) is the Fisher information, V_HC0 = sum(b_tilde^2 * e_star^2) is the HC0 sandwich variance, be = b_tilde * e_star, and h_ii are the hat matrix diagonals from the null model.

### 4.3 How beta is computed: score vs wald

The **wald beta** is the full maximum likelihood estimate (MLE). The GLM is fitted with IRLS (iteratively reweighted least squares), which iterates 10-25 times until convergence, fully optimizing the log-likelihood.

The **score beta** is a one-step approximation. Starting from the null model (beta_peak = 0), ASCENT computes a single Newton-Raphson update -- equivalent to one IRLS iteration with the peak variable added to the design matrix. This is much cheaper: one matrix solve instead of 10-25 iterative solves.

**Why this works:** Near the null, the log-likelihood is approximately quadratic, so one Newton step lands close to the true MLE. For the vast majority of peak-gene pairs (small-to-moderate effects), the score beta and wald beta are nearly identical.

**Where they diverge:** For large effect sizes, the log-likelihood surface is non-quadratic far from the null. A single Newton step from beta = 0 uses the curvature *at the null*, which overestimates the curvature at the true optimum. The result is that score beta **undershoots** the wald beta: |score beta| < |wald beta|.

### 4.4 Null model amortization

When multiple peaks are tested against the same gene, the null model (gene ~ covariates) is **fitted once and reused** for all peaks linked to that gene. The C++ implementation groups pairs by gene and parallelizes at the gene level with OpenMP. In a typical genome-wide analysis, this reduces the number of null model fits from number of peak-gene pairs to number of genes.

### 4.5 Why the score test is fast

The wald test cost per pair is dominated by bootstrap: up to 50,000 IRLS fits x 2-5 iterations each. The score test replaces this with:

- 1 null model fit (~10-25 IRLS iterations), amortized across all peaks for the same gene
- 1 matrix solve per peak (no iteration)
- 1 sandwich variance computation per peak

This reduces the total IRLS iterations from tens of millions (wald) to tens of thousands (score) -- a reduction of three to four orders of magnitude.

---

## 5. Empirical Calibration: Permutation Analysis

### 5.1 Null data generation

Two null generation strategies were used to evaluate Type I error calibration on PBMC CD14 Monocytes (2,286 cells, 1,000 peak-gene pairs, 100 null replicates):

- **rpois**: Normalize -> permute -> Poisson resample (scMultiMap method). Generates Poisson-like counts; only 767/1000 pairs pass the 5% nonzero filter because Poisson resampling adds zeros.
- **shuffle**: Single column permutation of the ATAC matrix (SCENT method). Preserves real count distribution (overdispersion, zero-inflation); all 1000 pairs pass filtering.

Shuffle is the harder, more realistic test because it preserves the exact count distribution that GLMs must handle.

### 5.2 Main finding: NegBin Wald (no bootstrap) is best calibrated

Across both null types, the ranking from best to worst calibrated:

1. **NegBin Wald** -- closest to the diagonal on QQ plots and calibration curves at all alpha levels
2. **NegBin Wald+Boot** -- slightly conservative, limited by bootstrap resolution floor
3. **Poisson Wald+Boot** -- similar to NegBin Boot
4. **NegBin Score (model-based)** -- slightly inflated
5. **NegBin Score (HC0)** -- slightly more inflated than model-based
6. **Poisson Score (HC0)** -- moderately inflated
7. **Poisson Wald** -- inflated
8. **Poisson Score (model-based)** -- most inflated

### 5.3 HC0 helps Poisson but hurts NegBin

A counterintuitive finding: HC0 sandwich correction moves p-values in opposite directions depending on the model family.

**Poisson (misspecified model):** The model-based variance I_info = sum(b_tilde^2 * W) uses W = mu, which assumes Var(Y) = mu. Real data is overdispersed, so residuals e_star^2 are much larger than mu. HC0 uses V_HC0 = sum(b_tilde^2 * e_star^2), which captures the true (larger) residual variance. V_HC0 >> I_info, giving a bigger denominator, smaller test statistic, and less inflation. HC0 partially corrects for the Poisson misspecification.

**NegBin (correctly specified model):** The model-based variance uses W = mu*theta/(theta+mu), which correctly accounts for overdispersion. When the model is correct, E[e_star^2] approximately equals W, so V_HC0 approximately equals I_info. However, HC0 has a known finite-sample downward bias (divides by n instead of n-p). This makes V_HC0 slightly smaller than I_info, yielding a slightly larger test statistic and slightly more inflation. HC0 adds noise when the model is already correct.

### 5.4 HC1/HC2/HC3 make no practical difference

All four HC variants were implemented and benchmarked:

| HC variant | Correction | Max |p_HC0 - p_HCx| |
|------------|-----------|---------------------|
| HC1 | V * n/(n-p) = V * 1.0009 | 0.0002 |
| HC2 | sum(e^2 / (1-h_ii)) | 0.007 |
| HC3 | sum(e^2 / (1-h_ii)^2) | 0.015 |

With ~2,300 cells and 2 parameters, the corrections are negligible. Correlations between all HC variants exceed 0.99999. At most 1-2 pairs flip significance at alpha=0.05 across 1,000 pairs. The calibration curves for HC0/HC1/HC2/HC3 are visually indistinguishable.

**Conclusion:** The Score test inflation relative to NegBin Wald is not caused by finite-sample sandwich bias (which HC1/2/3 fix). It comes from the score test approximation itself -- evaluating the likelihood at null parameters rather than the MLE.

### 5.5 Ranking comparison: Wald vs Score

On 1,000 real peak-gene pairs from PBMC CD14 Mono:

**Within the same model family, Wald and Score (model-based) give nearly identical rankings:**

| Comparison | Spearman rho | Top-10 overlap | Top-100 overlap |
|---|---|---|---|
| Poisson Wald vs Poisson Score | 1.0000 | 10/10 | 100/100 |
| NegBin Wald vs NegBin Score | 0.9996 | 10/10 | 97/100 |

The one-step score approximation is essentially exact for Poisson and nearly exact for NegBin. Wald and Score within the same family will produce identical enrichment curves.

**Across model families, rankings differ substantially at the top:**

| Comparison | Spearman rho | Top-10 overlap | Top-50 overlap |
|---|---|---|---|
| Poisson Wald vs NegBin Wald | 0.972 | 4/10 | 45/50 |

The model choice (Poisson vs NegBin) is the dominant factor in ranking differences. Poisson calls 129 pairs significant at alpha=0.05 vs NegBin's 86. All 86 NegBin hits are a subset of Poisson's 129 -- NegBin never finds something Poisson misses. The 43 Poisson-only pairs are overdispersion-inflated false positives.

**HC0 sandwich substantially disrupts rankings:**

| Comparison | Top-10 overlap | Max rank swap in top 100 |
|---|---|---|
| NegBin Wald vs NegBin Score (HC0) | 5/10 | 222 positions |
| NegBin Score vs NegBin Score (HC0) | -- | 133 positions |

HC0 reshuffles pairs based on per-pair residual patterns rather than true signal strength.

### 5.6 Implications for enrichment analysis

For enrichment-recall curves on orthogonal data (HiC, eQTL, CRISPRi):

- **Poisson vs NegBin will show different enrichment curves**, especially in the top 10-50 pairs. The 43 extra Poisson calls are likely false positives that dilute enrichment.
- **Wald vs Score (model-based) within the same family will show identical curves** -- they produce the same rankings.
- **HC0 vs model-based will show different curves** due to rank disruption, even though the calibration is similar. Model-based rankings are more stable.
- **Fair comparison across methods requires matching effective FDR**, not nominal alpha, since methods have different calibration.

---

## 6. Speed Reference

| Scenario | Score | NegBin Wald (no boot) | NegBin Wald+Boot |
|----------|-------|----------------------|------------------|
| 2k cells, 1.8k pairs | 9 sec | ~8 sec | 13 min |
| 35k cells, 200k pairs | ~4 hours | ~20 hours | ~15 days |
| 500k cells, 200k pairs | ~2.4 days | ~7 months | ~214 days |
