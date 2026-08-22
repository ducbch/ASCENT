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

### 4.1 Motivation: an analytical alternative to bootstrap

Recall the core problem (Section 1): model-based SEs are unreliable when the assumed variance function doesn't match real single-cell data. SCENT's solution is bootstrap -- resample the data thousands of times to empirically measure the true variability of beta, bypassing the model-based SE entirely.

The score test replaces the bootstrap with a cheap, one-step analytical test around the null model (details below). A sandwich (HC) standard error was originally added as an analytical stand-in for the bootstrap's robustness to variance misspecification. However, the permutation calibration study (Section 5) showed the sandwich correction does not improve calibration over the model-based score test on this data, so it was removed and the shipped score test is model-based only. The historical sandwich analysis is retained in Section 5 as the rationale.

### 4.2 Score test output columns

The score test outputs a single, model-based p-value per pair:

| Column | Formula | Description |
|--------|---------|-------------|
| `score_beta` | one-step then joint IRLS | Coefficient estimate: one-step Newton, refined to the MLE for significant pairs; see Section 4.3 |
| `score_se` | 1 / sqrt(I_info) | Model-based standard error |
| `score_z` | U / sqrt(I_info) | z-statistic |
| `score_p` | U^2 / I_info -> chi-sq(1) | Model-based analytic p-value |
| `boot_p` | adaptive bootstrap | Robust p-value for significant pairs (`\|z\| > score_boot_z`); NA otherwise. See Section 4.6 |
| `score_U` | sum(b_tilde * e_star) | Score statistic |
| `score_V` | sum(b_tilde^2 * W) | Fisher information (I_info) |
| `score_stat` | U^2 / I_info | Chi-squared(1) test statistic |

Where U is the score statistic and I_info = sum(b_tilde^2 * W) is the Fisher information.

> **Note on HC sandwich variants.** Earlier versions also emitted HC0/HC1/HC2/HC3 sandwich p-values (`score_p_HC0`, etc.). The permutation calibration study (Section 5) found that HC1/HC2/HC3 are indistinguishable from HC0 at n >> p, and that the sandwich correction does not fix the residual score-test inflation (which comes from the one-step approximation, not finite-sample variance bias). All HC variants have been **removed** from the package; only the model-based score test remains.

### 4.3 How beta is computed: score vs wald

The **wald beta** is the full maximum likelihood estimate (MLE). The GLM is fitted with IRLS (iteratively reweighted least squares), which iterates 10-25 times until convergence, fully optimizing the log-likelihood.

The **score beta** is a one-step approximation. Starting from the null model (beta_peak = 0), ASCENT computes a single Newton-Raphson update -- equivalent to one IRLS iteration with the peak variable added to the design matrix. This is much cheaper: one matrix solve instead of 10-25 iterative solves.

**Why this works:** Near the null, the log-likelihood is approximately quadratic, so one Newton step lands close to the true MLE. For the vast majority of peak-gene pairs (small-to-moderate effects), the score beta and wald beta are nearly identical.

**Where they diverge:** For large effect sizes, the log-likelihood surface is non-quadratic far from the null. A single Newton step from beta = 0 uses the curvature *at the null*, which overestimates the curvature at the true optimum. The result is that score beta **undershoots** the wald beta: |score beta| < |wald beta|.

**How `score_beta` avoids the undershoot:** Since the undershoot only matters for large effects -- which are exactly the pairs that reach significance -- ASCENT refines the estimate *only for significant pairs* (|z| > 2). For those, it runs a few (3) additional joint IRLS iterations of the full model [covariates | peak], warm-started at the one-step estimate, holding theta fixed for NegBin. Each iteration solves the (p_null + 1) system via a Schur complement so only the peak coefficient is profiled out, reusing the already-fitted null design. This is the same technique used by the `fasthurdle` score test (`refine_beta_joint_ztnb`), adapted to a single-part log-link GLM. Non-significant pairs keep the one-step value, which already matches the MLE for small effects.

Empirically, three iterations recover the full-model MLE to within ~1%. On real PBMC CD14 Monocyte data (1,000 pairs, one peak per gene), across the significant pairs the reported `score_beta` matches the full-GLM Wald MLE to a mean absolute error of 0.000 (Poisson) and 0.001 (NegBin, where ASCENT holds theta fixed while `glm.nb` re-estimates it jointly), versus a mean error of 0.05-0.06 and a maximum of ~0.6 for the un-refined one-step value. The overall correlation with the Wald MLE rises from 0.980 (one-step) to 0.994.

### 4.4 Null model amortization

When multiple peaks are tested against the same gene, the null model (gene ~ covariates) is **fitted once and reused** for all peaks linked to that gene. The C++ implementation groups pairs by gene and parallelizes at the gene level with OpenMP. In a typical genome-wide analysis, this reduces the number of null model fits from number of peak-gene pairs to number of genes.

### 4.5 Why the score test is fast

The wald test cost per pair is dominated by bootstrap: up to 50,000 IRLS fits x 2-5 iterations each. The score test replaces this with:

- 1 null model fit (~10-25 IRLS iterations), amortized across all peaks for the same gene
- 1 matrix solve per peak (no iteration)
- 1 model-based variance computation per peak

This reduces the total IRLS iterations from tens of millions (wald) to tens of thousands (score) -- a reduction of three to four orders of magnitude.

### 4.6 Selective bootstrap: screen with the score test, confirm the hits

The score test's one remaining weakness (Section 5) is a mild p-value inflation. That inflation only matters near the significance threshold -- exactly the pairs with a large `|score_z|`. ASCENT offers an **opt-in screen-then-confirm** stage controlled by the `score_boot_z` parameter (**default `Inf` = off**):

- **All pairs** get the cheap analytic `score_p` and a refined coefficient estimate.
- **When `score_boot_z` is finite**, pairs with `|score_z| > score_boot_z` are additionally fitted with the full GLM and run through the **same adaptive bootstrap as the Wald test** (100 -> 500 -> ... -> 50,000 replicates), producing a robust `boot_p`. For these pairs `score_beta` is the exact full-model MLE.

This gates the bootstrap at the *pair* level -- the same idea SCENT uses *within* a pair (escalating replicates for promising pairs), extended across pairs. It is **off by default**, for a measured reason: the bootstrap cost is not spread evenly across pairs but concentrates in the significant ones (the escalating pairs), which the gate keeps. On a benchmark of 1,800 pairs at 8 cores, enabling `score_boot_z = 2` bootstrapped ~15% of pairs but ran *as long as or longer than* the Wald+bootstrap over all 1,800 pairs -- because (a) both methods spend nearly all their time escalating the same significant pairs (the null pairs bail out at stage 0 almost for free), and (b) the score engine parallelizes over *genes*, so the bootstrap-heavy pairs cluster on a few threads and load-balance worse than the Wald engine's per-pair parallelism. The takeaway: the score test's speed advantage comes from being **bootstrap-free**; adding the selective bootstrap trades that away. It remains available for users who want resampled p-values on the hits, and is implemented in the rcpp engine only.

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
