# ASCENT: Technical Deep Dive

ASCENT (Accelerated Single-Cell ENhancer Target) is a high-performance reimplementation of the [SCENT](https://github.com/immunogenomics/SCENT) algorithm in Rcpp/C++ with OpenMP parallelism. This document describes the key design decisions and the score test methodology introduced in ASCENT.

---

## 1. The Problem: Model-Based P-values Are Unreliable

SCENT tests for statistical association between chromatin accessibility (ATAC-seq peaks) and gene expression (RNA) across single cells, using a GLM (Poisson or Negative Binomial) for each peak-gene pair.

A standard Wald test computes a p-value from the model-based standard error (derived from the Fisher information matrix). However, single-cell count data frequently violates the assumed variance function:

- **Poisson** assumes Var(y) = μ, but real scRNA-seq data is overdispersed (Var(y) >> μ) due to biological heterogeneity, dropout, and technical noise.
- **Negative Binomial** accounts for overdispersion via a dispersion parameter θ, but the assumed variance Var(y) = μ + μ²/θ may still not match the true data-generating process.

When the variance is underestimated, model-based SEs are too small, z-statistics are inflated, and **p-values are anti-conservative** (too many false positives). This is a well-known problem in count regression for genomics.

---

## 2. SCENT's Solution: Bootstrap

To obtain reliable p-values without trusting the model-based SE, SCENT uses an **adaptive bootstrap**:

1. **Fit the GLM** to get the observed coefficient β̂.
2. **Resample cells with replacement** and re-fit the GLM on each bootstrap sample to build an empirical distribution of β̂.
3. **Compute an empirical p-value** from the bootstrap distribution using the basic bootstrap method.

The adaptive schedule escalates the number of replicates for promising pairs (100 → 500 → 2,500 → 25,000 → 50,000), concentrating computation where it matters.

This approach is statistically robust — the bootstrap distribution captures the true sampling variability regardless of model misspecification. But it is computationally expensive: a single pair may require up to 50,000 GLM fits. Genome-wide analyses involve hundreds of thousands of pairs.

---

## 3. ASCENT Wald Test: Accelerating the Bootstrap

ASCENT preserves SCENT's bootstrap approach but replaces the computational engine. The key design changes:

### OpenMP parallelism across pairs

SCENT parallelizes *within* each pair — it uses `boot::boot(parallel = "multicore")` to distribute bootstrap replicates across forked R processes. Each fork duplicates the entire R session (copy-on-write), so memory scales linearly with core count.

ASCENT instead parallelizes *across* pairs using OpenMP threads with dynamic scheduling. All threads share the same memory space — the sparse matrices, covariate matrix, and metadata are read-only shared data. Only small per-thread working vectors (design matrix, weights, IRLS temporaries) are thread-local. Memory stays constant regardless of core count.

### Fixed-theta Negative Binomial bootstrap

For NegBin regression, SCENT uses `MASS::glm.nb()` which jointly re-estimates both beta and theta (dispersion) on every bootstrap replicate. This is computationally expensive (~50–100 IRLS iterations per replicate due to alternating beta-theta optimization) and statistically unnecessary — theta is a nuisance parameter that characterizes the gene's overdispersion, not the peak-gene association being tested.

ASCENT estimates theta once from the initial fit, then holds it fixed for all bootstrap replicates. This is the same approach used by DESeq2 and edgeR. Combined with warm-starting IRLS from the initial fit's converged coefficients, bootstrap replicates converge in 2–3 iterations instead of 50–100. This is the largest single source of speedup for NegBin regression (82x vs 22x for Poisson).

---

## 4. ASCENT Score Test: Eliminating the Bootstrap

### 4.1 Why sandwich SE replaces bootstrap

Recall the core problem (Section 1): model-based SEs are unreliable because the assumed variance function (Var(y) = μ for Poisson, Var(y) = μ + μ²/θ for NegBin) doesn't match real single-cell data. SCENT's solution is bootstrap — resample the data thousands of times to empirically measure the true variability of β̂, bypassing the model-based SE entirely.

The **HC0 sandwich standard error** solves the same problem analytically. Instead of assuming a variance function, it uses the observed squared residuals to estimate the true variability of each observation. Where the model-based SE asks "how variable would β̂ be *if the model were correct*?", the sandwich SE asks "how variable is β̂ *given what the residuals actually look like*?" — the same question the bootstrap answers, but without resampling.

This is a well-established approach in econometrics (White, 1980) and increasingly used in genomics. The sandwich SE is robust to the same misspecification that makes model-based SEs anti-conservative, providing reliable p-values in a single pass through the data.

### 4.2 Why score p-values are comparable to bootstrap p-values

Both the bootstrap and the sandwich SE are correcting for the same underlying problem — variance misspecification — so they arrive at similar answers. The bootstrap measures the true sampling distribution empirically; the sandwich SE estimates it analytically from the residuals. In benchmarks, the score test achieves 97–99% concordance with the wald+bootstrap test on significance calls at α = 0.05.

The score test is slightly **conservative** relative to bootstrap: sandwich SEs tend to be a bit larger than the empirical bootstrap SD, producing slightly higher p-values. This means the score test is less likely to produce false positives — a safe direction of disagreement.

One advantage of the analytic p-value: it has no resolution floor. Bootstrap p-values are bounded by 1/(B+1) — at 50,000 replicates, the minimum is 2×10⁻⁵. The score test can resolve arbitrarily small p-values.

### 4.3 How beta is computed: score vs wald

The **wald beta** is the full maximum likelihood estimate (MLE). The GLM is fitted with IRLS (iteratively reweighted least squares), which iterates 10–25 times until convergence, fully optimizing the log-likelihood.

The **score beta** is a one-step approximation. Starting from the null model (β_peak = 0), ASCENT computes a single Newton-Raphson update — equivalent to one IRLS iteration with the peak variable added to the design matrix. This is much cheaper: one matrix solve instead of 10–25 iterative solves.

**Why this works:** Near the null, the log-likelihood is approximately quadratic, so one Newton step lands close to the true MLE. For the vast majority of peak-gene pairs (small-to-moderate effects), the score beta and wald beta are nearly identical (Spearman ρ > 0.999).

**Where they diverge:** For large effect sizes, the log-likelihood surface is non-quadratic far from the null — the exponential link function means the curvature changes as β moves away from zero. A single Newton step from β = 0 uses the curvature *at the null*, which overestimates the curvature at the true optimum. The result is that score beta **undershoots** the wald beta: |score beta| < |wald beta|.

This undershoot matters little in practice — pairs with large effects are already highly significant under both tests. The undershoot affects the point estimate but not the significance call.

### 4.4 Combined effect on p-values

The two differences compound in the same direction:
- Slightly smaller |beta| (one-step undershoot) in the numerator of z = beta / SE
- Slightly larger SE (sandwich vs model-based) in the denominator

The result is attenuated |z| and more conservative p-values. This is why the score test is a safe drop-in replacement for wald+bootstrap — when they disagree, the score test errs on the side of caution.

The SE gap is **larger for Poisson** (where Var(y) = μ substantially underestimates true variance) and **smaller for Negative Binomial** (where the model already captures overdispersion via θ, leaving less for the sandwich correction to add).

### 4.5 Null model amortization

When multiple peaks are tested against the same gene, the null model (gene ~ covariates) is **fitted once and reused** for all peaks linked to that gene. The C++ implementation groups pairs by gene and parallelizes at the gene level with OpenMP. In a typical genome-wide analysis, this reduces the number of null model fits from number of peak-gene pairs to number of genes.

### 4.6 Why the score test is fast

The wald test cost per pair is dominated by bootstrap: up to 50,000 IRLS fits × 2–5 iterations each. The score test replaces this with:

- 1 null model fit (~10–25 IRLS iterations), amortized across all peaks for the same gene
- 1 matrix solve per peak (no iteration)
- 1 sandwich variance computation per peak

This reduces the total IRLS iterations from tens of millions (wald) to tens of thousands (score) — a reduction of three to four orders of magnitude.

---

## 5. Summary

| | SCENT (R) | ASCENT Wald (Rcpp) | ASCENT Score (Rcpp) |
|---|-----------|---------------|---------------------|
| Parallelism | Fork-based within bootstrap | OpenMP across pairs | OpenMP across genes |
| Memory | N forked processes (copy-on-write) | Shared address space | Shared address space |
| Bootstrap | Up to 50,000 replicates per pair | Up to 50,000 replicates per pair | **Eliminated** |
| NB theta | Re-estimated every replicate | Fixed from initial fit | Estimated in null model |
| Null model | Fitted per pair | Fitted per pair | **Fitted per gene** |
