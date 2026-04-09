# ASCENT

**Accelerated Single-Cell ENhancer Target gene linking**

ASCENT is a high-performance R package for linking cis-regulatory elements (ATAC-seq peaks) to their target genes using single-cell multimodal (RNA + ATAC) data. It is a reimplementation of the [SCENT](https://github.com/immunogenomics/SCENT) algorithm in Rcpp/C++ with OpenMP parallelism, delivering **~20x faster runtime** while producing numerically identical results. Unlike the R implementations which `fork()` the entire R session for each bootstrap core (memory scales linearly with core count), ASCENT uses shared-memory OpenMP threads — memory stays constant regardless of how many cores are used.

---

## Why ASCENT?

SCENT is a powerful statistical framework for identifying enhancer-gene links, but its native R implementation becomes prohibitively slow at genome-wide scale. A typical analysis involves tens of thousands of peak-gene pairs, each requiring up to 50,000 bootstrap GLM fits.

ASCENT provides three computational backends — the original SCENT algorithm (`glm`), an intermediate R-based optimisation that swaps `glm()` for the faster [`fastglm`](https://cran.r-project.org/package=fastglm) library (`fastglm`), and the full C++/OpenMP engine (`rcpp`):

| | glm (original SCENT) | fastglm | rcpp (ASCENT) |
|---|---|---|---|
| **Poisson (1,800 pairs, 8 cores)** | 12,918 sec | 7,200 sec | **597 sec** |
| **Poisson speedup** | 1.0x (baseline) | 1.8x | **21.6x** |
| **NegBin (32 pairs, 8 cores)** | 6,505 sec | 341 sec | **82 sec** |
| **NegBin speedup** | 1.0x (baseline) | 19.1x | **79.3x** |
| **Numerical agreement** | — | r = 1.000000 | r = 1.000000 |

*Benchmark on 2,286 CD14 Mono cells. The `glm` method uses base R `glm()` / `MASS::glm.nb()` and is identical to the original SCENT implementation. The `fastglm` method is an R-based optimisation developed as part of ASCENT that replaces the GLM solver with `fastglm::fastglm()`. The `rcpp` method is the full C++/OpenMP engine — the primary contribution of ASCENT.*

ASCENT is a **drop-in replacement** for SCENT. The same statistical model (Poisson/Negative Binomial GLM with adaptive bootstrap p-values) is preserved. Only the computational engine changes.

For a detailed explanation of the optimizations, see [DEEP_README.md](DEEP_README.md).

---

## Installation

### Prerequisites

- R >= 3.5.0
- A C++ compiler with C++14 and OpenMP support (GCC >= 5, Clang >= 3.7)
- [RcppArmadillo](https://cran.r-project.org/package=RcppArmadillo) and [Rcpp](https://cran.r-project.org/package=Rcpp)

### Install from source

```r
# Install dependencies
install.packages(c("Rcpp", "RcppArmadillo", "Matrix", "data.table", "stringr", "Hmisc"))

# Install ASCENT from local source
install.packages("/path/to/ASCENT", repos = NULL, type = "source")

# Or using devtools from a git repository
# devtools::install_github("username/ASCENT")
```

### Optional dependencies

These are only needed for the fallback methods (`method = "fastglm"` or `method = "glm"`):

```r
install.packages(c("fastglm", "boot", "MASS"))
```

### Verify installation

```r
library(ASCENT)
?ASCENT_algorithm
```

---

## Quick Start

### 1. Prepare your data

ASCENT requires:
- **RNA count matrix** (`dgCMatrix`): genes x cells
- **ATAC count matrix** (`dgCMatrix`): peaks x cells
- **Cell metadata** (`data.frame`): must include a `cell` column matching the matrix column names, covariates, and a cell-type column
- **Peak-gene pairs** (`data.frame`): two columns — `gene` (matching RNA rownames) and `peak` (matching ATAC rownames)

```r
library(ASCENT)
library(Matrix)

# Example: extract from a Seurat multiome object
rna_mat  <- GetAssayData(sobj, assay = "RNA",  layer = "counts")
atac_mat <- GetAssayData(sobj, assay = "ATAC", layer = "counts")

meta <- data.frame(
  cell      = colnames(rna_mat),
  log_nUMI  = log(sobj$nCount_RNA),
  log_nATAC = log(sobj$nCount_ATAC),
  celltype  = sobj$celltype
)

peak_info <- data.frame(
  gene = c("CD14", "CD14", "LYZ"),
  peak = c("chr5-140588413-140588992", "chr5-140595262-140596150", "chr12-69343211-69344057")
)
```

### 2. Create an ASCENT object

```r
obj <- CreateASCENTObj(
  rna        = rna_mat,
  atac       = atac_mat,
  meta.data  = meta,
  peak.info  = peak_info,
  covariates = c("log_nUMI", "log_nATAC"),
  celltypes  = "celltype"
)
```

### 3. Run the algorithm

```r
result <- ASCENT_algorithm(
  obj,
  celltype = "CD14 Mono",
  ncores   = 8L,          # OpenMP threads
  regr     = "poisson",   # or "negbin"
  bin      = TRUE,         # binarize ATAC counts
  method   = "rcpp"        # "rcpp" (default), "fastglm", or "glm"
)

# View results (sorted by bootstrap p-value)
res <- result@ASCENT.result
res[order(res$boot_basic_p), ]
#   gene                      peak        beta         se         z            p boot_basic_p
#   CD14  chr5-140588413-140588992  0.20862497 0.04127407  5.054626 4.312349e-07      0.00004
#   CD14  chr5-140595262-140596150  0.10122057 0.02974222  3.403262 6.658641e-04      0.00632
#    LYZ  chr12-69343211-69344057  0.08044277 0.01179825  6.818193 9.219273e-12      0.00688
```

### 4. Generate peak-gene pairs from genome annotation (optional)

If you don't have pre-defined pairs, use `CreatePeakToGeneList` with a gene body BED file:

```r
obj <- CreatePeakToGeneList(
  obj,
  genebed = "/path/to/GeneBody_500kb_margin.bed",
  nbatch  = 10
)
# Pairs are stored in obj@peak.info.list (batched for external parallelisation)
```

---

## Method Options

The `method` parameter in `ASCENT_algorithm` controls the computational engine:

| Method | Description | When to use |
|--------|-------------|-------------|
| `"rcpp"` | Full C++/OpenMP engine developed in ASCENT. Parallelizes across pairs. | Default. Production use. |
| `"fastglm"` | R-based optimisation developed in ASCENT. Replaces `glm()` with [`fastglm::fastglm()`](https://cran.r-project.org/package=fastglm) for ~2x speedup. Bootstrap via `boot::boot(parallel = "multicore")`. | Validation or when C++ compilation is unavailable. |
| `"glm"` | Original SCENT implementation using base `glm()` / `MASS::glm.nb()`. | Ground-truth reference. |

All three methods produce numerically equivalent results (beta correlation = 1.000000).

---

## Output

`ASCENT_algorithm` returns an ASCENT S4 object with the `@ASCENT.result` slot populated as a `data.frame`:

| Column | Description |
|--------|-------------|
| `gene` | Gene name |
| `peak` | Peak name |
| `beta` | ATAC coefficient from the GLM |
| `se` | Standard error of the ATAC coefficient |
| `z` | z-statistic (`beta / se`) |
| `p` | Wald p-value (two-sided) |
| `boot_basic_p` | Empirical p-value from adaptive bootstrap (primary result) |

Pairs that fail the quality filter (<=5% nonzero in either modality) are excluded from the output.

---

## Credits

ASCENT is developed by **Duc Nguyen** at the Lee Lab, Boston Children's Hospital.

The statistical framework (Poisson/NB GLM with adaptive bootstrap) is from the original SCENT package:

> **Sakaue, S.**, Weinand, K., Isaac, S., Dey, K.K. et al. (2024). Tissue-specific enhancer-gene maps from multimodal single-cell data identify causal disease alleles. *Nature Genetics*. [doi:10.1038/s41588-024-01682-1](https://doi.org/10.1038/s41588-024-01682-1)

SCENT source code: [https://github.com/immunogenomics/SCENT](https://github.com/immunogenomics/SCENT)

---

## License

MIT
