#!/usr/bin/env Rscript
## Generate figures for README.md from benchmark results
## Usage: Rscript scripts/generate_readme_figures.R

library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) script_path <- "scripts/generate_readme_figures.R"
basedir <- normalizePath(file.path(dirname(script_path), ".."))
figdir  <- file.path(basedir, "fig")
dir.create(figdir, showWarnings = FALSE)

t2dir <- file.path(basedir, "benchmark_t2_output")
t3dir <- file.path(basedir, "benchmark_t3_output")

## ── Load results ──────────────────────────────────────────────────────────────
# T2 Poisson
r2  <- readRDS(file.path(t2dir, "result_rcpp.rds"))
g2  <- readRDS(file.path(t2dir, "result_glm.rds"))
r2s <- readRDS(file.path(t2dir, "result_rcpp_score.rds"))

# T3 NegBin
r3  <- readRDS(file.path(t3dir, "result_rcpp.rds"))
g3  <- readRDS(file.path(t3dir, "result_glm.rds"))
r3s <- readRDS(file.path(t3dir, "result_rcpp_score.rds"))

# Timing — extract elapsed (3rd element) as plain numeric
load_elapsed <- function(path) as.numeric(readRDS(path)[3])

t2_glm_wald   <- load_elapsed(file.path(t2dir, "timing_glm.rds"))
t2_fglm_wald  <- load_elapsed(file.path(t2dir, "timing_fastglm.rds"))
t2_rcpp_wald  <- load_elapsed(file.path(t2dir, "timing_rcpp.rds"))
t2_glm_score  <- load_elapsed(file.path(t2dir, "timing_glm_score.rds"))
t2_fglm_score <- load_elapsed(file.path(t2dir, "timing_fastglm_score.rds"))
t2_rcpp_score <- load_elapsed(file.path(t2dir, "timing_rcpp_score.rds"))

t3_glm_wald   <- load_elapsed(file.path(t3dir, "timing_glm.rds"))
t3_fglm_wald  <- load_elapsed(file.path(t3dir, "timing_fastglm.rds"))
t3_rcpp_wald  <- load_elapsed(file.path(t3dir, "timing_rcpp.rds"))
t3_glm_score  <- load_elapsed(file.path(t3dir, "timing_glm_score.rds"))
t3_fglm_score <- load_elapsed(file.path(t3dir, "timing_fastglm_score.rds"))
t3_rcpp_score <- load_elapsed(file.path(t3dir, "timing_rcpp_score.rds"))

## ── Figure 1: Wald accuracy (rcpp vs glm) — Poisson ──────────────────────────
m2 <- merge(
  r2[, c("gene","peak","beta","se","boot_basic_p")],
  g2[, c("gene","peak","beta","se","boot_basic_p")],
  by = c("gene","peak"), suffixes = c(".rcpp", ".glm")
)

p1a <- ggplot(m2, aes(x = beta.glm, y = beta.rcpp)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Beta (glm)", y = "Beta (rcpp)",
       title = sprintf("Beta  (ρ = %.6f)", cor(m2$beta.glm, m2$beta.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p1b <- ggplot(m2, aes(x = se.glm, y = se.rcpp)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SE (glm)", y = "SE (rcpp)",
       title = sprintf("SE  (ρ = %.6f)", cor(m2$se.glm, m2$se.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p1c <- ggplot(m2, aes(x = -log10(boot_basic_p.glm), y = -log10(boot_basic_p.rcpp))) +
  geom_point(alpha = 0.3, size = 0.8, color = "#2166AC") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "-log10(boot_p) glm", y = "-log10(boot_p) rcpp",
       title = sprintf("Bootstrap p  (ρ = %.4f)", cor(m2$boot_basic_p.glm, m2$boot_basic_p.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

fig1 <- p1a + p1b + p1c +
  plot_annotation(
    title = "rcpp (ASCENT) vs glm (SCENT) — Wald test, Poisson",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(figdir, "ascent_accuracy.png"), fig1, width = 15, height = 5.5, dpi = 150)
cat("Saved ascent_accuracy.png\n")

## ── Figure 1b: Wald accuracy (rcpp vs glm) — NegBin ─────────────────────────
g3_complete <- g3[complete.cases(g3[, c("beta","se","boot_basic_p")]), ]

m3_wald <- merge(
  r3[, c("gene","peak","beta","se","boot_basic_p")],
  g3_complete[, c("gene","peak","beta","se","boot_basic_p")],
  by = c("gene","peak"), suffixes = c(".rcpp", ".glm")
)

p1d <- ggplot(m3_wald, aes(x = beta.glm, y = beta.rcpp)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Beta (glm)", y = "Beta (rcpp)",
       title = sprintf("Beta  (ρ = %.6f)", cor(m3_wald$beta.glm, m3_wald$beta.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p1e <- ggplot(m3_wald, aes(x = se.glm, y = se.rcpp)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SE (glm)", y = "SE (rcpp)",
       title = sprintf("SE  (ρ = %.6f)", cor(m3_wald$se.glm, m3_wald$se.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p1f <- ggplot(m3_wald, aes(x = -log10(boot_basic_p.glm), y = -log10(boot_basic_p.rcpp))) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "-log10(boot_p) glm", y = "-log10(boot_p) rcpp",
       title = sprintf("Bootstrap p  (ρ = %.4f)", cor(m3_wald$boot_basic_p.glm, m3_wald$boot_basic_p.rcpp, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

fig1b <- p1d + p1e + p1f +
  plot_annotation(
    title = "rcpp (ASCENT) vs glm (SCENT) — Wald test, Negative Binomial",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(figdir, "ascent_accuracy_negbin.png"), fig1b, width = 15, height = 5.5, dpi = 150)
cat("Saved ascent_accuracy_negbin.png\n")

## ── Figure 2: Score vs Wald accuracy (rcpp) — Poisson ─────────────────────────
m2s <- merge(
  r2s[, c("gene","peak","beta","se","score_p")],
  r2[, c("gene","peak","beta","se","p","boot_basic_p")],
  by = c("gene","peak"), suffixes = c(".score", ".wald")
)

p2a <- ggplot(m2s, aes(x = beta.wald, y = beta.score)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#4393C3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Beta (wald)", y = "Beta (score)",
       title = sprintf("Beta  (ρ = %.4f)", cor(m2s$beta.wald, m2s$beta.score, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p2b <- ggplot(m2s, aes(x = se.wald, y = se.score)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#4393C3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SE (wald)", y = "SE (score)",
       title = sprintf("SE  (ρ = %.4f)", cor(m2s$se.wald, m2s$se.score, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p2_pmax <- max(c(-log10(m2s$boot_basic_p), -log10(m2s$score_p)), na.rm = TRUE)
p2c <- ggplot(m2s, aes(x = -log10(boot_basic_p), y = -log10(score_p))) +
  geom_point(alpha = 0.3, size = 0.8, color = "#4393C3") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "-log10(boot_p) wald", y = "-log10(score_p) score",
       title = sprintf("score_p vs boot_p  (ρ = %.4f)", cor(m2s$boot_basic_p, m2s$score_p, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed(xlim = c(0, p2_pmax), ylim = c(0, p2_pmax))

fig2 <- p2a + p2b + p2c +
  plot_annotation(
    title = "Score test vs Wald test (rcpp) — Poisson",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(figdir, "ascent_score_accuracy.png"), fig2, width = 15, height = 5.5, dpi = 150)
cat("Saved ascent_score_accuracy.png\n")

## ── Figure 3: Score vs Wald accuracy (rcpp) — NegBin ──────────────────────────
m3s <- merge(
  r3s[, c("gene","peak","beta","se","score_p")],
  r3[, c("gene","peak","beta","se","p","boot_basic_p")],
  by = c("gene","peak"), suffixes = c(".score", ".wald")
)

p3a <- ggplot(m3s, aes(x = beta.wald, y = beta.score)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Beta (wald)", y = "Beta (score)",
       title = sprintf("Beta  (ρ = %.4f)", cor(m3s$beta.wald, m3s$beta.score, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p3b <- ggplot(m3s, aes(x = se.wald, y = se.score)) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "SE (wald)", y = "SE (score)",
       title = sprintf("SE  (ρ = %.4f)", cor(m3s$se.wald, m3s$se.score, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed()

p3_pmax <- max(c(-log10(m3s$boot_basic_p), -log10(m3s$score_p)), na.rm = TRUE)
p3c <- ggplot(m3s, aes(x = -log10(boot_basic_p), y = -log10(score_p))) +
  geom_point(alpha = 0.3, size = 0.8, color = "#7FBC41") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "-log10(boot_p) wald", y = "-log10(score_p) score",
       title = sprintf("score_p vs boot_p  (ρ = %.4f)", cor(m3s$boot_basic_p, m3s$score_p, method = "spearman"))) +
  theme_bw(base_size = 11) +
  coord_fixed(xlim = c(0, p3_pmax), ylim = c(0, p3_pmax))

fig3 <- p3a + p3b + p3c +
  plot_annotation(
    title = "Score test vs Wald test (rcpp) — Negative Binomial",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(file.path(figdir, "ascent_score_accuracy_negbin.png"), fig3, width = 15, height = 5.5, dpi = 150)
cat("Saved ascent_score_accuracy_negbin.png\n")

## ── Figure 4: Runtime barplot (Poisson, all 6 methods) ───────────────────────
method_levels <- c("rcpp (score)", "fastglm (score)", "glm (score)",
                   "rcpp (wald)", "fastglm (wald)", "glm (wald)")

method_colors <- c(
  "rcpp (score)"   = "#4393C3", "fastglm (score)" = "#FDDBC7", "glm (score)"    = "#D6604D",
  "rcpp (wald)"    = "#2166AC", "fastglm (wald)"  = "#F4A582", "glm (wald)"     = "#B2182B"
)

bench_t2 <- data.frame(
  Method   = factor(c("glm (wald)", "fastglm (wald)", "rcpp (wald)",
                       "glm (score)", "fastglm (score)", "rcpp (score)"),
                    levels = method_levels),
  Time_sec = c(t2_glm_wald, t2_fglm_wald, t2_rcpp_wald,
               t2_glm_score, t2_fglm_score, t2_rcpp_score),
  stringsAsFactors = FALSE
)
bench_t2$Time_min <- bench_t2$Time_sec / 60
ref2 <- bench_t2$Time_sec[bench_t2$Method == "glm (wald)"]
bench_t2$Speedup <- ref2 / bench_t2$Time_sec
bench_t2$Time_label <- ifelse(bench_t2$Time_sec < 60,
                              sprintf("%.1f s", bench_t2$Time_sec),
                       ifelse(bench_t2$Time_sec < 3600,
                              sprintf("%.1f min", bench_t2$Time_min),
                              sprintf("%.1f hr", bench_t2$Time_sec / 3600)))

t2_use_hr <- max(bench_t2$Time_sec) >= 3600
bench_t2$Time_y <- if (t2_use_hr) bench_t2$Time_sec / 3600 else bench_t2$Time_min

p4 <- ggplot(bench_t2, aes(x = Method, y = Time_y, fill = Method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.0fx)", Time_label, Speedup)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = method_colors) +
  labs(title = sprintf("Poisson — %d pairs, 8 cores", nrow(r2)),
       y = if (t2_use_hr) "Time (hours)" else "Time (minutes)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  expand_limits(y = max(bench_t2$Time_y) * 1.3)

ggsave(file.path(figdir, "ascent_benchmark.png"), p4, width = 8, height = 5, dpi = 150)
cat("Saved ascent_benchmark.png\n")

## ── Figure 5: Runtime barplot (NegBin, all 6 methods) ────────────────────────
bench_t3 <- data.frame(
  Method   = factor(c("glm (wald)", "fastglm (wald)", "rcpp (wald)",
                       "glm (score)", "fastglm (score)", "rcpp (score)"),
                    levels = method_levels),
  Time_sec = c(t3_glm_wald, t3_fglm_wald, t3_rcpp_wald,
               t3_glm_score, t3_fglm_score, t3_rcpp_score),
  stringsAsFactors = FALSE
)
bench_t3$Time_min <- bench_t3$Time_sec / 60
ref3 <- bench_t3$Time_sec[bench_t3$Method == "glm (wald)"]
bench_t3$Speedup <- ref3 / bench_t3$Time_sec
bench_t3$Time_label <- ifelse(bench_t3$Time_sec < 60,
                              sprintf("%.1f s", bench_t3$Time_sec),
                       ifelse(bench_t3$Time_sec < 3600,
                              sprintf("%.1f min", bench_t3$Time_min),
                              sprintf("%.1f hr", bench_t3$Time_sec / 3600)))

t3_use_hr <- max(bench_t3$Time_sec) >= 3600
bench_t3$Time_y <- if (t3_use_hr) bench_t3$Time_sec / 3600 else bench_t3$Time_min

p5 <- ggplot(bench_t3, aes(x = Method, y = Time_y, fill = Method)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%s\n(%.0fx)", Time_label, Speedup)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = method_colors) +
  labs(title = sprintf("Negative Binomial — %d pairs, 8 cores", nrow(r3)),
       y = if (t3_use_hr) "Time (hours)" else "Time (minutes)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1)) +
  expand_limits(y = max(bench_t3$Time_y) * 1.3)

ggsave(file.path(figdir, "ascent_benchmark_negbin.png"), p5, width = 8, height = 5, dpi = 150)
cat("Saved ascent_benchmark_negbin.png\n")

cat("\nAll figures saved to:", figdir, "\n")
