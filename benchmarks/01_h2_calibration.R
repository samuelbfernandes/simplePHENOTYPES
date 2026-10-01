#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark 01 -- per-gene h2 and cis-fraction calibration
# ------------------------------------------------------------------------------
# Question: does simulate_transcriptome() realize the per-gene expression
# heritability (h2) and cis fraction that were requested?
#
# We sweep a range of TARGET h2 values (with cis_fraction held fixed) and a range
# of TARGET cis-fraction values (with h2 held fixed) on the real maize panel, and
# compare the REALIZED per-gene values reported in `$genes` with the targets. We
# report bias (mean realized - target) and RMSE, save a target-vs-realized
# scatter, and FAIL (stop) if a calibration tolerance is exceeded.
#
# Which columns are informative (audit TX-F1/TX-F2):
#   * `h2_allocated` = Var(G)/(Var(G)+Var(R)) is bounded in [0, 1] but is a
#     variance ALLOCATION, not a heritability: it equals the target by
#     construction on the reference population, so it is reported only as a
#     sanity column. The informative one is the realized heritability
#     `h2_realized` = Var(G)/Var(P) = Var(G)/(Var(G)+Var(R)+gr_cov), which is
#     h2/(1 + gr_cov) when Var(G)+Var(R) = 1 (also stored as `h2_var_ratio`), which carries the finite-sample genetic-residual
#     covariance.
#   * `cis_fraction_realized` = v_cis/Var(G) (covariance included). It scatters
#     around the target by the cis-trans covariance, so this panel can fail.
#
# Run:  Rscript benchmarks/01_h2_calibration.R
# ==============================================================================

# --- bootstrap: source shared helpers relative to this script -----------------
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- if (length(f)) dirname(normalizePath(f[1])) else file.path(getwd(), "benchmarks")
  source(file.path(d, "_common.R"))
})

data("SNP55K_maize282_maf04")
g <- SNP55K_maize282_maf04

N_GENES <- 200L         # genes per target (modest; whole script < 1 min)

# ------------------------------------------------------------------------------
# 1. h2 calibration: sweep target h2, fixed cis_fraction
# ------------------------------------------------------------------------------
h2_targets <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
h2_rows <- lapply(seq_along(h2_targets), function(i) {
  h <- h2_targets[i]
  tx <- simulate_transcriptome(g, n_genes = N_GENES, seed = 100 + i,
                               h2 = h, cis_fraction = 0.5)
  real <- tx$genes$h2_realized              # realized Var(G)/Var(P)
  data.frame(target = h, realized = real,
             bounded = tx$genes$h2_allocated,   # allocation, == target by construction
             gene = tx$genes$gene_id, stringsAsFactors = FALSE)
})
h2_df <- do.call(rbind, h2_rows)

h2_summary <- do.call(rbind, lapply(split(h2_df, h2_df$target), function(d) {
  data.frame(
    target      = d$target[1],
    n_genes     = nrow(d),
    mean_real   = mean(d$realized),
    bias        = mean(d$realized - d$target[1]),
    rmse        = sqrt(mean((d$realized - d$target[1])^2)),
    frac_zero   = mean(d$realized == 0),   # degenerate (no realizable genetic var)
    allocated_h2_max_dev = max(abs(d$bounded - d$target[1])),  # ~0 by construction
    stringsAsFactors = FALSE
  )
}))
rownames(h2_summary) <- NULL

# ------------------------------------------------------------------------------
# 2. cis-fraction calibration: sweep target cis_fraction, fixed (high) h2
#    Reported only over genes that actually carry a cis eQTL (n_cis > 0); a gene
#    with no cis marker cannot realize a cis fraction.
# ------------------------------------------------------------------------------
cf_targets <- c(0.1, 0.25, 0.5, 0.75, 0.9)
cf_rows <- lapply(seq_along(cf_targets), function(i) {
  cf <- cf_targets[i]
  tx <- simulate_transcriptome(g, n_genes = N_GENES, seed = 200 + i,
                               h2 = 0.7, cis_fraction = cf)
  keep <- tx$genes$n_cis > 0 & tx$genes$trans_scale > 0 &
    tx$genes$h2_realized > 0            # both a cis and a trans part exist
  data.frame(target = cf,
             realized = tx$genes$cis_fraction_realized[keep],
             stringsAsFactors = FALSE)
})
cf_df <- do.call(rbind, cf_rows)

cf_summary <- do.call(rbind, lapply(split(cf_df, cf_df$target), function(d) {
  data.frame(
    target    = d$target[1],
    n_genes   = nrow(d),
    mean_real = mean(d$realized),
    bias      = mean(d$realized - d$target[1]),
    rmse      = sqrt(mean((d$realized - d$target[1])^2)),
    stringsAsFactors = FALSE
  )
}))
rownames(cf_summary) <- NULL

# ------------------------------------------------------------------------------
# report
# ------------------------------------------------------------------------------
cat("\n=== Benchmark 01: h2 & cis-fraction calibration ===\n")
cat("\n-- per-gene h2 calibration: realized Var(G)/Var(P) (cis_fraction = 0.5,",
    N_GENES, "genes/target) --\n")
print(format(h2_summary, digits = 3))
cat(sprintf("\n  overall h2 bias = %+.4f   overall h2 RMSE = %.4f\n",
            mean(h2_df$realized - h2_df$target),
            sqrt(mean((h2_df$realized - h2_df$target)^2))))

cat("\n-- per-gene cis-fraction calibration: v_cis/Var(G) (h2 = 0.7, genes with a cis and a trans part) --\n")
print(format(cf_summary, digits = 3))
cat(sprintf("\n  overall cis-fraction bias = %+.4f   RMSE = %.4f\n",
            mean(cf_df$realized - cf_df$target),
            sqrt(mean((cf_df$realized - cf_df$target)^2))))

# ---- pass/fail (the panels can fail: neither is a restated target) -----------
h2_bias <- mean(h2_df$realized - h2_df$target)
h2_rmse <- sqrt(mean((h2_df$realized - h2_df$target)^2))
cf_bias <- mean(cf_df$realized - cf_df$target)
cf_rmse <- sqrt(mean((cf_df$realized - cf_df$target)^2))
CHECKS <- c(h2_bias = abs(h2_bias) < 0.01, h2_rmse = h2_rmse < 0.05,
            cis_bias = abs(cf_bias) < 0.02, cis_rmse = cf_rmse < 0.06,
            cis_not_tautological = cf_rmse > 1e-6)
cat("\n-- checks --\n"); print(CHECKS)
if (!all(CHECKS)) stop("Benchmark 01: calibration check(s) failed: ",
                       paste(names(CHECKS)[!CHECKS], collapse = ", "))

bench_write_csv(h2_summary, "01_h2_calibration_summary.csv")
bench_write_csv(cf_summary, "01_cisfraction_calibration_summary.csv")

# ------------------------------------------------------------------------------
# scatter: target vs realized (per gene), with y = x reference
# ------------------------------------------------------------------------------
bench_png("01_h2_calibration_scatter.png", {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  jx <- function(x) x + stats::runif(length(x), -0.012, 0.012)
  graphics::plot(jx(h2_df$target), h2_df$realized,
                 pch = 16, col = grDevices::adjustcolor("#2c7fb8", 0.28),
                 xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "target h2", ylab = "realized h2 (Var(G)/Var(P))",
                 main = "Expression heritability")
  graphics::abline(0, 1, col = "grey40", lwd = 2, lty = 2)
  graphics::points(h2_summary$target, h2_summary$mean_real,
                   pch = 23, bg = "#d95f0e", cex = 1.4)
  graphics::legend("topleft", bty = "n", cex = 0.9,
                   pt.cex = c(1, 1.3), pch = c(16, 23),
                   col = c("#2c7fb8", "black"), pt.bg = c(NA, "#d95f0e"),
                   legend = c("per gene", "target mean"))

  graphics::plot(jx(cf_df$target), cf_df$realized,
                 pch = 16, col = grDevices::adjustcolor("#31a354", 0.28),
                 xlim = c(0, 1), ylim = c(0, 1),
                 xlab = "target cis fraction", ylab = "realized v_cis/Var(G)",
                 main = "cis fraction (genes with cis eQTL)")
  graphics::abline(0, 1, col = "grey40", lwd = 2, lty = 2)
  graphics::points(cf_summary$target, cf_summary$mean_real,
                   pch = 23, bg = "#d95f0e", cex = 1.4)
})

cat("\nBenchmark 01 complete.\n")
