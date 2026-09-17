#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark 04 -- mediation split recovery & the derived/real H2 asymmetry
# ------------------------------------------------------------------------------
# A derived (genome-traced) transcriptome() layer splits the expression-mediated
# phenotype component into a genetic-mediated part Tx_g (counts toward H2) and an
# environmental part Tx_e. We check two things:
#
#   (1) mediation_split()'s genetic-mediated share recovers the intended fraction.
#       If every causal gene has expression heritability h2*, then a share
#       ~ prop * h2* of the phenotype should be genetic-mediated. We sweep h2* and
#       compare mediation_split() to prop * h2*.
#
#   (2) The genetic-mediated part ENTERS realized H2, whereas the SAME expression
#       matrix fed as a REAL/observed source (expression = ) does not: its genetic
#       content is not asserted, so genetic_values() excludes it and H2 ~ 0.
#
# derived vs real use identical causal genes and slopes, so the ONLY difference is
# the provenance of the expression (genome-derived vs observed).
#
# Run:  Rscript benchmarks/04_mediation_recovery.R
# ==============================================================================

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- if (length(f)) dirname(normalizePath(f[1])) else file.path(getwd(), "benchmarks")
  source(file.path(d, "_common.R"))
})

data("SNP55K_maize282_maf04")
g <- SNP55K_maize282_maf04

N_GENES   <- 120L
N_CAUSAL  <- 20L
PROP      <- 0.5                 # target share of V_P for the expression component
H2_STAR   <- c(0.2, 0.5, 0.8)    # per-gene expression heritability to sweep

set.seed(99)
causal_idx <- sort(sample.int(N_GENES, N_CAUSAL))
slopes     <- stats::rnorm(N_CAUSAL)     # fixed slopes shared by derived & real

realized_h2 <- function(sim) {
  gv <- genetic_values(sim)[, 1]
  y  <- phenotypes_wide(sim)$Trait_1
  vy <- stats::var(y)
  if (!is.finite(vy) || vy <= 0) return(NA_real_)
  stats::var(gv) / vy
}

rows <- lapply(seq_along(H2_STAR), function(i) {
  h2s <- H2_STAR[i]
  tx <- simulate_transcriptome(g, n_genes = N_GENES, seed = 300 + i,
                               h2 = h2s, cis_fraction = 0.5)
  causal_h2 <- mean(tx$genes$h2_realized[causal_idx])

  # (1) derived: genome-traced expression -> genetic/environmental split
  ph_d <- simulate_phenotype(g, seed = 400 + i, transcriptome = tx) |>
    transcriptome(prop = PROP, genes = causal_idx, slopes = slopes)
  ms <- mediation_split(ph_d)

  # (2) real: the SAME matrix as an observed source -> no asserted genetics
  ph_r <- simulate_phenotype(g, seed = 400 + i, expression = tx$expression) |>
    transcriptome(prop = PROP, genes = causal_idx, slopes = slopes)

  data.frame(
    h2_star            = h2s,
    causal_h2_realized = causal_h2,
    prop               = PROP,
    intended_gm_target = PROP * h2s,          # prop * h2*
    intended_gm_real   = PROP * causal_h2,    # prop * realized causal h2
    med_genetic        = ms$genetic_mediated,
    med_env            = ms$env_mediated,
    med_cov            = ms$covariance,
    med_total          = ms$genetic_mediated + ms$env_mediated + ms$covariance,
    H2_derived         = realized_h2(ph_d),
    H2_real_source     = realized_h2(ph_r),
    real_med_is_null   = is.null(mediation_split(ph_r)),
    stringsAsFactors = FALSE
  )
})
res <- do.call(rbind, rows)
rownames(res) <- NULL

cat("\n=== Benchmark 04: mediation recovery & derived/real H2 asymmetry ===\n")
cat(sprintf("  %d causal genes, prop = %.2f (target expression-mediated share)\n",
            N_CAUSAL, PROP))
cat("\n-- (1) genetic-mediated share vs intended (prop * h2*) --\n")
print(format(res[, c("h2_star", "causal_h2_realized", "intended_gm_target",
                     "intended_gm_real", "med_genetic", "med_env", "med_total")],
             digits = 3))
cat(sprintf("\n  genetic-mediated tracks prop*h2*: corr(intended, realized) = %.3f\n",
            stats::cor(res$intended_gm_target, res$med_genetic)))
cat(sprintf("  mean |realized - prop*realized_causal_h2| = %.3f\n",
            mean(abs(res$med_genetic - res$intended_gm_real))))

cat("\n-- (2) the SAME expression matrix: derived enters H2, real does not --\n")
print(format(res[, c("h2_star", "med_genetic", "H2_derived", "H2_real_source",
                     "real_med_is_null")], digits = 3))
cat("\n  Derived H2 ~ the genetic-mediated share; a real/observed source of the\n")
cat("  identical matrix yields H2 ~ 0 and NULL mediation (genetics not asserted).\n")

bench_write_csv(res, "04_mediation_recovery.csv")

# plot: intended vs realized genetic-mediated; derived vs real H2 bars.
bench_png("04_mediation_recovery.png", {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  lim <- c(0, max(res$intended_gm_target, res$med_genetic) * 1.1)
  graphics::plot(res$intended_gm_target, res$med_genetic, pch = 19, cex = 1.4,
                 col = "#d95f0e", xlim = lim, ylim = lim,
                 xlab = "intended genetic-mediated (prop * h2*)",
                 ylab = "realized (mediation_split)",
                 main = "mediation recovery")
  graphics::abline(0, 1, col = "grey40", lwd = 2, lty = 2)
  graphics::text(res$intended_gm_target, res$med_genetic,
                 labels = sprintf("h2*=%.1f", res$h2_star), pos = 4, cex = 0.8)

  bh <- t(as.matrix(res[, c("H2_derived", "H2_real_source")]))
  colnames(bh) <- sprintf("h2*=%.1f", res$h2_star)
  graphics::barplot(bh, beside = TRUE, ylim = c(0, max(bh, na.rm = TRUE) * 1.25),
                    col = c("#31a354", "#bdbdbd"), ylab = "realized H2",
                    main = "expression enters H2\nonly when genome-derived")
  graphics::legend("topleft", bty = "n", fill = c("#31a354", "#bdbdbd"),
                   legend = c("derived (transcriptome=)", "real (expression=)"),
                   cex = 0.9)
})

cat("\nBenchmark 04 complete.\n")
