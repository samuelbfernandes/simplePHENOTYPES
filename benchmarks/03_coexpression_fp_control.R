#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark 03 -- co-expression WITHOUT a genetic basis  (HEADLINE)
# ------------------------------------------------------------------------------
# The pitch: simulate_transcriptome(geno = NULL) generates a purely NON-GENETIC
# transcriptome -- real co-expression modules, but every gene h2 = 0. This is a
# ground truth no genotype->trait simulator provides: co-expression structure
# that is provably NOT genetic.
#
# We show (a) strong within-module vs between-module co-expression, and (b) that
# a naive "genetic co-expression" test -- the everyday inference that a tight
# co-expression module reflects shared genetic regulation -- FALSE-POSITIVES on
# every module here, because the truth is h2 = 0. As a foil we generate a
# genotype-driven transcriptome (real trans-hotspots, h2 > 0) whose co-expression
# looks the same, so co-expression alone cannot tell the two apart.
#
# Run:  Rscript benchmarks/03_coexpression_fp_control.R
# ==============================================================================

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- if (length(f)) dirname(normalizePath(f[1])) else file.path(getwd(), "benchmarks")
  source(file.path(d, "_common.R"))
})

N_IND    <- 280L
N_GENES  <- 150L
N_FACTOR <- 8L
KAPPA    <- 0.6      # strong shared module factor -> tight co-expression

# ---- (1) genotype-free transcriptome: real modules, every gene h2 = 0 --------
txf <- simulate_transcriptome(geno = NULL, n_ind = N_IND, n_genes = N_GENES,
                              n_factors = N_FACTOR,
                              residual_module_fraction = KAPPA, seed = 11)
stopifnot(all(txf$genes$h2_realized == 0))     # ground truth: NOT genetic
mod <- txf$genes$module

# pairwise gene-gene correlations, split into within-module and between-module
cormat <- stats::cor(t(txf$expression))
ut <- upper.tri(cormat)
same <- outer(mod, mod, `==`)
within_r  <- cormat[ut & same]
between_r <- cormat[ut & !same]

# ---- naive "genetic co-expression" test --------------------------------------
# A common (wrong) inference: a module whose members co-express far more tightly
# than background is "co-regulated" -> attributed to shared genetics. We formalize
# it as: flag module m as genetically co-regulated if its within-module |r| is
# stochastically greater than background |r| (one-sided Wilcoxon, p < 0.05).
# Because every gene here has h2 = 0, EVERY flag is a false positive.
naive_flags <- do.call(rbind, lapply(sort(unique(mod)), function(m) {
  idx <- which(mod == m)
  if (length(idx) < 3L) return(NULL)
  sub <- cormat[idx, idx]
  w <- abs(sub[upper.tri(sub)])
  b <- abs(between_r)
  p <- tryCatch(stats::wilcox.test(w, b, alternative = "greater")$p.value,
                error = function(e) NA_real_)
  data.frame(module = m, n_genes = length(idx),
             mean_within_absr = mean(w), mean_between_absr = mean(b),
             p_value = p, flagged = isTRUE(p < 0.05),
             stringsAsFactors = FALSE)
}))
rownames(naive_flags) <- NULL

fpr <- mean(naive_flags$flagged)     # every flagged module is a false positive

# ---- (2) foil: genotype-driven transcriptome (h2 > 0), same module strength --
data("SNP55K_maize282_maf04")
g <- SNP55K_maize282_maf04
txg <- simulate_transcriptome(g, n_genes = N_GENES, n_factors = N_FACTOR,
                              residual_module_fraction = KAPPA,
                              h2 = 0.5, cis_fraction = 0.2, seed = 12)
modg <- txg$genes$module
cormatg <- stats::cor(t(txg$expression))
utg <- upper.tri(cormatg); sameg <- outer(modg, modg, `==`)
within_rg  <- cormatg[utg & sameg]
between_rg <- cormatg[utg & !sameg]

summary_tbl <- data.frame(
  dataset            = c("genotype-free (h2=0)", "genotype-driven (h2>0)"),
  mean_h2_realized   = c(mean(txf$genes$h2_realized), mean(txg$genes$h2_realized)),
  mean_within_absr   = c(mean(abs(within_r)),  mean(abs(within_rg))),
  mean_between_absr  = c(mean(abs(between_r)), mean(abs(between_rg))),
  within_over_between = c(mean(abs(within_r)) / mean(abs(between_r)),
                          mean(abs(within_rg)) / mean(abs(between_rg))),
  stringsAsFactors = FALSE
)

cat("\n=== Benchmark 03: co-expression without a genetic basis (HEADLINE) ===\n")
cat(sprintf("\n  genotype-free transcriptome: %d genes x %d individuals, %d modules\n",
            N_GENES, N_IND, length(unique(mod))))
cat(sprintf("  every gene h2_realized = 0  (range %.3f .. %.3f)\n",
            min(txf$genes$h2_realized), max(txf$genes$h2_realized)))
cat(sprintf("  within-module mean |r| = %.3f   between-module mean |r| = %.3f  (%.1fx)\n",
            mean(abs(within_r)), mean(abs(between_r)),
            mean(abs(within_r)) / mean(abs(between_r))))

cat("\n-- naive 'genetic co-expression' test on the h2 = 0 data --\n")
print(format(naive_flags, digits = 3))
cat(sprintf(
  "\n  %d of %d modules flagged as 'genetically co-regulated' -> FALSE-POSITIVE\n",
  sum(naive_flags$flagged), nrow(naive_flags)))
cat(sprintf("  rate = %.0f%%, yet the true genetic co-regulation is ZERO (h2 = 0).\n",
            100 * fpr))

cat("\n-- foil: co-expression looks the same with vs without a genetic basis --\n")
print(format(summary_tbl, digits = 3))

bench_write_csv(naive_flags, "03_naive_test_false_positives.csv")
bench_write_csv(summary_tbl, "03_coexpression_summary.csv")

# ---- plots: block-ordered correlation heatmap + within/between distribution ---
bench_png("03_coexpression_fp_control.png", {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 2))
  on.exit(graphics::par(op), add = TRUE)

  ord <- order(mod)
  M <- cormat[ord, ord]
  pal <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(64)
  graphics::image(seq_len(nrow(M)), seq_len(ncol(M)), M, zlim = c(-1, 1),
                  col = pal, xlab = "gene (ordered by module)",
                  ylab = "gene (ordered by module)",
                  main = "h2 = 0 co-expression\n(module blocks on the diagonal)")

  db <- stats::density(abs(between_r)); dw <- stats::density(abs(within_r))
  graphics::plot(dw, col = "#d95f0e", lwd = 2,
                 xlim = c(0, max(dw$x, db$x)),
                 ylim = c(0, max(dw$y, db$y)),
                 xlab = "|gene-gene correlation|", main = "within vs between module")
  graphics::lines(db, col = "#2c7fb8", lwd = 2)
  graphics::legend("topright", bty = "n", lwd = 2,
                   col = c("#d95f0e", "#2c7fb8"),
                   legend = c("within module", "between module"))
})

cat("\nBenchmark 03 complete.\n")
