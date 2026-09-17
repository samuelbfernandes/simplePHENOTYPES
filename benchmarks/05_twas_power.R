#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark 05 -- TWAS-style power for a transcriptome-mediated phenotype
# ------------------------------------------------------------------------------
# For a derived transcriptome-mediated phenotype we ask how well a TWAS-style
# gene->phenotype association recovers the causal genes, and how that depends on
# how much of the phenotype the transcriptome explains (prop).
#
#   A) Observed-expression TWAS: correlate each gene's OBSERVED expression with
#      the phenotype. Power (fraction of causal genes passing a Bonferroni
#      threshold) rises with prop. Both cis- and trans-driven causal genes are
#      detectable because expression is observed directly.
#
#   B) cis-predicted TWAS (the twas_sim setting): correlate each gene's
#      CIS-PREDICTED expression (its true cis-eQTL dosages x effects) with the
#      phenotype -- the cis-only imputation a single-gene cis-TWAS builds.
#      cis-DRIVEN causal genes stay associated; TRANS-driven causal genes lose
#      association, because their expression is not a cis function of the genome.
#
# To make the cis vs trans contrast clean we plant a controlled architecture:
# half the genes cis-driven (cis fraction 0.9), half purely trans (cis fraction 0,
# so a cis-only model can reconstruct nothing for them).
#
# Run:  Rscript benchmarks/05_twas_power.R
# ==============================================================================

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", grep("^--file=", a, value = TRUE))
  d <- if (length(f)) dirname(normalizePath(f[1])) else file.path(getwd(), "benchmarks")
  source(file.path(d, "_common.R"))
})

data("SNP55K_maize282_maf04")
g <- SNP55K_maize282_maf04
snp_index <- stats::setNames(seq_len(nrow(g)), as.character(g$snp))

N_GENES     <- 120L
N_CIS_GENES <- 60L                # genes 1..60 cis-driven, 61..120 trans-driven
N_CAUSAL_EA <- 12L                # causal genes drawn per group (cis / trans)
PROPS       <- c(0.05, 0.1, 0.2, 0.4, 0.6)

# controlled cis fraction: strongly cis-driven vs purely trans-driven genes.
cis_frac_vec <- c(rep(0.9, N_CIS_GENES), rep(0.0, N_GENES - N_CIS_GENES))
tx <- simulate_transcriptome(g, n_genes = N_GENES, seed = 21,
                             h2 = 0.8, cis_fraction = cis_frac_vec)
E  <- tx$expression
ce <- tx$cis_eqtl
cf_real <- tx$genes$cis_fraction_realized

set.seed(55)
# cis-driven causal genes need a cis eQTL (so cis-prediction is defined); purely
# trans causal genes need a realized genetic variance (from the trans hub).
cis_pool   <- (seq_len(N_GENES) <= N_CIS_GENES) & tx$genes$n_cis > 0
trans_pool <- (seq_len(N_GENES) >  N_CIS_GENES) & tx$genes$h2_realized > 0
causal_cis   <- sort(sample(which(cis_pool),   N_CAUSAL_EA))
causal_trans <- sort(sample(which(trans_pool), N_CAUSAL_EA))
causal_idx   <- sort(c(causal_cis, causal_trans))
causal_ids   <- rownames(E)[causal_idx]
is_cis       <- causal_idx %in% causal_cis
slopes       <- stats::rnorm(length(causal_idx))

cat(sprintf("Causal genes: %d cis-driven (realized cis frac %.2f), %d trans-driven (%.2f)\n",
            sum(is_cis),  mean(cf_real[causal_cis]),
            sum(!is_cis), mean(cf_real[causal_trans])))

# cis-predicted expression for EVERY gene (true cis-eQTL dosages x effects).
cis_pred <- matrix(0, N_GENES, ncol(E), dimnames = dimnames(E))
for (gid in rownames(E)) {
  cr <- ce[ce$gene_id == gid, , drop = FALSE]
  if (!nrow(cr)) next
  D  <- bench_dosage(g, cr$snp, snp_index)
  Zc <- scale(D, center = TRUE, scale = FALSE)
  cis_pred[gid, ] <- as.numeric(Zc %*% cr$effect)
}

# per-gene correlation + p-value of a predictor matrix (genes x ind) vs y.
cor_p <- function(mat, y) {
  n <- length(y); yc <- y - mean(y); sy <- sqrt(sum(yc^2))
  r <- apply(mat, 1L, function(row) {
    xc <- row - mean(row); sx <- sqrt(sum(xc^2))
    if (sx == 0 || sy == 0) return(0)
    sum(xc * yc) / (sx * sy)
  })
  tt <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  list(absr = abs(r), p = 2 * stats::pt(-abs(tt), df = n - 2))
}

bonf <- 0.05 / N_GENES     # genome-wide threshold for the observed-expression TWAS

rows <- lapply(seq_along(PROPS), function(i) {
  prop <- PROPS[i]
  ph <- simulate_phenotype(g, seed = 500 + i, transcriptome = tx) |>
    transcriptome(prop = prop, genes = causal_idx, slopes = slopes)
  y <- phenotypes_wide(ph)$Trait_1

  obs  <- cor_p(E, y)
  cisp <- cor_p(cis_pred, y)
  ci <- causal_idx[is_cis]; tr <- causal_idx[!is_cis]
  noncausal <- setdiff(seq_len(N_GENES), causal_idx)

  data.frame(
    prop = prop,
    # A) observed-expression TWAS: power over ALL causal genes (Bonferroni).
    # Observed expression detects a causal gene regardless of cis/trans origin.
    power_obs_all = mean(obs$p[causal_idx] < bonf),
    fpr_obs       = mean(obs$p[noncausal] < bonf),
    # B) cis-predicted TWAS: mean association + nominal-alpha power, cis vs trans.
    # Purely trans genes have no cis prediction, so their association is ~0.
    cispred_meanr_cis   = mean(cisp$absr[ci]),
    cispred_meanr_trans = mean(cisp$absr[tr]),
    cispred_pow_cis     = mean(cisp$p[ci] < 0.05),
    cispred_pow_trans   = mean(cisp$p[tr] < 0.05),
    stringsAsFactors = FALSE
  )
})
res <- do.call(rbind, rows)
rownames(res) <- NULL

cat("\n=== Benchmark 05: TWAS-style power vs transcriptome prop ===\n")
cat(sprintf("  %d causal genes (%d cis-driven, %d purely trans); %d genes tested\n",
            length(causal_idx), sum(is_cis), sum(!is_cis), N_GENES))
cat("\n-- (A) observed-expression TWAS: power over causal genes rises with prop --\n")
cat(sprintf("       (Bonferroni p < %.2e; fpr_obs = non-causal false positives)\n", bonf))
print(format(res[, c("prop", "power_obs_all", "fpr_obs")], digits = 3))
cat("\n-- (B) cis-predicted TWAS: cis-driven detectable, purely trans invisible --\n")
cat("       (mean |cor| of cis-predicted expression vs phenotype; power at p<0.05)\n")
print(format(res[, c("prop", "cispred_meanr_cis", "cispred_meanr_trans",
                     "cispred_pow_cis", "cispred_pow_trans")], digits = 3))
cat(sprintf(
  "\n  At prop = %.2f a cis-only TWAS recovers cis-driven causal genes (power %.2f,\n",
  max(PROPS), res$cispred_pow_cis[nrow(res)]))
cat(sprintf(
  "  mean |r| %.3f) but is blind to purely trans-driven ones (power %.2f, |r| %.3f)\n",
  res$cispred_meanr_cis[nrow(res)], res$cispred_pow_trans[nrow(res)],
  res$cispred_meanr_trans[nrow(res)]))
cat("  -- exactly the cis-only limit a single-gene cis-TWAS (e.g. twas_sim) sees.\n")

bench_write_csv(res, "05_twas_power.csv")

# plot: observed-expression power vs prop; cis-predicted association cis vs trans.
bench_png("05_twas_power.png", {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.4, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(res$prop, res$power_obs_all, type = "b", pch = 19, lwd = 2,
                 col = "#1c9099", ylim = c(0, 1),
                 xlab = "transcriptome prop", ylab = "detection rate",
                 main = "observed-expression TWAS")
  graphics::lines(res$prop, res$fpr_obs, type = "b", pch = 1, lwd = 2, lty = 2,
                  col = "grey45")
  graphics::legend("topleft", bty = "n", lwd = 2, lty = c(1, 2), pch = c(19, 1),
                   col = c("#1c9099", "grey45"),
                   legend = c("power (causal)", "FPR (non-causal)"), cex = 0.9)

  graphics::plot(res$prop, res$cispred_meanr_cis, type = "b", pch = 17, lwd = 2,
                 col = "#d95f0e",
                 ylim = c(0, max(res$cispred_meanr_cis) * 1.15),
                 xlab = "transcriptome prop", ylab = "mean |cor| (cis-predicted)",
                 main = "cis-predicted TWAS\n(twas_sim setting)")
  graphics::lines(res$prop, res$cispred_meanr_trans, type = "b", pch = 15, lwd = 2,
                  col = "#756bb1")
  graphics::legend("topleft", bty = "n", lwd = 2, pch = c(17, 15),
                   col = c("#d95f0e", "#756bb1"),
                   legend = c("cis-driven", "purely trans"), cex = 0.9)
})

cat("\nBenchmark 05 complete.\n")
