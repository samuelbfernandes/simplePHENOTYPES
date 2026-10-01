#!/usr/bin/env Rscript
# ==============================================================================
# Benchmark 02 -- cis-eQTL recovery by a marginal eQTL scan
# ------------------------------------------------------------------------------
# Question: are the true cis-eQTL that simulate_transcriptome() plants actually
# recoverable by a naive per-gene association scan?
#
# We simulate a cis-heavy transcriptome (high h2, cis_fraction ~ 1) so genes have
# a strong cis signal, then for a modest subset of genes run a marginal lm of the
# gene's expression on each marker in a small candidate set = {the gene's true
# cis-eQTL} + {random decoy markers}. We rank markers by association strength and
# report where the true cis-eQTL land (rank 1 = top hit) and the detection rate
# (true eQTL in the top-K).
#
# Two gene sets are scanned so the result is not selection-biased (audit TX-F7):
#   * `random`   -- N_TEST genes drawn at random from ALL genes with a cis eQTL
#                   (the headline number);
#   * `strongest`-- the N_TEST genes with the largest |cis effect| (the easiest
#                   subset; shown for comparison only).
# The chance comparator is EXACT per gene: with m true markers among m + D
# candidates, P(at least one true marker in the top K under random ranking) =
# 1 - choose(D, K) / choose(D + m, K); it is averaged over the scanned genes.
#
# Run:  Rscript benchmarks/02_eqtl_recovery.R
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

N_GENES  <- 300L    # genes to simulate
N_TEST   <- 40L     # genes scanned per set (random / strongest)
N_DECOY  <- 200L    # random decoy markers added to each gene's candidate set
TOPK     <- 5L      # "detected" if a true cis-eQTL ranks within the top-K

set.seed(42)

# cis-driven transcriptome: high heritability, (near) all-cis genetic variance.
tx <- simulate_transcriptome(g, n_genes = N_GENES, seed = 7,
                             h2 = 0.8, cis_fraction = 0.95)
E  <- tx$expression                 # genes x individuals
ce <- tx$cis_eqtl                   # true cis-eQTL truth table
stopifnot(!is.null(ce), nrow(ce) > 0)

# candidate genes: every gene with a cis eQTL. `random` = a random N_TEST of them;
# `strongest` = the N_TEST with the largest |effect| (selection-biased).
eff_by_gene <- tapply(abs(ce$effect), ce$gene_id, max)
cand_genes  <- intersect(names(sort(eff_by_gene, decreasing = TRUE)), rownames(E))
strongest_genes <- utils::head(cand_genes, N_TEST)
random_genes    <- sample(cand_genes, N_TEST)

all_snps <- as.character(g$snp)

scan_one <- function(gene) {
  true_snps <- unique(ce$snp[ce$gene_id == gene])
  true_snps <- true_snps[!is.na(true_snps)]
  if (!length(true_snps)) return(NULL)
  decoys <- sample(setdiff(all_snps, true_snps), N_DECOY)
  cand   <- c(true_snps, decoys)
  D      <- bench_dosage(g, cand, snp_index)          # ind x markers
  y      <- as.numeric(E[gene, ])
  # marginal association: |t| statistic of a simple lm per marker (vectorized by
  # correlation, which is monotone in |t| for a single predictor).
  yc <- y - mean(y)
  Dc <- scale(D, center = TRUE, scale = FALSE)
  sx <- sqrt(colSums(Dc^2))
  sx[sx == 0] <- NA_real_
  r  <- as.numeric(crossprod(Dc, yc)) / (sx * sqrt(sum(yc^2)))
  stat <- abs(r)
  ord  <- order(stat, decreasing = TRUE)
  ranks <- match(seq_along(cand), ord)               # rank of each candidate
  names(ranks) <- cand
  best_true_rank <- min(ranks[true_snps], na.rm = TRUE)
  data.frame(
    gene           = gene,
    n_true         = length(true_snps),
    chance_topK    = 1 - choose(N_DECOY, TOPK) / choose(N_DECOY + length(true_snps), TOPK),
    n_candidates   = length(cand),
    best_true_rank = best_true_rank,
    top1           = best_true_rank == 1L,
    topK           = best_true_rank <= TOPK,
    max_true_absr  = max(stat[match(true_snps, cand)], na.rm = TRUE),
    median_decoy_absr = stats::median(stat[match(decoys, cand)], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

res_random    <- do.call(rbind, lapply(random_genes, scan_one))
res_strongest <- do.call(rbind, lapply(strongest_genes, scan_one))
res_random$set <- "random"; res_strongest$set <- "strongest"
res <- rbind(res_random, res_strongest)
rownames(res) <- NULL

summarize <- function(r, label) data.frame(
  gene_set         = label,
  n_genes_tested   = nrow(r),
  n_decoys_each    = N_DECOY,
  topK             = TOPK,
  detection_top1   = mean(r$top1),
  detection_topK   = mean(r$topK),
  chance_topK      = mean(r$chance_topK),   # exact random-ranking comparator
  median_best_rank = stats::median(r$best_true_rank),
  mean_best_rank   = mean(r$best_true_rank),
  mean_n_true      = mean(r$n_true),
  stringsAsFactors = FALSE
)
detect <- rbind(summarize(res_random, "random"),
                summarize(res_strongest, "strongest"))

cat("\n=== Benchmark 02: cis-eQTL recovery (marginal scan) ===\n")
cat(sprintf("  %d random + %d strongest genes scanned; candidate set = true cis-eQTL + %d decoys\n",
            nrow(res_random), nrow(res_strongest), N_DECOY))
cat("\n-- per-gene (first 10 random) --\n")
print(format(utils::head(res_random[, c("gene", "n_true", "best_true_rank",
                                        "top1", "topK", "max_true_absr",
                                        "median_decoy_absr")], 10), digits = 3))
cat("\n-- detection summary --\n")
print(format(detect, digits = 3))
r <- detect[detect$gene_set == "random", ]
cat(sprintf(
  "\n  RANDOM genes: a true cis-eQTL is the single top hit for %.0f%% and within the\n",
  100 * r$detection_top1))
cat(sprintf(
  "  top-%d for %.0f%% (exact random-ranking chance, averaged over genes: %.2f%%).\n",
  TOPK, 100 * r$detection_topK, 100 * r$chance_topK))

bench_write_csv(res, "02_eqtl_recovery_pergene.csv")
bench_write_csv(detect, "02_eqtl_recovery_summary.csv")
res_plot <- res_random

# plot: distribution of the true cis-eQTL rank, and true vs decoy association.
bench_png("02_eqtl_recovery.png", {
  op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  graphics::hist(res_plot$best_true_rank, breaks = seq(0.5, max(res_plot$best_true_rank) + 0.5, 1),
                 col = "#2c7fb8", border = "white",
                 xlab = "rank of best true cis-eQTL", main = "eQTL rank")
  graphics::abline(v = TOPK + 0.5, col = "#d95f0e", lwd = 2, lty = 2)
  graphics::legend("topright", bty = "n", lty = 2, lwd = 2, col = "#d95f0e",
                   legend = sprintf("top-%d cutoff", TOPK), cex = 0.9)

  yl <- range(c(res_plot$max_true_absr, res_plot$median_decoy_absr), na.rm = TRUE)
  graphics::plot(seq_len(nrow(res_plot)), res_plot$max_true_absr, pch = 16, col = "#d95f0e",
                 ylim = yl, xlab = "gene (scanned)", ylab = "|correlation|",
                 main = "true cis-eQTL vs decoys")
  graphics::points(seq_len(nrow(res_plot)), res_plot$median_decoy_absr, pch = 1,
                   col = "grey40")
  graphics::legend("topright", bty = "n", pch = c(16, 1),
                   col = c("#d95f0e", "grey40"), cex = 0.9,
                   legend = c("max true eQTL", "median decoy"))
})

cat("\nBenchmark 02 complete.\n")
