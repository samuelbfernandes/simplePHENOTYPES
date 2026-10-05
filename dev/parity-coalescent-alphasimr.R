# Cross-engine evidence for the coalescent founders (docs/SPEC-coalescent.md, phase b):
# simplePHENOTYPES' SMC' core vs AlphaSimR 2.1.0 runMacs() (MaCS with its default
# 1-base history window, i.e. effectively SMC') on the GENERIC preset, read from
# the runMacs() source: "1E8 -t 1E-5 -r 4E-6" -> theta = 1000, rho = 400 per
# chromosome, -eN history below, map length 1 Morgan. runMacs() is not seed
# reproducible, so the comparison is between distributions over replicates.
# Needs AlphaSimR (not a package dependency). Run from the repo root:
#   Rscript dev/parity-coalescent-alphasimr.R
suppressMessages({devtools::load_all(quiet = TRUE); library(AlphaSimR)})

n_ind <- 100L          # diploid, non-inbred -> 200 haplotypes
reps <- as.integer(Sys.getenv("REPS", "20"))
set.seed(as.integer(Sys.getenv("SEED", "20261004")))   # one seed for the whole run
hist <- data.frame(time = c(0.25, 2.5, 25, 250, 2500), size = c(5, 15, 60, 120, 1000))
bins <- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)   # distance, fraction of chromosome

summarise <- function(h, pos) {
  # h: sites x haplotypes 0/1; pos in [0, 1)
  p <- rowMeans(h)
  maf <- pmin(p, 1 - p)
  # r^2 for a random sample of site pairs, binned by distance. The pairs are drawn
  # from the ambient stream (no set.seed here: resetting the global RNG would make
  # every later replicate's founder seed depend on the previous one).
  m <- nrow(h)
  i <- sample.int(m, 20000, replace = TRUE)
  j <- sample.int(m, 20000, replace = TRUE)
  keep <- i != j
  i <- i[keep]; j <- j[keep]
  d <- abs(pos[i] - pos[j])
  hi <- h[i, , drop = FALSE]; hj <- h[j, , drop = FALSE]
  pi <- rowMeans(hi); pj <- rowMeans(hj)
  dd <- rowMeans(hi * hj) - pi * pj
  r2 <- dd^2 / (pi * (1 - pi) * pj * (1 - pj))
  b <- cut(d, bins, include.lowest = TRUE)
  c(S = m, mean_maf = mean(maf),
    sfs = tabulate(cut(maf, seq(0, 0.5, 0.1), include.lowest = TRUE, labels = FALSE), 5) / m,
    r2 = tapply(r2, b, mean))
}

sp <- replicate(reps, {
  x <- .coalescent_chromosome(2L * n_ind, 1000, 400, history = hist)
  summarise(x$hap, x$pos)
})
asr <- replicate(reps, {
  fp <- runMacs(nInd = n_ind, nChr = 1, species = "GENERIC")
  h <- pullSegSiteHaplo(fp, chr = 1)          # haplotypes x sites
  pos <- fp@genMap[[1]] / max(fp@genMap[[1]])   # map is linear in position
  summarise(t(h), as.numeric(pos))
})

out <- data.frame(stat = rownames(sp),
                  simplePHENOTYPES = rowMeans(sp), sp_se = apply(sp, 1, sd) / sqrt(reps),
                  AlphaSimR = rowMeans(asr), asr_se = apply(asr, 1, sd) / sqrt(reps))
out$z <- (out$simplePHENOTYPES - out$AlphaSimR) / sqrt(out$sp_se^2 + out$asr_se^2)
base::print.data.frame(format(out, digits = 4), row.names = FALSE)
