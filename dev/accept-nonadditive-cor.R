# dev/accept-nonadditive-cor.R -- DECISION-023 acceptance run (SPEC-nonadditive-
# correlation.md section 6, criteria 1-2). Not a unit test: 6 models x 3 targets x
# 30 seeds on a simulated outbred HWE panel (600 individuals, 1500 markers), so it
# takes several minutes. The committed tests/testthat/test-nonadditive-cor.R is a
# lighter 4-seed regression of the same behaviour.
#
# Usage (from the package root):  Rscript dev/accept-nonadditive-cor.R
# Pass: in every cell, the mean realized TOTAL correlation and every component's
# mean correlation are each within max(0.06, 3 se) of cor (criteria 1 and 2).
#
suppressMessages(devtools::load_all(".", quiet = TRUE))
set.seed(123)
N <- 600; M <- 1500
p <- runif(M, 0.15, 0.5)
geno <- t(vapply(p, function(pp) stats::rbinom(N, 2, pp) - 1L, integer(N)))
df <- data.frame(snp = paste0("m", seq_len(M)), allele = "A/G",
                 chr = rep(1:10, each = M / 10),
                 pos = rep(seq_len(M / 10), 10) * 1e5, cm = NA_real_,
                 stringsAsFactors = FALSE)
colnames(geno) <- paste0("i", seq_len(N))
df <- cbind(df, as.data.frame(geno))

K <- 60
models <- list(
  AD      = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, seed = s) |>
              additive(prop = 0.3, n_qtn = K) |> dominance(prop = 0.2),
  AD_new  = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, seed = s) |>
              additive(prop = 0.3, n_qtn = K) |>
              dominance(prop = 0.2, same_as_add = FALSE, n_qtn = K),
  AE      = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, seed = s) |>
              additive(prop = 0.3, n_qtn = K) |> epistasis(prop = 0.2, n_pairs = K),
  AExd    = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, seed = s) |>
              additive(prop = 0.3, n_qtn = K) |>
              epistasis(prop = 0.2, n_pairs = K, interaction_type = c("a", "d")),
  ADE_pi  = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, pi = 0.7, seed = s) |>
              additive(prop = 0.25, n_qtn = K) |> dominance(prop = 0.15) |>
              epistasis(prop = 0.15, n_pairs = K),
  oneAD   = function(r, s) simulate_phenotype(df, architecture = "pleiotropy",
              n_traits = 2, cor = r, h2 = 0.5, n_qtn = K, model = "AD", seed = s)
)
rs <- c(-0.5, 0, 0.5); seeds <- 1:30
comp_cor <- function(sim) {
  g <- .genetic_value_matrix(sim, 1L)
  c(total = stats::cor(g[, 1], g[, 2]),
    vapply(sim$layers, function(ly) stats::cor(.component_raw(ly, sim, 1L, 1L),
                                               .component_raw(ly, sim, 2L, 1L)),
           numeric(1)))
}
cat(sprintf("%-7s %5s | %-17s | %s\n", "model", "cor", "TOTAL mean(se)", "per-component means"))
fails <- 0L
for (m in names(models)) for (r in rs) {
  res <- t(vapply(seeds, function(s)
    suppressWarnings(suppressMessages(comp_cor(models[[m]](r, s)))),
    numeric(1 + length(suppressWarnings(suppressMessages(models[[m]](r, 1))$layers)))))
  within <- function(x) {
    abs(mean(x) - r) < max(0.06, 3 * stats::sd(x) / sqrt(length(x)))
  }
  tot <- res[, 1]; mu <- mean(tot); se <- stats::sd(tot) / sqrt(length(tot))
  ok <- within(tot) && all(apply(res[, -1, drop = FALSE], 2L, within))
  if (!ok) fails <- fails + 1L
  cat(sprintf("%-7s %5.1f | %6.3f (%.3f) %s | %s\n", m, r, mu, se,
              if (ok) "PASS" else "FAIL",
              paste(sprintf("%.3f", colMeans(res[, -1, drop = FALSE])), collapse = " ")))
}
cat("\nFAILS:", fails, "of", length(models) * length(rs), "\n")
