# test-fix2-grammar.R
#
# Round-2 regression tests for the independent review of the grammar/effects
# fixes (FIX_LIST B1, B12 V2 half, C1, C2, A3 V2 side).

data("SNP55K_maize282_maf04")

.fx2_matrix <- function(n = 200, m = 300, p_lo = 0.15, p_hi = 0.5, seed = 99) {
  set.seed(seed)
  p <- stats::runif(m, p_lo, p_hi)
  g <- vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1, numeric(n))
  dimnames(g) <- list(paste0("i", seq_len(n)), paste0("m", seq_len(m)))
  g
}

# ---------------------------------------------------------------------------
# B1: vqtl(same_as_add = TRUE) after PleioArch zero-effect pruning
# ---------------------------------------------------------------------------
test_that("vqtl(same_as_add = TRUE) sizes effects per retained locus after pleiotropy pruning (B1)", {
  pl <- suppressMessages(suppressWarnings(
    simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2,
                       architecture = "pleiotropy", pi = c(1, 0.5), cor = 0.3,
                       seed = 1) |>
      additive(prop = 0.3, n_qtn = 20)))
  add <- pl$layers[[1]]
  expect_false(length(unique(lengths(add$qtn))) == 1L)   # unequal retained loci
  v <- suppressMessages(vqtl(pl, prop = 0.1))            # used to error
  ly <- v$layers[[2]]
  expect_identical(lengths(ly$effect), lengths(ly$qtn))
  expect_identical(ly$qtn, add$qtn)
  expect_true(all(is.finite(v$pheno$value)))
  # the full chain (dominance, vQTL, epistasis, table, long output) is coherent
  full <- suppressMessages(suppressWarnings(
    pl |> dominance(prop = 0.1) |> vqtl(prop = 0.1) |>
      epistasis(prop = 0.1, n_pairs = 4)))
  for (ly in full$layers) {
    if (ly$type %in% c("additive", "dominance", "vqtl")) {
      expect_identical(lengths(ly$effect), lengths(ly$qtn))
    }
  }
  expect_gt(nrow(qtn_table(full)), 0L)
  # print names the retained per-trait counts instead of the requested 20
  out <- paste(utils::capture.output(print(full)), collapse = "\n")
  expect_match(out, "retained per trait: 15/20", fixed = TRUE)
})

test_that("vqtl(same_as_add = TRUE) works with vary_qtn replications after pruning (B1)", {
  pl <- suppressMessages(suppressWarnings(
    simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, n_reps = 3,
                       vary_qtn = TRUE, architecture = "pleiotropy",
                       pi = c(1, 0.5), cor = 0.3, seed = 1) |>
      additive(prop = 0.3, n_qtn = 20) |> vqtl(prop = 0.1)))
  ly <- pl$layers[[2]]
  for (r in seq_len(3)) {
    expect_identical(lengths(ly$effect_reps[[r]]), lengths(ly$qtn_reps[[r]]))
  }
})

test_that("dominance reusing pruned pleiotropy loci keeps the shared set when a trait has pi = 0 (B1 audit)", {
  # trait 1 has no shared variance, so it no longer holds the shared loci and
  # the intersection of the retained sets is empty; the shared loci must come
  # from the additive draw, not from the intersection.
  pl <- suppressMessages(suppressWarnings(
    simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2,
                       architecture = "pleiotropy", pi = c(0, 1), cor = 0,
                       seed = 1) |>
      additive(prop = 0.3, n_qtn = 20)))
  expect_false(is.null(attr(pl$layers[[1]]$qtn, "pleio_shared")))
  x <- suppressMessages(suppressWarnings(dominance(pl, prop = 0.1)))
  dom <- x$layers[[2]]
  expect_identical(lengths(dom$effect), lengths(dom$qtn))
  expect_true(all(vapply(dom$effect, function(e) any(e != 0), logical(1))))
  expect_true(all(is.finite(x$pheno$value)))
})

test_that("a balanced pleiotropy draw is unchanged by the pruning bookkeeping (B1)", {
  a <- suppressMessages(suppressWarnings(
    simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2,
                       architecture = "pleiotropy", pi = 0.5, cor = 0.3,
                       seed = 1) |> additive(prop = 0.3, n_qtn = 20)))
  expect_identical(unique(lengths(a$layers[[1]]$qtn)), 20L)
  out <- paste(utils::capture.output(print(a)), collapse = "\n")
  expect_no_match(out, "retained per trait", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# B12 (V2 half): orthogonal guard message for an all-heterozygous locus
# ---------------------------------------------------------------------------
test_that("the orthogonal d guard names the all-heterozygous cause too (B12)", {
  M <- .fx2_matrix(n = 60, m = 40)
  M[, 1] <- 0                                   # heterozygous in every individual
  expect_error(
    simulate_phenotype(M, seed = 1) |>
      additive(prop = 0.3, orthogonal = TRUE, a = 0.5, d = 0.5, qtn = c(1, 2)),
    "heterozygous in every individual")
  expect_error(
    simulate_phenotype(M, seed = 1) |>
      additive(prop = 0.3, orthogonal = TRUE, a = 0.5, d = 0.5, qtn = c(1, 2)),
    "no heterozygous individuals")
  # d = 0 on that locus is inert and stays valid
  ok <- simulate_phenotype(M, seed = 1) |>
    additive(prop = 0.3, orthogonal = TRUE, a = 0.5, d = c(0, 0.5), qtn = c(1, 2))
  expect_true(all(is.finite(ok$pheno$value)))
})

# ---------------------------------------------------------------------------
# C1: exact A+D identities; multi-layer; sign of the cross term
# ---------------------------------------------------------------------------
test_that("Cov(dosage, het) changes sign across p = 0.5, so the cross term is not one-signed (C1)", {
  hwe <- function(p, n = 4000) {
    q <- 1 - p
    x <- rep(c(-1, 0, 1), times = round(n * c(q^2, 2 * p * q, p^2)))
    h <- as.numeric(x == 0)
    c(cov = stats::cov(x, h), theory = -(2 * p - 1) * 2 * p * q)
  }
  lo <- hwe(0.2); hi <- hwe(0.8)
  expect_gt(lo[["cov"]], 0)
  expect_lt(hi[["cov"]], 0)
  expect_equal(lo[["cov"]], -hi[["cov"]], tolerance = 1e-3)
  expect_equal(lo[["cov"]], lo[["theory"]], tolerance = 2e-3)
})

.fx2_ad <- function(seed = 5, nA = 1L) {
  M <- .fx2_matrix(n = 150, m = 60, seed = seed)
  sim <- simulate_phenotype(M, h2 = 0.5, seed = seed)
  sim <- additive(sim, prop = 0.2, qtn = 1:8)
  if (nA == 2L) sim <- additive(sim, prop = 0.2, qtn = 1:8)
  dominance(sim, prop = 0.1)
}

test_that("ad_report closes on the aggregate component variances, one layer each (C1)", {
  s <- .fx2_ad(nA = 1L)
  ar <- s$ad_report
  expect_true(all(c("var_cA", "var_cD", "cov2_comp") %in% names(ar)))
  # exact identities
  expect_equal(ar$var_A + ar$var_D + ar$cov2_AD, ar$realized, tolerance = 1e-8)
  expect_equal(ar$var_cA + ar$var_cD + ar$cov2_comp, ar$realized,
               tolerance = 1e-8)
  # one layer of each type: Var(c_A) = prop_A, Var(c_D) = prop_D (over V_P), and
  # the gap to the request is 2Cov / V_P, NOT realized - requested
  vp <- stats::var(s$pheno$value[s$pheno$trait == "Trait_1" & s$pheno$rep == 1])
  expect_equal(ar$var_cA, 0.2 / vp, tolerance = 1e-8)
  expect_equal(ar$var_cD, 0.1 / vp, tolerance = 1e-8)
  expect_equal(ar$realized - ar$requested / vp, ar$cov2_comp, tolerance = 1e-8)
})

test_that("with two additive layers the request-based right-hand side misses the within-additive covariance (C1)", {
  s <- .fx2_ad(nA = 2L)
  ar <- s$ad_report
  vp <- stats::var(s$pheno$value[s$pheno$trait == "Trait_1" & s$pheno$rep == 1])
  expect_equal(ar$requested, 0.5)
  # identical additive layers add perfectly: Var(c_A) = (2 * sqrt(.2))^2 = 0.8
  expect_equal(ar$var_cA * vp, 0.8, tolerance = 1e-8)
  expect_equal(ar$var_cD * vp, 0.1, tolerance = 1e-8)
  # the aggregate identity closes ...
  expect_equal(ar$var_cA + ar$var_cD + ar$cov2_comp, ar$realized,
               tolerance = 1e-8)
  expect_equal(ar$var_A + ar$var_D + ar$cov2_AD, ar$realized, tolerance = 1e-8)
  # ... the request-based one does not, by exactly the within-additive term
  within_A <- (ar$var_cA - 0.4 / vp)
  expect_gt(abs(within_A), 0.1)
  expect_equal(ar$realized - ar$requested / vp,
               within_A + (ar$var_cD - 0.1 / vp) + ar$cov2_comp,
               tolerance = 1e-8)
  expect_gt(abs(ar$realized - ar$requested / vp - ar$cov2_comp), 0.1)
})

test_that("print states the requested/V_P identity and no longer claims one-signed effects (C1)", {
  out <- paste(utils::capture.output(print(.fx2_ad())), collapse = "\n")
  expect_match(out, "requested/V_P", fixed = TRUE)
  expect_match(out, "Var(c_A)", fixed = TRUE)
  expect_match(out, "2Cov(c_A,c_D)", fixed = TRUE)
  expect_match(out, "orthogonal = TRUE", fixed = TRUE)
  expect_no_match(out, "one-signed", fixed = TRUE)
})

test_that("the sources no longer assert a one-signed cross term (C1)", {
  root <- testthat::test_path("..", "..")
  f <- file.path(root, "R", c("grammar_realize.R", "grammar_simulate_phenotype.R",
                              "grammar_layers.R"))
  skip_if_not(all(file.exists(f)), "R sources not available (installed package)")
  txt <- unlist(lapply(f, readLines))
  claims <- grep("one-signed across loci|the fixed-sign geometric effect series is one-signed",
                 txt, value = TRUE)
  # the only mention left is the conditional one ("ONLY when ...")
  expect_true(all(grepl("ONLY|only", claims)) || length(claims) == 0L)
})

# ---------------------------------------------------------------------------
# C2: 31-bit sub-seed is collision-resistant, not injective
# ---------------------------------------------------------------------------
test_that(".layer_seed has a documented cross-family collision far outside ordinary ranges (C2)", {
  expect_identical(.layer_seed(123, "transcriptome_rep106", 0L), 2116039371L)
  expect_identical(.layer_seed(123, "residual_t40160", 3L), 2116039371L)
  # ...while ordinary ranges stay collision-free (regression of O1)
  expect_false(.layer_seed(123, "additive_rep12", 0L) ==
                 .layer_seed(123, "additive_rep21", 0L))
  seeds <- vapply(paste0("additive_rep", 1:2000), .layer_seed, 0L, seed = 123)
  expect_identical(anyDuplicated(seeds), 0L)
})

# ---------------------------------------------------------------------------
# A3 (V2 side): duplicate chr_pos is accepted (no error)
# ---------------------------------------------------------------------------
test_that("duplicate chr/pos markers are accepted by the grammar (A3)", {
  M <- .fx2_matrix(n = 40, m = 30)
  df <- cbind(data.frame(snp = colnames(M), allele = "A/G", chr = 1L,
                         pos = rep(1:15, each = 2) * 1e4, cm = NA_real_,
                         stringsAsFactors = FALSE),
              as.data.frame(t(M)))
  expect_false(anyDuplicated(df[, c("chr", "pos")]) == 0L)
  s <- simulate_phenotype(df, h2 = 0.4, n_qtn = 5, seed = 1)
  expect_true(all(is.finite(s$pheno$value)))
  s2 <- simulate_phenotype(df, seed = 1) |> additive(prop = 0.3, n_qtn = 5)
  expect_true(all(is.finite(s2$pheno$value)))
})
