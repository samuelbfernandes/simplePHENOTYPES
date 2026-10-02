# test-audit-grammar.R
#
# Regression tests for the 2026-09-29 audit of the v2 grammar and the
# effects/architecture engines (GRAM-*, EFF-*). One block per finding; the
# finding id is in each test name.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

# individuals x markers -1/0/1 matrix with HWE genotypes at frequencies `p`
.ag_matrix <- function(n = 200, m = 300, p_lo = 0.15, p_hi = 0.5, seed = 99) {
  set.seed(seed)
  p <- stats::runif(m, p_lo, p_hi)
  g <- vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1, numeric(n))
  dimnames(g) <- list(paste0("i", seq_len(n)), paste0("m", seq_len(m)))
  g
}

# numeric-format data frame (chr/pos map) from a matrix
.ag_df <- function(M, n_chr = 10L) {
  m <- ncol(M)
  cbind(data.frame(snp = colnames(M), allele = "A/G",
                   chr = rep(seq_len(n_chr), length.out = m),
                   pos = seq_len(m) * 1e4, cm = NA_real_,
                   stringsAsFactors = FALSE),
        as.data.frame(t(M)))
}

# ---------------------------------------------------------------------------
# O1 / EFF-F1 (D2): position-sensitive layer sub-seed
# ---------------------------------------------------------------------------
test_that(".layer_seed is position sensitive: permuted digits no longer collide (O1)", {
  expect_false(.layer_seed(123, "additive_rep12", 0L) ==
                 .layer_seed(123, "additive_rep21", 0L))
  expect_false(.layer_seed(123, "residual_t12", 0L) ==
                 .layer_seed(123, "residual_t21", 0L))
  # anagram layer types no longer alias either
  expect_false(.layer_seed(1, "abc", 0L) == .layer_seed(1, "cba", 0L))
})

test_that(".layer_seed is injective on the production labels and layer x occurrence grid (O1)", {
  types <- c("additive", "dominance", "epistasis", "vqtl", "transcriptome")
  for (s in c(0L, 1L, 123L, 2147483000L)) {
    for (ty in types) {
      seeds <- vapply(paste0(ty, "_rep", 1:500), .layer_seed, 0L, seed = s)
      expect_identical(anyDuplicated(seeds), 0L)
    }
    for (fam in c("residual_t", "vqtl_residual_t", "complex_resid_t")) {
      seeds <- vapply(paste0(fam, 1:100), .layer_seed, 0L, seed = s)
      expect_identical(anyDuplicated(seeds), 0L)
    }
    grid <- expand.grid(type = types, occ = 0:20, stringsAsFactors = FALSE)
    seeds <- mapply(function(ty, o) .layer_seed(s, ty, o), grid$type, grid$occ)
    expect_identical(anyDuplicated(seeds), 0L)
  }
})

test_that(".layer_seed keeps NULL passthrough and is safe for the largest seed", {
  expect_null(.layer_seed(NULL, "additive", 0L))
  big <- .layer_seed(.Machine$integer.max, "residual_t99", 20L)
  expect_true(is.integer(big) && !is.na(big) && big >= 0L)
})

test_that("vary_qtn replications 12 and 21 draw different QTNs and effects (EFF-F1)", {
  M <- .ag_matrix()
  ph <- simulate_phenotype(M, n_reps = 21, vary_qtn = TRUE, seed = 1) |>
    additive(prop = 0.3, n_qtn = 5)
  ly <- ph$layers[[1]]
  key <- vapply(seq_len(21), function(r) {
    paste(ly$qtn_reps[[r]][[1]], collapse = ",")
  }, "")
  expect_identical(anyDuplicated(key), 0L)
  expect_false(identical(ly$qtn_reps[[12]], ly$qtn_reps[[21]]))
  # pleiotropy engine: same guard
  pl <- simulate_phenotype(M, architecture = "pleiotropy", n_traits = 2,
                           cor = 0.3, n_reps = 21, vary_qtn = TRUE, seed = 1) |>
    additive(prop = 0.3, n_qtn = 10)
  ply <- pl$layers[[1]]
  pkey <- vapply(seq_len(21), function(r) {
    paste(ply$qtn_reps[[r]][[1]], collapse = ",")
  }, "")
  expect_identical(anyDuplicated(pkey), 0L)
})

test_that("traits 12 and 21 get different residuals (O1)", {
  M <- .ag_matrix(n = 60, m = 50)
  ph <- simulate_phenotype(M, n_traits = 21, seed = 123)
  w <- phenotypes_wide(ph)
  expect_lt(abs(stats::cor(w$Trait_12, w$Trait_21)), 0.9)
  tr <- w[, paste0("Trait_", 1:21)]
  expect_identical(anyDuplicated(as.list(tr)), 0L)
})

test_that("the sub-seed rule is (seed, layer type, occurrence of that type) (O5/F6)", {
  M <- .ag_matrix()
  base <- simulate_phenotype(M, seed = 77)
  one <- additive(base, prop = 0.2, n_qtn = 5)
  expect_identical(
    one$layers[[1]]$qtn[[1]],
    .draw_qtn(base, 5L, .layer_seed(77L, "additive", 0L))[[1]])
  # a different-type layer inserted before it changes nothing (reordering
  # invariance across types) ...
  with_dom <- base |> dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 3) |>
    additive(prop = 0.2, n_qtn = 5)
  expect_identical(with_dom$layers[[2]]$qtn[[1]], one$layers[[1]]$qtn[[1]])
  # ... but the second additive layer is keyed by occurrence 1, so inserting a
  # same-type layer shifts the later same-type layers' draws (documented).
  two <- base |> additive(prop = 0.1, n_qtn = 3) |> additive(prop = 0.2, n_qtn = 5)
  expect_identical(
    two$layers[[2]]$qtn[[1]],
    .draw_qtn(base, 5L, .layer_seed(77L, "additive", 1L))[[1]])
  expect_false(identical(two$layers[[2]]$qtn[[1]], one$layers[[1]]$qtn[[1]]))
})

# ---------------------------------------------------------------------------
# F1 (D3): A + D on shared loci -- realized Var(A), Var(D), 2Cov(A,D)
# ---------------------------------------------------------------------------
.ag_ad_pair <- function(seed) {
  M <- .ag_matrix(n = 1000, m = 400, p_lo = 0.10, p_hi = 0.20, seed = 5)
  df <- .ag_df(M)
  flip <- df
  flip[, -(1:5)] <- -df[, -(1:5)]                 # same loci, opposite +1 allele
  mk <- function(x) simulate_phenotype(x, h2 = 0.5, seed = seed) |>
    additive(prop = 0.4, n_qtn = 20) |> dominance(prop = 0.1)
  list(orig = mk(df), flipped = mk(flip))
}

test_that("A+D on shared loci reports the realized partition and its coding dependence (F1)", {
  res <- lapply(1:3, .ag_ad_pair)
  for (r in res) {
    expect_false(is.null(r$orig$ad_report))
    ar <- r$orig$ad_report
    expect_true(all(c("trait", "requested", "realized", "var_A", "var_D",
                      "cov2_AD", "cov2_comp") %in% names(ar)))
    expect_equal(ar$requested, 0.5)
    # the partition closes on the realized genetic share of V_P
    expect_equal(ar$var_A + ar$var_D + ar$cov2_AD, ar$realized, tolerance = 1e-8)
    # the requested budget stays honest (unchanged), the realized one is reported
    expect_equal(r$orig$var_budget$prop, r$flipped$var_budget$prop)
  }
  # flipping the +1 allele on the same loci/seeds moves the realized quantities
  h2 <- vapply(res, function(r) c(r$orig$ad_report$realized,
                                  r$flipped$ad_report$realized), numeric(2))
  expect_gt(mean(h2[1, ]) - mean(h2[2, ]), 0.2)      # ~0.63 vs ~0.25 (audit)
  vA <- vapply(res, function(r) c(r$orig$ad_report$var_A,
                                  r$flipped$ad_report$var_A), numeric(2))
  expect_gt(mean(abs(vA[1, ] - vA[2, ])), 0.05)
  c2 <- vapply(res, function(r) c(r$orig$ad_report$cov2_comp,
                                  r$flipped$ad_report$cov2_comp), numeric(2))
  expect_gt(mean(c2[1, ]) - mean(c2[2, ]), 0.2)      # opposite-signed cross term
  # the reported realized share is the print's realized H2
  r1 <- res[[1]]$orig
  expect_equal(r1$ad_report$realized, .realized_h2(r1), tolerance = 1e-8)
})

test_that("print recommends the orthogonal model when A+D share loci (F1)", {
  r <- .ag_ad_pair(1)$orig
  out <- paste(utils::capture.output(print(r)), collapse = "\n")
  expect_match(out, "orthogonal = TRUE", fixed = TRUE)
  expect_match(out, "Var(A)", fixed = TRUE)
  expect_match(out, "2Cov(A,D)", fixed = TRUE)
  # a plain additive model prints no such note
  a <- simulate_phenotype(.ag_matrix(), seed = 1) |> additive(prop = 0.3, n_qtn = 4)
  expect_null(a$ad_report)
  expect_no_match(paste(utils::capture.output(print(a)), collapse = "\n"),
                  "orthogonal = TRUE", fixed = TRUE)
})

test_that("A+D on disjoint loci reports nothing; the orthogonal model reports its own split (F1)", {
  M <- .ag_matrix(n = 300)
  d <- simulate_phenotype(M, seed = 3) |> additive(prop = 0.3, qtn = 1:4) |>
    dominance(prop = 0.1, qtn = 11:14)
  expect_null(d$ad_report)
  o <- simulate_phenotype(M, seed = 3) |>
    additive(prop = 0.4, orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 5)
  expect_null(o$ad_report)
  expect_true("add_dom_cov" %in% o$var_budget$component)
})

test_that("the documentation no longer calls the A+D bias finite-sample (F1)", {
  skip_if_no_source("R", "grammar_simulate_phenotype.R")
  txt <- paste(readLines(testthat::test_path("..", "..", "R",
                                             "grammar_simulate_phenotype.R")),
               collapse = "\n")
  skip_if(nchar(txt) == 0L)
  expect_no_match(txt, "in a finite sample, especially", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# O2 / O6 / F3 / O9: complex_phenotypes()
# ---------------------------------------------------------------------------
.ag_complex <- function(nt = 2L, h2 = 0.4) {
  M <- .ag_matrix(n = 150, m = 200)
  a <- simulate_phenotype(M, n_traits = nt, seed = 10) |>
    additive(prop = h2, n_qtn = 4)
  b <- simulate_phenotype(M, n_traits = nt, seed = 11) |>
    additive(prop = h2, n_qtn = 4)
  list(a = a, b = b,
       z = suppressWarnings(complex_phenotypes(a, b, h2 = h2)))   # differing seeds
}

test_that("a layer added to a complex model errors instead of being ignored (O2)", {
  for (nt in 1:2) {
    z <- .ag_complex(nt)$z
    expect_error(additive(z, prop = 0.1, n_qtn = 3), "complex")
    expect_error(dominance(z, prop = 0.1), "complex")
    expect_error(epistasis(z, prop = 0.1, n_pairs = 2), "complex")
    expect_error(vqtl(z, prop = 0.1), "complex")
  }
})

test_that("complex_phenotypes() rejects h2-incomplete inputs like every accessor (O6)", {
  M <- .ag_matrix(n = 100, m = 100)
  s1 <- simulate_phenotype(M, h2 = 0.8, seed = 1) |> additive(prop = 0.1, n_qtn = 3)
  s2 <- simulate_phenotype(M, h2 = 0.8, seed = 2) |> additive(prop = 0.8, n_qtn = 3)
  expect_error(phenotypes_long(s1), "Incomplete h2 allocation")
  expect_error(complex_phenotypes(s1, s2, h2 = 0.5), "Incomplete h2 allocation")
})

test_that("a complex object carries no stale model-1 state (F3)", {
  z <- .ag_complex(2L)$z
  expect_null(z$mediation)
  expect_null(mediation_split(z))
  expect_false(isTRUE(z$one_call))
  expect_null(z$ad_report)
  expect_null(z$expression)
  # one-call input: the exhausted-budget hint must not leak
  M <- .ag_matrix(n = 150, m = 200)
  o1 <- simulate_phenotype(M, h2 = 0.4, n_qtn = 3, seed = 4)
  o2 <- simulate_phenotype(M, h2 = 0.4, n_qtn = 3, seed = 5)
  zz <- suppressWarnings(complex_phenotypes(o1, o2, h2 = 0.4))
  expect_false(isTRUE(zz$one_call))
  expect_identical(zz$n_qtn, 0L)
})

test_that("a complex model with a derived transcriptome input reports no mediation split (F3)", {
  tx <- simulate_transcriptome(G, n_genes = 60, seed = 1)
  m1 <- simulate_phenotype(G, h2 = 0.3, seed = 2, transcriptome = tx) |>
    transcriptome(prop = 0.4, n_genes = 15) |> additive(prop = 0.3, n_qtn = 3)
  expect_false(is.null(mediation_split(m1)))
  m2 <- simulate_phenotype(G, h2 = 0.3, seed = 3) |> additive(n_qtn = 3)
  zz <- suppressWarnings(complex_phenotypes(m1, m2, h2 = 0.5))
  expect_null(mediation_split(zz))
})

test_that("a one-trait complex model is a well-formed matrix and plots its own target (O9)", {
  z1 <- .ag_complex(1L, h2 = 0.5)$z
  gm <- .genetic_matrix(z1, 1L)
  expect_true(is.matrix(gm))
  expect_identical(dim(gm), c(z1$n_ind, 1L))
  gv <- genetic_values(z1)
  expect_identical(dim(gv), c(z1$n_ind, 1L))
  expect_equal(.plot_genetic_target(z1), 0.5)
  expect_equal(.plot_genetic_target(simulate_phenotype(.ag_matrix(), seed = 1) |>
                                      additive(prop = 0.3, n_qtn = 3)), 0.3)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp)
  expect_no_error(plot(z1, which = "cor"))
  grDevices::dev.off()
  unlink(tmp)
})

test_that("complex genetic value is sqrt(h2) * scale(G1 + G2) (Codex proposal 10)", {
  parts <- .ag_complex(2L, h2 = 0.4)
  z <- parts$z
  raw <- genetic_values(parts$a) + genetic_values(parts$b)
  for (t in 1:2) {
    expect_equal(as.numeric(genetic_values(z)[, t]),
                 as.numeric((raw[, t] - mean(raw[, t])) / stats::sd(raw[, t]) *
                              sqrt(0.4)), tolerance = 1e-10)
    expect_equal(stats::var(genetic_values(z)[, t]), 0.4, tolerance = 1e-10)
  }
})

# ---------------------------------------------------------------------------
# O3: minimum number of individuals
# ---------------------------------------------------------------------------
test_that("fewer than three individuals is an error, never NaN / 1e31 (O3)", {
  M2 <- matrix(c(-1, 1), 2, 1, dimnames = list(c("a", "b"), "m1"))
  for (s in 1:10) {
    expect_error(simulate_phenotype(M2, seed = s) |>
                   additive(prop = 0.5, qtn = "m1"), "at least three")
  }
  M <- .ag_matrix(n = 20, m = 30)
  expect_error(simulate_phenotype(M, individuals = 1:2, seed = 1), "at least three",
               ignore.case = TRUE)
  # three individuals is the documented minimum and runs
  ph <- simulate_phenotype(M, individuals = 1:3, seed = 1) |>
    additive(prop = 0.5, n_qtn = 3)
  expect_true(all(is.finite(ph$pheno$value)))
  expect_true(is.finite(.realized_h2(ph)))
})

test_that(".draw_residual rejects n < 2 and hits the exact variance at n = 2 (EFF-O6)", {
  expect_error(.draw_residual(0L, 0.5), "n >= 2")
  expect_error(.draw_residual(1L, 0.5), "n >= 2")
  e <- .draw_residual(2L, 0.5)
  expect_equal(mean(e), 0)
  expect_equal(stats::var(e), 0.5)
  set.seed(9); before <- .Random.seed
  expect_equal(.draw_residual(10L, 0), rep(0, 10))
  expect_identical(before, .Random.seed)                # no RNG at resid_var = 0
})

# ---------------------------------------------------------------------------
# O4 / F7: epistasis(qtn = list(...)) is one element per trait
# ---------------------------------------------------------------------------
test_that("epistasis(qtn = list) gives one set matrix per trait, never a replicated cor = 1 (O4/F7)", {
  M <- .ag_matrix(n = 200, m = 60)
  base <- simulate_phenotype(M, n_traits = 2, seed = 1)
  per_trait <- list(matrix(c(1, 2, 3, 4), 2, byrow = TRUE),
                    matrix(c(5, 6, 7, 8), 2, byrow = TRUE))
  ep <- epistasis(base, prop = 0.3, qtn = per_trait)
  expect_equal(ep$layers[[1]]$qtn[[1]], per_trait[[1]])
  expect_equal(ep$layers[[1]]$qtn[[2]], per_trait[[2]])
  g <- genetic_values(ep)
  expect_lt(abs(stats::cor(g[, 1], g[, 2])), 0.9)
  # the old reading (a list of sets) is now refused, not silently replicated
  expect_error(epistasis(base, prop = 0.3, qtn = list(c(1, 2), c(3, 4), c(5, 6))),
               "one element per trait")
  # a single matrix is shared by every trait (documented)
  sh <- epistasis(base, prop = 0.3, qtn = matrix(c(1, 2, 3, 4), 2, byrow = TRUE))
  expect_identical(sh$layers[[1]]$qtn[[1]], sh$layers[[1]]$qtn[[2]])
  # per-trait sets must have the same number of sets
  expect_error(epistasis(base, prop = 0.3,
                         qtn = list(matrix(1:4, 2), c(5, 6))), "same number")
})

# ---------------------------------------------------------------------------
# F2: constant-dosage (all-heterozygous) loci are not QTN candidates
# ---------------------------------------------------------------------------
test_that("all-heterozygous loci are excluded from the candidate pool (F2)", {
  M <- .ag_matrix(n = 50, m = 20)
  M[, 1] <- 0                                       # every individual heterozygous
  sim <- simulate_phenotype(M, seed = 1)
  expect_equal(sim$maf[1], 0.5)
  expect_false(1L %in% .candidate_markers(sim))
  expect_true(all(2:20 %in% .candidate_markers(sim)))
  for (s in 1:15) {
    ph <- simulate_phenotype(M, seed = s) |> additive(prop = 0.3, n_qtn = 5)
    expect_false(1L %in% ph$layers[[1]]$qtn[[1]])
  }
  # an F1-like panel (all heterozygous): accurate message, not "polymorphic"
  F1 <- matrix(0, 30, 10, dimnames = list(paste0("i", 1:30), paste0("m", 1:10)))
  expect_error(simulate_phenotype(F1, seed = 1) |> additive(prop = 0.3, n_qtn = 2),
               "constant|heterozygous in every")
  # dominance on a constant het indicator is refused with the specific message
  M2 <- .ag_matrix(n = 50, m = 20)
  M2[, 3] <- 0
  expect_error(suppressWarnings(simulate_phenotype(M2, seed = 1) |>
                 additive(prop = 0.3, qtn = 4) |> dominance(prop = 0.1, qtn = 3)),
               "no heterozygous|every individual")
})

test_that("polymorphic panels keep the same candidate set (F2 behaviour-preserving)", {
  M <- .ag_matrix(n = 100, m = 40)
  sim <- simulate_phenotype(M, seed = 1)
  expect_identical(.candidate_markers(sim), which(sim$maf > 0))
})

# ---------------------------------------------------------------------------
# F4 / F8 / O7: qtn= overrides same_as_add; ignored arguments warn
# ---------------------------------------------------------------------------
test_that("qtn= wins over same_as_add and the layer says so (F4/O7)", {
  M <- .ag_matrix(n = 300, m = 60)
  sim <- simulate_phenotype(M, seed = 1)
  d <- dominance(sim, prop = 0.1, qtn = c("m1", "m2"))
  expect_false(isTRUE(d$layers[[1]]$same_as_add))
  out <- paste(utils::capture.output(print(d)), collapse = "\n")
  expect_no_match(out, "same QTNs as additive", fixed = TRUE)
  v <- suppressMessages(vqtl(sim, prop = 0.1, qtn = c("m1", "m2")))
  expect_false(isTRUE(v$layers[[1]]$same_as_add))
  a <- additive(sim, prop = 0.2, qtn = 5:8)
  d2 <- dominance(a, prop = 0.1, qtn = 1:3)
  expect_identical(d2$layers[[2]]$qtn[[1]], 1:3)
  expect_false(isTRUE(d2$layers[[2]]$same_as_add))
  # reusing the additive loci still says so
  d3 <- dominance(a, prop = 0.1)
  expect_true(isTRUE(d3$layers[[2]]$same_as_add))
  expect_match(paste(utils::capture.output(print(d3)), collapse = "\n"),
               "same QTNs as additive", fixed = TRUE)
})

test_that("an n_qtn that is ignored produces a warning (F8)", {
  M <- .ag_matrix(n = 300, m = 60)
  sim <- simulate_phenotype(M, seed = 1)
  expect_warning(additive(sim, prop = 0.2, n_qtn = 10, qtn = 1:3), "ignored")
  a <- additive(sim, prop = 0.2, qtn = 1:3)
  expect_warning(dominance(a, prop = 0.1, same_as_add = TRUE, n_qtn = 9), "ignored")
  expect_warning(suppressMessages(vqtl(a, prop = 0.1, n_qtn = 9)), "ignored")
  # no spurious warning when n_qtn is actually used
  expect_no_warning(additive(sim, prop = 0.2, n_qtn = 3))
})

# ---------------------------------------------------------------------------
# F5: partially spent multi-trait budget
# ---------------------------------------------------------------------------
test_that("prop = NULL fills only the unspent traits when the budget is partial (F5)", {
  M <- .ag_matrix(n = 400, m = 100)
  ph <- simulate_phenotype(M, n_traits = 2, h2 = c(0.5, 0.3), seed = 1) |>
    additive(prop = c(0.2, 0.3), n_qtn = 4) |>
    dominance(same_as_add = TRUE)
  expect_equal(ph$layers[[2]]$prop, c(0.3, 0))
  # fully spent for every trait still errors
  expect_error(
    simulate_phenotype(M, n_traits = 2, h2 = c(0.5, 0.3), seed = 1) |>
      additive(prop = c(0.5, 0.3), n_qtn = 4) |> dominance(),
    "No heritability budget")
})

# ---------------------------------------------------------------------------
# O8: prop = 0 layers on hetless loci are clean no-ops
# ---------------------------------------------------------------------------
test_that("prop = 0 layers on hetless loci do not error or warn (O8)", {
  M <- .ag_matrix(n = 100, m = 20)
  M[, 1:10] <- ifelse(M[, 1:10] == 0, 1, M[, 1:10])   # markers 1-10: no hets
  sim <- simulate_phenotype(M, seed = 1)
  expect_no_error(d <- dominance(sim, prop = 0, qtn = "m1"))
  expect_no_error(o <- additive(sim, prop = 0, orthogonal = TRUE, a = 1, d = 1,
                                qtn = "m1"))
  expect_no_warning(e <- epistasis(sim, prop = 0, qtn = matrix(c(1, 2), 1),
                                   interaction_type = "d"))
  expect_equal(phenotypes_long(d)$value, phenotypes_long(sim)$value)
  # a positive prop on hetless loci still errors as before
  expect_error(dominance(sim, prop = 0.1, qtn = "m1"), "no heterozygous")
})

# ---------------------------------------------------------------------------
# F9: degenerate-locus messages
# ---------------------------------------------------------------------------
test_that("degenerate loci give accurate messages (F9)", {
  M <- .ag_matrix(n = 100, m = 20)
  sim <- simulate_phenotype(M, seed = 1)
  err <- tryCatch(additive(sim, prop = 0.3, n_qtn = 3, effect = c(0, 0, 0)),
                  error = function(e) conditionMessage(e))
  expect_match(err, "zero effects|constant", ignore.case = TRUE)
  expect_no_match(err, "Choose polymorphic loci", fixed = TRUE)
  err0 <- tryCatch(additive(sim, prop = 0.3, n_qtn = 3, effect = 0),
                   error = function(e) conditionMessage(e))
  expect_match(err0, "every effect is zero", fixed = TRUE)
  # epistasis with every "d" position hetless: specific error, not "absorb"
  Mh <- M
  Mh[Mh == 0] <- 1
  sim2 <- simulate_phenotype(Mh, seed = 1)
  expect_error(epistasis(sim2, prop = 0.1, n_pairs = 2, interaction_type = "d"),
               "no heterozygous")
})

# ---------------------------------------------------------------------------
# EFF-F6 / O2 / F5: effect series domain
# ---------------------------------------------------------------------------
test_that("effect series overflow, underflow and a zero base fail with a series message (EFF-F6/O2)", {
  expect_error(.effect_series(1100L, effect = 2), "overflow|not finite")
  expect_error(.effect_series(1100L, effect = 0.5), "underflow")
  # a zero base is the documented all-zero series (orthogonal a = 0 = pure
  # dominance); it is rejected downstream only where it cannot realize `prop`
  expect_equal(.effect_series(3L, effect = 0), c(0, 0, 0))
  ok <- .effect_series(50L, effect = 0.5)
  expect_true(all(is.finite(ok)) && all(ok != 0))
  err <- tryCatch(
    additive(simulate_phenotype(.ag_matrix(n = 60, m = 1200), seed = 1),
             prop = 0.4, n_qtn = 1100, effect = 2),
    error = function(e) conditionMessage(e))
  expect_no_match(err, "polymorphic", fixed = TRUE)
})

test_that("an invalid dist is rejected even when a custom series is supplied (EFF-F5)", {
  expect_error(.effect_series(3L, dist = "uniform", effect = 1:3), "geometric")
  expect_error(.effect_series(3L, dist = "uniform"), "geometric")
  expect_equal(.effect_series(3L, effect = c(1, 2, 3)), c(1, 2, 3))
})

# ---------------------------------------------------------------------------
# EFF-F3 / F2: the "ld" partner rule and r2 window guards
# ---------------------------------------------------------------------------
test_that("ld window guards: r2_max = 1 never pairs identical columns, degenerate windows error (EFF-F3)", {
  for (s in 1:6) {
    ph <- simulate_phenotype(G, architecture = "ld", n_traits = 2,
                             ld_type = "direct", r2_min = 0.95, r2_max = 1,
                             n_qtn = 3, seed = s) |> additive(prop = 0.4)
    ld <- ph$layers[[1]]$ld
    expect_true(all(ld$r2 < 1 - 1e-12))
    for (i in seq_len(nrow(ld))) {
      x <- .geno_cols(ph, c(ld$qtn_t1[i], ld$qtn_t2[i]))
      expect_false(identical(unname(x[, 1]), unname(x[, 2])))
    }
  }
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 2,
                                  r2_min = 0, r2_max = 0), "strictly between")
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 2,
                                  r2_min = 1, r2_max = 1), "strictly between")
})

test_that("ld partner = 'random' is opt-in; the default is the strongest in-window partner (EFF-F2)", {
  expect_error(simulate_phenotype(G, architecture = "ld", n_traits = 2,
                                  partner = "best"), "partner")
  strongest <- function(s, partner) {
    args <- list(G, architecture = "ld", n_traits = 2, ld_type = "direct",
                 r2_min = 0.2, r2_max = 0.8, n_qtn = 6, seed = s)
    if (!is.null(partner)) args$partner <- partner
    do.call(simulate_phenotype, args) |> additive(prop = 0.4)
  }
  d <- strongest(3, NULL)
  s <- strongest(3, "strongest")
  expect_identical(d$layers[[1]]$qtn, s$layers[[1]]$qtn)
  # every default pick equals the focal's best in-window partner
  ph <- d
  ld <- ph$layers[[1]]$ld
  chr <- ph$map$chr
  for (i in seq_len(nrow(ld))) {
    f <- ld$qtn_t1[i]
    on <- setdiff(which(chr == chr[f] & ph$maf > 0), f)
    r2v <- as.numeric(suppressWarnings(stats::cor(.geno_cols(ph, f),
                                                  .geno_cols(ph, on))))^2
    ok <- is.finite(r2v) & r2v >= 0.2 & r2v <= 0.8
    expect_equal(ld$r2[i], max(r2v[ok & !on %in% ld$qtn_t1[seq_len(i - 1)] &
                                     !on %in% ld$qtn_t2[seq_len(i - 1)]]))
  }
  r <- strongest(3, "random")
  expect_true(all(r$layers[[1]]$ld$r2 >= 0.2 & r$layers[[1]]$ld$r2 <= 0.8))
  mean_r2 <- function(partner) mean(vapply(1:6, function(s) {
    mean(strongest(s, partner)$layers[[1]]$ld$r2)
  }, 0))
  expect_lt(mean_r2("random"), mean_r2("strongest"))
})

# ---------------------------------------------------------------------------
# EFF-F4: pi_t = 1 / 0 traits carry no zero-effect QTNs
# ---------------------------------------------------------------------------
test_that("traits with no trait-specific variance report no zero-effect QTNs (EFF-F4)", {
  M <- .ag_matrix(n = 300, m = 200)
  ph <- simulate_phenotype(M, architecture = "pleiotropy", n_traits = 2,
                           cor = 0.3, pi = c(1, 0.5), seed = 2) |>
    additive(prop = 0.4, n_qtn = 20)
  ly <- ph$layers[[1]]
  expect_true(all(ly$effect[[1]] != 0))
  expect_true(all(ly$effect[[2]] != 0))
  expect_true(all(vapply(seq_len(2), function(t) length(ly$qtn[[t]]) ==
                           length(ly$effect[[t]]), TRUE)))
  tab <- qtn_table(ph)
  expect_true(all(tab$effect[tab$layer == "additive"] != 0))
  expect_true(all(is.finite(genetic_values(ph))))
})

# ---------------------------------------------------------------------------
# GRAM-O11 / EFF-O1: the corrected identity (prose, not code)
# ---------------------------------------------------------------------------
test_that("the non-additive raw component variance targets V_t, not Sigma_tt (EFF-O1)", {
  M <- .ag_matrix(n = 400, m = 500, p_lo = 0.3, p_hi = 0.5)
  vs <- vapply(1:20, function(s) {
    sim <- simulate_phenotype(M, architecture = "pleiotropy", n_traits = 2,
                              cor = 0.3, pi = 0.5, seed = s) |>
      additive(prop = 0.3, n_qtn = 40) |>
      dominance(prop = 0.2, same_as_add = FALSE, n_qtn = 40)
    ly <- sim$layers[[2]]
    stats::var(.component_raw(ly, sim, 1L, 1L))
  }, 0)
  # raw (unscaled) component variance is on the design scale; it is not the
  # requested Sigma_tt = pi * V = 0.1 (the source comment used to say so)
  expect_gt(mean(vs), 0.15)
})

# ---------------------------------------------------------------------------
# Codex proposals: candidate pool edge cases, .gabriel_blocks contract
# ---------------------------------------------------------------------------
test_that(".candidate_markers drops NA / NaN / Inf / zero-MAF markers (Codex proposal 11)", {
  sim <- list(maf = c(0.1, NA, NaN, Inf, 0, 0.5, 0.3),
              all_het = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE))
  expect_identical(.candidate_markers(sim), c(1L, 7L))
  # a sim built without the all_het field (older objects) still works
  sim$all_het <- NULL
  expect_identical(.candidate_markers(sim), c(1L, 6L, 7L))
})

test_that(".gabriel_blocks returns a clean logical for a huge block span and leaves rare markers alone", {
  d <- G[G$chr == 1, ][1:300, ]
  dose <- as.matrix(d[, -(1:5)]) + 1L
  p <- rowSums(dose) / (2 * ncol(dose))
  maf <- pmin(p, 1 - p)
  keep <- rep(TRUE, nrow(d))
  out <- .gabriel_blocks(dose, d$chr, d$pos, keep, maf, max_kb = 1e12)
  expect_type(out, "logical")
  expect_false(anyNA(out))
  expect_length(out, nrow(d))
  expect_true(all(out[maf < 0.05]))          # markers below the Haploview floor are untouched
})
