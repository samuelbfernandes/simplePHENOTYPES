# test-audit-v1-pleio.R
#
# Regression tests for the v1 pleiotropy / genetic-effect / vQTL / genotype
# loader audit (group v1-pleio). The frozen engine is not "fixed" numerically:
# these tests lock (a) the exact v1 conventions, (b) informative errors for
# inputs that used to give silently wrong results or a cryptic crash, and
# (c) file-content corrections.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
.v1_run <- function(...) {
  tmp <- tempfile("v1pleio_")
  dir.create(tmp)
  args <- c(list(geno_obj = get("SNP55K_maize282_maf04"), to_r = TRUE,
                 output_format = "long", home_dir = tmp, verbose = FALSE),
            list(...))
  out <- suppressWarnings(suppressMessages(do.call(create_phenotypes, args)))
  attr(out, "dir") <- tmp
  out
}

.v1_err <- function(...) {
  tryCatch({
    suppressWarnings(suppressMessages(.v1_run(...)))
    NA_character_
  }, error = function(e) conditionMessage(e))
}

# create_phenotypes() reports failures with message() and returns NULL (v1
# behaviour, top-level handler owned by the core group); this returns the text
# that reached the message stream, or NA when the run completed.
.v1_msg <- function(...) {
  msgs <- character(0)
  res <- withCallingHandlers(
    tryCatch(suppressWarnings(.v1_run(...)), error = function(e) {
      msgs <<- c(msgs, conditionMessage(e)); NULL
    }),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  list(result = res, text = paste(msgs, collapse = "\n"))
}

.v1_file <- function(out, name) {
  hits <- list.files(attr(out, "dir"), pattern = paste0("^", name, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0L) return(NULL)
  data.table::fread(hits[[1L]], data.table = FALSE)
}

# synthetic marker panel: 5 meta columns + samples coded -1/0/1
.v1_panel <- function(m = 40, n = 24, het_markers = integer(0), seed = 1,
                      pos = NULL) {
  set.seed(seed)
  X <- matrix(sample(c(-1L, 1L), m * n, replace = TRUE), m, n)
  for (k in het_markers) X[k, sample(n, 4)] <- 0L
  colnames(X) <- paste0("s", seq_len(n))
  data.frame(snp = paste0("m", seq_len(m)), allele = "A/G", chr = 1L,
             pos = if (is.null(pos)) seq_len(m) * 100L else pos, cm = 0,
             X, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# V1P-F1 make_pd
# ---------------------------------------------------------------------------
test_that("make_pd: frozen eigenvalue clamp is documented arithmetic, not a correlation matrix", {
  m <- matrix(c(1, .9, -.9, .9, 1, .9, -.9, .9, 1), 3)
  out <- suppressWarnings(make_pd(m, verbose = FALSE))
  expect_equal(unname(diag(out)), rep(1.27, 3))
  expect_equal(out[1, 2], 0.63)
  expect_equal(out[1, 3], -0.63)
  expect_equal(min(eigen(out)$values), 0.01, tolerance = 1e-8)
  # realized correlation is not the requested 0.9
  expect_equal(abs(out[1, 2]) / out[1, 1], 0.63 / 1.27)
  expect_no_error(chol(out))
})

test_that("make_pd: positive-definite input is returned bit-identical", {
  m <- matrix(c(1, .5, .5, 1), 2)
  expect_identical(make_pd(m, verbose = FALSE), m)
})

test_that("make_pd: singular input fails with an informative error, not a chol() abort", {
  m <- matrix(1, 2, 2)
  expect_error(suppressWarnings(make_pd(m, verbose = FALSE)),
               "not positive definite", fixed = FALSE)
  expect_error(suppressWarnings(make_pd(m, verbose = FALSE,
                                        what = "sample covariance of the standardized genetic values")),
               "perfectly collinear")
})

test_that("make_pd: whatever it returns can be Cholesky-factorized (200 random matrices)", {
  set.seed(42)
  bad_msg <- 0L
  for (it in 1:200) {
    A <- matrix(runif(9, -1, 1), 3)
    m <- cov2cor(crossprod(A) + diag(0.001, 3))
    m[] <- pmin(pmax(m + matrix(rnorm(9, 0, .3), 3) %*% diag(3), -1), 1)
    m <- (m + t(m)) / 2
    diag(m) <- 1
    res <- tryCatch(suppressWarnings(make_pd(m, verbose = FALSE)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      expect_match(conditionMessage(res), "not positive definite")
      bad_msg <- bad_msg + 1L
    } else {
      expect_no_error(chol(res))
    }
  }
  expect_gt(bad_msg, 0L)  # the failure mode really occurs and is reported
})

test_that("make_pd: asymmetric input is rejected", {
  m <- matrix(c(1, .5, .2, 1), 2)
  expect_error(make_pd(m, verbose = FALSE), "not symmetric")
})

# ---------------------------------------------------------------------------
# V1P-F2 cor imposed on genetic values (documented behaviour is locked)
# ---------------------------------------------------------------------------
test_that("cor: exact realization and variance preservation on genetic values", {
  set.seed(1)
  n <- 200
  X <- matrix(sample(c(-1, 0, 1), n * 6, replace = TRUE), n, 6,
              dimnames = list(paste0("i", 1:n), NULL))
  colnames(X) <- paste0("Chr_1_", 1:6)
  add_obj <- list(X)
  eff <- list(c(.2, .04, .008, .5, .3, .1), c(.1, .3, .05, .02, .4, .2))
  cr <- matrix(c(1, .5, .5, 1), 2)
  res <- base_line_multi_traits(
    add_obj = add_obj, dom_obj = list(NULL), epi_obj = list(NULL),
    add_effect = eff, dom_effect = list(NULL, NULL),
    epi_effect = list(NULL, NULL), epi_interaction = 2, ntraits = 2,
    cor = cr, architecture = "pleiotropic", rep = 1, rep_by = "experiment",
    add = TRUE, dom = FALSE, epi = FALSE, sim_method = "custom",
    verbose = FALSE)
  g <- res[[1]]$base_line
  expect_equal(cor(g)[1, 2], 0.5, tolerance = 1e-10)
  raw1 <- scale(X %*% eff[[1]], scale = FALSE)
  raw2 <- scale(X %*% eff[[2]], scale = FALSE)
  expect_equal(unname(g[, 1]), unname(raw1[, 1]), tolerance = 1e-10)
  expect_equal(var(g[, 2]), var(raw2[, 1]), tolerance = 1e-10)
  # trait 2 is NOT its own QTN model any more (documented leakage)
  expect_false(isTRUE(all.equal(unname(g[, 2]), unname(raw2[, 1]))))
})

test_that("cor: size mismatch and zero-variance traits give informative errors", {
  set.seed(2)
  n <- 30
  X <- matrix(sample(c(-1, 0, 1), n * 3, replace = TRUE), n, 3)
  colnames(X) <- paste0("Chr_1_", 1:3)
  args <- list(add_obj = list(X), dom_obj = list(NULL), epi_obj = list(NULL),
               dom_effect = list(NULL, NULL), epi_effect = list(NULL, NULL),
               epi_interaction = 2, ntraits = 2, architecture = "pleiotropic",
               rep = 1, rep_by = "experiment", add = TRUE, dom = FALSE,
               epi = FALSE, sim_method = "custom", verbose = FALSE)
  expect_error(do.call(base_line_multi_traits,
                       c(args, list(add_effect = list(c(.1, .2, .3), c(.3, .2, .1)),
                                    cor = 0.5))),
               "ntraits|2 x 2", ignore.case = TRUE)
  expect_error(do.call(base_line_multi_traits,
                       c(args, list(add_effect = list(c(.1, .2, .3), c(0, 0, 0)),
                                    cor = matrix(c(1, .5, .5, 1), 2)))),
               "zero .* genetic variance")
  # identical effect series: perfectly collinear traits -> informative error
  expect_error(suppressWarnings(
    do.call(base_line_multi_traits,
            c(args, list(add_effect = list(c(.1, .2, .3), c(.1, .2, .3)),
                         cor = matrix(c(1, .5, .5, 1), 2))))),
    "perfectly collinear")
})

# ---------------------------------------------------------------------------
# V1P-F3 / F4 effect vectors; V1P-F14 conventions; O7 positional epistasis
# ---------------------------------------------------------------------------
test_that("genetic_effect: hand-computed A + D + E with the v1 conventions", {
  Xa <- cbind(a1 = c(-1, 0, 1, 1, 0, -1, 1, 0, 1, -1, 0, 1),
              a2 = c(1, 1, 0, -1, 0, 0, -1, 1, 1, 0, -1, 1))
  Xd <- cbind(d1 = c(0, 1, -1, 0, 1, 1, -1, 0, 0, 1, 1, -1),
              d2 = c(1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1))  # hetless
  Xe <- cbind(e1 = c(1, 0, -1, 1, 1, -1, 0, 1, -1, 1, 0, 1),
              e2 = c(-1, 1, 1, 0, 1, -1, 1, -1, 1, 0, 1, 1))
  rownames(Xa) <- rownames(Xd) <- rownames(Xe) <- paste0("i", 1:12)
  a <- c(.4, .16); d <- c(.3, .2); e <- .25
  res <- suppressWarnings(genetic_effect(
    add_obj = Xa, dom_obj = Xd, epi_obj = Xe, add_effect = a,
    dom_effect = d, epi_effect = e, epi_interaction = 2,
    add = TRUE, dom = TRUE, epi = TRUE))
  manual <- Xa %*% a + (Xd[, 1] == 0) * d[1] + (Xd[, 2] == 0) * d[2] +
    Xe[, 1] * Xe[, 2] * e            # uncentred product
  manual <- manual - mean(manual)
  expect_equal(unname(unlist(res$base_line)), as.numeric(manual),
               tolerance = 1e-12)
  expect_equal(res$VA, var(as.numeric(Xa %*% a)), tolerance = 1e-12)
  expect_equal(res$var_dom[2], 0)                       # hetless dominance QTN
  expect_equal(res$var_epi, var(Xe[, 1] * Xe[, 2] * e), tolerance = 1e-12)
})

test_that("genetic_effect: effect vectors are never recycled", {
  X <- cbind(c(-1, 0, 1, 1), c(1, 1, 0, -1), c(0, 1, -1, 1))
  rownames(X) <- letters[1:4]
  expect_error(genetic_effect(add_obj = X, add_effect = c(.2, .09, .008, .3, .04, .027),
                              add = TRUE, dom = FALSE, epi = FALSE),
               "additive effect vector has 6")
  expect_error(genetic_effect(add_obj = X, add_effect = c(.2, .3),
                              add = TRUE, dom = FALSE, epi = FALSE),
               "effect vectors are never recycled")
  expect_error(genetic_effect(dom_obj = X, dom_effect = .2,
                              add = FALSE, dom = TRUE, epi = FALSE),
               "dominance effect vector has 1")
  expect_error(genetic_effect(epi_obj = X[, 1:2], epi_effect = c(.1, .2),
                              epi_interaction = 2, add = FALSE, dom = FALSE,
                              epi = TRUE),
               "epistatic effect vector has 2")
})

test_that("create_phenotypes: a length-2 effect for 3 QTNs is rejected (was silently recycled)", {
  r <- .v1_msg(add_QTN_num = 3, add_effect = c(.2, .3), h2 = 1, model = "A",
               rep = 1, seed = 5)
  expect_null(r$result)
  # rejected either up front (argument validation) or by genetic_effect()
  expect_match(r$text, "effect")
})

test_that("genetic_effect: epistasis groups columns by position, duplicate names are harmless (O7)", {
  X <- cbind(c(-1, 1, 1, -1), c(1, -1, 1, -1))
  colnames(X) <- c("Chr_1_10", "Chr_1_10")             # duplicated names
  rownames(X) <- letters[1:4]
  res <- genetic_effect(epi_obj = X, epi_effect = 1, epi_interaction = 2,
                        add = FALSE, dom = FALSE, epi = TRUE)
  # product of the two POSITIONAL columns, centred: (-1, -1, 1, 1); the
  # name-indexed version multiplied column 1 by itself (all ones -> 0)
  expect_equal(unname(unlist(res$base_line)), c(-1, -1, 1, 1))
})

test_that("genetic_effect: epi_interaction is required for an epistatic component", {
  X <- cbind(c(-1, 1, 1, -1), c(1, -1, 1, -1))
  rownames(X) <- letters[1:4]
  expect_error(genetic_effect(epi_obj = X, epi_effect = 1, epi_interaction = NULL,
                              add = FALSE, dom = FALSE, epi = TRUE),
               "epi_interaction")
})

test_that("genetic_effect: NA genotypes at QTNs are rejected instead of giving NA genetic values", {
  X <- cbind(c(-1, NA, 1, 1))
  rownames(X) <- letters[1:4]
  expect_error(genetic_effect(add_obj = X, add_effect = .2, add = TRUE,
                              dom = FALSE, epi = FALSE), "missing values")
})

test_that("genetic_effect: A + D with hetless dominance loci warns (dominance is silently zero otherwise)", {
  Xa <- cbind(c(-1, 0, 1, 1, 0, -1))
  Xd <- cbind(c(-1, 1, 1, -1, 1, -1))
  rownames(Xa) <- rownames(Xd) <- letters[1:6]
  expect_warning(genetic_effect(add_obj = Xa, dom_obj = Xd, add_effect = .3,
                                dom_effect = .2, add = TRUE, dom = TRUE,
                                epi = FALSE), "heterozyg")
  # zero-count dummy (effect 0) does not warn
  expect_no_warning(genetic_effect(add_obj = Xa, dom_obj = Xd, add_effect = .3,
                                   dom_effect = 0, add = TRUE, dom = TRUE,
                                   epi = FALSE))
})

test_that("epi_QTN_num = 0 with an additive model runs (dummy epistatic locus is zeroed)", {
  out <- .v1_run(add_QTN_num = 2, add_effect = .2, epi_QTN_num = 0,
                 epi_effect = .2, h2 = .5, model = "AE", rep = 1, seed = 3)
  expect_true(is.data.frame(out))
})

# ---------------------------------------------------------------------------
# V1P-F6 partial pleiotropy: overlapping trait-specific sets are detected
# ---------------------------------------------------------------------------
.partial_call <- function(G, seed, same = TRUE, spec = c(2, 2, 2)) {
  nt <- length(spec)
  eff <- lapply(seq_len(nt), function(t) rep(.2, 1 + spec[t]))
  withr::with_tempdir(qtn_partially_pleiotropic(
    genotypes = G, seed = seed, pleio_a = 1, pleio_d = 1,
    trait_spec_a_QTN_num = spec, trait_spec_d_QTN_num = spec,
    add_effect = eff, dom_effect = eff, epi_effect = NULL, epi_interaction = 2,
    ntraits = nt, rep = 1, rep_by = "experiment", export_gt = FALSE,
    same_add_dom_QTN = same, add = TRUE, dom = TRUE, epi = FALSE,
    verbose = FALSE))
}

test_that("partial pleiotropy: het re-sampling either yields disjoint sets or a clear error", {
  G <- .v1_panel(m = 40, n = 24, het_markers = c(5, 17, 33))
  ok <- 0L; err <- 0L
  for (sd in 1:40) {
    r <- tryCatch(.partial_call(G, sd), error = function(e) e)
    if (inherits(r, "error")) {
      expect_match(conditionMessage(r), "more than one trait|not be partially pleiotropic")
      err <- err + 1L
    } else {
      ok <- ok + 1L
      for (rp in seq_along(r$add_ef_trait_obj)) {
        cols <- unlist(lapply(r$add_ef_trait_obj[[rp]], function(m) {
          colnames(m)
        }))
        # pleiotropic marker appears once per trait; trait-specific markers once
        tab <- table(cols)
        expect_equal(sum(tab > 1), 1L)          # only the shared marker repeats
        expect_equal(max(tab), 3L)              # ... in each of the 3 traits
      }
    }
  }
  expect_gt(err, 0L)   # the overlap does occur under the frozen sampling
  expect_gt(ok, 0L)
})

test_that("partial pleiotropy: additive-only (no dominance) selection is unchanged and disjoint", {
  G <- .v1_panel(m = 40, n = 24, het_markers = c(5, 17, 33))
  withr::local_dir(withr::local_tempdir())
  eff <- lapply(1:3, function(t) rep(.2, 3))
  r <- qtn_partially_pleiotropic(
    genotypes = G, seed = 3, pleio_a = 1, trait_spec_a_QTN_num = c(2, 2, 2),
    add_effect = eff, epi_interaction = 2, ntraits = 3, rep = 1,
    rep_by = "experiment", export_gt = FALSE, same_add_dom_QTN = FALSE,
    add = TRUE, dom = FALSE, epi = FALSE, verbose = FALSE)
  cols <- unlist(lapply(r$add_ef_trait_obj[[1]], colnames))
  expect_equal(sum(table(cols) > 1), 1L)
})

# ---------------------------------------------------------------------------
# V1P-F7 partial Epistatic_QTNs.txt effect column
# ---------------------------------------------------------------------------
test_that("partial: Epistatic_QTNs.txt reports each trait's own effects", {
  G <- .v1_panel(m = 60, n = 24)
  withr::local_dir(withr::local_tempdir())
  ee <- list(c(.1, .2, .3), c(.4, .5, .6))
  qtn_partially_pleiotropic(
    genotypes = G, seed = 4, pleio_e = 1, trait_spec_e_QTN_num = c(2, 2),
    epi_effect = ee, epi_interaction = 2, ntraits = 2, rep = 1,
    rep_by = "experiment", export_gt = FALSE, same_add_dom_QTN = FALSE,
    add = FALSE, dom = FALSE, epi = TRUE, verbose = FALSE)
  f <- data.table::fread("Epistatic_QTNs.txt", data.table = FALSE)
  t1 <- f[f$trait == "trait_1", "epistatic_effect"]
  t2 <- f[f$trait == "trait_2", "epistatic_effect"]
  expect_equal(t1, rep(c(.1, .2, .3), each = 2))
  expect_equal(t2, rep(c(.4, .5, .6), each = 2))
})

# ---------------------------------------------------------------------------
# V1P-F11 partial dominance seed files
# ---------------------------------------------------------------------------
test_that("partial: dominance seed files carry the right name and the seeds actually used", {
  G <- .v1_panel(m = 60, n = 24, het_markers = c(2, 9, 20, 31, 44, 50))
  withr::local_dir(withr::local_tempdir())
  eff <- lapply(1:2, function(t) rep(.2, 4))
  eff2 <- lapply(1:2, function(t) rep(.2, 3))
  suppressMessages(qtn_partially_pleiotropic(
    genotypes = G, seed = 5, pleio_a = 2, pleio_d = 3,
    trait_spec_a_QTN_num = c(2, 2), trait_spec_d_QTN_num = c(1, 1),
    add_effect = eff, dom_effect = list(rep(.2, 4), rep(.2, 4)),
    epi_interaction = 2, ntraits = 2, rep = 2, rep_by = "QTN",
    export_gt = FALSE, same_add_dom_QTN = FALSE, add = TRUE, dom = TRUE,
    epi = FALSE, verbose = TRUE))
  expect_true(file.exists("Seed_number_for_3Pleiotropic_Dom_QTN.txt"))
  expect_true(file.exists("Seed_number_for_1_1Trait_specific_Dom_QTN.txt"))
  expect_equal(scan("Seed_number_for_3Pleiotropic_Dom_QTN.txt", quiet = TRUE),
               5 + 1:2 + 2)                      # seed + j + rep
  # every (replicate, trait) seed is recorded, replicate-major: seed + i + j + rep
  expect_equal(scan("Seed_number_for_1_1Trait_specific_Dom_QTN.txt", quiet = TRUE),
               c(5 + 1 + 1 + 2, 5 + 2 + 1 + 2, 5 + 1 + 2 + 2, 5 + 2 + 2 + 2))
  expect_equal(scan("Seed_number_for_2_2Trait_specific_Add_QTN.txt", quiet = TRUE),
               c(5 + 1 + 1, 5 + 2 + 1, 5 + 1 + 2, 5 + 2 + 2))
})

# ---------------------------------------------------------------------------
# V1P-F13 duplicated chr_pos
# ---------------------------------------------------------------------------
test_that("partial and fully pleiotropic selections treat duplicated chr_pos the same way", {
  G <- .v1_panel(m = 30, n = 20, pos = rep(100L, 30))   # every marker same position
  withr::local_dir(withr::local_tempdir())
  eff <- list(rep(.2, 3), rep(.2, 3))
  expect_no_error(qtn_partially_pleiotropic(
    genotypes = G, seed = 3, pleio_a = 1, trait_spec_a_QTN_num = c(2, 2),
    add_effect = eff, epi_interaction = 2, ntraits = 2, rep = 1,
    rep_by = "experiment", export_gt = FALSE, same_add_dom_QTN = FALSE,
    add = TRUE, dom = FALSE, epi = FALSE, verbose = FALSE))
  expect_no_error(qtn_pleiotropic(
    genotypes = G, seed = 3, ntraits = 2, same_add_dom_QTN = FALSE,
    same_mv_QTN = FALSE, add_QTN_num = 3, add_effect = eff,
    epi_interaction = 2, rep = 1, rep_by = "experiment", export_gt = FALSE,
    add = TRUE, dom = FALSE, epi = FALSE, var = FALSE, verbose = FALSE))
})

# ---------------------------------------------------------------------------
# O4 zero trait-specific count for one trait
# ---------------------------------------------------------------------------
test_that("partial: a zero trait-specific count for one trait is allowed when pleio > 0 (O4)", {
  out <- .v1_run(ntraits = 2, pleio_a = 1, trait_spec_a_QTN_num = c(0, 1),
                 add_effect = list(.2, c(.3, .09)), h2 = c(.5, .5),
                 model = "A", architecture = "partially", seed = 2, rep = 1)
  expect_true(is.data.frame(out))
  qf <- .v1_file(out, "Additive_QTNs.txt")
  expect_equal(as.vector(table(qf$trait)), c(1L, 2L))
  # no trait-specific QTN at all: pure shared architecture
  out0 <- .v1_run(ntraits = 2, pleio_a = 2, trait_spec_a_QTN_num = c(0, 0),
                  add_effect = list(c(.2, .04), c(.3, .09)), h2 = c(.5, .5),
                  model = "A", architecture = "partially", seed = 2, rep = 1)
  expect_true(is.data.frame(out0))
})

test_that("partial: a trait with no QTN of a simulated class is rejected with a clear message", {
  r <- .v1_msg(ntraits = 2, pleio_a = 0, trait_spec_a_QTN_num = c(0, 1),
               add_effect = list(.2, .3), h2 = c(.5, .5), model = "A",
               architecture = "partially", seed = 2, rep = 1)
  expect_null(r$result)
  expect_match(r$text, "at least one a QTN per trait|at least one")
})

# ---------------------------------------------------------------------------
# O9 preflight arithmetic
# ---------------------------------------------------------------------------
test_that("partial preflight demands pleio + sum(trait-specific), not sum(pleio + specific) (O9)", {
  G <- .v1_panel(m = 20, n = 20)
  withr::local_dir(withr::local_tempdir())
  eff <- list(rep(.2, 7), rep(.2, 13), rep(.2, 4))
  # 20 eligible markers, need 3 + 4 + 10 + 1 = 18 (the old check demanded 24)
  expect_no_error(qtn_partially_pleiotropic(
    genotypes = G, seed = 3, pleio_a = 3, trait_spec_a_QTN_num = c(4, 10, 1),
    add_effect = eff, epi_interaction = 2, ntraits = 3, rep = 1,
    rep_by = "experiment", export_gt = FALSE, same_add_dom_QTN = FALSE,
    constraints = list(maf_above = 0.01, maf_below = NULL, hets = NULL),
    add = TRUE, dom = FALSE, epi = FALSE, verbose = FALSE))
  # genuinely too few markers -> the informative message
  expect_error(qtn_partially_pleiotropic(
    genotypes = G, seed = 3, pleio_a = 3, trait_spec_a_QTN_num = c(4, 10, 5),
    add_effect = list(rep(.2, 7), rep(.2, 13), rep(.2, 8)), epi_interaction = 2,
    ntraits = 3, rep = 1, rep_by = "experiment", export_gt = FALSE,
    same_add_dom_QTN = FALSE,
    constraints = list(maf_above = 0.01, maf_below = NULL, hets = NULL),
    add = TRUE, dom = FALSE, epi = FALSE, verbose = FALSE),
    "Not enough SNP left")
})

test_that("epistatic preflight counts physical loci (epi_QTN_num * epi_interaction)", {
  G <- .v1_panel(m = 3, n = 20)
  withr::local_dir(withr::local_tempdir())
  expect_error(qtn_pleiotropic(
    genotypes = G, seed = 3, ntraits = 1, same_add_dom_QTN = FALSE,
    same_mv_QTN = FALSE, epi_QTN_num = 2, epi_effect = list(c(.2, .04)),
    epi_interaction = 2, rep = 1, rep_by = "experiment", export_gt = FALSE,
    add = FALSE, dom = FALSE, epi = TRUE, var = FALSE, verbose = FALSE),
    "distinct markers")
})

# ---------------------------------------------------------------------------
# O13 sample_cor replicate-local
# ---------------------------------------------------------------------------
test_that("sample_cor does not leak from one replicate to the next (O13)", {
  set.seed(3)
  n <- 20
  X1 <- matrix(sample(c(-1, 0, 1), n * 3, replace = TRUE), n, 3,
               dimnames = list(letters[1:n], paste0("Chr_1_", 1:3)))
  X2 <- X1 * 0
  res <- base_line_multi_traits(
    add_obj = list(X1, X2), dom_obj = list(NULL, NULL), epi_obj = list(NULL, NULL),
    add_effect = list(c(.2, .3, .1), c(.1, .2, .3)),
    dom_effect = list(NULL, NULL), epi_effect = list(NULL, NULL),
    epi_interaction = 2, ntraits = 2, cor = NULL, architecture = "pleiotropic",
    rep = 2, rep_by = "QTN", add = TRUE, dom = FALSE, epi = FALSE,
    sim_method = "custom", verbose = FALSE)
  expect_false(is.null(res[[1]]$sample_cor))
  expect_null(res[[2]]$sample_cor)
})

# ---------------------------------------------------------------------------
# V1P-F10 vQTL
# ---------------------------------------------------------------------------
.vqtl_call <- function(h2, var_effect = .5, rep = 2, QTN = NULL, base = NULL,
                       var_QTN_num = length(var_effect)) {
  n <- 60
  if (is.null(QTN) || is.null(base)) set.seed(9)
  if (is.null(QTN)) {
    QTN <- matrix(sample(c(-1, 0, 1), n * var_QTN_num, replace = TRUE), n,
                  var_QTN_num, dimnames = list(paste0("i", 1:n), NULL))
  }
  if (is.null(base)) base <- rnorm(n)
  withr::with_tempdir({
    out <- capture.output(res <- suppressMessages(vQTL(
      QTN = QTN, base_line_trait = base, var_QTN_num = var_QTN_num,
      var_effect = var_effect, h2 = h2, rep = rep, seed = 11, mean = 0,
      output_format = "long", fam = NULL, to_r = TRUE)))
  })
  list(out = out, res = res)
}

test_that("vQTL: h2 = 0, negative var_effect, multi-trait input give informative errors", {
  expect_error(.vqtl_call(matrix(0)), "h2. must be in")
  expect_error(.vqtl_call(matrix(.5), var_effect = -.6),
               "negative")
  expect_error(.vqtl_call(matrix(.5), base = cbind(rnorm(60), rnorm(60))),
               "single trait")
  expect_error(.vqtl_call(matrix(.5), var_effect = c(.2, .3), var_QTN_num = 1),
               "one finite value per vQTN")
})

test_that("vQTL: a single replicate runs (was 'dim(X) must have a positive length')", {
  r <- .vqtl_call(matrix(.5), rep = 1)
  expect_equal(ncol(r$res$simulated_data), 2L)
})

test_that("vQTL: sample heritability is reported per row of h2", {
  r <- .vqtl_call(matrix(c(.2, .8)), rep = 3)
  i <- grep("Sample Heritability", r$out)
  vals <- scan(text = sub("^\\[1\\]", "", r$out[i + 2L]), quiet = TRUE)
  # 1 / var(phenotype): 0.2 row must be clearly smaller than 0.8 row
  expect_length(vals, 2L)
  expect_lt(vals[1], vals[2])
  expect_false(isTRUE(all.equal(vals[1], vals[2])))
})

test_that("vQTL: within-individual variance ratio follows sigma = 1 + var_effect * (dosage + 1)", {
  n <- 40
  QTN <- matrix(rep(c(-1, 1), each = n / 2), n, 1,
                dimnames = list(paste0("i", 1:n), NULL))
  r <- .vqtl_call(matrix(.5), var_effect = .5, rep = 3000, QTN = QTN,
                  base = rnorm(n))
  P <- as.matrix(r$res$simulated_data[, -1])
  v_lo <- mean(apply(P[QTN[, 1] == -1, ], 1, var))
  v_hi <- mean(apply(P[QTN[, 1] == 1, ], 1, var))
  expect_gt(v_hi / v_lo, 3.6)     # sigma ratio (1 + .5 * 2) / (1 + 0) = 2 -> variance ratio 4
  expect_lt(v_hi / v_lo, 4.4)
})

# ---------------------------------------------------------------------------
# V1P-F12 genotype loader
# ---------------------------------------------------------------------------
.num_panel <- function(vals) {
  data.frame(snp = paste0("m", seq_len(nrow(vals))), allele = "A/G", chr = 1L,
             pos = seq_len(nrow(vals)) * 10L, cm = 0, vals,
             check.names = FALSE)
}

test_that("genotypes(): 0/1/2 panels are shifted, -1/0/1 kept, 0/1-only rejected", {
  p012 <- .num_panel(cbind(s1 = c(0, 1), s2 = c(1, 2), s3 = c(2, 2)))
  g <- genotypes(geno_obj = p012, verbose = FALSE)$geno_obj
  expect_equal(unname(as.matrix(g[, -(1:5)])), rbind(c(-1, 0, 1), c(0, 1, 1)))
  p101 <- .num_panel(cbind(s1 = c(-1, 0), s2 = c(1, 1), s3 = c(0, -1)))
  g <- genotypes(geno_obj = p101, verbose = FALSE)$geno_obj
  expect_equal(unname(as.matrix(g[, -(1:5)])), unname(as.matrix(p101[, -(1:5)])))
  p01 <- .num_panel(cbind(s1 = c(0, 1), s2 = c(1, 1), s3 = c(0, 0)))
  expect_error(genotypes(geno_obj = p01, verbose = FALSE), "0/1")
})

test_that("genotypes(): out-of-range numeric codes are rejected after normalization", {
  bad <- .num_panel(cbind(s1 = c(-1, 0), s2 = c(0, 7), s3 = c(1, 0)))
  expect_error(genotypes(geno_obj = bad, verbose = FALSE), "outside")
  mix <- .num_panel(cbind(s1 = c(-1, 0), s2 = c(2, 1), s3 = c(1, 0)))
  expect_error(genotypes(geno_obj = mix, verbose = FALSE), "outside")
})

test_that("genotypes(): imputation rule, 'None' leaves NA, all-missing marker and maf_cutoff", {
  p <- .num_panel(cbind(s1 = c(-1, NA), s2 = c(1, NA), s3 = c(NA, NA), s4 = c(0, NA)))
  val <- function(imp) unname(as.matrix(genotypes(geno_obj = p, SNP_impute = imp,
                                                  verbose = FALSE)$geno_obj[, -(1:5)]))
  expect_equal(val("Middle")[2, ], rep(0, 4))
  expect_equal(val("Minor")[2, ], rep(-1, 4))
  expect_equal(val("Major")[2, ], rep(1, 4))
  expect_true(all(is.na(val("None")[2, ])))          # not silently 0
  expect_true(is.na(val("None")[1, 3]))
  # maf_cutoff is applied after imputation (documented): "None" -> undefined MAF
  kept <- genotypes(geno_obj = p, SNP_impute = "None", maf_cutoff = 0.05,
                    verbose = FALSE)$geno_obj
  expect_equal(kept$snp, "m1")
  kept <- genotypes(geno_obj = p, SNP_impute = "Minor", maf_cutoff = 0.05,
                    verbose = FALSE)$geno_obj
  expect_equal(kept$snp, "m1")   # constant -1 -> MAF 0
})

test_that("genotypes(): tied MAF = 0.5 goes to the first-listed allele (documented rule)", {
  hm <- function(alleles, calls) {
    d <- data.frame(`rs#` = paste0("m", seq_along(alleles)), alleles = alleles,
                    chrom = 1, pos = seq_along(alleles) * 100, strand = "+",
                    `assembly#` = NA, center = NA, protLSID = NA,
                    assayLSID = NA, panelLSID = NA, QCcode = NA,
                    check.names = FALSE)
    for (k in seq_along(calls[[1]])) {
      d[[paste0("s", k)]] <- vapply(calls, function(x) x[k], "")
    }
    d
  }
  h <- hm(c("A/G", "G/A"), list(c("AA", "GG", "AA", "GG"),
                                c("GG", "AA", "GG", "AA")))
  g <- genotypes(geno_obj = h, verbose = FALSE)$geno_obj
  expect_equal(unname(as.matrix(g[1, -(1:5)])), rbind(c(1, -1, 1, -1)))
  expect_equal(unname(as.matrix(g[2, -(1:5)])), rbind(c(1, -1, 1, -1)))
})

# ---------------------------------------------------------------------------
# RNG hygiene: seeded helpers leave the caller's RNG stream untouched
# ---------------------------------------------------------------------------
test_that("QTN selection and vQTL restore the caller's RNG state", {
  G <- .v1_panel(m = 40, n = 24, het_markers = c(5, 17, 33))
  set.seed(101)
  s0 <- .Random.seed
  withr::with_tempdir({
    qtn_pleiotropic(genotypes = G, seed = 3, ntraits = 1,
                    same_add_dom_QTN = FALSE, same_mv_QTN = FALSE,
                    add_QTN_num = 3, add_effect = list(rep(.2, 3)),
                    epi_interaction = 2, rep = 1, rep_by = "experiment",
                    export_gt = FALSE, add = TRUE, dom = FALSE, epi = FALSE,
                    var = FALSE, verbose = FALSE)
  })
  expect_identical(.Random.seed, s0)
  withr::with_tempdir({
    qtn_partially_pleiotropic(
      genotypes = G, seed = 3, pleio_a = 1, trait_spec_a_QTN_num = c(2, 2),
      add_effect = list(rep(.2, 3), rep(.2, 3)), epi_interaction = 2,
      ntraits = 2, rep = 1, rep_by = "experiment", export_gt = FALSE,
      same_add_dom_QTN = FALSE, add = TRUE, dom = FALSE, epi = FALSE,
      verbose = FALSE)
  })
  expect_identical(.Random.seed, s0)
  Q <- matrix(sample(c(-1, 0, 1), 60, replace = TRUE), 60, 1,
              dimnames = list(paste0("i", 1:60), NULL))
  b <- rnorm(60)
  set.seed(101)
  s0 <- .Random.seed
  invisible(.vqtl_call(matrix(.5), rep = 2, QTN = Q, base = b))
  expect_identical(.Random.seed, s0)
})
