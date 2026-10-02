# test-v1ld-search.R
#
# Round 9, group v1ld: bounded retry of the direct-LD marker search of the
# frozen v1 engine, create_phenotypes(architecture = "LD", type_of_ld = "direct").
#
# Before: the dominance walks (model "D", "AD" with and without
# same_add_dom_QTN) kept the neighbour pointers of the old marker after a
# re-draw (and the same_add_dom branch reported the original draw instead of
# the re-drawn marker), so only a minority of seeds met the LD contract; the
# rest stopped with the "LD contract" error.
# Now: the first attempt is the frozen search, unchanged; only when it fails
# the replicate is searched again (<= 50 attempts, derived seeds). Every
# accepted pair is verified (same chromosome, distinct, |LD| inside
# [ld_min, ld_max], reported LD == recomputed LD).

.ls_geno <- function() {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  e$SNP55K_maize282_maf04
}

.ls_home <- function(env = parent.frame()) {
  d <- tempfile("ls_")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}

.ls_read <- function(home, file) {
  f <- file.path(home, file)
  if (!file.exists(f)) return(NULL)
  utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
}

# run create_phenotypes() in a fresh home_dir; returns list(res, home, msgs)
.ls_cp <- function(G, ..., seed, home = NULL, verbose = TRUE,
                   ld_min = 0.2, ld_max = 0.8) {
  if (is.null(home)) home <- .ls_home(parent.frame())
  msgs <- character()
  res <- tryCatch(
    withCallingHandlers(
      suppressWarnings(
        create_phenotypes(geno_obj = G, home_dir = home, output_dir = "",
                          to_r = TRUE, verbose = verbose, architecture = "LD",
                          type_of_ld = "direct", ld_method = "corr",
                          ld_min = ld_min, ld_max = ld_max, seed = seed, ...)),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      }),
    error = function(e) conditionMessage(e))
  list(res = res, home = home, msgs = msgs)
}

.ls_ld <- function(G, a, b, method = "corr") {
  g <- function(s) as.numeric(G[match(s, G$snp), -(1:5)]) + 1
  abs(SNPRelate::snpgdsLDpair(g(a), g(b), method = method))[1]
}

# independent verification of an LD summary table against the genotypes
.ls_expect_contract <- function(G, ld, lo = 0.2, hi = 0.8, method = "corr") {
  a <- ld[["QTN_for_trait_1"]]
  b <- ld[["QTN_for_trait_2"]]
  chr <- function(s) G$chr[match(s, G$snp)]
  rep_col <- ld[["rep"]]
  for (r in unique(rep_col)) {
    k <- rep_col == r
    expect_false(any(a[k] == b[k]))
    expect_false(anyDuplicated(a[k]) > 0)
    expect_false(anyDuplicated(b[k]) > 0)
    expect_length(intersect(a[k], b[k]), 0)
  }
  expect_equal(chr(a), chr(b))
  act <- vapply(seq_along(a), function(i) .ls_ld(G, a[i], b[i], method),
                numeric(1))
  expect_true(all(act >= lo - 1e-9 & act <= hi + 1e-9))
  expect_equal(ld[[grep("^Actual_LD", names(ld))[1]]], act, tolerance = 1e-6)
  invisible(act)
}

cfg_d   <- list(dom_QTN_num = 3, dom_effect = c(0.5, 0.3), h2 = c(0.5, 0.5),
                model = "D")
cfg_ads <- list(add_QTN_num = 3, add_effect = c(0.5, 0.3), h2 = c(0.5, 0.5),
                model = "AD", same_add_dom_QTN = TRUE)
cfg_adp <- list(add_QTN_num = 3, dom_QTN_num = 3, add_effect = c(0.5, 0.3),
                dom_effect = c(0.5, 0.3), h2 = c(0.5, 0.5), model = "AD",
                same_add_dom_QTN = FALSE)
cfg_a   <- list(add_QTN_num = 3, add_effect = c(0.5, 0.3), h2 = c(0.5, 0.5),
                model = "A")

.ls_run <- function(cfg, seed, ..., rep = 1) {
  do.call(.ls_cp, c(list(G = .ls_geno(), seed = seed, rep = rep), cfg, list(...)),
          envir = parent.frame())
}

# ---------------------------------------------------------------------------
# 1. bit-identity: calls that succeeded before keep their exact output
# ---------------------------------------------------------------------------
test_that("v1ld: seeds that met the LD contract on the first attempt are unchanged", {
  skip_if_not_installed("SNPRelate")
  # QTN pairs and LD recorded from the pre-change engine (HEAD 7db73fc)
  ref <- list(
    list(cfg = cfg_d, seed = 12, file = "LD_Summary_Dominance.txt",
         t1 = c("ss196520033", "ss196472773", "ss196508280"),
         t2 = c("ss196441177", "ss196475012", "ss196497037"),
         ld = c(0.202648, 0.20921, 0.795529)),
    list(cfg = cfg_d, seed = 29, file = "LD_Summary_Dominance.txt",
         t1 = c("ss196503498", "ss196496839", "ss196454358"),
         t2 = c("ss196503502", "ss196496813", "ss196454360"),
         ld = c(0.725133, 0.601196, 0.743601)),
    list(cfg = cfg_ads, seed = 30, file = "LD_Summary.txt",
         t1 = c("ss196503498", "ss196496839", "ss196454358"),
         t2 = c("ss196503502", "ss196496813", "ss196454360"),
         ld = c(0.725133, 0.601196, 0.743601)),
    list(cfg = cfg_ads, seed = 40, file = "LD_Summary.txt",
         t1 = c("ss196438301", "ss196510530", "ss196467481"),
         t2 = c("ss196438305", "ss196498144", "ss196467483"),
         ld = c(0.757532, 0.719199, 0.586575)),
    list(cfg = cfg_adp, seed = 12, file = "LD_Summary_Additive.txt",
         t1 = c("ss196520033", "ss196468948", "ss196451156"),
         t2 = c("ss196441177", "ss196468967", "ss196451154"),
         ld = c(0.202648, 0.417754, 0.778482)),
    list(cfg = cfg_a, seed = 7, file = "LD_Summary_Additive.txt",
         t1 = c("ss196509786", "ss196437647", "ss196523370"),
         t2 = c("ss196509794", "ss196437621", "ss196442460"),
         ld = c(0.339382, 0.438891, 0.302941))
  )
  for (r in ref) {
    o <- .ls_run(r$cfg, r$seed)
    expect_false(is.character(o$res), info = paste(r$seed, o$res))
    ld <- .ls_read(o$home, r$file)
    expect_equal(ld[["QTN_for_trait_1"]], r$t1, info = r$file)
    expect_equal(ld[["QTN_for_trait_2"]], r$t2, info = r$file)
    expect_equal(ld[[grep("^Actual_LD", names(ld))[1]]], r$ld, tolerance = 1e-5)
    # no retry happened
    expect_false(any(grepl("contract met on attempt", o$msgs)))
  }
})

# ---------------------------------------------------------------------------
# 2. seeds that used to stop with the contract error now succeed, verified
# ---------------------------------------------------------------------------
test_that("v1ld: former 'LD contract' failures of the dominance branches now return valid pairs", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  G <- .ls_geno()
  # at the pre-change engine every one of these seeds stopped with the
  # "LD contract" error (pair spanning two chromosomes / outside the window /
  # reported LD differing from the actual one)
  cases <- list(
    list(cfg = cfg_ads, seeds = c(1, 2, 3), file = "LD_Summary.txt"),
    list(cfg = cfg_adp, seeds = c(1, 2),    file = "LD_Summary_Additive.txt"),
    list(cfg = cfg_adp, seeds = c(1, 2),    file = "LD_Summary_Dominance.txt")
  )
  n_retry <- 0
  for (cs in cases) {
    for (sd in cs$seeds) {
      o <- .ls_run(cs$cfg, sd)
      expect_false(is.character(o$res), info = paste(sd, o$res))
      if (is.character(o$res)) next
      ld <- .ls_read(o$home, cs$file)
      .ls_expect_contract(G, ld)
      n_retry <- n_retry + any(grepl("contract met on attempt", o$msgs))
    }
  }
  expect_gt(n_retry, 0)
  # dominance only: contract met; the only remaining stop (if any) is the
  # separate "homozygote" guard on the selected dominance markers
  for (sd in 2:5) {
    o <- .ls_run(cfg_d, sd)
    if (is.character(o$res)) {
      expect_match(o$res, "All individuals are homozygote")
    } else {
      .ls_expect_contract(G, .ls_read(o$home, "LD_Summary_Dominance.txt"))
    }
  }
})

test_that("v1ld: with several QTN replicates every replicate is searched and verified on its own", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  G <- .ls_geno()
  o <- .ls_run(cfg_ads, 5, rep = 3, vary_QTN = TRUE)
  expect_false(is.character(o$res), info = o$res)
  ld <- .ls_read(o$home, "LD_Summary.txt")
  expect_equal(sort(unique(ld[["rep"]])), 1:3)
  .ls_expect_contract(G, ld)
})

test_that("v1ld: additive-only direct LD with several QTN replicates closes one GDS handle per replicate", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  G <- .ls_geno()
  # before: the second replicate failed with "The file ... has been created or
  # opened" because only the last GDS handle was closed
  o <- .ls_run(cfg_a, 3, rep = 3, vary_QTN = TRUE)
  expect_false(is.character(o$res), info = o$res)
  ld <- .ls_read(o$home, "LD_Summary_Additive.txt")
  expect_equal(sort(unique(ld[["rep"]])), 1:3)
  .ls_expect_contract(G, ld)
})

test_that("v1ld: the window is inclusive and never relaxed (narrow window still verified)", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  G <- .ls_geno()
  o <- .ls_run(cfg_ads, 4, ld_min = 0.6, ld_max = 0.7)
  if (is.character(o$res)) {
    expect_match(o$res, "LD contract|None of the selected SNPs")
  } else {
    .ls_expect_contract(G, .ls_read(o$home, "LD_Summary.txt"), lo = 0.6, hi = 0.7)
  }
})

test_that("v1ld: a successful retry names its attempt and derived seed", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  o <- .ls_run(cfg_ads, 1)
  expect_false(is.character(o$res), info = o$res)
  m <- grep("contract met on attempt", o$msgs, value = TRUE)
  expect_length(m, 1)
  att <- as.integer(sub(".*attempt ([0-9]+) .*", "\\1", m))
  expect_gt(att, 1L)
  derived <- as.numeric(sub(".*derived seed (-?[0-9]+)\\).*", "\\1", m))
  expect_equal(derived, simplePHENOTYPES:::.ld_attempt_seed(1, att))
})

# ---------------------------------------------------------------------------
# 3. the retry rule (unit level)
# ---------------------------------------------------------------------------
test_that("v1ld: .ld_attempt_seed keeps attempt 1 and walks towards zero", {
  f <- simplePHENOTYPES:::.ld_attempt_seed
  stride <- simplePHENOTYPES:::.ld_retry_stride
  maxa <- simplePHENOTYPES:::.ld_max_attempts()
  expect_identical(maxa, 50L)
  expect_identical(f(123, 1L), 123)
  expect_identical(f(-7, 1L), -7)
  expect_equal(f(123, 2L), 123 - stride)
  expect_equal(f(-7, 2L), -7 + stride)
  expect_equal(f(0, 3L), 2 * stride)
  for (sd in c(1, 200, 2e6, 214748364, -214748364, 0)) {
    a <- vapply(1:maxa, function(k) f(sd, k), numeric(1))
    expect_true(all(abs(a) <= max(abs(sd), (maxa - 1) * stride)))
    expect_equal(anyDuplicated(a), 0)   # distinct attempts get distinct seeds
  }
})

test_that("v1ld: seed validator stays exact for LD and covers the retry seeds", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  f <- simplePHENOTYPES:::.ld_attempt_seed
  maxa <- simplePHENOTYPES:::.ld_max_attempts()
  h2 <- matrix(0.5)
  M <- .Machine$integer.max
  # the accepted interval of an LD call is unchanged: floor((M - extra) / 10),
  # extra = 2 * rep + n_qtn
  expect_silent(v(214748364, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3))
  expect_error(v(214748365, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3), "too large")
  expect_silent(v(-214748364, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3))
  # every derived retry seed of an accepted call stays inside the integer
  # range after the multiplier 10 and the replicate/QTN offsets
  for (sd in c(1, 5000, 214748364, -214748364)) {
    a <- vapply(1:maxa, function(k) f(sd, k), numeric(1))
    expect_lte(10 * max(abs(a)) + 2 * 1 + 3 + 1, M)
  }
  # small seeds (below the retry span) are accepted: the span itself fits
  expect_silent(v(1, 1, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3))
  # a replicate count that leaves no room for the retry span is refused
  expect_error(v(1, 8e8, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3), "too large|No seed")
  # non-LD calls are not affected
  expect_silent(v(214748364, 0, matrix(1), FALSE, FALSE))
})

test_that("v1ld: a seed at the validator bound runs without integer overflow", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  o <- .ls_run(cfg_ads, 214748300)
  # the call either succeeds or stops with an informative LD error; derived
  # seeds never overflow (that would end in "NAs produced by integer overflow")
  if (is.character(o$res)) {
    expect_match(o$res, "LD contract|None of the selected SNPs")
  } else {
    expect_s3_class(o$res, "data.frame")
  }
})

test_that("v1ld: .ld_run_attempt returns search failures, holds retry warnings, passes other errors", {
  run_attempt <- simplePHENOTYPES:::.ld_run_attempt
  stop_ <- simplePHENOTYPES:::.ld_search_stop
  # a frozen "None of the selected SNPs" stop is returned, not raised
  r <- run_attempt(stop_("None of the selected SNPs met the minimum LD threshold."))
  expect_s3_class(r$failure, "ld_search_failed")
  expect_match(conditionMessage(r$failure), "None of the selected SNPs")
  # a read past the last marker is a search failure too
  r <- run_attempt(stop("'start' is invalid"))
  expect_match(conditionMessage(r$failure), "'start' is invalid")
  # the frozen walk's exhaustion of every candidate marker (sample.int() on an
  # empty set) is a search failure with an informative message
  r <- run_attempt(sample(setdiff(1:5, 1:5), 1))
  expect_s3_class(r$failure, "ld_search_exhausted")
  expect_match(conditionMessage(r$failure), "used up every candidate marker")
  # an unrelated "invalid first argument" error is not mistaken for it
  expect_error(run_attempt(stop("invalid first argument")), "invalid first argument")
  # other errors (e.g. "Monomorphic SNPs") propagate unchanged
  expect_error(run_attempt(stop("Monomorphic SNPs are not accepted", call. = FALSE)),
               "Monomorphic")
  # assignments of the evaluated expression are visible to the caller
  run_attempt({ visible <- 42 })
  expect_identical(visible, 42)
  # attempt 1 lets warnings through, retries hold them back
  expect_warning(run_attempt(warning("w1"), quiet = FALSE), "w1")
  expect_silent(h <- run_attempt(warning("w2"), quiet = TRUE))
  expect_length(h$warnings, 1)
  expect_warning(simplePHENOTYPES:::.ld_emit_warnings(h$warnings), "w2")
})

test_that("v1ld: the retry loop stops after the bound, or at once when a retry exhausts the markers", {
  stop_now <- simplePHENOTYPES:::.ld_give_up_now
  ex <- structure(class = c("ld_search_exhausted", "ld_search_failed", "error", "condition"),
                  list(message = "x", call = NULL))
  other <- "a pair spans two chromosomes"
  expect_false(stop_now(other, 1L))
  expect_false(stop_now(other, 49L))
  expect_true(stop_now(other, 50L))
  expect_false(stop_now(ex, 1L))   # the frozen walk may exhaust because of its pointer defect
  expect_true(stop_now(ex, 2L))    # the repaired walk exhausting the markers is conclusive
})

test_that("v1ld: giving up re-raises the search error or the contract error with the attempt count", {
  giveup <- simplePHENOTYPES:::.ld_search_giveup
  f <- simplePHENOTYPES:::.ld_search_stop
  e <- tryCatch(f("None of the selected SNPs met the maximum LD threshold. Try another seed number."),
                error = function(e) e)
  expect_error(giveup(e, 1, 50L), "None of the selected SNPs met the maximum")
  expect_error(giveup("a pair spans two chromosomes", 2, 50L),
               "LD contract.*direct LD, replicate 2.*a pair spans two chromosomes.*49 more time")
})

test_that("v1ld: the violation helper keeps the frozen checks (same messages through .ld_check_direct)", {
  skip_if_not_installed("SNPRelate")
  G <- .ls_geno()
  viol <- simplePHENOTYPES:::.ld_direct_violation
  chk <- simplePHENOTYPES:::.ld_check_direct
  expect_match(viol(G, 5, 5, 0.5, "corr", 0.2, 0.8), "paired with itself")
  expect_match(viol(G, c(5, 5), c(6, 7), c(0.5, 0.5), "corr", 0.2, 0.8), "duplicated")
  expect_error(chk(G, 5, 5, 0.5, "corr", 0.2, 0.8, 1), "LD contract.*paired with itself")
  i <- which(G$chr[-1] != G$chr[-nrow(G)])[1]   # last marker of the first chromosome
  expect_match(viol(G, i, i + 1, 0.5, "corr", 0.2, 0.8), "two chromosomes")
  # a pair that meets the contract returns NULL / TRUE
  ld <- .ls_ld(G, G$snp[2], G$snp[3])
  expect_null(viol(G, 2, 3, ld, "corr", 0, 1))
  expect_true(chk(G, 2, 3, ld, "corr", 0, 1, 1))
  # a reported LD that differs from the recomputed one is rejected
  expect_match(viol(G, 2, 3, ld + 0.1, "corr", 0, 1), "differs")
})

# ---------------------------------------------------------------------------
# 4. a search that never meets the window ends in the informative error
# ---------------------------------------------------------------------------
test_that("v1ld: an unattainable window gives up with the contract error after the bounded retries", {
  skip_if_not_installed("SNPRelate")
  skip_on_cran()
  # two attempts only, to keep the test short; the window [0.9999, 0.99995]
  # holds no pair of distinct markers
  testthat::local_mocked_bindings(.ld_max_attempts = function() 2L,
                                  .package = "simplePHENOTYPES")
  o <- .ls_run(cfg_ads, 1, ld_min = 0.9999, ld_max = 0.99995)
  expect_true(is.character(o$res))
  expect_match(o$res, "LD contract.*direct LD, replicate 1.*repeated 1 more time")
})

test_that("v1ld: a window with no pair in a small panel ends in an informative error, not 'invalid first argument'", {
  skip_if_not_installed("SNPRelate")
  G <- .ls_geno()[1:60, ]
  for (cfg in list(cfg_d, cfg_a)) {
    o <- do.call(.ls_cp, c(list(G = G, seed = 1, rep = 1, ld_min = 0.9999,
                                ld_max = 0.99995), cfg))
    expect_true(is.character(o$res))
    expect_false(grepl("invalid first argument", o$res, fixed = TRUE))
    expect_match(o$res, "used up every candidate marker|LD contract|ran past")
  }
})
