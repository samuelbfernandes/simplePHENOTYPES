# test-adopt-v1-core.R
#
# Round-9 adoption of the independent audit's remaining test proposals for the
# frozen v1 engine, group "v1 core + linkage" (reconciliation/v1-core-linkage.md
# section 6, Fable T1-T24 and Codex rows). Proposals that were already adopted
# in test-audit-v1-core.R / test-fix*-v1.R / test-v130-parity.R are not repeated
# (see the coverage table in the round-9 handoff). Nothing here edits the
# frozen RDS references or the create_phenotypes() signature; the references
# are only READ.
#
# Deterministic seeds, tempdir() writes only.

# ---------------------------------------------------------------------------
# helpers (own prefix: test files share a package namespace, not a scope)
# ---------------------------------------------------------------------------
.ad1_geno <- function() {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  e$SNP55K_maize282_maf04
}

# every run gets its own sub-folder of one per-file root under tempdir()
.ad1_root <- tempfile("ad1root_")
dir.create(.ad1_root)
.ad1_n <- 0L
.ad1_home <- function() {
  .ad1_n <<- .ad1_n + 1L
  d <- file.path(.ad1_root, sprintf("run%04d", .ad1_n))
  dir.create(d)
  d
}

# roxygen text of an R source file with the "#'" prefixes and line breaks removed
.ad1_roxy <- function(file) {
  l <- readLines(testthat::test_path("..", "..", "R", file), warn = FALSE)
  gsub("[[:space:]]+", " ", paste(sub("^#' ?", "", l), collapse = " "))
}

# create_phenotypes() into a fresh home_dir; returns the data (or the error
# message, as a character) together with the folder
.ad1_run <- function(G, ..., home = .ad1_home()) {
  r <- tryCatch(
    suppressMessages(suppressWarnings(
      create_phenotypes(geno_obj = G, home_dir = home, output_dir = "",
                        to_r = TRUE, verbose = FALSE, ...))),
    error = function(e) conditionMessage(e))
  list(r = r, home = home)
}

.ad1_read <- function(home, file) {
  f <- file.path(home, file)
  if (!file.exists(f)) return(NULL)
  as.data.frame(data.table::fread(f))
}

.ad1_ld <- function(G, s1, s2) {
  g <- function(s) as.numeric(G[match(s, G$snp), -(1:5)]) + 1
  abs(SNPRelate::snpgdsLDpair(g(s1), g(s2), method = "corr"))[1]
}

.ad1_maf <- function(G) {
  ns <- ncol(G) - 5
  m <- apply(G[, -(1:5)], 1, function(x) {
    p <- ((sum(x) + ns) / ns * 0.5)
    min(p, 1 - p)
  })
  stats::setNames(m, G$snp)
}

.ad1_ref_dir <- function() {
  installed <- system.file("extdata", "v1_3_0_reference", package = "simplePHENOTYPES")
  if (nzchar(installed) && dir.exists(installed)) installed
  else testthat::test_path("..", "..", "inst", "extdata", "v1_3_0_reference")
}

# ---------------------------------------------------------------------------
# T24 / T19 / Codex row 11: exact frozen-reference parity for the partial
# architecture (test-v130-parity.R compares only the additive count to +-3)
# ---------------------------------------------------------------------------
test_that("T24/Codex-11: partial pleiotropy reproduces the frozen QTN identities and every trait column", {
  f <- file.path(.ad1_ref_dir(), "partial_pleiotropy.rds")
  skip_if_not(file.exists(f), "frozen partial_pleiotropy reference not available")
  ref <- readRDS(f)
  G <- .ad1_geno()
  o <- do.call(.ad1_run, c(list(G), ref$call, list(output_format = "long")))
  expect_s3_class(o$r, "data.frame")

  add <- .ad1_read(o$home, "Additive_QTNs.txt")
  epi <- .ad1_read(o$home, "Epistatic_QTNs.txt")
  # exact marker identity AND order, trait labels and types (rep 1)
  expect_identical(add$snp, ref$qtns$add$snp)
  expect_identical(add$trait, ref$qtns$add$trait)
  expect_identical(add$type, ref$qtns$add$type)
  expect_identical(epi$snp, ref$qtns$epi$snp)
  expect_identical(epi$trait, ref$qtns$epi$trait)

  r_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  g_ord <- o$r[order(o$r[[1L]]), ]
  expect_identical(names(g_ord), names(r_ord))
  for (cn in setdiff(names(r_ord), c("<Trait>", "Rep"))) {
    expect_equal(g_ord[[cn]], r_ord[[cn]], tolerance = 1e-10, label = cn)
  }
})

test_that("T24: the single-trait reference is reproduced for every column and the QTN effect", {
  f <- file.path(.ad1_ref_dir(), "single_trait.rds")
  skip_if_not(file.exists(f), "frozen single_trait reference not available")
  ref <- readRDS(f)
  G <- .ad1_geno()
  o <- do.call(.ad1_run, c(list(G), ref$call, list(output_format = "long")))
  add <- .ad1_read(o$home, "Additive_QTNs.txt")
  expect_identical(add$snp[add$rep == 1], ref$qtns$add$snp[ref$qtns$add$rep == 1])
  r_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  g_ord <- o$r[order(o$r[[1L]]), ]
  expect_identical(g_ord[[1L]], r_ord[[1L]])          # same individuals, same order
  expect_equal(ncol(g_ord), ncol(r_ord))
  expect_equal(g_ord[[2L]], r_ord[[2L]], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# T3 / Codex row 4 (frozen quirk): residual seeds collide for heritabilities
# that round to the same multiple of 0.1. The audit proposed "distinct"; the
# owner decision kept the frozen seed formula and DOCUMENTED the consequence
# (roxygen of create_phenotypes(), `seed`), so the test locks both the
# behaviour and its documentation.
# ---------------------------------------------------------------------------
test_that("T3: h2 values rounding to the same tenth share a residual stream; others do not (documented)", {
  G <- .ad1_geno()
  run <- function(h2) {
    o <- .ad1_run(G, add_QTN_num = 3, add_effect = 0.2, h2 = h2, rep = 1,
                  model = "A", seed = 7)
    gv <- .ad1_read(o$home, "Genetic_values.txt")[[1]]
    lapply(o$r, function(d) d[[2]] - gv)
  }
  same <- run(c(0.26, 0.34))           # round(10 * h2) is 3 for both
  expect_equal(stats::cor(same[[1]], same[[2]]), 1, tolerance = 1e-10)
  # ... the residual scale still follows h2 (only the standardized draw is shared)
  expect_gt(stats::sd(same[[1]]) / stats::sd(same[[2]]), 1)
  diff <- run(c(0.26, 0.45))           # 3 vs 4 (R rounds 4.5 to 4)
  expect_lt(abs(stats::cor(diff[[1]], diff[[2]])), 0.5)

  skip_if_no_source("R", "legacy_create_phenotypes.R")
  src <- .ad1_roxy("legacy_create_phenotypes.R")
  expect_match(src, "0.26 and 0.34", fixed = TRUE)
  expect_match(src, "perfectly correlated", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# T23: lock the direct additive LD branch that was already correct
# ---------------------------------------------------------------------------
test_that("T23: direct additive LD pairs are distinct, same-chromosome, inside the window and reported correctly (6 seeds)", {
  G <- .ad1_geno()
  for (s in c(1, 2, 3, 5, 7, 200)) {
    o <- .ad1_run(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
                  rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
                  ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = s)
    expect_s3_class(o$r, "data.frame")
    ld <- .ad1_read(o$home, "LD_Summary_Additive.txt")
    a <- ld[["QTN_for_trait_1"]]
    b <- ld[["QTN_for_trait_2"]]
    expect_length(a, 3)
    expect_false(anyNA(match(c(a, b), G$snp)))
    expect_false(any(a == b))
    expect_equal(G$chr[match(a, G$snp)], G$chr[match(b, G$snp)])
    act <- vapply(seq_along(a), function(k) .ad1_ld(G, a[k], b[k]), numeric(1))
    expect_true(all(act >= 0.2 - 1e-9 & act <= 0.8 + 1e-9), info = paste("seed", s))
    expect_equal(ld[["Actual_LD"]], act, tolerance = 1e-6)
    # the QTN table agrees with the summary (same markers per trait)
    q <- .ad1_read(o$home, "Additive_QTNs.txt")
    expect_setequal(c(a, b), q$snp)
  }
})

# ---------------------------------------------------------------------------
# T6 / T24: the indirect-LD seeds that aborted with cryptic errors
# ---------------------------------------------------------------------------
test_that("T6/T24: indirect LD seeds 3-7 and 200 end in a valid result or the LD-contract error, never a cryptic abort; seed 200 (the parity seed) runs", {
  G <- .ad1_geno()
  go <- function(s) {
    .ad1_run(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
             rep = 1, model = "A", architecture = "LD", type_of_ld = "indirect",
             ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = s)
  }
  n_ok <- 0
  for (s in c(3:7, 200)) {
    o <- go(s)
    if (is.character(o$r)) {
      expect_match(o$r, "LD contract", info = paste("seed", s))
      expect_false(grepl("duplicate 'row.names'|undefined columns|subscript out of bounds",
                         o$r), info = paste("seed", s))
      next
    }
    n_ok <- n_ok + 1
    expect_s3_class(o$r, "data.frame")
    ld <- .ad1_read(o$home, "LD_Summary_Additive.txt")
    cause <- ld[["SNP_causing_LD"]]
    t1 <- ld[["QTN_for_trait_1"]]
    t2 <- ld[["QTN_for_trait_2"]]
    # three distinct markers per row; no marker is causal for both traits
    expect_false(any(cause == t1 | cause == t2 | t1 == t2), info = paste("seed", s))
    expect_length(intersect(t1, t2), 0)
    expect_false(anyDuplicated(t1) > 0)
    expect_false(anyDuplicated(t2) > 0)
    for (k in seq_along(cause)) {
      for (tt in list(t1, t2)) {
        ldk <- .ad1_ld(G, cause[k], tt[k])
        expect_gte(ldk, 0.2 - 1e-9)
        expect_lte(ldk, 0.8 + 1e-9)
      }
    }
  }
  expect_gt(n_ok, 0)                        # seed 200 is a frozen-reference seed
  expect_s3_class(go(200)$r, "data.frame")
})

# ---------------------------------------------------------------------------
# T8 (AD half): indirect LD with dominance QTNs
# ---------------------------------------------------------------------------
test_that("T8: indirect LD with additive + dominance QTNs either runs inside the window or stops with a controlled message", {
  G <- .ad1_geno()
  n_err <- 0
  for (s in 1:6) {
    for (cons in list(NULL, list(maf_above = NULL, maf_below = NULL, hets = "include"))) {
      args <- list(add_QTN_num = 3, add_effect = c(0.5, 0.3), dom_QTN_num = 3,
                   dom_effect = c(0.5, 0.3), h2 = c(0.5, 0.5), rep = 1,
                   model = "AD", architecture = "LD", type_of_ld = "indirect",
                   ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = s)
      if (!is.null(cons)) args$constraints <- cons
      o <- do.call(.ad1_run, c(list(G), args))
      if (is.character(o$r)) {
        n_err <- n_err + 1
        expect_match(o$r, "LD contract|re-sample an intermediate marker")
        # an up-front stop must not leave a half-written output folder behind
        expect_false(file.exists(file.path(o$home, "Dominance_QTNs.txt")))
      } else {
        d <- .ad1_read(o$home, "Dominance_QTNs.txt")
        expect_s3_class(o$r, "data.frame")
        expect_true(all(d$snp %in% G$snp))
      }
    }
  }
  # on the heterozygote-poor maize panel most runs are (informatively) rejected
  expect_gt(n_err, 0)
})

# ---------------------------------------------------------------------------
# Codex row 6: constraint eligibility. The documented contract (roxygen of
# `constraints`, DOC-ONLY decision for Codex F14) is that constraints filter the
# randomly drawn anchors, NOT the partner markers; lock exactly that.
# ---------------------------------------------------------------------------
test_that("Codex-6: MAF and heterozygote constraints hold for every sampled direct-LD anchor (partners are documented as unfiltered)", {
  G <- .ad1_geno()
  maf <- .ad1_maf(G)
  hets <- stats::setNames(apply(G[, -(1:5)], 1, function(x) sum(x == 0)), G$snp)
  run <- function(s, cons) {
    .ad1_run(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
             rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
             ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = s,
             constraints = cons)
  }
  for (s in c(1:4, 200)) {
    o <- run(s, list(maf_above = 0.47, maf_below = NULL, hets = NULL))
    expect_s3_class(o$r, "data.frame")
    q <- .ad1_read(o$home, "Additive_QTNs.txt")
    anchors <- q$snp[q$type == "QTN_selected"]
    expect_length(anchors, 3)
    expect_true(all(maf[anchors] > 0.47), info = paste("seed", s))
  }
  for (s in 1:4) {
    o <- run(s, list(maf_above = NULL, maf_below = NULL, hets = "remove"))
    q <- .ad1_read(o$home, "Additive_QTNs.txt")
    anchors <- q$snp[q$type == "QTN_selected"]
    expect_true(all(hets[anchors] == 0), info = paste("seed", s))
  }
  # an unknown option value is a clear error, not a silent no-op
  expect_match(run(1, list(maf_above = NULL, maf_below = NULL, hets = "exclude"))$r,
               "hets option")

  skip_if_no_source("R", "legacy_create_phenotypes.R")
  src <- .ad1_roxy("legacy_create_phenotypes.R")
  expect_match(src, "filter only the randomly drawn QTNs, not the partner markers",
               fixed = TRUE)
})

# ---------------------------------------------------------------------------
# Codex row 5: first / last row of the marker set
# ---------------------------------------------------------------------------
test_that("Codex-5: direct LD at the first and last markers ends in a valid pair or a precise no-partner error, never the GDS 'start is invalid' error", {
  G <- .ad1_geno()
  go <- function(GG, s) {
    .ad1_run(GG, add_QTN_num = 1, add_effect = c(0.2, 0.3), h2 = c(0.5, 0.5),
             rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
             ld_min = 0.05, ld_max = 0.8, ld_method = "corr", seed = s)
  }
  n_err <- 0
  n_ok <- 0
  for (s in 1:3) {
    for (GG in list(G[1:2, ], G[(nrow(G) - 2):nrow(G), ],
                    G[1:6, ], G[(nrow(G) - 5):nrow(G), ])) {
      o <- go(GG, s)
      if (is.character(o$r)) {
        n_err <- n_err + 1
        expect_false(grepl("'start' is invalid", o$r, fixed = TRUE))
        expect_match(o$r, "ran past the first or last marker|LD contract")
      } else {
        n_ok <- n_ok + 1
        q <- .ad1_read(o$home, "Additive_QTNs.txt")
        expect_length(unique(q$snp), 2)
        expect_true(all(q$snp %in% GG$snp))
        expect_gte(.ad1_ld(GG, q$snp[1], q$snp[2]), 0.05 - 1e-9)
        expect_lte(.ad1_ld(GG, q$snp[1], q$snp[2]), 0.8 + 1e-9)
      }
    }
  }
  expect_gt(n_err, 0)     # a 2-marker panel cannot host a pair away from the border
  expect_gt(n_ok, 0)      # a 6-marker panel can, from either end
})

# ---------------------------------------------------------------------------
# Codex row 13: missing marker ids stop before anything is written
# ---------------------------------------------------------------------------
test_that("Codex-13: a QTN_list with an unknown marker names it and leaves no output behind", {
  G <- .ad1_geno()
  ql <- list(add = list(c("ss196521092", "not_a_marker"), "ss196469503"))
  o <- .ad1_run(G, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
                h2 = c(0.5, 0.5), rep = 1, model = "A", seed = 1)
  expect_type(o$r, "character")
  expect_match(o$r, "not_a_marker")
  expect_false(grepl("Chr_NA", o$r, fixed = TRUE))
  expect_length(list.files(o$home, recursive = TRUE, all.files = TRUE), 0)
})

# ---------------------------------------------------------------------------
# Codex rows 15-16 / T16: user epistatic list with more groups than the width
# ---------------------------------------------------------------------------
test_that("Codex-15/T16: six user epistatic markers, arity 2, three effects give QTN ids 1,1,2,2,3,3 and effects per group", {
  G <- .ad1_geno()
  m1 <- G$snp[c(10, 50, 90, 130, 170, 210)]
  m2 <- G$snp[c(11, 51, 91, 131, 171, 211)]
  o <- .ad1_run(G, QTN_list = list(epi = list(m1, m2)), ntraits = 2,
                epi_effect = list(c(0.1, 0.2, 0.3), c(0.3, 0.2, 0.1)),
                epi_interaction = 2, h2 = c(0.5, 0.5), rep = 1, model = "E", seed = 1)
  expect_s3_class(o$r, "data.frame")
  e <- .ad1_read(o$home, "Epistatic_QTNs.txt")
  t1 <- e[e$trait == "trait_1", ]
  t2 <- e[e$trait == "trait_2", ]
  expect_equal(nrow(t1), 6)
  expect_equal(t1$QTN, c(1, 1, 2, 2, 3, 3))
  expect_equal(t1$epistatic_effect, c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3))
  expect_equal(t2$epistatic_effect, c(0.3, 0.3, 0.2, 0.2, 0.1, 0.1))
  expect_identical(t1$snp, m1)
  expect_identical(t2$snp, m2)
})
