# test-audit-v1-core.R
#
# Regression tests for the v1 (frozen legacy) engine audit, group "v1-core"
# (reconciliation files v1-core-linkage.md and v1-pleio-effects.md).
#
# Owner decision D1: bad inputs are rejected up front with a message that names
# the argument and the remedy; every currently-correct output stays
# bit-identical (test-v130-parity.R guards that). These tests therefore mostly
# assert *informative errors*, plus the few places where diagnostics were wrong
# (labels, per-replicate PVE) and are now fixed.
#
# create_phenotypes() re-signals its inner errors (it used to swallow them and
# return NULL), so `expect_error()` is a meaningful assertion here.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
.a1_geno <- function() {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  e$SNP55K_maize282_maf04
}

.a1_home <- function(env = parent.frame()) {
  d <- tempfile("a1_")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}

# run create_phenotypes() into a fresh home_dir (no run sub-folder)
.a1_cp <- function(G, ..., home = NULL, to_r = TRUE, verbose = TRUE) {
  if (is.null(home)) home <- .a1_home(parent.frame())
  suppressMessages(suppressWarnings(
    create_phenotypes(geno_obj = G, home_dir = home, output_dir = "",
                      to_r = to_r, verbose = verbose, ...)
  ))
}

.a1_read <- function(home, file) {
  f <- file.path(home, file)
  if (!file.exists(f)) return(NULL)
  utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
}

.a1_ld <- function(G, s1, s2) {
  g <- function(s) as.numeric(G[match(s, G$snp), -(1:5)]) + 1
  abs(SNPRelate::snpgdsLDpair(g(s1), g(s2), method = "corr"))[1]
}

# ---------------------------------------------------------------------------
# V1C-F1: single-trait QTN_list
# ---------------------------------------------------------------------------
test_that("V1C-F1: single-trait QTN_list is an explicit, informative error", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, QTN_list = list(add = list("ss196521092")), add_effect = 0.2,
           h2 = 0.5, rep = 1, model = "A", seed = 1),
    "QTN_list.*ntraits = 1|ntraits = 1.*QTN_list"
  )
})

test_that("multi-trait QTN_list still runs (positive control)", {
  G <- .a1_geno()
  ql <- list(add = list(c("ss196521092", "ss196469503"),
                        c("ss196521092", "ss196520297")))
  ph <- .a1_cp(G, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
               h2 = c(0.5, 0.5), rep = 1, model = "A", seed = 1)
  expect_s3_class(ph, "data.frame")
  expect_true(all(is.finite(as.matrix(ph[, 2:3]))))
})

# ---------------------------------------------------------------------------
# V1P-O1: single-trait DE / ADE ; V1P-F5: single-trait same_add_dom_QTN
# ---------------------------------------------------------------------------
test_that("V1P-O1: single-trait models with dominance and epistasis are rejected", {
  G <- .a1_geno()
  for (m in c("ADE", "DE")) {
    expect_error(
      .a1_cp(G, add_QTN_num = 2, dom_QTN_num = 2, epi_QTN_num = 1,
             add_effect = 0.3, dom_effect = 0.3, epi_effect = 0.3,
             epi_interaction = 2, h2 = 0.5, rep = 1, model = m, seed = 17),
      "single trait.*not supported|not supported.*single trait"
    )
  }
  # controls: AE and AD for one trait still run
  expect_s3_class(
    .a1_cp(G, add_QTN_num = 2, epi_QTN_num = 1, add_effect = 0.3,
           epi_effect = 0.3, h2 = 0.5, rep = 1, model = "AE", seed = 17),
    "data.frame")
})

test_that("V1P-F5: single-trait same_add_dom_QTN is rejected explicitly", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 2, add_effect = 0.3, same_add_dom_QTN = TRUE,
           degree_of_dom = 1, h2 = 0.5, rep = 1, model = "AD", seed = 3),
    "same_add_dom_QTN.*single trait|single trait.*same_add_dom_QTN"
  )
})

# ---------------------------------------------------------------------------
# V1C-F5b: model D + indirect LD ; RX-1 / ld window
# ---------------------------------------------------------------------------
test_that("V1C-F5b: model D with indirect LD is rejected up front", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, dom_QTN_num = 1, dom_effect = c(0.2, 0.3), h2 = c(0.5, 0.5),
           rep = 1, model = "D", architecture = "LD", type_of_ld = "indirect",
           seed = 200),
    "indirect"
  )
})

test_that("RX-1: ld_max >= 1 and ld_min > ld_max are rejected (direct and indirect)", {
  G <- .a1_geno()
  for (ty in c("direct", "indirect")) {
    expect_error(
      .a1_cp(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
             rep = 1, model = "A", architecture = "LD", type_of_ld = ty,
             ld_min = 0.2, ld_max = 1, seed = 200),
      "ld_max"
    )
  }
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
           rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
           ld_min = 0.9, ld_max = 0.3, seed = 200),
    "ld_min"
  )
})

test_that("direct LD pair invariants hold at ld_max = 0.99 (control for RX-1)", {
  G <- .a1_geno()
  home <- .a1_home()
  .a1_cp(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
         rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
         ld_min = 0.2, ld_max = 0.99, ld_method = "corr", seed = 200, home = home)
  ld <- .a1_read(home, "LD_Summary_Additive.txt")
  a <- ld[["QTN_for_trait_1"]]
  b <- ld[["QTN_for_trait_2"]]
  expect_false(any(a == b))
  chr <- function(s) G$chr[match(s, G$snp)]
  expect_equal(chr(a), chr(b))
  act <- vapply(seq_along(a), function(k) .a1_ld(G, a[k], b[k]), numeric(1))
  expect_true(all(act >= 0.2 & act <= 0.99))
})

# ---------------------------------------------------------------------------
# V1C-F4: indirect LD -- informative error or valid distinct pairs
# ---------------------------------------------------------------------------
test_that("V1C-F4: indirect LD never returns duplicate/shared QTNs and aborts informatively", {
  G <- .a1_geno()
  outcome <- function(seed) {
    home <- .a1_home()
    r <- tryCatch(
      .a1_cp(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
             rep = 1, model = "A", architecture = "LD", type_of_ld = "indirect",
             ld_min = 0.2, ld_max = 0.8, seed = seed, home = home),
      error = function(e) conditionMessage(e))
    list(r = r, home = home)
  }
  msgs <- character()
  for (s in c(3, 4, 5, 6, 200)) {
    o <- outcome(s)
    if (is.character(o$r)) {
      msgs <- c(msgs, o$r)
      expect_match(o$r, "LD contract")
    } else {
      q <- .a1_read(o$home, "Additive_QTNs.txt")
      t1 <- q$snp[q$trait == "trait_1"]
      t2 <- q$snp[q$trait == "trait_2"]
      expect_length(intersect(t1, t2), 0)
      expect_false(anyDuplicated(t1) > 0)
      expect_false(anyDuplicated(t2) > 0)
    }
  }
  # seeds 3, 4 (duplicate row.names) and 5 (marker causal for both traits)
  # are the audit's failing seeds: they must now produce the LD-contract error.
  expect_gte(length(msgs), 3)
  # the parity seed still runs
  expect_false(is.character(outcome(200)$r))
})

# ---------------------------------------------------------------------------
# V1C-F5a: direct LD dominance / same_add_dom branches: window is enforced
# ---------------------------------------------------------------------------
test_that("V1C-F5a: direct LD same_add_dom pairs are inside [ld_min, ld_max] or rejected", {
  G <- .a1_geno()
  n_err <- 0
  n_ok <- 0
  for (s in c(1:6, 30)) {
    home <- .a1_home()
    r <- tryCatch(
      .a1_cp(G, add_QTN_num = 3, add_effect = c(0.5, 0.3), h2 = c(0.5, 0.5),
             rep = 1, model = "AD", same_add_dom_QTN = TRUE,
             architecture = "LD", type_of_ld = "direct", ld_method = "corr",
             ld_min = 0.2, ld_max = 0.8, seed = s, home = home),
      error = function(e) conditionMessage(e))
    if (is.character(r)) {
      n_err <- n_err + 1
      expect_match(r, "LD contract")
    } else {
      n_ok <- n_ok + 1
      ld <- .a1_read(home, "LD_Summary.txt")
      a <- ld[["QTN_for_trait_1"]]; b <- ld[["QTN_for_trait_2"]]
      act <- vapply(seq_along(a), function(k) .a1_ld(G, a[k], b[k]), numeric(1))
      expect_true(all(act >= 0.2 - 1e-9 & act <= 0.8 + 1e-9))
      expect_equal(ld[[grep("^Actual_LD", names(ld))[1]]], act, tolerance = 1e-6)
    }
  }
  # before the fix the old branch silently returned out-of-window pairs for
  # ~90% of seeds; those runs must now be rejected, and seed 30 still runs
  expect_gt(n_err, 0)
  expect_gt(n_ok, 0)
})

test_that("V1C-F5a: direct LD dominance-only pairs are linked (same chromosome, in window) or rejected", {
  G <- .a1_geno()
  n_err <- 0
  n_ok <- 0
  for (s in c(1, 2, 12)) {
    home <- .a1_home()
    r <- tryCatch(
      .a1_cp(G, dom_QTN_num = 3, dom_effect = c(0.5, 0.3), h2 = c(0.5, 0.5),
             rep = 1, model = "D", architecture = "LD", type_of_ld = "direct",
             ld_method = "corr", ld_min = 0.2, ld_max = 0.8, seed = s,
             home = home),
      error = function(e) conditionMessage(e))
    if (is.character(r)) {
      n_err <- n_err + 1
      expect_match(r, "LD contract")
    } else {
      n_ok <- n_ok + 1
      ld <- .a1_read(home, "LD_Summary_Dominance.txt")
      a <- ld[["QTN_for_trait_1"]]; b <- ld[["QTN_for_trait_2"]]
      expect_equal(G$chr[match(a, G$snp)], G$chr[match(b, G$snp)])
      act <- vapply(seq_along(a), function(k) .a1_ld(G, a[k], b[k]), numeric(1))
      expect_true(all(act >= 0.2 - 1e-9 & act <= 0.8 + 1e-9))
    }
  }
  # the frozen dominance walk does not reset its neighbour pointers after a
  # re-sample, so pairs span chromosomes for most seeds: those are rejected
  expect_gt(n_err, 0)
  expect_gt(n_ok, 0)
})

# ---------------------------------------------------------------------------
# V1C-F6: LD diagnostic trait labels
# ---------------------------------------------------------------------------
test_that("V1C-F6: direct LD trait labels match the markers that generate each trait", {
  G <- .a1_geno()
  home <- .a1_home()
  .a1_cp(G, add_QTN_num = 2, add_effect = c(0.9, 0.5), h2 = c(0.5, 0.5),
         rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
         ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = 200,
         home = home)
  q <- .a1_read(home, "Additive_QTNs.txt")
  gv <- .a1_read(home, "Genetic_values.txt")
  for (tr in 1:2) {
    rows <- q[q$trait == paste0("trait_", tr), ]
    eff <- rows$additive_effect
    X <- as.matrix(t(G[match(rows$snp, G$snp), -(1:5)]))
    val <- as.numeric(X %*% eff)
    expect_gt(abs(stats::cor(val, gv[[paste0("Trait_", tr)]])), 0.999)
  }
  ld <- .a1_read(home, "LD_Summary_Additive.txt")
  expect_equal(sort(ld[["QTN_for_trait_1"]]), sort(q$snp[q$trait == "trait_1"]))
  expect_equal(sort(ld[["QTN_for_trait_2"]]), sort(q$snp[q$trait == "trait_2"]))
})

test_that("V1C-F6: indirect LD summary reports trait-1 LD/markers in the trait-1 columns", {
  G <- .a1_geno()
  home <- .a1_home()
  .a1_cp(G, add_QTN_num = 3, add_effect = c(0.02, 0.05), h2 = c(0.2, 0.4),
         rep = 1, model = "A", architecture = "LD", type_of_ld = "indirect",
         ld_min = 0.2, ld_max = 0.8, ld_method = "corr", seed = 200, home = home)
  q <- .a1_read(home, "Additive_QTNs.txt")
  ld <- .a1_read(home, "LD_Summary_Additive.txt")
  cause <- ld[["SNP_causing_LD"]]
  t1 <- ld[["QTN_for_trait_1"]]
  t2 <- ld[["QTN_for_trait_2"]]
  expect_equal(t1, q$snp[q$trait == "trait_1"])
  expect_equal(t2, q$snp[q$trait == "trait_2"])
  act1 <- vapply(seq_along(cause), function(k) .a1_ld(G, cause[k], t1[k]), numeric(1))
  act2 <- vapply(seq_along(cause), function(k) .a1_ld(G, cause[k], t2[k]), numeric(1))
  expect_equal(ld[["Actual_LD_with_QTN_of_Trait_1"]], act1, tolerance = 1e-6)
  expect_equal(ld[["Actual_LD_with_QTN_of_Trait_2"]], act2, tolerance = 1e-6)
  # trait 1 markers sit downstream (higher marker index) of their cause
  expect_true(all(match(t1, G$snp) > match(cause, G$snp)))
  expect_true(all(q$type[q$trait == "trait_1"] == "QTN_downstream"))
})

# ---------------------------------------------------------------------------
# V1C-F2: residual seed
# ---------------------------------------------------------------------------
test_that("V1C-F2: h2 < 0.05 with rep > 1 (identical replicates) is rejected", {
  G <- .a1_geno()
  for (h in c(0.04, 0.01)) {
    expect_error(
      .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 3, h2 = h,
             model = "A", seed = 7, output_format = "wide"),
      "identical replicates"
    )
  }
  # a single replicate, or h2 >= 0.05, is fine
  expect_s3_class(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.04,
           model = "A", seed = 7), "data.frame")
  ph <- .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 2, h2 = 0.06,
               model = "A", seed = 7, output_format = "wide")
  expect_false(isTRUE(all.equal(ph[[2]], ph[[3]])))
})

test_that("V1C-F2 / Codex F07: the documented residual seed formula is the one in the code", {
  G <- .a1_geno()
  home <- .a1_home()
  .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 3, h2 = 0.5, model = "A",
         seed = 7, output_format = "wide", home = home)
  f <- list.files(home, "^Seed_number_for_3_Reps", full.names = TRUE)
  used <- scan(f, quiet = TRUE)
  expect_equal(used, (7 + 1:3) * round(10 * 0.5))
})

test_that("Codex F07: an overflowing seed is rejected before any draw", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.7,
           model = "A", seed = .Machine$integer.max),
    "seed"
  )
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.7,
           model = "A", seed = "abc"),
    "seed"
  )
})

# ---------------------------------------------------------------------------
# V1C-F7, V1C-F8, V1C-F9, V1C-F17
# ---------------------------------------------------------------------------
test_that("V1C-F7: QTN_list together with same_add_dom_QTN is rejected", {
  G <- .a1_geno()
  ql <- list(add = list(c("ss196521092"), c("ss196469503")))
  expect_error(
    .a1_cp(G, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
           same_add_dom_QTN = TRUE, degree_of_dom = 1, h2 = c(0.5, 0.5),
           rep = 1, model = "AD", seed = 1),
    "same_add_dom_QTN.*QTN_list|QTN_list.*same_add_dom_QTN"
  )
})

test_that("V1C-F8: h2 outside [0, 1] or not finite is rejected", {
  G <- .a1_geno()
  for (h in list(1.5, -0.1, NA_real_, Inf)) {
    expect_error(
      .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = h,
             model = "A", seed = 1),
      "h2"
    )
  }
  # h2 = 1 (no residual) is valid
  expect_s3_class(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 1,
           model = "A", seed = 1), "data.frame")
})

test_that("V1C-F9: to_r with several h2 rows: multi-trait and vary_QTN are rejected, single-trait works", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 2, add_effect = c(0.3, 0.2), ntraits = 2,
           h2 = rbind(c(0.2, 0.3), c(0.5, 0.6)), rep = 1, model = "A", seed = 1),
    "to_r"
  )
  expect_error(
    .a1_cp(G, add_QTN_num = 2, add_effect = 0.3, h2 = c(0.3, 0.6), rep = 2,
           vary_QTN = TRUE, model = "A", seed = 1),
    "to_r"
  )
  res <- .a1_cp(G, add_QTN_num = 2, add_effect = 0.3, h2 = c(0.3, 0.6), rep = 2,
                model = "A", seed = 1)
  expect_type(res, "list")
  expect_named(res, c("h2_0.3", "h2_0.6"))
})

test_that("V1C-F17: ntraits that disagrees with QTN_list is rejected (no silent NA traits)", {
  G <- .a1_geno()
  ql <- list(add = list(c("ss196521092", "ss196469503"),
                        c("ss196521092", "ss196520297")))
  expect_error(
    .a1_cp(G, QTN_list = ql, add_effect = c(0.3, 0.2), h2 = c(0.5, 0.5),
           rep = 1, model = "A", seed = 1),
    "ntraits"
  )
  expect_error(
    .a1_cp(G, QTN_list = ql, ntraits = 3, add_effect = c(0.3, 0.2, 0.1),
           h2 = c(0.5, 0.5, 0.5), rep = 1, model = "A", seed = 1),
    "ntraits.*match|match.*ntraits"
  )
})

# ---------------------------------------------------------------------------
# V1C-F10: residual correlation label
# ---------------------------------------------------------------------------
test_that("V1C-F10: the residual correlation is labelled as the input, not a sample estimate", {
  G <- .a1_geno()
  home <- .a1_home()
  cr <- matrix(c(1, 0.3, 0.3, 1), 2)
  .a1_cp(G, add_QTN_num = 3, add_effect = c(0.3, 0.2), ntraits = 2,
         h2 = c(0.5, 0.5), cor_res = cr, rep = 2, model = "A", seed = 1,
         home = home)
  log <- readLines(file.path(home, "Log_Sim.txt"))
  expect_false(any(grepl("Sample Residual Correlation", log, fixed = TRUE)))
  expect_true(any(grepl("as specified by `cor_res`", log, fixed = TRUE)))
})

# ---------------------------------------------------------------------------
# V1C-F11 / V1P-F9: heterozygote guard
# ---------------------------------------------------------------------------
test_that("V1C-F11: dominance on data without heterozygotes is caught (A+D warns, D stops)", {
  G <- .a1_geno()
  G0 <- G
  v <- as.matrix(G0[, -(1:5)])
  v[v == 0] <- 1
  G0[, -(1:5)] <- v
  expect_warning(
    create_phenotypes(geno_obj = G0, add_QTN_num = 2, dom_QTN_num = 2,
                      add_effect = 0.2, dom_effect = 0.3, h2 = 0.5, rep = 1,
                      model = "AD", seed = 5, to_r = TRUE, verbose = FALSE,
                      home_dir = .a1_home(), output_dir = ""),
    "selected dominance QTNs"
  )
  expect_error(
    .a1_cp(G0, dom_QTN_num = 2, dom_effect = 0.3, h2 = 0.5, rep = 1,
           model = "D", seed = 5),
    "homozygote for the selected"
  )
})

test_that("V1C-F11: a locus set with heterozygotes but no AA homozygote is not rejected", {
  G <- .a1_geno()
  G1 <- G
  v <- as.matrix(G1[, -(1:5)])
  v[v == 1] <- -1                       # only -1 / 0 remain
  G1[, -(1:5)] <- v
  ph <- .a1_cp(G1, dom_QTN_num = 2, dom_effect = 0.5, h2 = 0.5, rep = 1,
               model = "D", seed = 3)
  expect_s3_class(ph, "data.frame")
})

# ---------------------------------------------------------------------------
# V1C-F12: vary_QTN + QTN_variance
# ---------------------------------------------------------------------------
test_that("V1C-F12: vary_QTN + QTN_variance reports each replicate's own QTN variances", {
  G <- .a1_geno()
  home <- .a1_home()
  ph <- .a1_cp(G, add_QTN_num = 2, add_effect = c(0.5, 0.3), ntraits = 2,
               rep = 2, h2 = c(0.5, 0.5), model = "A",
               architecture = "pleiotropic", vary_QTN = TRUE,
               QTN_variance = TRUE, seed = 3, home = home)
  q <- .a1_read(home, "Additive_QTNs.txt")
  pve <- .a1_read(home, "PVE_of_ADD_QTNs_trait_1.txt")
  for (p in 1:2) {
    vp <- stats::var(ph[ph$Rep == p, 2])
    rows <- q[q$rep == p, ]
    for (k in seq_len(nrow(rows))) {
      g <- as.numeric(G[match(rows$snp[k], G$snp), -(1:5)])
      expect_equal(pve[[paste0("QTN_", k)]][p],
                   stats::var(rows$add_eff_t1[k] * g) / vp, tolerance = 1e-8)
    }
  }
})

# ---------------------------------------------------------------------------
# V1C-F13 / Codex F08: user-specified QTN lists
# ---------------------------------------------------------------------------
test_that("V1C-F13: user-specified epistatic lists with more groups than interaction width work", {
  G <- .a1_geno()
  ql <- list(epi = list(c("ss196521092", "ss196469503", "ss196520297", "ss196433798"),
                        c("ss196469557", "ss196527787", "ss196521092", "ss196469503")))
  home <- .a1_home()
  ph <- .a1_cp(G, QTN_list = ql, ntraits = 2, epi_effect = list(c(0.3, 0.2), c(0.1, 0.4)),
               epi_interaction = 2, h2 = c(0.5, 0.5), rep = 1, model = "E",
               seed = 1, home = home)
  expect_s3_class(ph, "data.frame")
  e <- .a1_read(home, "Epistatic_QTNs.txt")
  expect_equal(e$epistatic_effect[e$trait == "trait_1"], c(0.3, 0.3, 0.2, 0.2))
})

test_that("Codex F08: a partially missing marker in QTN_list is named in the error", {
  G <- .a1_geno()
  ql <- list(add = list(c("ss196521092", "not_a_marker"), "ss196469503"))
  expect_error(
    .a1_cp(G, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
           h2 = c(0.5, 0.5), rep = 1, model = "A", seed = 1),
    "not_a_marker"
  )
})

test_that("V1C-F13: qtn_from_user assembles the variance-QTN table from the whole marker vector", {
  G <- .a1_geno()
  old <- setwd(.a1_home())
  on.exit(setwd(old), add = TRUE)
  res <- simplePHENOTYPES:::qtn_from_user(
    genotypes = G,
    QTN_list = list(var = list(c("ss196521092", "ss196469503"))),
    export_gt = FALSE, architecture = "User-Defined",
    same_add_dom_QTN = FALSE, same_mv_QTN = FALSE,
    add = FALSE, dom = FALSE, epi = FALSE, var = TRUE,
    var_effect = list(c(0.1, 0.2)), ntraits = 1, verbose = FALSE)
  expect_equal(ncol(res$var_ef_trait_obj[[1]]), 2)
  vt <- utils::read.delim("Variance_QTNs.txt")
  expect_equal(nrow(vt), 2)
  expect_equal(vt$variance_effect, c(0.1, 0.2))
})

# ---------------------------------------------------------------------------
# V1C-F14: SNP_effect
# ---------------------------------------------------------------------------
test_that("V1C-F14: SNP_effect on numeric geno_obj is rejected; it reaches the loader for file input", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
           model = "A", seed = 1, SNP_effect = "Dom"),
    "SNP_effect"
  )
  hmp <- normalizePath(file.path(testthat::test_path(), "..", "test.hmp.txt"),
                       mustWork = FALSE)
  skip_if_not(file.exists(hmp))
  run <- function(eff) {
    home <- .a1_home()
    suppressMessages(suppressWarnings(create_phenotypes(
      geno_file = hmp, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.7,
      model = "A", seed = 1, SNP_effect = eff, to_r = TRUE, verbose = FALSE,
      home_dir = home, output_dir = "")))
  }
  expect_false(isTRUE(all.equal(run("Add"), run("Dom"))))
})

# ---------------------------------------------------------------------------
# V1C-F15: gds collision
# ---------------------------------------------------------------------------
test_that("V1C-F15: out_geno = 'gds' into a folder that already has the file is an up-front error", {
  G <- .a1_geno()
  home <- .a1_home()
  args <- list(add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
               model = "A", seed = 1, out_geno = "gds", home = home)
  do.call(.a1_cp, c(list(G), args))
  expect_true(any(grepl("\\.gds$", list.files(home))))
  expect_error(do.call(.a1_cp, c(list(G), args)), "already exists")
})

# ---------------------------------------------------------------------------
# V1C-F16 + Codex F15: errors propagate, cleanup is scoped
# ---------------------------------------------------------------------------
test_that("V1C-F16: an inner error is re-signalled (not swallowed into NULL)", {
  G <- .a1_geno()
  testthat::local_mocked_bindings(
    qtn_pleiotropic = function(...) stop("boom from the QTN step", call. = FALSE),
    .package = "simplePHENOTYPES")
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
           model = "A", seed = 1),
    "boom from the QTN step"
  )
})

test_that("Codex F15: the error cleanup removes only the run folder, not unrelated files", {
  G <- .a1_geno()
  home <- .a1_home()
  testthat::local_mocked_bindings(
    qtn_pleiotropic = function(...) {
      dir.create(file.path(home, "unrelated"))
      writeLines("keep me", file.path(home, "unrelated", "important.txt"))
      stop("boom from the QTN step", call. = FALSE)
    },
    .package = "simplePHENOTYPES")
  expect_error(
    suppressMessages(create_phenotypes(
      geno_obj = G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
      model = "A", seed = 1, home_dir = home, verbose = FALSE)),
    "boom")
  expect_true(file.exists(file.path(home, "unrelated", "important.txt")))
  expect_false(dir.exists(file.path(home, "simplePHENOTYPES_output")))
})

# ---------------------------------------------------------------------------
# V1C-F18 / F19 / F20: argument validation
# ---------------------------------------------------------------------------
test_that("V1C-F18/F19 + V1P-F1: cor and cor_res are validated; non-PD is rejected, not repaired", {
  G <- .a1_geno()
  base <- list(add_QTN_num = 3, add_effect = c(0.3, 0.2, 0.1), ntraits = 3,
               h2 = c(0.5, 0.5, 0.5), rep = 1, model = "A", seed = 1)
  bad_pd <- matrix(c(1, .9, -.9, .9, 1, .9, -.9, .9, 1), 3)
  asym   <- matrix(c(1, .3, .1, .5, 1, .2, .1, .2, 1), 3)
  expect_error(do.call(.a1_cp, c(list(G), base, list(cor = bad_pd))), "positive definite")
  expect_error(do.call(.a1_cp, c(list(G), base, list(cor = asym))), "symmetric")
  expect_error(do.call(.a1_cp, c(list(G), base, list(cor = diag(2)))), "3 x 3")
  expect_error(do.call(.a1_cp, c(list(G), base, list(cor_res = bad_pd))), "positive semi-definite")
  expect_error(do.call(.a1_cp, c(list(G), base, list(cor_res = asym))), "symmetric")
  ok <- matrix(c(1, .3, .1, .3, 1, .2, .1, .2, 1), 3)
  ph <- do.call(.a1_cp, c(list(G), base, list(cor = ok)))
  expect_s3_class(ph, "data.frame")
})

test_that("V1C-F19: output_format is validated; wide with ntraits > 1 needs rep >= 2", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
           model = "A", seed = 1, output_format = "foo"),
    "output_format")
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = c(0.3, 0.2), ntraits = 2, rep = 1,
           h2 = c(0.5, 0.5), model = "A", seed = 1, output_format = "wide"),
    "wide.*rep")
})

test_that("V1C-F20: cryptic failures for rep/h2/model/seed/LD arguments become clear errors", {
  G <- .a1_geno()
  ok <- list(add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
             model = "A", seed = 1)
  mod <- function(...) { a <- ok; a[names(list(...))] <- list(...); a }
  expect_error(do.call(.a1_cp, c(list(G), mod(rep = NULL))), "`rep`")
  expect_error(do.call(.a1_cp, c(list(G), mod(rep = 0))), "`rep`")
  expect_error(do.call(.a1_cp, c(list(G), mod(h2 = NULL))), "h2")
  expect_error(do.call(.a1_cp, c(list(G), mod(model = "a"))), "model")
  expect_error(do.call(.a1_cp, c(list(G), mod(add_QTN_num = 1e6))), "Not enough markers")
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = c(0.3, 0.2, 0.1), ntraits = 3,
           h2 = c(0.5, 0.5, 0.5), rep = 1, model = "A", seed = 1,
           architecture = "LD"),
    "ntraits = 2")
  expect_error(do.call(.a1_cp, c(list(G), mod(seed = "abc"))), "seed")
})

test_that("V1C-F20 / Codex F03: QTN_list cannot be combined with architecture = 'LD'", {
  G <- .a1_geno()
  ql <- list(add = list("ss196521092", "ss196469503"))
  expect_error(
    .a1_cp(G, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
           h2 = c(0.5, 0.5), rep = 1, model = "A", architecture = "LD",
           type_of_ld = "direct", ld_min = 0.2, ld_max = 0.9, seed = 22),
    "architecture = \"LD\""
  )
})

# ---------------------------------------------------------------------------
# V1C-F21 / F22 / F23
# ---------------------------------------------------------------------------
test_that("V1C-F21: out_name survives do.call() (no deparse of the whole data set)", {
  G <- .a1_geno()
  home <- .a1_home()
  args <- list(geno_obj = G, add_QTN_num = 3, add_effect = 0.2, rep = 1,
               h2 = 0.5, model = "A", seed = 1, out_geno = "numeric",
               home_dir = home, output_dir = "", verbose = FALSE, to_r = TRUE)
  suppressMessages(do.call(create_phenotypes, args))
  expect_true(file.exists(file.path(home, "geno_obj_numeric.txt")))
})

test_that("V1C-F22: in-memory dose detection: 0/1-only errors clearly, NA in the first column works", {
  G <- .a1_geno()
  G01 <- G[1:200, ]
  v <- as.matrix(G01[, -(1:5)])
  v[v == 0] <- 1; v[v == -1] <- 0
  G01[, -(1:5)] <- v
  expect_error(
    .a1_cp(G01, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
           model = "A", seed = 1),
    "only the values 0 and 1"
  )
  # a 0/1/2 dosage panel whose first sample column contains an NA and no 2/-1
  Gd <- G[1:300, ]
  v <- as.matrix(Gd[, -(1:5)]) + 1
  v[, 1] <- ifelse(v[, 1] == 2, 1, v[, 1])
  v[1, 1] <- NA
  Gd[, -(1:5)] <- v
  ph <- .a1_cp(Gd, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
               model = "A", seed = 1)
  expect_s3_class(ph, "data.frame")
  expect_true(all(is.finite(ph[[2]])))
})

test_that("V1C-F23: constraint() with no constraint returns all rows; MAF bounds are strict", {
  G <- .a1_geno()[1:20, ]
  expect_equal(simplePHENOTYPES:::constraint(G, NULL, NULL, NULL, verbose = FALSE),
               seq_len(20))
  ns <- ncol(G) - 5
  maf <- apply(G[, -(1:5)], 1, function(x) {
    sumx <- ((sum(x) + ns) / ns * 0.5); min(sumx, 1 - sumx)
  })
  m <- maf[3]
  above <- simplePHENOTYPES:::constraint(G, maf_above = m, verbose = FALSE)
  below <- simplePHENOTYPES:::constraint(G, maf_below = m, verbose = FALSE)
  expect_false(3 %in% above)   # strictly greater than maf_above
  expect_false(3 %in% below)   # strictly smaller than maf_below
})

# ---------------------------------------------------------------------------
# Codex F12: RNG kind / state hygiene
# ---------------------------------------------------------------------------
test_that("Codex F12: the caller's RNG kind and state are restored (explicit seed)", {
  G <- .a1_geno()
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(99)
  before_kind <- RNGkind()
  before_state <- .Random.seed
  .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
         model = "A", seed = 1)
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_state)
  # ... and after an error
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 5,
           model = "A", seed = 1))
  expect_identical(RNGkind(), before_kind)
  expect_identical(.Random.seed, before_state)
})

test_that("Codex F12: seed = NULL draws under the caller's stream and successive calls differ", {
  G <- .a1_geno()
  old <- RNGkind()
  on.exit(RNGkind(old[1], old[2], old[3]), add = TRUE)
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(99)
  before_kind <- RNGkind()
  a <- .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5, model = "A")
  expect_identical(RNGkind(), before_kind)
  b <- .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5, model = "A")
  expect_false(isTRUE(all.equal(a[[2]], b[[2]])))
  # the same caller state gives the same master seed (reproducible NULL-seed runs)
  set.seed(99)
  a2 <- .a1_cp(G, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5, model = "A")
  expect_equal(a2[[2]], a[[2]])
})

# ---------------------------------------------------------------------------
# V1P-F3 / V1P-F4: effect specification
# ---------------------------------------------------------------------------
test_that("V1P-F3: an add_effect vector that is neither 1 nor one-per-QTN is rejected", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 3, add_effect = c(0.2, 0.3), h2 = 1, rep = 1,
           model = "A", seed = 1),
    "add_effect.*3 QTN"
  )
})

test_that("V1P-F4: mixing custom and geometric effect specifications is rejected", {
  G <- .a1_geno()
  expect_error(
    .a1_cp(G, add_QTN_num = 2, dom_QTN_num = 2, add_effect = c(0.1, 0.2),
           dom_effect = 0.3, h2 = 0.5, rep = 1, model = "AD", seed = 3),
    "inconsistently"
  )
  # all custom, and all geometric, are fine
  expect_s3_class(
    .a1_cp(G, add_QTN_num = 2, dom_QTN_num = 2, add_effect = c(0.1, 0.2),
           dom_effect = c(0.3, 0.1), h2 = 0.5, rep = 1, model = "AD", seed = 3),
    "data.frame")
  expect_s3_class(
    .a1_cp(G, add_QTN_num = 2, dom_QTN_num = 2, add_effect = 0.1,
           dom_effect = 0.3, h2 = 0.5, rep = 1, model = "AD", seed = 3),
    "data.frame")
})

# ---------------------------------------------------------------------------
# V1P-F10: vQTL argument validation
# ---------------------------------------------------------------------------
test_that("V1P-F10: vQTL with h2 = 0, negative variance effects or several traits is rejected", {
  G <- .a1_geno()
  ok <- list(add_QTN_num = 2, var_QTN_num = 2, add_effect = 0.3,
             var_effect = 0.2, h2 = 0.5, rep = 1, model = "AV", seed = 1,
             same_mv_QTN = FALSE)
  expect_error(do.call(.a1_cp, c(list(G), modifyList(ok, list(h2 = 0)))), "h2")
  expect_error(do.call(.a1_cp, c(list(G), modifyList(ok, list(var_effect = -3)))),
               "negative")
  expect_error(
    .a1_cp(G, add_QTN_num = 2, var_QTN_num = 2, add_effect = c(0.3, 0.2),
           var_effect = 0.2, ntraits = 2, h2 = c(0.5, 0.5), rep = 1,
           model = "AV", seed = 1),
    "single trait")
  ph <- do.call(.a1_cp, c(list(G), ok))
  expect_s3_class(ph, "data.frame")
})

# ---------------------------------------------------------------------------
# tiny marker sets: informative error instead of "'start' is invalid"
# ---------------------------------------------------------------------------
test_that("Codex F04: a marker set too small for the LD search fails informatively", {
  G <- .a1_geno()[1:3, ]
  for (s in 1:4) {
    r <- tryCatch(
      .a1_cp(G, add_QTN_num = 1, add_effect = c(0.2, 0.3), h2 = c(0.5, 0.5),
             rep = 1, model = "A", architecture = "LD", type_of_ld = "direct",
             ld_min = 0.05, ld_max = 0.8, ld_method = "corr", seed = s),
      error = function(e) conditionMessage(e))
    expect_true(is.character(r))
    expect_false(grepl("'start' is invalid", r, fixed = TRUE))
    expect_match(r, "ran past|LD contract|Monomorphic")
  }
})

# ---------------------------------------------------------------------------
# strengthened smoke tests (test-shim.R uses expect_no_error() around a
# function that used to swallow errors)
# ---------------------------------------------------------------------------
test_that("shim (strengthened): AD, pleiotropic 3 traits and residual-correlated runs return finite data", {
  G <- .a1_geno()
  n <- ncol(G) - 5L
  ad <- .a1_cp(G, add_QTN_num = 2, dom_QTN_num = 2, add_effect = 0.3,
               dom_effect = 0.2, rep = 1, h2 = 0.5, model = "AD", seed = 2)
  expect_equal(nrow(ad), n)
  expect_true(all(is.finite(ad[[2]])))
  p3 <- .a1_cp(G, add_QTN_num = 3, add_effect = c(0.04, 0.2, 0.1), ntraits = 3,
               rep = 1, h2 = c(0.2, 0.4, 0.4), architecture = "pleiotropic",
               model = "A", seed = 10)
  expect_equal(dim(p3), c(n, 5))
  expect_true(all(is.finite(as.matrix(p3[, 2:4]))))
})

test_that("Codex F13: in-memory numeric genotypes with values outside -1, 0, 1 are rejected", {
  G <- .a1_geno()
  G7 <- G[1:200, ]
  G7[3, 10] <- 7
  expect_error(
    .a1_cp(G7, add_QTN_num = 3, add_effect = 0.2, rep = 1, h2 = 0.5,
           model = "A", seed = 1),
    "outside -1, 0, 1")
})

test_that("Codex F08: duplicated marker names in the genotype data are rejected for QTN_list markers", {
  G <- .a1_geno()
  Gd <- rbind(G[1:50, ], G[3, ])
  ql <- list(add = list(c(G$snp[3], G$snp[5]), G$snp[7]))
  expect_error(
    .a1_cp(Gd, QTN_list = ql, ntraits = 2, add_effect = c(0.3, 0.2),
           h2 = c(0.5, 0.5), rep = 1, model = "A", seed = 1),
    "duplicated marker names")
})

test_that("V1C-F3: multi-trait same_add_dom_QTN with custom effect lists runs; big effect + full-length vector is rejected", {
  G <- .a1_geno()
  home <- .a1_home()
  ph <- .a1_cp(G, add_QTN_num = 2, add_effect = list(c(0.3, 0.1), c(0.2, 0.05)),
               same_add_dom_QTN = TRUE, degree_of_dom = 0.5, ntraits = 2,
               h2 = c(0.5, 0.5), rep = 1, model = "AD", sim_method = "custom",
               seed = 3, home = home)
  expect_s3_class(ph, "data.frame")
  q <- .a1_read(home, "Additive_QTNs.txt")
  expect_equal(q$add_eff_t1, c(0.3, 0.1))
  # big_add_QTN_effect fills the first QTN, so add_effect must hold 2 values for 3 QTNs
  expect_error(
    .a1_cp(G, add_QTN_num = 3, big_add_QTN_effect = 0.9,
           add_effect = c(0.1, 0.2, 0.3), h2 = 0.5, rep = 1, model = "A",
           seed = 3, sim_method = "geometric"),
    "add_effect")
})
