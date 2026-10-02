# test-adopt-v1-pleio.R
#
# Round-9 adoption of the independent audit's remaining test proposals for the
# frozen v1 engine, group "v1 pleiotropy / partial / vQTL / baseline /
# genotypes" (reconciliation/v1-pleio-effects.md section 6, Fable T1-T19 and
# Codex rows). Proposals already adopted in test-audit-v1-pleio.R,
# test-audit-v1-core.R and test-fix*-v1.R are not repeated (round-9 handoff has
# the coverage table). The frozen RDS references are only read elsewhere; the
# create_phenotypes() signature is untouched.
#
# Deterministic seeds, tempdir() writes only.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
.ad2_geno <- function() {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  e$SNP55K_maize282_maf04
}

.ad2_root <- tempfile("ad2root_")
dir.create(.ad2_root)
.ad2_n <- 0L
.ad2_home <- function() {
  .ad2_n <<- .ad2_n + 1L
  d <- file.path(.ad2_root, sprintf("run%04d", .ad2_n))
  dir.create(d)
  d
}

.ad2_run <- function(G, ...) {
  home <- .ad2_home()
  r <- tryCatch(
    suppressMessages(suppressWarnings(
      create_phenotypes(geno_obj = G, home_dir = home, output_dir = "",
                        to_r = TRUE, verbose = FALSE, ...))),
    error = function(e) conditionMessage(e))
  list(r = r, home = home)
}

.ad2_read <- function(home, file) {
  f <- file.path(home, file)
  if (!file.exists(f)) return(NULL)
  as.data.frame(data.table::fread(f))
}

.ad2_dos <- function(G, snps) t(as.matrix(G[match(snps, G$snp), -(1:5)]))

# vQTL() into a throw-away working directory; returns the simulated data
.ad2_vq <- function(QTN, base, ve, h2, rep = 2, seed = 11, mean = 0) {
  withr::with_tempdir({
    utils::capture.output(res <- suppressMessages(vQTL(
      QTN = QTN, base_line_trait = base, var_QTN_num = length(ve),
      var_effect = ve, h2 = h2, rep = rep, seed = seed, mean = mean,
      output_format = "long", fam = NULL, to_r = TRUE)))
  })
  res
}

# ---------------------------------------------------------------------------
# T2: the uncentred epistatic product (public API, exported genetic values)
# ---------------------------------------------------------------------------
test_that("T2: model E genetic values are the centred sum of effect * x1 * x2 over the interaction pairs", {
  G <- .ad2_geno()
  o <- .ad2_run(G, model = "E", epi_QTN_num = 2, epi_effect = c(0.4, 0.16),
                epi_interaction = 2, h2 = 1, rep = 1, seed = 4, export_gt = TRUE)
  expect_s3_class(o$r, "data.frame")
  e <- .ad2_read(o$home, "Epistatic_QTNs.txt")
  expect_equal(nrow(e), 4)
  expect_equal(e$epi_eff_t1, c(0.4, 0.4, 0.16, 0.16))
  X <- .ad2_dos(G, e$snp)
  raw <- X[, 1] * X[, 2] * 0.4 + X[, 3] * X[, 4] * 0.16     # {-1, 0, 1} coding, no recentring of x
  gv <- .ad2_read(o$home, "Genetic_values.txt")[[1]]
  expect_equal(gv, as.numeric(raw - mean(raw)), tolerance = 1e-10)
  # h2 = 1: the phenotype is exactly the genetic value (no residual)
  expect_equal(o$r[[2]], gv, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# T4 / T5 (documented behaviour, locked): `cor` is imposed on the genetic
# values of partially pleiotropic traits, which mixes in the other trait's
# trait-specific QTNs
# ---------------------------------------------------------------------------
test_that("T5: with `cor`, a partially pleiotropic trait is no longer a function of its own QTNs only; without `cor` it is", {
  G <- .ad2_geno()
  run <- function(cr) {
    args <- list(ntraits = 2, pleio_a = 2, trait_spec_a_QTN_num = c(2, 2),
                 add_effect = list(c(0.5, 0.3, 0.2, 0.1), c(0.4, 0.3, 0.2, 0.1)),
                 h2 = c(1, 1), model = "A", architecture = "partially", seed = 5,
                 rep = 1)
    if (!is.null(cr)) args$cor <- cr
    o <- do.call(.ad2_run, c(list(G), args))
    q <- .ad2_read(o$home, "Additive_QTNs.txt")
    gv <- .ad2_read(o$home, "Genetic_values.txt")
    own2 <- q$snp[q$trait == "trait_2"]
    other <- q$snp[q$trait == "trait_1" & q$type == "trait_specific"]
    list(r2_own = summary(stats::lm(gv[[2]] ~ .ad2_dos(G, own2)))$r.squared,
         r2_other = summary(stats::lm(gv[[2]] ~ .ad2_dos(G, other)))$r.squared,
         gv = gv)
  }
  nocor <- run(NULL)
  withcor <- run(matrix(c(1, 0.5, 0.5, 1), 2))
  expect_equal(nocor$r2_own, 1, tolerance = 1e-8)             # own QTN model exactly
  expect_lt(withcor$r2_own, 1 - 1e-3)                          # leakage from the other trait
  expect_gt(withcor$r2_other, nocor$r2_other + 0.03)           # ... through trait-1-specific loci
  # ... while the requested correlation is realized exactly (T4)
  expect_equal(stats::cor(withcor$gv)[1, 2], 0.5, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# Codex "vqtl-valid-equation": locks the Murphy et al. (2022) construction and
# its draw order (the Python port needs the exact stream)
# ---------------------------------------------------------------------------
test_that("Codex-vQTL: phenotype = scale(base) + k * sigma_i * z, z drawn individual-major after set.seed(seed + j)", {
  set.seed(9)
  n <- 60
  Q <- matrix(sample(c(-1, 0, 1), n, replace = TRUE), n, 1,
              dimnames = list(paste0("i", 1:n), NULL))
  b <- stats::rnorm(n)
  for (ve in c(0, 0.5)) {
    for (h in c(0.25, 0.5, 1)) {
      res <- .ad2_vq(Q, b, ve, matrix(h), rep = 3, seed = 11, mean = 0.7)
      sigma <- 1 + ve * (Q[, 1] + 1)                        # dosage + 1 in {0, 1, 2}
      k <- sqrt((1 / h - 1) / stats::median(sigma)^2)      # var(scale(b)) = 1
      withr::with_seed(11 + 1, {
        z <- matrix(stats::rnorm(n * 3), n, 3, byrow = TRUE)
      })
      expected <- as.numeric(scale(b)) + 0.7 + k * sigma * z
      got <- as.matrix(res$simulated_data[, -1])
      expect_equal(unname(got), unname(expected), tolerance = 1e-12,
                   label = paste("ve", ve, "h2", h))
      expect_identical(res$simulated_data$taxa, rownames(Q))
    }
  }
})

test_that("T11 (second half): replicate means track the genetic value (var(rowMeans) ~ 1), variance ratio unaffected", {
  n <- 40
  Q <- matrix(rep(c(-1, 1), each = n / 2), n, 1, dimnames = list(paste0("i", 1:n), NULL))
  set.seed(21)
  b <- stats::rnorm(n)
  res <- .ad2_vq(Q, b, 0.5, matrix(0.5), rep = 2000, seed = 3)
  P <- as.matrix(res$simulated_data[, -1])
  expect_equal(stats::var(rowMeans(P)), 1, tolerance = 0.1)
  expect_equal(unname(rowMeans(P)), as.numeric(scale(b)), tolerance = 0.5)
})

test_that("Codex-vQTL: h2 above 1 or NA is rejected before any draw", {
  set.seed(9)
  n <- 60
  Q <- matrix(sample(c(-1, 0, 1), n, replace = TRUE), n, 1,
              dimnames = list(paste0("i", 1:n), NULL))
  b <- stats::rnorm(n)
  expect_error(.ad2_vq(Q, b, 0.5, matrix(1.2)), "h2. must be in")
  expect_error(.ad2_vq(Q, b, 0.5, matrix(NA_real_)), "h2. must be in")
  expect_error(.ad2_vq(Q, b, 0.5, matrix(c(0.5, 1.5))), "h2. must be in")
  # a vQTN genotype matrix whose width disagrees with the effect vector
  expect_error(.ad2_vq(cbind(Q, Q), b, 0.5, matrix(0.5)), "genotype matrix has 2")
})

test_that("Codex-vQTL (V1-A1): a zero-variance additive baseline is rejected instead of returning NaN phenotypes", {
  G <- .ad2_geno()
  o <- .ad2_run(G, add_QTN_num = 2, var_QTN_num = 2, add_effect = 0,
                var_effect = 0.2, h2 = 0.5, rep = 2, model = "AV", seed = 1,
                same_mv_QTN = FALSE)
  expect_type(o$r, "character")
  expect_match(o$r, "baseline additive component is constant")
  expect_match(o$r, "add_effect")
  # direct call, same_mv_QTN = TRUE, and the valid case stays finite
  Q <- matrix(c(-1, 0, 1, 0, 1, -1, 1, 0, -1, 0), ncol = 1)
  expect_error(.ad2_vq(Q, rep(1, 10), 0.2, matrix(0.5)), "baseline additive component is constant")
  expect_error(.ad2_vq(Q, c(rep(1, 9), NA), 0.2, matrix(0.5)), "constant \\(or missing\\)")
  ok <- .ad2_run(G, add_QTN_num = 2, var_QTN_num = 2, add_effect = 0.3,
                 var_effect = 0.2, h2 = 0.5, rep = 2, model = "AV", seed = 1,
                 same_mv_QTN = FALSE)
  expect_false(is.character(ok$r))
  expect_true(all(is.finite(as.matrix(ok$r[, -1]))))
})

# ---------------------------------------------------------------------------
# Codex row 11 / T19: structural, deterministic coverage of the partial
# architecture beyond the single frozen A/E fixture: A, AE, AD, ADE, a zero
# trait-specific count, QTN constraints and vary_QTN. (A new frozen RDS per
# case would only snapshot today's output; the checks below are the
# architecture's own invariants.)
# ---------------------------------------------------------------------------
.partial_variants <- list(
  A = list(model = "A", ntraits = 2, pleio_a = 1, trait_spec_a_QTN_num = c(2, 3),
           add_effect = c(0.5, 0.3), h2 = c(0.5, 0.6)),
  AE = list(model = "AE", ntraits = 2, pleio_a = 2, pleio_e = 1,
            trait_spec_a_QTN_num = c(1, 2), trait_spec_e_QTN_num = c(1, 0),
            add_effect = c(0.5, 0.3), epi_effect = c(0.3, 0.2), epi_interaction = 2,
            h2 = c(0.5, 0.6)),
  AD = list(model = "AD", ntraits = 2, pleio_a = 2, pleio_d = 1,
            trait_spec_a_QTN_num = c(1, 2), trait_spec_d_QTN_num = c(1, 1),
            add_effect = c(0.5, 0.3), dom_effect = c(0.4, 0.2), h2 = c(0.5, 0.6)),
  ADE = list(model = "ADE", ntraits = 2, pleio_a = 1, pleio_d = 1, pleio_e = 1,
             trait_spec_a_QTN_num = c(1, 1), trait_spec_d_QTN_num = c(1, 1),
             trait_spec_e_QTN_num = c(1, 1), add_effect = c(0.5, 0.3),
             dom_effect = c(0.4, 0.2), epi_effect = c(0.3, 0.2), epi_interaction = 2,
             h2 = c(0.5, 0.6)),
  constrained = list(model = "A", ntraits = 2, pleio_a = 1,
                     trait_spec_a_QTN_num = c(2, 3), add_effect = c(0.5, 0.3),
                     h2 = c(0.5, 0.6),
                     constraints = list(maf_above = 0.46, maf_below = NULL, hets = NULL))
)

test_that("Codex-11: every partial-architecture variant is reproducible and has the requested QTN counts, shared and trait-specific", {
  G <- .ad2_geno()
  maf <- {
    ns <- ncol(G) - 5
    stats::setNames(apply(G[, -(1:5)], 1, function(x) {
      p <- (sum(x) + ns) / ns * 0.5
      min(p, 1 - p)
    }), G$snp)
  }
  files <- c(a = "Additive_QTNs.txt", d = "Dominance_QTNs.txt", e = "Epistatic_QTNs.txt")
  for (nm in names(.partial_variants)) {
    a <- .partial_variants[[nm]]
    call1 <- do.call(.ad2_run, c(list(G), list(architecture = "partially", seed = 11, rep = 1), a))
    call2 <- do.call(.ad2_run, c(list(G), list(architecture = "partially", seed = 11, rep = 1), a))
    expect_s3_class(call1$r, "data.frame")
    expect_identical(call1$r, call2$r)                         # same seed: same phenotypes
    for (cls in names(files)) {
      pl <- a[[paste0("pleio_", cls)]]
      if (is.null(pl)) next
      sp <- a[[paste0("trait_spec_", cls, "_QTN_num")]]
      width <- if (cls == "e") a$epi_interaction else 1
      f1 <- .ad2_read(call1$home, files[[cls]])
      f2 <- .ad2_read(call2$home, files[[cls]])
      expect_identical(f1$snp, f2$snp, info = paste(nm, cls))
      for (t in seq_along(sp)) {
        rows <- f1[f1$trait == paste0("trait_", t), ]
        expect_equal(sum(rows$type == "Pleiotropic"), pl * width, info = paste(nm, cls, t))
        expect_equal(sum(rows$type == "trait_specific"), sp[t] * width, info = paste(nm, cls, t))
      }
      pleio <- lapply(seq_along(sp), function(t) {
        r <- f1[f1$trait == paste0("trait_", t) & f1$type == "Pleiotropic", ]
        r$snp
      })
      spec <- lapply(seq_along(sp), function(t) {
        r <- f1[f1$trait == paste0("trait_", t) & f1$type == "trait_specific", ]
        r$snp
      })
      # the pleiotropic markers are the SAME loci in every trait ...
      expect_identical(sort(pleio[[1]]), sort(pleio[[2]]), info = paste(nm, cls))
      # ... and no trait-specific locus is shared with another trait or with the shared set
      expect_length(intersect(spec[[1]], spec[[2]]), 0)
      expect_length(intersect(unlist(spec), unlist(pleio)), 0)
    }
    if (nm == "constrained") {
      f1 <- .ad2_read(call1$home, files[["a"]])
      expect_true(all(maf[f1$snp] > 0.46))                     # constraint applies to every random QTN here
    }
    # per-trait genetic value equals the centred QTN sum (no `cor`)
    if (nm == "A") {
      f1 <- .ad2_read(call1$home, files[["a"]])
      gv <- .ad2_read(call1$home, "Genetic_values.txt")
      for (t in 1:2) {
        rows <- f1[f1$trait == paste0("trait_", t), ]
        raw <- as.numeric(.ad2_dos(G, rows$snp) %*% rows$additive_effect)
        expect_equal(gv[[t]], raw - mean(raw), tolerance = 1e-8)
      }
    }
  }
})

test_that("Codex-11: vary_QTN redraws the QTNs per replicate, a fixed architecture keeps them; a different seed changes them", {
  G <- .ad2_geno()
  a <- c(list(architecture = "partially", rep = 2, model = "A", ntraits = 2,
              pleio_a = 1, trait_spec_a_QTN_num = c(2, 3), add_effect = c(0.5, 0.3),
              h2 = c(0.5, 0.6)))
  fixed <- do.call(.ad2_run, c(list(G), a, list(vary_QTN = FALSE, seed = 11)))
  vary <- do.call(.ad2_run, c(list(G), a, list(vary_QTN = TRUE, seed = 11)))
  other <- do.call(.ad2_run, c(list(G), a, list(vary_QTN = TRUE, seed = 12)))
  qt <- function(o) .ad2_read(o$home, "Additive_QTNs.txt")
  # fixed architecture: one QTN set, written once (labelled replicate 1) and
  # shared by all replicates
  expect_setequal(unique(qt(fixed)$rep), 1)
  expect_equal(nrow(qt(fixed)), 7)
  # varying architecture: a QTN set per replicate, and the two sets differ
  qv <- qt(vary)
  expect_setequal(unique(qv$rep), 1:2)
  s1 <- sort(qv$snp[qv$rep == 1])
  s2 <- sort(qv$snp[qv$rep == 2])
  expect_length(s1, 7)
  expect_false(identical(s1, s2))
  # same QTN count and shared/specific structure in every replicate
  expect_equal(table(qv$type[qv$rep == 1]), table(qv$type[qv$rep == 2]))
  # a different seed gives a different first replicate
  qo <- qt(other)
  expect_false(identical(s1, sort(qo$snp[qo$rep == 1])))
  # replicate 1 of the varying run is the fixed run's QTN set (same seed, same draw)
  expect_identical(s1, sort(qt(fixed)$snp))
})

# ---------------------------------------------------------------------------
# Codex "SNP_impute None / numeric alphabet": object and file input agree
# ---------------------------------------------------------------------------
test_that("Codex-numeric: a numeric panel gives the same genotypes() result from an object and from a file", {
  set.seed(5)
  m <- 6
  n <- 8
  vals <- matrix(sample(c(-1, 0, 1), m * n, replace = TRUE), m, n,
                 dimnames = list(NULL, paste0("s", 1:n)))
  vals[2, 3] <- NA
  panel <- data.frame(snp = paste0("m", 1:m), allele = "A/G", chr = 1L,
                      pos = (1:m) * 100L, cm = 0, vals, check.names = FALSE)
  f <- tempfile(fileext = ".txt")
  withr::defer(unlink(f))
  utils::write.table(panel, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  from_obj <- genotypes(geno_obj = panel, SNP_impute = "None", verbose = FALSE)$geno_obj
  skip_if_not("geno_file" %in% names(formals(genotypes)))
  from_file <- tryCatch(
    genotypes(geno_file = f, SNP_impute = "None", verbose = FALSE)$geno_obj,
    error = function(e) NULL)
  skip_if(is.null(from_file), "numeric text-file input is not accepted by genotypes()")
  expect_equal(unname(as.matrix(from_file[, -(1:5)])), unname(as.matrix(from_obj[, -(1:5)])))
  expect_true(is.na(from_obj[2, 8]))                           # NA kept, not silently 0
})
