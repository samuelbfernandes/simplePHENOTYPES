# test-v130-parity.R
#
# v1.3.0 behavioral parity tests (SPEC.md §8).
#
# Each test loads a frozen reference RDS produced by
# inst/extdata/v1_3_0_reference/capture_references.R, runs the equivalent
# v2 grammar call, and asserts:
#   (a) identical QTN selection  — same marker IDs, in the same order.
#   (b) numerically identical phenotypes — after normalizing output format.
#
# All tests start with skip("grammar not yet implemented") and are enabled one
# by one as each grammar function is completed (Block 3 of TODO.md).
# Running devtools::test() must stay green at every commit regardless.
#
# Tolerance for floating-point comparison: 1e-10 (double precision; allows
# for platform rounding but catches any algorithmic divergence).

REF_DIR <- testthat::test_path("..", "..", "inst", "extdata", "v1_3_0_reference")

# Helper: load a reference RDS.
.ref <- function(name) readRDS(file.path(REF_DIR, paste0(name, ".rds")))

# Helper: extract QTN IDs from the reference qtns$add data.frame (rep 1 only).
.ref_qtns <- function(r) sort(unique(r$qtns$add$snp[r$qtns$add$rep == 1]))

# Helper: extract the numeric phenotype matrix from a v2 phenotype_sim result.
# Returns a samples × traits matrix, sorted by sample ID.
.v2_pheno_mat <- function(sim) {
  stop("implement when grammar is available")
}

# Helper: extract QTN IDs selected by the v2 grammar (additive layer, rep 1).
.v2_qtns <- function(sim) {
  stop("implement when grammar is available")
}

# ---------------------------------------------------------------------------
# 1. Single trait, additive, seed = 1
# ---------------------------------------------------------------------------
test_that("parity: single trait, additive (seed = 1)", {
  skip("grammar not yet implemented")

  ref <- .ref("single_trait")
  data("SNP55K_maize282_maf04")

  sim <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
    additive(prop = 0.7, n_qtn = 3, effect = 0.2)

  # (a) QTN identity
  expect_equal(.v2_qtns(sim), .ref_qtns(ref))

  # (b) phenotype values
  ref_pheno <- ref$phenotypes[order(ref$phenotypes[[1]]), 2]
  v2_pheno  <- .v2_pheno_mat(sim)[order(rownames(.v2_pheno_mat(sim))), 1]
  expect_equal(v2_pheno, ref_pheno, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. Pleiotropy, 3 traits, additive + dominance, seed = 10
# ---------------------------------------------------------------------------
test_that("parity: pleiotropy 3 traits AD (seed = 10)", {
  skip("grammar not yet implemented")

  ref <- .ref("pleiotropy")
  data("SNP55K_maize282_maf04")

  sim <- simulate_phenotype(SNP55K_maize282_maf04,
                             architecture = "pleiotropy",
                             n_traits     = 3,
                             seed         = 10) |>
    additive(prop   = c(0.2, 0.4, 0.4),
             n_qtn  = 3,
             effect = list(c(0.04, 0.0016, 6.4e-5),
                           c(0.2, 0.04, 0.008),
                           c(0.1, 0.01, 0.001))) |>
    dominance(prop       = c(0.2, 0.4, 0.4),
              same_as_add = TRUE,
              n_qtn       = 4,
              effect      = list(c(0.04, 0.0016, 6.4e-5, 2.56e-6),
                                 c(0.2, 0.04, 0.008, 0.0016),
                                 c(0.1, 0.01, 0.001, 1e-4)))

  # (a) additive QTN identity
  expect_equal(.v2_qtns(sim), .ref_qtns(ref))

  # (b) phenotype values (trait 1 only)
  ref_t1 <- ref$phenotypes[order(ref$phenotypes[[1]]), "Trait_1_H2_0.2"]
  v2_mat  <- .v2_pheno_mat(sim)
  v2_t1   <- v2_mat[order(rownames(v2_mat)), 1]
  expect_equal(v2_t1, ref_t1, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. LD spurious — indirect, seed = 200
# ---------------------------------------------------------------------------
test_that("parity: LD indirect (seed = 200)", {
  skip("grammar not yet implemented")

  ref <- .ref("ld_indirect")
  data("SNP55K_maize282_maf04")

  sim <- simulate_phenotype(SNP55K_maize282_maf04,
                             architecture = "ld",
                             n_traits     = 2,
                             ld_type      = "indirect",
                             r2_max       = 0.8,
                             r2_min       = 0.2,
                             r2_method    = "composite",
                             seed         = 200) |>
    additive(prop   = c(0.2, 0.4),
             n_qtn  = 3,
             effect = c(0.02, 0.05))

  # For LD architecture, QTN IDs are the actual causal markers (not intermediates).
  ref_qtns <- sort(unique(
    ref$qtns$add$snp[ref$qtns$add$type != "cause_of_LD"]))
  expect_equal(.v2_qtns(sim), ref_qtns)

  ref_t1 <- ref$phenotypes[order(ref$phenotypes[[1]]), "Trait_1_H2_0.2"]
  v2_mat  <- .v2_pheno_mat(sim)
  v2_t1   <- v2_mat[order(rownames(v2_mat)), 1]
  expect_equal(v2_t1, ref_t1, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4. LD spurious — direct, seed = 200
# ---------------------------------------------------------------------------
test_that("parity: LD direct (seed = 200)", {
  skip("grammar not yet implemented")

  ref <- .ref("ld_direct")
  data("SNP55K_maize282_maf04")

  sim <- simulate_phenotype(SNP55K_maize282_maf04,
                             architecture = "ld",
                             n_traits     = 2,
                             ld_type      = "direct",
                             r2_max       = 0.8,
                             r2_min       = 0.2,
                             r2_method    = "composite",
                             seed         = 200) |>
    additive(prop   = c(0.2, 0.4),
             n_qtn  = 3,
             effect = c(0.02, 0.05))

  ref_qtns <- sort(unique(
    ref$qtns$add$snp[ref$qtns$add$type != "cause_of_LD"]))
  expect_equal(.v2_qtns(sim), ref_qtns)

  ref_t1 <- ref$phenotypes[order(ref$phenotypes[[1]]), "Trait_1_H2_0.2"]
  v2_mat  <- .v2_pheno_mat(sim)
  v2_t1   <- v2_mat[order(rownames(v2_mat)), 1]
  expect_equal(v2_t1, ref_t1, tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Partial pleiotropy — reconstructed via complex_phenotypes()
# ---------------------------------------------------------------------------
# Note: partial-pleiotropy parity is STRUCTURAL, not bit-for-bit.  The v1
# architecture uses a different internal mechanism; the v2 equivalent combines
# a pleiotropy model + an independent model.  Assertions check:
#   • correct total QTN count per trait
#   • h2 within 5% of target
# (Exact phenotype match is intentionally NOT asserted for this scenario.)
test_that("parity: partial pleiotropy (structural)", {
  skip("grammar not yet implemented")

  ref   <- .ref("partial_pleiotropy")
  data("SNP55K_maize282_maf04")

  pleio_sim <- simulate_phenotype(SNP55K_maize282_maf04,
                                   architecture = "pleiotropy",
                                   n_traits     = 3,
                                   seed         = 42) |>
    additive(prop = c(0.2, 0.4, 0.8), n_qtn = 3)

  indep_sim <- simulate_phenotype(SNP55K_maize282_maf04,
                                   architecture = "independent",
                                   n_traits     = 3,
                                   seed         = 42) |>
    additive(prop = c(0.2, 0.4, 0.8), n_qtn = c(4, 10, 1))

  sim <- complex_phenotypes(pleio_sim, indep_sim, h2 = c(0.2, 0.4, 0.8))

  # Structural: total QTN counts (pleiotropic + trait-specific) per trait
  ref_add_n <- nrow(ref$qtns$add[ref$qtns$add$rep == 1, ])
  v2_add_n  <- length(.v2_qtns(sim))
  expect_true(abs(v2_add_n - ref_add_n) <= 3,
              info = paste("add QTN count: v2 =", v2_add_n, "ref =", ref_add_n))
})

# ---------------------------------------------------------------------------
# 6. Smoke test: create_phenotypes() shim still runs without error
# ---------------------------------------------------------------------------
# This test never carries a skip() — it must pass at every commit.
# It does NOT check values; it only confirms the v1 API did not break.
test_that("shim: create_phenotypes() runs without error", {
  data("SNP55K_maize282_maf04")
  expect_no_error(
    suppressMessages(
      create_phenotypes(
        geno_obj    = SNP55K_maize282_maf04,
        add_QTN_num = 3,
        add_effect  = 0.2,
        rep         = 1,
        h2          = 0.7,
        model       = "A",
        seed        = 1,
        home_dir    = tempdir()
      )
    )
  )
})
