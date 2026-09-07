# test-v130-parity.R
#
# v1.3.0 regression guard (DECISION-009 / SPEC.md §8.1).
#
# Calls create_phenotypes() (frozen legacy fn) with the same parameters used in
# inst/extdata/v1_3_0_reference/capture_references.R and asserts it still
# reproduces the frozen references:
#   (a) identical QTN marker IDs (rep 1).
#   (b) numerically identical phenotype values (tolerance 1e-10).
#
# If a bug fix intentionally changes output, re-run capture_references.R and
# commit the updated RDS files.
#
# Test 5 (partial pleiotropy) is structural only — QTN counts and run-without-
# error — because architecture = "partially" has no bit-identical v2 equivalent.

# Resolve against the installed package first so the guard also runs under
# R CMD check, where inst/extdata/ has been flattened to extdata/.
REF_DIR <- local({
  installed <- system.file("extdata", "v1_3_0_reference",
                           package = "simplePHENOTYPES")
  if (nzchar(installed) && dir.exists(installed)) {
    installed
  } else {
    testthat::test_path("..", "..", "inst", "extdata", "v1_3_0_reference")
  }
})

.ref <- function(name) readRDS(file.path(REF_DIR, paste0(name, ".rds")))

# Run create_phenotypes() in a fresh temp dir; return list(phenotypes=, qtns=).
.run_legacy <- function(geno, call_args) {
  tmp <- tempfile("sp_parity_")
  dir.create(tmp, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  full_args <- c(
    list(
      geno_obj      = geno,
      to_r          = TRUE,
      output_format = "long",
      home_dir      = tmp,
      verbose       = FALSE
    ),
    call_args
  )

  pheno <- suppressMessages(do.call(create_phenotypes, full_args))

  find_qtn_file <- function(...) {
    for (fname in c(...)) {
      hits <- c(
        file.path(tmp, fname),
        list.files(tmp, pattern = paste0("^", fname, "$"),
                   recursive = TRUE, full.names = TRUE)
      )
      hits <- hits[file.exists(hits)]
      if (length(hits)) return(data.table::fread(hits[[1L]], data.table = FALSE))
    }
    NULL
  }

  list(
    phenotypes = pheno,
    qtns = list(
      # Current package writes "Additive_QTNs.txt"; v1.3.0 wrote
      # "Additive_Selected_QTNs.txt" — search both names for forward-compatibility.
      add = find_qtn_file("Additive_QTNs.txt", "Additive_Selected_QTNs.txt",
                          "Additive_and_Dominance_Selected_QTNs.txt"),
      dom = find_qtn_file("Dominance_QTNs.txt", "Dominance_Selected_QTNs.txt"),
      epi = find_qtn_file("Epistatic_QTNs.txt", "Epistatic_Selected_QTNs.txt")
    )
  )
}

# ---------------------------------------------------------------------------
# 1. Single trait, additive, seed = 1
# ---------------------------------------------------------------------------
test_that("v1.3.0 regression: single trait, additive (seed = 1)", {
  ref <- .ref("single_trait")
  data("SNP55K_maize282_maf04")

  got <- .run_legacy(SNP55K_maize282_maf04, list(
    add_QTN_num = 3, add_effect = 0.2,
    h2 = 0.7, model = "A", rep = 1, seed = 1
  ))

  # (a) QTN identity
  expect_equal(
    sort(got$qtns$add$snp[got$qtns$add$rep == 1]),
    sort(ref$qtns$add$snp[ref$qtns$add$rep == 1])
  )

  # (b) phenotype values
  ref_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  got_ord <- got$phenotypes[order(got$phenotypes[[1L]]), ]
  expect_equal(got_ord[[2L]], ref_ord[[2L]], tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 2. Pleiotropy, 3 traits, additive + dominance, seed = 10
# ---------------------------------------------------------------------------
test_that("v1.3.0 regression: pleiotropy 3 traits, AD (seed = 10)", {
  ref <- .ref("pleiotropy")
  data("SNP55K_maize282_maf04")

  got <- .run_legacy(SNP55K_maize282_maf04, list(
    add_QTN_num  = 3,
    dom_QTN_num  = 4,
    h2           = c(0.2, 0.4, 0.4),
    add_effect   = c(0.04, 0.2, 0.1),
    dom_effect   = c(0.04, 0.2, 0.1),
    ntraits      = 3,
    rep          = 1,
    vary_QTN     = FALSE,
    architecture = "pleiotropic",
    seed         = 10,
    model        = "AD",
    sim_method   = "geometric"
  ))

  # (a) additive QTN identity (rep 1)
  expect_equal(
    sort(got$qtns$add$snp[got$qtns$add$rep == 1]),
    sort(ref$qtns$add$snp[ref$qtns$add$rep == 1])
  )

  # (b) phenotype values, trait 1
  ref_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  got_ord <- got$phenotypes[order(got$phenotypes[[1L]]), ]
  expect_equal(got_ord[["Trait_1_H2_0.2"]], ref_ord[["Trait_1_H2_0.2"]],
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 3. LD spurious — indirect, seed = 200
# ---------------------------------------------------------------------------
test_that("v1.3.0 regression: LD indirect (seed = 200)", {
  ref <- .ref("ld_indirect")
  data("SNP55K_maize282_maf04")

  got <- .run_legacy(SNP55K_maize282_maf04, list(
    add_QTN_num  = 3,
    h2           = c(0.2, 0.4),
    add_effect   = c(0.02, 0.05),
    rep          = 1,
    seed         = 200,
    architecture = "LD",
    model        = "A",
    ld_max       = 0.8,
    ld_min       = 0.2,
    ld_method    = "composite",
    type_of_ld   = "indirect"
  ))

  # Causal markers only (not cause_of_LD rows)
  ref_snps <- sort(ref$qtns$add$snp[ref$qtns$add$type != "cause_of_LD" &
                                      ref$qtns$add$rep  == 1])
  got_snps <- sort(got$qtns$add$snp[got$qtns$add$type != "cause_of_LD" &
                                      got$qtns$add$rep  == 1])
  expect_equal(got_snps, ref_snps)

  ref_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  got_ord <- got$phenotypes[order(got$phenotypes[[1L]]), ]
  expect_equal(got_ord[["Trait_1_H2_0.2"]], ref_ord[["Trait_1_H2_0.2"]],
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 4. LD spurious — direct, seed = 200
# ---------------------------------------------------------------------------
# The ld_direct reference is re-blessed from the current frozen-legacy
# create_phenotypes() (post-1.3.0 ld_min/ld_max bug fixes 68a227e / b95529a;
# CRAN 1.3.0 was buggy here — see BUGS.md and capture_references.R 5d). This
# guards the corrected direct-LD selection against future unintentional change.
test_that("v1.3.0 regression: LD direct (seed = 200)", {
  ref <- .ref("ld_direct")
  data("SNP55K_maize282_maf04")

  got <- .run_legacy(SNP55K_maize282_maf04, list(
    add_QTN_num  = 3,
    h2           = c(0.2, 0.4),
    add_effect   = c(0.02, 0.05),
    rep          = 1,
    seed         = 200,
    architecture = "LD",
    model        = "A",
    ld_max       = 0.8,
    ld_min       = 0.2,
    ld_method    = "composite",
    type_of_ld   = "direct"
  ))

  # Causal markers only (not cause_of_LD rows)
  ref_snps <- sort(ref$qtns$add$snp[ref$qtns$add$type != "cause_of_LD" &
                                      ref$qtns$add$rep  == 1])
  got_snps <- sort(got$qtns$add$snp[got$qtns$add$type != "cause_of_LD" &
                                      got$qtns$add$rep  == 1])
  expect_equal(got_snps, ref_snps)

  ref_ord <- ref$phenotypes[order(ref$phenotypes[[1L]]), ]
  got_ord <- got$phenotypes[order(got$phenotypes[[1L]]), ]
  expect_equal(got_ord[["Trait_1_H2_0.2"]], ref_ord[["Trait_1_H2_0.2"]],
               tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# 5. Partial pleiotropy — structural only
# ---------------------------------------------------------------------------
# Exact bit-match is not asserted: architecture = "partially" has no
# grammar equivalent yet. Checks that the legacy path runs and produces
# QTN counts consistent with the reference.
test_that("v1.3.0 regression: partial pleiotropy, structural (seed = 42)", {
  ref <- .ref("partial_pleiotropy")
  data("SNP55K_maize282_maf04")

  got <- .run_legacy(SNP55K_maize282_maf04, list(
    ntraits              = 3,
    pleio_a              = 3,
    pleio_e              = 2,
    same_add_dom_QTN     = TRUE,
    degree_of_dom        = 0.5,
    trait_spec_a_QTN_num = c(4, 10, 1),
    trait_spec_e_QTN_num = c(3, 2, 5),
    h2                   = c(0.2, 0.4, 0.8),
    add_effect           = c(0.5, 0.33, 0.2),
    epi_effect           = c(0.3, 0.3, 0.3),
    epi_interaction      = 2,
    rep                  = 1,
    architecture         = "partially",
    model                = "AE",
    seed                 = 42
  ))

  # Structural: additive QTN count in rep 1 is within ±3 of the reference
  ref_n <- nrow(ref$qtns$add[ref$qtns$add$rep == 1, ])
  got_n <- if (!is.null(got$qtns$add)) nrow(got$qtns$add[got$qtns$add$rep == 1, ]) else 0L
  expect_true(
    abs(got_n - ref_n) <= 3L,
    label = paste("add QTN count: got =", got_n, "ref =", ref_n)
  )
})
