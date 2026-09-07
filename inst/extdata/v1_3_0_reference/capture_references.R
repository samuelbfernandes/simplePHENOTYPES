# capture_references.R
#
# Run this script ONCE before implementing the v2 grammar.
# It installs simplePHENOTYPES 1.3.0 into a versioned local library, runs each
# README/vignette scenario with big_add_QTN_effect removed, and saves QTN
# marker names + phenotype values as RDS files in this directory.
#
# These RDS files gate all v1.3.0 parity tests (test-v130-parity.R).
#
# Usage (from package root):
#   Rscript inst/extdata/v1_3_0_reference/capture_references.R
#
# Requirements: remotes, withr, data.table

stopifnot(
  requireNamespace("remotes",    quietly = TRUE),
  requireNamespace("withr",      quietly = TRUE),
  requireNamespace("data.table", quietly = TRUE)
)

# ---------------------------------------------------------------------------
# 1. Install v1.3.0 into a dedicated library
# ---------------------------------------------------------------------------

lib130 <- path.expand("~/.R/simplePHENOTYPES_v130_lib")
dir.create(lib130, recursive = TRUE, showWarnings = FALSE)

installed_ver <- tryCatch(
  packageDescription("simplePHENOTYPES", lib.loc = lib130)$Version,
  error = function(e) NA_character_
)

if (is.na(installed_ver) || installed_ver != "1.3.0") {
  message("Installing simplePHENOTYPES 1.3.0 into ", lib130, " ...")
  remotes::install_version(
    "simplePHENOTYPES",
    version = "1.3.0",
    lib     = lib130,
    repos   = "https://cloud.r-project.org",
    upgrade = "never"
  )
} else {
  message("simplePHENOTYPES 1.3.0 already present at ", lib130)
}

# Load the v1.3.0 namespace without attaching it — dev version stays active.
ns130 <- loadNamespace("simplePHENOTYPES", lib.loc = lib130)

# ---------------------------------------------------------------------------
# 2. Dataset
# ---------------------------------------------------------------------------

data("SNP55K_maize282_maf04", package = "simplePHENOTYPES",
     envir = environment())
geno <- SNP55K_maize282_maf04

# ---------------------------------------------------------------------------
# 3. Resolve output directory (where RDS files are saved)
# ---------------------------------------------------------------------------

this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(f) && nzchar(f)) normalizePath(f) else stop()
}, error = function(e) {
  # Fallback when sourced() interactively from the package root
  normalizePath("inst/extdata/v1_3_0_reference/capture_references.R")
})
out_dir <- dirname(this_file)

# ---------------------------------------------------------------------------
# 4. Helpers
# ---------------------------------------------------------------------------

# After create_phenotypes() returns, read whichever *_QTNs.txt files exist
# inside `dir`. No output_dir is passed to create_phenotypes, so v1.3.0
# writes QTN files directly into home_dir (= the temp working dir).
read_qtns <- function(dir) {
  # v1.3.0 writes *_Selected_QTNs.txt to home_dir. The additive file name
  # varies by architecture: "Additive_Selected_QTNs.txt" for pleiotropic/LD/
  # independent, and "Additive_and_Dominance_Selected_QTNs.txt" when
  # same_add_dom_QTN = TRUE in the partially-pleiotropic architecture.
  find_file <- function(...) {
    for (fname in c(...)) {
      paths <- c(
        file.path(dir, fname),
        list.files(dir, pattern = paste0("^", fname, "$"),
                   recursive = TRUE, full.names = TRUE)
      )
      hit <- paths[file.exists(paths)]
      if (length(hit)) return(data.table::fread(hit[[1]], data.table = FALSE))
    }
    NULL
  }

  # The current package writes "Additive_QTNs.txt" etc.; v1.3.0 wrote the
  # "*_Selected_QTNs.txt" names. Search both so this works whether a scenario is
  # captured from the v1.3.0 namespace or the current frozen-legacy one (5d).
  list(
    add = find_file("Additive_QTNs.txt", "Additive_Selected_QTNs.txt",
                    "Additive_and_Dominance_Selected_QTNs.txt"),
    dom = find_file("Dominance_QTNs.txt", "Dominance_Selected_QTNs.txt"),
    epi = find_file("Epistatic_QTNs.txt", "Epistatic_Selected_QTNs.txt")
  )
}

# Run one scenario inside a clean temp directory and return the reference list.
# output_dir is intentionally omitted so QTN files land in home_dir (= tmp).
#
# `cp` is the create_phenotypes() to drive: ns130$create_phenotypes (CRAN 1.3.0)
# for the bit-frozen scenarios, or the current package's create_phenotypes for
# scenarios that must reflect post-1.3.0 frozen-legacy bug fixes (see 5d).
# `version` labels the saved reference accordingly.
run_scenario <- function(scenario, call_args,
                         cp = ns130$create_phenotypes, version = "1.3.0") {
  tmp <- tempfile(pattern = paste0("sp130_", scenario, "_"))
  dir.create(tmp)
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

  message("\n--- Scenario: ", scenario, " ---")
  pheno <- withr::with_dir(tmp, {
    do.call(cp, full_args)
  })

  qtns <- read_qtns(tmp)

  list(
    version    = version,
    captured   = Sys.time(),
    scenario   = scenario,
    call       = call_args,
    qtns       = qtns,
    phenotypes = pheno
  )
}

# ---------------------------------------------------------------------------
# 5. Reference scenarios
# ---------------------------------------------------------------------------

# -- 5a. Single trait --------------------------------------------------------
# Model A, 1 trait. big_add_QTN_effect omitted (SPEC §8.1).
# seed = 1 assigned (README has no seed; 1 is the canonical reference).
ref <- run_scenario(
  "single_trait",
  list(
    add_QTN_num = 3,
    add_effect  = 0.2,
    h2          = 0.7,
    model       = "A",
    rep         = 1,
    seed        = 1
  )
)
saveRDS(ref, file.path(out_dir, "single_trait.rds"))
message("Saved single_trait.rds")

# -- 5b. Pleiotropy ----------------------------------------------------------
# Model AD, 3 traits, architecture = "pleiotropic", seed = 10.
# big_add_QTN_effect omitted (SPEC §8.1). rep = 1 locks the RNG stream.
ref <- run_scenario(
  "pleiotropy",
  list(
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
  )
)
saveRDS(ref, file.path(out_dir, "pleiotropy.rds"))
message("Saved pleiotropy.rds")

# -- 5c. LD indirect ---------------------------------------------------------
# Model A, 2 traits, architecture = "LD", type_of_ld = "indirect", seed = 200.
# README has no big_add_QTN_effect for this scenario. rep = 1.
ref <- run_scenario(
  "ld_indirect",
  list(
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
  )
)
saveRDS(ref, file.path(out_dir, "ld_indirect.rds"))
message("Saved ld_indirect.rds")

# -- 5d. LD direct -----------------------------------------------------------
# Identical to 5c except type_of_ld = "direct".
#
# NOTE (re-blessed): unlike the other scenarios, this one is captured from the
# CURRENT frozen-legacy create_phenotypes(), NOT CRAN 1.3.0. CRAN 1.3.0 had a
# bug in the direct-LD ld_min/ld_max candidate-acceptance loop, fixed after
# release by commits 68a227e and b95529a (both "fixed bug on ld_min/ld_max").
# The frozen-legacy engine includes those fixes, so its correct direct-LD output
# differs from 1.3.0 at the 3rd causal marker. Per DECISION-009 this reference is
# deliberately re-blessed to the post-fix output. `create_phenotypes` here is the
# current package's (the dev namespace stays active alongside ns130).
ref <- run_scenario(
  "ld_direct",
  list(
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
  ),
  cp      = create_phenotypes,
  version = "1.3.0+ldfix (re-blessed from current frozen create_phenotypes)"
)
saveRDS(ref, file.path(out_dir, "ld_direct.rds"))
message("Saved ld_direct.rds (re-blessed from current package)")

# -- 5e. Partial pleiotropy --------------------------------------------------
# architecture = "partially"; seed = 42 (assigned — README has no seed).
# Omits cor / cor_res: v2 uses complex_phenotypes() which does not replicate
# the v1 genetic-correlation mechanism numerically. Parity for this scenario
# is structural (correct QTN counts and h2 range), not bit-for-bit.
ref <- run_scenario(
  "partial_pleiotropy",
  list(
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
  )
)
saveRDS(ref, file.path(out_dir, "partial_pleiotropy.rds"))
message("Saved partial_pleiotropy.rds")

# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

rds_files <- sort(list.files(out_dir, pattern = "\\.rds$", full.names = TRUE))
message("\nReference files written to: ", out_dir)
for (f in rds_files) {
  r        <- readRDS(f)
  add_rows <- if (!is.null(r$qtns$add)) nrow(r$qtns$add) else 0L
  phe_rows <- if (!is.null(r$phenotypes)) nrow(r$phenotypes) else 0L
  message(sprintf("  %-28s  add QTN rows: %3d  phenotype rows: %5d",
                  basename(f), add_rows, phe_rows))
}
message("\nDone. Commit the .rds files to inst/extdata/v1_3_0_reference/")
message("Do NOT re-run this script after v2 implementation begins.")
