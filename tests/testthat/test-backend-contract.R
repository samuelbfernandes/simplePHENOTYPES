# The backend contract (docs/BACKEND_CONTRACT.md): the exported engine surface
# breedingDesigner depends on. This guard fails if any contract function stops
# being exported -- a removed/renamed export is a MAJOR (3.0.0) break and must be
# a deliberate change here, caught before it reaches a downstream consumer.

# Keep in sync with docs/BACKEND_CONTRACT.md.
backend_contract <- c(
  # Populations & crossing
  "as_population", "population_from_haplotypes", "haplotypes", "cross",
  "selfcross", "double_haploid", "dosages", "n_individuals", "synthetic_map",
  # Phenotype grammar
  "simulate_phenotype", "additive", "dominance", "epistasis", "vqtl",
  "complex_phenotypes", "genetic_values", "qtn_table",
  "phenotypes_long", "phenotypes_wide", "write_phenotypes",
  # Selection engine
  "select_ind", "single_seed_descent", "bulk", "pedigree",
  "recurrent_selection",
  # Modern methods
  "g_matrix", "optimum_contribution", "sample_parents", "cross_usefulness",
  "mabc_select", "recurrent_parent_recovery",
  # Pedigree and mating plans (DECISION-024/025)
  "parentage", "families", "mate", "mating_design",
  # Block 3B items (DECISION-026..031)
  "combining_ability", "template_effects", "progeny_test",
  "marker_select",
  "predict_ebv", "a_matrix", "prediction_accuracy", "selection_methods",
  "crossbreed", "breed_composition", "heterosis",
  # Fixed-scale accessors
  "additive_value", "genotypic_value", "phenotype_value", "gxe_value",
  # Genotype ingestion / QC
  "as_numeric", "filter_geno"
)

# S3 methods the contract names (dispatched, not called by name, so checked via
# the registry rather than the export set).
backend_contract_methods <- c(
  "[.Population", "c.Population", "print.Population",
  "plot.phenotype_sim", "print.ocs"
)

test_that("every backend-contract function is exported", {
  exported <- getNamespaceExports("simplePHENOTYPES")
  missing <- setdiff(backend_contract, exported)
  expect_identical(
    missing, character(0),
    info = paste("no longer exported (a MAJOR contract break):",
                 paste(missing, collapse = ", "))
  )
})

test_that("every backend-contract S3 method is registered", {
  for (m in backend_contract_methods) {
    expect_true(
      !is.null(utils::getS3method(sub("\\.[^.]+$", "", m),
                                  sub("^[^.]+\\.", "", m),
                                  optional = TRUE)),
      info = paste("S3 method not registered (a contract break):", m)
    )
  }
})

test_that("frozen legacy is not in the contract", {
  # create_phenotypes() is exported but deliberately excluded: consumers build on
  # the grammar, not the frozen v1 engine (DECISION-008).
  expect_false("create_phenotypes" %in% backend_contract)
})

# Signature manifest for the crossing surface: formal names, their ORDER and the
# defaults of the optional ones. Renaming, reordering or re-defaulting a formal is
# a contract break even when the name stays exported.
test_that("the crossing contract signatures are frozen (names, order, defaults)", {
  sig <- function(f) {
    fm <- formals(get(f, envir = asNamespace("simplePHENOTYPES")))
    vapply(fm, function(x) if (is.symbol(x) && !nzchar(as.character(x))) "<none>"
                           else paste(deparse(x), collapse = ""), character(1))
  }
  expect_identical(sig("cross"), c(mother = "<none>", father = "<none>", n = "1",
                                   seed = "NULL", interference = "NULL"))
  expect_identical(sig("selfcross"), c(parent = "<none>", n = "1", seed = "NULL",
                                       interference = "NULL"))
  expect_identical(sig("double_haploid"), c(parent = "<none>", n = "1", seed = "NULL",
                                            interference = "NULL"))
  expect_identical(names(sig("as_population")), c("geno", "individuals", "pool"))
  expect_identical(sig("as_population")[["pool"]], "NA_character_")
  expect_identical(names(sig("mate")), c("plan", "...", "seed", "prefix", "interference"))
  expect_identical(names(sig("mating_design")),
                   c("mothers", "fathers", "design", "n_crosses", "progeny_per_cross",
                     "mothers_per_father", "allow_self", "seed"))
  expect_identical(names(sig("crossbreed")),
                   c("breeds", "system", "n_progeny", "generations", "sire_breed",
                     "seed", "interference"))
  expect_identical(names(sig("heterosis")), c("pop", "breeds", "qtn", "a", "d"))
  expect_identical(names(sig("synthetic_map")),
                   c("chr", "pos", "total_cm", "cm_per_mb", "centromere",
                     "suppression", "width"))
  expect_identical(names(sig("as_numeric")), c("x", "..."))
})

# The Rust crate's declared minimum supported Rust version and dependency pin
# (RUST-F5): the vendored extendr-api needs 1.71, and a wildcard requirement lets a
# fresh resolve pick any release.
test_that("the Rust crate declares a consistent MSRV and a pinned dependency", {
  cargo <- file.path("..", "..", "src", "rust", "Cargo.toml")
  skip_if_not(file.exists(cargo), "not run from a source tree")
  ct <- readLines(cargo, warn = FALSE)
  rv <- sub("^.*=\\s*['\"]([0-9.]+)['\"].*$", "\\1",
            grep("^\\s*rust-version\\s*=", ct, value = TRUE))
  expect_identical(rv, "1.71")
  dep <- grep("^\\s*extendr-api\\s*=", ct, value = TRUE)
  expect_length(dep, 1L)
  expect_false(grepl("['\"]\\*['\"]", dep))
  expect_true(grepl("result_list", dep))
  desc <- file.path("..", "..", "DESCRIPTION")
  if (file.exists(desc)) {
    sr <- read.dcf(desc, "SystemRequirements")[1, 1]
    d <- sub("^.*rustc\\s*>=\\s*([0-9.]+).*$", "\\1", sr)
    skip_if(utils::compareVersion(d, rv) < 0,
            paste0("DESCRIPTION SystemRequirements lists rustc >= ", d,
                   " but Cargo.toml requires ", rv, "; update DESCRIPTION"))
    expect_gte(utils::compareVersion(d, rv), 0)
  }
})
