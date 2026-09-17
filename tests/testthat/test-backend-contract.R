# The backend contract (docs/BACKEND_CONTRACT.md): the exported engine surface
# breedingDesigner depends on. This guard fails if any contract function stops
# being exported -- a removed/renamed export is a MAJOR (3.0.0) break and must be
# a deliberate change here, caught before it reaches a downstream consumer.

# Keep in sync with docs/BACKEND_CONTRACT.md.
backend_contract <- c(
  # Populations & crossing
  "as_population", "cross", "selfcross", "double_haploid", "dosages",
  "n_individuals", "synthetic_map",
  # Phenotype grammar
  "simulate_phenotype", "additive", "dominance", "epistasis", "vqtl",
  "complex_phenotypes", "genetic_values", "qtn_table",
  "phenotypes_long", "phenotypes_wide", "write_phenotypes",
  # Selection engine
  "select_ind", "single_seed_descent", "bulk", "pedigree",
  "recurrent_selection",
  # Modern methods
  "g_matrix", "optimum_contribution", "sample_parents", "cross_usefulness",
  # Fixed-scale accessors
  "additive_value", "genotypic_value", "phenotype_value",
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
