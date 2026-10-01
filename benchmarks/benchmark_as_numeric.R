# Dev-only script: microbenchmark is not a package dependency.
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("install.packages('microbenchmark') to run this benchmark", call. = FALSE)
}
if (!requireNamespace("withr", quietly = TRUE)) {
  stop("install.packages('withr') to run this benchmark", call. = FALSE)
}
library(microbenchmark)
library(simplePHENOTYPES)
library(data.table)

# create_phenotypes() writes output files to the working directory. The whole benchmark
# runs inside withr::with_dir() on a fresh temp dir, so no files are created in the
# project root and the working directory is restored even if a step fails. (A top-level
# on.exit() in a script is a no-op, which is what this replaces.)
tmp      <- tempfile("bench_as_numeric_")
dir.create(tmp)
vcf_file <- file.path(getwd(), "tests", "test.vcf")
stopifnot(file.exists(vcf_file))

run_benchmark <- function(vcf_file, tmp) {
  # ---------------------------------------------------------------------------
  # 1. Timing comparison
  # ---------------------------------------------------------------------------
  # as_numeric() uses impute = "Middle" to match the create_phenotypes() default.
  # Note: create_phenotypes() converts VCF with method = "biallelic.only" while
  # as_numeric() uses method = "copy.num.of.ref"; see parity section for impact.

  mbm <- microbenchmark(

    as_numeric = suppressMessages(
      as_numeric(vcf_file, to_r = TRUE, impute = "Middle", verbose = FALSE)
    ),

    create_phenotypes = suppressMessages(
      create_phenotypes(
        geno_file   = vcf_file,
        add_QTN_num = 3,
        add_effect  = 0.2,
        out_geno    = "numeric",
        rep         = 10,
        h2          = 0.7,
        model       = "A",
        quiet       = TRUE,
        home_dir    = tmp
      )
    ),

    times = 5L
  )

  print(mbm)

  # ---------------------------------------------------------------------------
  # 2. Parity check
  # ---------------------------------------------------------------------------
  # create_phenotypes() writes <home_dir>/simplePHENOTYPES_output/test_numeric.txt.
  # as_numeric() returns the same encoding in memory.

  num_file <- file.path(tmp, "simplePHENOTYPES_output", "test_numeric.txt")

  if (!file.exists(num_file)) {
    suppressMessages(
      create_phenotypes(
        geno_file   = vcf_file,
        add_QTN_num = 3,
        add_effect  = 0.2,
        out_geno    = "numeric",
        rep         = 1,
        h2          = 0.7,
        model       = "A",
        home_dir    = tmp
      )
    )
  }

  cp_geno <- data.table::fread(num_file, data.table = FALSE)
  an_geno <- suppressMessages(
    as_numeric(vcf_file, to_r = TRUE, impute = "Middle", verbose = FALSE)
  )

  common_snps  <- intersect(an_geno$snp, cp_geno[[1]])
  common_samps <- intersect(
    colnames(an_geno)[6:ncol(an_geno)],
    colnames(cp_geno)[6:ncol(cp_geno)]
  )

  an_mat <- as.matrix(an_geno[match(common_snps, an_geno$snp),  common_samps])
  cp_mat <- as.matrix(cp_geno[match(common_snps, cp_geno[[1]]), common_samps])

  n_cells <- length(an_mat)
  n_match <- sum(an_mat == cp_mat, na.rm = TRUE)
  n_diff  <- sum(an_mat != cp_mat, na.rm = TRUE)

  cat("\n--- Parity check ---\n")
  cat("  SNPs compared:           ", length(common_snps),  "\n")
  cat("  Samples compared:        ", length(common_samps), "\n")
  cat("  Matching cells:          ", n_match, "of", n_cells,
      sprintf("(%.2f%%)\n", 100 * n_match / n_cells))
  cat("  Differing cells:         ", n_diff, "\n")
  cat("  as_numeric range:        ", range(an_mat, na.rm = TRUE), "\n")
  cat("  create_phenotypes range: ", range(cp_mat, na.rm = TRUE), "\n")
}

withr::with_dir(tmp, run_benchmark(vcf_file, tmp))
unlink(tmp, recursive = TRUE)
