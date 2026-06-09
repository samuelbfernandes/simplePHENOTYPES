# benchmark_as_numeric.R
#
# Compares three numericalization strategies on HapMap data:
#
#   OLD-file-loader:   apply(chars, 1, numericalization())
#                      — the per-SNP R loop that file_loader used before the
#                        restructure branch.
#
#   OLD-table-to-num:  table_to_numeric() → apply(xx, 1, make_numeric())
#                      — the per-SNP R loop that handle_hapmap / as_numeric
#                        used before the restructure branch.
#
#   NEW-kernel:        parse_hapmap_chars_to_raw() + numericalize_core()
#                      — the unified Rust kernel used by both as_numeric()
#                        and create_phenotypes() after the restructure.
#
# Each strategy is benchmarked at two levels:
#   (a) core numericalization only  (data already in memory, I/O excluded)
#   (b) full end-to-end call        (includes file read + numericalization)
#
# Requirements: microbenchmark, simplePHENOTYPES (installed from source)
#
# Usage:
#   Rscript benchmarks/benchmark_as_numeric.R
#   OR source("benchmarks/benchmark_as_numeric.R") from the package root.

library(microbenchmark)
library(simplePHENOTYPES)
library(data.table)

# ---------------------------------------------------------------------------
# Paths — adjust if running from outside the package root
# ---------------------------------------------------------------------------

PKG_ROOT  <- here::here()   # requires 'here' package; or set manually
HMP_SMALL <- file.path(PKG_ROOT, "tests", "test.hmp.txt")        # ~1.4 MB
HMP_LARGE <- file.path(PKG_ROOT, "tests",
                        "SNP55K_maize282_AGPv2_20100513_1.hmp.txt") # ~45 MB

# Fall back to relative paths if 'here' is not available
if (!requireNamespace("here", quietly = TRUE)) {
  HMP_SMALL <- "tests/test.hmp.txt"
  HMP_LARGE <- "tests/SNP55K_maize282_AGPv2_20100513_1.hmp.txt"
}

stopifnot(file.exists(HMP_SMALL))
has_large <- file.exists(HMP_LARGE)

TIMES <- 20L   # microbenchmark repetitions

# ---------------------------------------------------------------------------
# Pre-load data once (exclude I/O from core benchmarks)
# ---------------------------------------------------------------------------

G_small    <- data.table::fread(HMP_SMALL, data.table = FALSE)
chars_small <- as.matrix(G_small[, -(1:11)])
bit_small   <- nchar(as.character(G_small[2, 12]))  # 1 or 2 char genotypes

cat(sprintf("Small dataset: %d SNPs x %d samples\n",
            nrow(chars_small), ncol(chars_small)))

if (has_large) {
  G_large    <- data.table::fread(HMP_LARGE, data.table = FALSE)
  chars_large <- as.matrix(G_large[, -(1:11)])
  bit_large   <- nchar(as.character(G_large[2, 12]))
  cat(sprintf("Large dataset: %d SNPs x %d samples\n",
              nrow(chars_large), ncol(chars_large)))
}

# Pre-parse raw matrix for kernel-only benchmark (shared setup)
raw_small <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(chars_small)
flip_small <- simplePHENOTYPES:::compute_flip(raw_small)

if (has_large) {
  raw_large  <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(chars_large)
  flip_large <- simplePHENOTYPES:::compute_flip(raw_large)
}

# ---------------------------------------------------------------------------
# Helper: print a compact summary table
# ---------------------------------------------------------------------------

.print_summary <- function(mbm, dataset_label) {
  cat("\n==========================================================\n")
  cat(" Dataset:", dataset_label, "\n")
  cat("==========================================================\n")
  s <- summary(mbm)
  s$mean_s <- s$mean / 1e9
  cat(sprintf("  %-30s %8s %8s %8s\n", "Benchmark", "median(s)", "mean(s)", "neval"))
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-30s %8.3f %8.3f %8d\n",
                s$expr[i], s$median[i] / 1e9, s$mean_s[i], s$neval[i]))
  }
  # Speedup relative to slowest
  base_med <- max(s$median)
  cat("\n  Speedup vs slowest:\n")
  for (i in seq_len(nrow(s))) {
    cat(sprintf("  %-30s %.1fx\n", s$expr[i], base_med / s$median[i]))
  }
  invisible(mbm)
}

# ---------------------------------------------------------------------------
# Benchmark 1 — core numericalization (in-memory, no file I/O)
# ---------------------------------------------------------------------------

cat("\n--- Benchmark 1: core numericalization (in-memory) ---\n")

mbm_core_small <- microbenchmark(

  # Old create_phenotypes path: per-SNP apply using numericalization()
  `OLD file_loader (apply+numericalization)` = {
    apply(chars_small, 1, function(one)
      simplePHENOTYPES:::numericalization(
        one, bit = bit_small, effect = "Add", impute = "None"))
  },

  # Old as_numeric path: per-SNP apply using make_numeric() via table_to_numeric()
  `OLD as_numeric (table_to_numeric)` = {
    suppressMessages(
      simplePHENOTYPES:::table_to_numeric(
        G_small[, -(1:11)],
        code_as  = "-101",
        model    = "Add",
        impute   = "None",
        method   = "frequency",
        verbose  = FALSE
      )
    )
  },

  # New path: char parsing in R + Rust coding kernel
  `NEW parse_to_raw + numericalize_core` = {
    raw  <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(chars_small)
    flip <- simplePHENOTYPES:::compute_flip(raw)
    simplePHENOTYPES:::numericalize_core(
      raw_dosage = as.integer(raw),
      n_snp      = nrow(raw),
      n_samp     = ncol(raw),
      flip       = flip,
      code_as    = "-101",
      model      = "Add",
      impute     = "None"
    )
  },

  # New path: Rust coding only (char parsing already done)
  `NEW numericalize_core only (Rust)` = {
    simplePHENOTYPES:::numericalize_core(
      raw_dosage = as.integer(raw_small),
      n_snp      = nrow(raw_small),
      n_samp     = ncol(raw_small),
      flip       = flip_small,
      code_as    = "-101",
      model      = "Add",
      impute     = "None"
    )
  },

  times = TIMES
)

.print_summary(mbm_core_small, paste("small:", nrow(chars_small), "SNPs"))
print(mbm_core_small)

# ---------------------------------------------------------------------------
# Benchmark 2 — full end-to-end pipeline (includes file read + numericalization)
# ---------------------------------------------------------------------------

cat("\n--- Benchmark 2: end-to-end pipeline (file read included) ---\n")

mbm_e2e_small <- microbenchmark(

  # New as_numeric() public API
  `NEW as_numeric(file)` = {
    suppressMessages(
      as_numeric(HMP_SMALL, to_r = TRUE, verbose = FALSE)
    )
  },

  # New create_phenotypes numericalization path (genotypes() + file_loader())
  `NEW genotypes/file_loader (create_phen. path)` = {
    suppressMessages(
      simplePHENOTYPES:::genotypes(geno_file = HMP_SMALL,
                                   SNP_impute = "None",
                                   verbose    = FALSE)
    )
  },

  times = TIMES
)

.print_summary(mbm_e2e_small, paste("small:", nrow(chars_small), "SNPs"))
print(mbm_e2e_small)

# ---------------------------------------------------------------------------
# Benchmark 3 — large dataset (skip if file not present)
# ---------------------------------------------------------------------------

if (has_large) {
  cat("\n--- Benchmark 3: core numericalization, large dataset ---\n")

  mbm_core_large <- microbenchmark(

    `OLD file_loader (apply+numericalization)` = {
      apply(chars_large, 1, function(one)
        simplePHENOTYPES:::numericalization(
          one, bit = bit_large, effect = "Add", impute = "None"))
    },

    `OLD as_numeric (table_to_numeric)` = {
      suppressMessages(
        simplePHENOTYPES:::table_to_numeric(
          G_large[, -(1:11)],
          code_as = "-101", model = "Add",
          impute  = "None", method = "frequency", verbose = FALSE
        )
      )
    },

    `NEW parse_to_raw + numericalize_core` = {
      raw  <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(chars_large)
      flip <- simplePHENOTYPES:::compute_flip(raw)
      simplePHENOTYPES:::numericalize_core(
        raw_dosage = as.integer(raw),
        n_snp      = nrow(raw),
        n_samp     = ncol(raw),
        flip       = flip,
        code_as    = "-101", model = "Add", impute = "None"
      )
    },

    `NEW numericalize_core only (Rust)` = {
      simplePHENOTYPES:::numericalize_core(
        raw_dosage = as.integer(raw_large),
        n_snp      = nrow(raw_large),
        n_samp     = ncol(raw_large),
        flip       = flip_large,
        code_as    = "-101", model = "Add", impute = "None"
      )
    },

    times = 5L   # fewer reps; each iteration is slow on 51K SNPs
  )

  .print_summary(mbm_core_large, paste("large:", nrow(chars_large), "SNPs"))
  print(mbm_core_large)

  cat("\n--- Benchmark 4: end-to-end pipeline, large dataset ---\n")

  mbm_e2e_large <- microbenchmark(

    `NEW as_numeric(file)` = {
      suppressMessages(as_numeric(HMP_LARGE, to_r = TRUE, verbose = FALSE))
    },

    `NEW genotypes/file_loader (create_phen. path)` = {
      suppressMessages(
        simplePHENOTYPES:::genotypes(geno_file = HMP_LARGE,
                                     SNP_impute = "None",
                                     verbose    = FALSE)
      )
    },

    times = 5L
  )

  .print_summary(mbm_e2e_large, paste("large:", nrow(chars_large), "SNPs"))
  print(mbm_e2e_large)

} else {
  cat("\nSkipping large-dataset benchmarks: SNP55K file not found.\n")
  cat("Place tests/SNP55K_maize282_AGPv2_20100513_1.hmp.txt to enable them.\n")
}

cat("\nDone.\n")
