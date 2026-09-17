# Shared helpers for the transcriptome benchmark suite.
# Sourced by each 0X_*.R script. Keeps every script standalone-runnable:
#   Rscript benchmarks/0X_name.R
# from the package root. It loads the package with devtools::load_all(".") and
# provides an output directory plus a plotting wrapper that degrades gracefully
# when no graphics device is available.

# --- locate the package root (works whether run from root or benchmarks/) -----
.bench_find_root <- function() {
  # Prefer the directory that contains DESCRIPTION, walking up from this file or
  # the working directory.
  cand <- c(getwd(), file.path(getwd(), ".."))
  # If sourced, try to use the script's own location.
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa)) {
    sp <- normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)
    cand <- c(dirname(dirname(sp)), dirname(sp), cand)
  }
  for (d in cand) {
    if (file.exists(file.path(d, "DESCRIPTION"))) {
      return(normalizePath(d))
    }
  }
  normalizePath(getwd())
}

BENCH_ROOT <- .bench_find_root()
BENCH_OUT  <- file.path(BENCH_ROOT, "benchmarks", "output")
dir.create(BENCH_OUT, showWarnings = FALSE, recursive = TRUE)

# --- load the package (rextendr) ----------------------------------------------
message("Loading simplePHENOTYPES from ", BENCH_ROOT, " ...")
ok <- suppressWarnings(tryCatch({
  suppressMessages(devtools::load_all(BENCH_ROOT, quiet = TRUE))
  TRUE
}, error = function(e) {
  message("!! devtools::load_all failed: ", conditionMessage(e))
  FALSE
}))
if (!ok) {
  stop("Could not load the package with devtools::load_all(). ",
       "See the error above.", call. = FALSE)
}

# --- helpers ------------------------------------------------------------------

# Save a CSV under benchmarks/output and echo the path.
bench_write_csv <- function(df, name) {
  path <- file.path(BENCH_OUT, name)
  utils::write.csv(df, path, row.names = FALSE)
  message("  wrote ", path)
  invisible(path)
}

# Run plotting code into a PNG, but never let a missing/broken graphics device
# crash the benchmark: on failure we warn and continue.
bench_png <- function(name, expr, width = 900, height = 700, res = 110) {
  path <- file.path(BENCH_OUT, name)
  opened <- tryCatch({
    grDevices::png(path, width = width, height = height, res = res)
    TRUE
  }, error = function(e) {
    message("  (skipping plot ", name, ": ", conditionMessage(e), ")")
    FALSE
  })
  if (!opened) return(invisible(NULL))
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) {
                   message("  (plot ", name, " failed: ", conditionMessage(e), ")")
                   FALSE
                 })
  try(grDevices::dev.off(), silent = TRUE)
  if (ok) message("  wrote ", path)
  invisible(if (ok) path else NULL)
}

# Pull a dosage matrix (individuals x markers, coded -1/0/1) for a set of marker
# names straight from a numeric-format genotype data frame, without touching any
# internal package function.
bench_dosage <- function(geno_df, snp_names, snp_index = NULL) {
  if (is.null(snp_index)) {
    snp_index <- stats::setNames(seq_len(nrow(geno_df)), as.character(geno_df$snp))
  }
  rows <- snp_index[as.character(snp_names)]
  if (anyNA(rows)) {
    stop("bench_dosage(): marker(s) not found: ",
         paste(utils::head(snp_names[is.na(rows)], 5), collapse = ", "))
  }
  D <- t(as.matrix(geno_df[rows, -(1:5), drop = FALSE]))
  storage.mode(D) <- "double"
  colnames(D) <- as.character(snp_names)
  D
}

set.seed(1)   # global default; each script also seeds its own draws
options(width = 100)
