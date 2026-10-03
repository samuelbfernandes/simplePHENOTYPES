#' Realized phenotypes in long format
#'
#' The canonical output of the v2 grammar: one row per individual x trait x rep
#' with columns `id`, `trait`, `rep`, `value`.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in long format.
#' @seealso [phenotypes_wide()], [write_phenotypes()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' long <- phenotypes_long(ph)
#' head(long)
#'
#' # One row per individual x trait x rep, which is the shape most
#' # mixed-model and plotting tools expect.
#' table(long$trait)
phenotypes_long <- function(sim) {
  .check_sim(sim)
  .check_h2_complete(sim)
  sim$pheno
}

#' Realized phenotypes in wide format
#'
#' One row per individual x rep; one column per trait.
#'
#' @param sim a `phenotype_sim`.
#' @return a data frame in wide format.
#' @seealso [phenotypes_long()], [write_phenotypes()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' wide <- phenotypes_wide(ph)
#' head(wide)
#'
#' # One column per trait, which makes cross-trait comparisons direct.
#' round(cor(wide$Trait_1, wide$Trait_2), 2)
phenotypes_wide <- function(sim) {
  .check_sim(sim)
  .check_h2_complete(sim)
  long <- sim$pheno
  wide <- stats::reshape(
    long[, c("id", "rep", "trait", "value")],
    idvar = c("id", "rep"), timevar = "trait", direction = "wide"
  )
  names(wide) <- sub("^value\\.", "", names(wide))
  rownames(wide) <- NULL
  wide
}

#' Write realized phenotypes to disk
#'
#' Writes the long (default) or wide table as a delimited text file or as JSON.
#' Specialized exporters (gemma / plink / multi-file) are out of scope for the
#' grammar core.
#'
#' Delimited files are written with `data.table::fwrite()`'s default of up to 15
#' significant digits, so reading the file back reproduces the phenotypes to
#' about 1e-14 (relative), not bit-for-bit.
#'
#' JSON (`file_type = "json"`) is an array with one object per row of the chosen
#' layout, e.g. `[{"id": "33-16", "trait": "Trait_1", "rep": 1, "value": 0.52}, ...]`
#' for `"long"`, written as UTF-8 with 17 significant digits so every value reads
#' back exactly (e.g. `jsonlite::read_json(file, simplifyVector = TRUE)`, or
#' `pandas.read_json()` in Python). It needs the \pkg{jsonlite} package.
#'
#' @section Companion files (QTN table and split marker data):
#' The phenotype file alone does not say which markers produced it. Two options
#' write that information next to it, in the same `file_type` (and `sep`):
#'
#' * `qtn_file`: the QTN table of [qtn_table()], every column, through
#'   [write_qtn_table()].
#' * `split_markers = TRUE`: the marker data are split into two files, one with
#'   the **causal** markers and one with every **non-causal** marker, so a
#'   downstream analysis can be run on the non-causal markers alone (or on both)
#'   with the truth kept separately. The QTN table is written too (to `qtn_file`,
#'   or to a default name) because it is the companion of the causal file.
#'
#' Default names are derived from `file`: `<stem>_qtn_table.<ext>`,
#' `<stem>_qtn_markers.<ext>` and `<stem>_noncausal_markers.<ext>`, where
#' `stem` is `file` without its extension (`ext`; `txt` or `json` when `file`
#' has none). `markers_files = c(causal = ..., noncausal = ...)` overrides the
#' two marker paths. All paths must be distinct.
#'
#' **Causal set.** A marker is causal when any layer (`additive()`,
#' `dominance()`, `epistasis()` -- every member of an interacting set --,
#' `vqtl()`, and under `architecture = "ld"` the loci of both traits) of any
#' trait uses it in any of the replications selected by `rep`: the union of the
#' `snp` values of [qtn_table()] over those replications, excluding the gene
#' rows of `transcriptome()` layers (genes are not markers). With
#' `vary_qtn = TRUE` the replications differ in their QTNs, so pass
#' `rep = "all"` (or the vector of replications of interest) to make the causal
#' file cover every replication of the phenotype file; a marker causal in one
#' selected replication but not another is in the causal file, never in the
#' non-causal one. `rep` only affects the QTN table and the causal set -- the
#' phenotype file always holds every replication, as before.
#'
#' **Text layout of the marker files** (`file_type = "text"`): the package's
#' numeric format -- columns `snp`, `allele`, `chr`, `pos`, `cm` (then `counted`
#' when the input carried it), then one column per simulated individual with
#' the -1/0/1 dosage -- so [as_numeric()] and [simulate_phenotype()] read each
#' file back. [as_population()] reads a file too, provided it is not empty and
#' its map is complete (`chr`, `pos` and `cm` without missing values): a
#' matrix-origin simulation has no such map, and a file with no markers (no
#' causal marker, or every marker causal) cannot found a population. Missing
#' metadata (e.g. no allele label for a matrix input) is written as `NA`.
#' Markers keep their map order. The causal text file has no embedded QTN
#' table; the QTN table text file is its companion.
#'
#' **JSON layout of the marker files** (`file_type = "json"`): one object
#' `{"individuals": [...], "markers": [...]}`. `individuals` lists the
#' simulated individuals' ids; each element of `markers` is an object with
#' `snp`, `allele`, `chr`, `pos`, `cm`, `maf` (and `counted` when available) and
#' `genotypes`, an array with one dosage per individual in `individuals` order
#' (`null` for a missing value). The causal file adds, at the top level,
#' `qtn_table` (the rows of [write_qtn_table()] for the selected replications,
#' with a `rep` field when several), and in every marker object a `causal_for`
#' array of `{trait, layer, set}` objects (plus `rep` when several
#' replications) naming the layers that use the marker, so a reader can join
#' without the table. Gene rows of `transcriptome()` layers appear in
#' `qtn_table` but have no marker object. Numbers carry 17 significant digits,
#' so `jsonlite::read_json(file, simplifyVector = TRUE)` returns the
#' `individuals` vector, a `markers` data frame whose `genotypes` column is a
#' list of vectors, and (causal file) the `qtn_table` data frame; in Python,
#' `json.load()` then `pandas.DataFrame(obj["markers"])`.
#'
#' Marker data can be large (tens of thousands of markers by hundreds of
#' individuals); both layouts are written in chunks of markers. For a
#' data-frame or matrix input the chunks stream from the user's object, so no
#' second copy of the genotypes is made. A `Population` stores haplotypes, not
#' dosages, so its one markers-by-individuals dosage matrix is built once and
#' held in memory for the duration of the export.
#'
#' **Safety.** Before anything is written, every path is checked: its
#' directory must exist and be writable, it must not be a directory, and no
#' two paths may resolve to the same file (relative and absolute spellings,
#' symlinked directories and symlinked files -- live or dangling -- are
#' resolved; on macOS the comparison also ignores case and Unicode
#' normalisation form, on Windows case; two *hard links* to one file have
#' unrelated names and are not detected). A file name longer than 255 bytes --
#' including a default companion name derived from a long `file` -- is
#' rejected before anything is written. With companions, all files are first
#' written into a private staging directory (`.<token>.stage`, created
#' exclusively in each destination directory) under short fixed names that
#' keep the destination's extension, so a `.gz` destination is still
#' compressed by the text writer and no staged file can share a name with a
#' requested output or with anything else in the directory. The set is then
#' moved into place as a group: an existing file keeps its permission mode and
#' is first moved into the staging directory as its backup, and if any step of
#' that move fails the backups are put back. Should a backup itself fail to go
#' back, its staging directory is kept and the error names it and the
#' destination, so the previous content is never silently lost. A destination
#' that is a symbolic link is written through (the link stays, its target is
#' updated or created); a link loop or a chain of more than 40 links is an
#' error. This is a best-effort guarantee built on `file.rename()` within one
#' directory; it has been exercised on POSIX file systems, not on Windows, and
#' a staging directory that cannot be removed (its parent made read-only
#' meanwhile) is named in a warning rather than silently left. The plain
#' one-file call writes directly, as earlier versions did. JSON numbers use a
#' period as decimal mark whatever the session's `LC_NUMERIC`.
#'
#' @param sim a `phenotype_sim`.
#' @param file output path.
#' @param format "long" (default) or "wide".
#' @param sep field separator (default tab); text files only.
#' @param file_type `"text"` (default, a delimited file) or `"json"`.
#' @param qtn_file optional path; when given, the QTN table of [qtn_table()] is
#'   also written there with [write_qtn_table()], in the same `file_type` and
#'   `sep`. Default `NULL` (not written, unless `split_markers = TRUE`).
#' @param split_markers `TRUE` also writes the marker data as two files, causal
#'   and non-causal markers (see Details), plus the QTN table. Default `FALSE`.
#'   Needs genotypes: a phenotype built from expression alone has no markers.
#' @param markers_files optional named character vector
#'   `c(causal = path, noncausal = path)` overriding the default marker file
#'   names used by `split_markers = TRUE`.
#' @param rep replication(s) whose QTN table and causal set are written: a
#'   vector of replication numbers or `"all"`. Default `1L`. Only matters with
#'   `vary_qtn = TRUE`; see Details. Does not affect the phenotype file.
#' @return `file`, invisibly, when no companion file is written; otherwise a
#'   named character vector of every path written (`phenotypes`, and whichever
#'   of `qtn_table`, `causal`, `noncausal` apply), invisibly.
#' @seealso [phenotypes_long()], [phenotypes_wide()], [write_qtn_table()],
#'   [qtn_table()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#'
#' # Written to a temporary directory here; use your own path in practice.
#' out <- file.path(tempdir(), "phenotypes.txt")
#' write_phenotypes(ph, file = out)
#' head(read.delim(out))
#'
#' # Wide layout, comma separated.
#' out_wide <- file.path(tempdir(), "phenotypes_wide.csv")
#' write_phenotypes(ph, file = out_wide, format = "wide", sep = ",")
#' head(read.csv(out_wide))
#'
#' # JSON, one object per row (needs the jsonlite package).
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   out_json <- file.path(tempdir(), "phenotypes.json")
#'   write_phenotypes(ph, file = out_json, file_type = "json")
#'   head(jsonlite::read_json(out_json, simplifyVector = TRUE))
#'   unlink(out_json)
#' }
#'
#' # Phenotypes plus the QTN table that produced them.
#' out_qtn <- file.path(tempdir(), "phenotypes_qtn.txt")
#' write_phenotypes(ph, file = out, qtn_file = out_qtn)
#' head(read.delim(out_qtn))
#'
#' # Phenotypes, QTN table, and the marker data split into causal and
#' # non-causal markers (default names derived from `file`).
#' files <- write_phenotypes(ph, file = out, split_markers = TRUE)
#' files
#' causal <- read.delim(files[["causal"]])
#' causal$snp                              # the QTNs, in map order
#' nrow(read.delim(files[["noncausal"]]))  # every other marker
#'
#' unlink(c(out, out_wide, out_qtn, files))
write_phenotypes <- function(sim, file, format = c("long", "wide"),
                             sep = "\t", file_type = c("text", "json"),
                             qtn_file = NULL, split_markers = FALSE,
                             markers_files = NULL, rep = 1L) {
  .check_sim(sim)
  format <- match.arg(format)
  file_type <- match.arg(file_type)
  .validate_flag(split_markers, "split_markers")
  .check_path_arg(file, "file")
  tab <- if (format == "long") phenotypes_long(sim) else phenotypes_wide(sim)
  if (file_type == "json") {
    .require_jsonlite("write_phenotypes")
  }
  # Resolve and validate every companion path before anything is written, so
  # a bad request leaves no partial output behind.
  companions <- .companion_paths(sim, file, file_type, qtn_file, split_markers,
                                 markers_files)
  write_pheno <- function(path) {
    if (file_type == "json") {
      .write_json_rows(tab, path)
    } else {
      data.table::fwrite(tab, file = path, sep = sep)
    }
  }
  if (is.null(companions)) {
    # the one-file call of earlier versions: written in place, as before
    if (nchar(basename(file), type = "bytes") > 255L) {
      stop("write_phenotypes(): the file name of ", sQuote(file), " is longer ",
           "than 255 bytes, which the file system does not allow. Give a ",
           "shorter name.", call. = FALSE)
    }
    write_pheno(file)
    return(invisible(file))
  }
  reps <- .resolve_reps(sim, rep)
  if (split_markers && (is.null(sim$n_markers) || sim$n_markers < 1L)) {
    stop("write_phenotypes(): `split_markers = TRUE` needs genotypes, but ",
         "this phenotype was built from expression alone (no `geno`), so ",
         "there are no markers to split. Only the QTN table can be written ",
         "(`qtn_file`).", call. = FALSE)
  }
  paths <- c(phenotypes = file, qtn_table = companions$qtn_table)
  if (split_markers) {
    paths <- c(paths, causal = companions$causal,
               noncausal = companions$noncausal)
  }
  # The QTN table is assembled once and shared by the table file and the
  # causal marker file; a Population's dosage matrix is built once for the
  # whole export, QTN table included (see .export_sim()).
  sim <- .export_sim(sim)
  tab_q <- .qtn_table_reps(sim, reps)
  # Every file is written to a temporary sibling and moved into place only
  # after all of them succeeded, so a failing companion (unwritable path, disk
  # full, ...) never replaces an existing phenotype file with a partial set.
  .staged_write(paths, "write_phenotypes", function(tmp) {
    write_pheno(tmp[["phenotypes"]])
    .write_qtn_table_to(sim, reps, tmp[["qtn_table"]], file_type, sep,
                        tab = tab_q)
    if (split_markers) {
      .write_split_markers(sim, .causal_markers(sim, reps, tab = tab_q),
                           tmp[["causal"]], tmp[["noncausal"]], file_type, sep)
    }
  })
  invisible(paths)
}

#' Write the QTN table of a simulation to disk
#'
#' Writes every column of [qtn_table()] -- the v2 equivalent of the v1
#' `Additive_QTNs.txt` / `QTN_effects_summary.txt` files -- as a delimited text
#' file or as JSON. One replication is written by default; several (or
#' `rep = "all"`) are stacked with a leading `rep` column, which is how the
#' per-replication architectures of a `vary_qtn = TRUE` simulation are kept
#' apart. Gene rows of `transcriptome()` layers are included (their `snp` is
#' the gene identifier; see [qtn_table()]).
#'
#' Text files are written with `data.table::fwrite()` (up to 15 significant
#' digits, so effects read back to about 1e-14 relative; missing values are
#' empty fields). JSON follows the conventions of [write_phenotypes()]: an
#' array with one object per row, UTF-8, 17 significant digits (every value
#' reads back exactly), `null` for a missing value; it needs the \pkg{jsonlite}
#' package. Read it back with `read.delim(file)` or
#' `jsonlite::read_json(file, simplifyVector = TRUE)`.
#'
#' @param sim a `phenotype_sim`.
#' @param file output path.
#' @param rep replication(s) whose QTN architecture to write: a vector of
#'   replication numbers, or `"all"`. Default `1L`. With more than one, the
#'   rows are stacked and a leading `rep` column says which replication each
#'   row belongs to.
#' @param file_type `"text"` (default, a delimited file) or `"json"`.
#' @param sep field separator (default tab); text files only.
#' @return `file`, invisibly.
#' @seealso [qtn_table()] for the columns, [write_phenotypes()] which can call
#'   this and also split the marker data into causal and non-causal files.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, seed = 1) |>
#'   additive(prop = 0.3, n_qtn = 3) |>
#'   epistasis(prop = 0.2, n_pairs = 2)
#'
#' # Written to a temporary directory here; use your own path in practice.
#' out <- file.path(tempdir(), "qtn_table.txt")
#' write_qtn_table(ph, file = out)
#' read.delim(out)
#'
#' # Every replication of a vary_qtn simulation, stacked with a `rep` column.
#' ph3 <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, n_reps = 3,
#'                           vary_qtn = TRUE, seed = 1) |>
#'   additive(n_qtn = 2)
#' out3 <- file.path(tempdir(), "qtn_table_all_reps.csv")
#' write_qtn_table(ph3, file = out3, rep = "all", sep = ",")
#' read.csv(out3)
#'
#' # JSON (needs the jsonlite package).
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   out_json <- file.path(tempdir(), "qtn_table.json")
#'   write_qtn_table(ph, file = out_json, file_type = "json")
#'   jsonlite::read_json(out_json, simplifyVector = TRUE)
#'   unlink(out_json)
#' }
#'
#' unlink(c(out, out3))
write_qtn_table <- function(sim, file, rep = 1L, file_type = c("text", "json"),
                            sep = "\t") {
  .check_sim(sim)
  .check_h2_complete(sim)
  file_type <- match.arg(file_type)
  .check_path_arg(file, "file")
  reps <- .resolve_reps(sim, rep)
  if (file_type == "json") .require_jsonlite("write_qtn_table")
  sim <- .export_sim(sim)        # a Population's dosages once, not per layer
  .staged_write(c(qtn_table = file), "write_qtn_table", function(tmp) {
    .write_qtn_table_to(sim, reps, tmp[["qtn_table"]], file_type, sep)
  })
  invisible(file)
}

#' Write the stacked QTN table of `reps` to `path` (text or JSON)
#' @keywords internal
#' @noRd
.write_qtn_table_to <- function(sim, reps, path, file_type, sep,
                                tab = .qtn_table_reps(sim, reps)) {
  if (file_type == "json") {
    .write_json_rows(tab, path)
  } else {
    data.table::fwrite(tab, file = path, sep = sep)
  }
  invisible(path)
}

# ---------------------------------------------------------------------------
# export helpers: replication sets, JSON rows, companion paths
# ---------------------------------------------------------------------------

#' Stop unless jsonlite is installed (it is in Suggests)
#' @keywords internal
#' @noRd
.require_jsonlite <- function(fn) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(fn, "(): `file_type = \"json\"` needs the jsonlite ",
         "package; install it with install.packages(\"jsonlite\").",
         call. = FALSE)
  }
}

#' A single, non-empty output path
#' @keywords internal
#' @noRd
.check_path_arg <- function(x, arg) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", arg, "` must be one non-empty file path.", call. = FALSE)
  }
}

#' Resolve a `rep` export argument to a vector of replication indices
#'
#' Accepts `"all"` or a vector of whole numbers in `1..n_reps`; duplicates are
#' dropped, the order given is kept (it is the stacking order of the table).
#' @keywords internal
#' @noRd
.resolve_reps <- function(sim, rep) {
  if (is.character(rep) && length(rep) == 1L && identical(rep, "all")) {
    return(seq_len(sim$n_reps))
  }
  if (!is.numeric(rep) || length(rep) == 0L || any(!is.finite(rep)) ||
      any(rep != floor(rep)) || any(rep < 1L)) {
    stop("`rep` must be \"all\" or a vector of positive whole numbers; got ",
         paste(utils::head(rep, 5), collapse = ", "), ".", call. = FALSE)
  }
  if (any(rep > sim$n_reps)) {
    stop("`rep` must be between 1 and n_reps (", sim$n_reps, "); got ",
         paste(rep[rep > sim$n_reps], collapse = ", "), ".", call. = FALSE)
  }
  unique(as.integer(rep))
}

#' The QTN table for one or several replications
#'
#' One replication returns exactly `qtn_table(sim, rep)`. Several are stacked
#' in the order given with a leading integer `rep` column (also when every
#' replication happens to share the same loci, so the layout only depends on
#' how many replications were asked for).
#' @keywords internal
#' @noRd
.qtn_table_reps <- function(sim, reps) {
  if (length(reps) == 1L) {
    return(qtn_table(sim, reps))
  }
  parts <- lapply(reps, function(r) {
    tab <- qtn_table(sim, r)
    cbind(rep = rep_len(as.integer(r), nrow(tab)), tab,
          stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Write a data frame as a JSON array of row objects
#'
#' The one JSON row convention of the package: UTF-8, 17 significant digits
#' (every double reads back exactly), `null` for a missing value of any type,
#' scalars unboxed. Numbers are formatted in the C numeric locale.
#' @keywords internal
#' @noRd
.write_json_rows <- function(tab, file) {
  .with_c_numeric(
    jsonlite::write_json(tab, path = file, dataframe = "rows",
                         digits = I(17), na = "null", auto_unbox = TRUE))
}

#' Evaluate `expr` with `LC_NUMERIC = "C"`
#'
#' JSON needs a period as the decimal mark whatever the session locale; both
#' `sprintf()` and jsonlite's number formatting follow `LC_NUMERIC`, so every
#' JSON writer runs inside this. A failing `Sys.setlocale()` is tolerated (the
#' expression then runs in the current locale, which is "C" in practice since R
#' itself warns against changing LC_NUMERIC). The previous value is restored on
#' exit.
#' @keywords internal
#' @noRd
.with_c_numeric <- function(expr) {
  current <- tryCatch(Sys.getlocale("LC_NUMERIC"), error = function(e) "")
  if (nzchar(current) && !identical(current, "C")) {
    set <- tryCatch(suppressWarnings(Sys.setlocale("LC_NUMERIC", "C")),
                    error = function(e) "")
    if (nzchar(set)) {
      on.exit(tryCatch(suppressWarnings(Sys.setlocale("LC_NUMERIC", current)),
                       error = function(e) NULL), add = TRUE)
    }
  }
  expr
}

#' `<stem>_<suffix>.<ext>` next to `file`
#'
#' `stem` is `file` without its extension; a `file` without one gets `.txt` or
#' `.json` according to `file_type`.
#' @keywords internal
#' @noRd
.sibling_path <- function(file, suffix, file_type) {
  base <- basename(file)
  has_ext <- grepl("^.+\\.[^.]+$", base)
  ext <- if (has_ext) sub("^.*\\.([^.]+)$", "\\1", base) else
    if (file_type == "json") "json" else "txt"
  stem <- if (has_ext) sub("\\.[^.]+$", "", file) else file
  paste0(stem, "_", suffix, ".", ext)
}

#' Resolve the companion paths of write_phenotypes(), or NULL when none
#'
#' Returns a list with `qtn_table`, `causal` and `noncausal` entries (the latter
#' two only with `split_markers`), after checking that every path, `file`
#' included, is distinct so no output overwrites another.
#' @keywords internal
#' @noRd
.companion_paths <- function(sim, file, file_type, qtn_file, split_markers,
                             markers_files) {
  if (!is.null(qtn_file)) .check_path_arg(qtn_file, "qtn_file")
  if (!split_markers && !is.null(markers_files)) {
    stop("`markers_files` is only used with `split_markers = TRUE`.",
         call. = FALSE)
  }
  if (is.null(qtn_file) && !split_markers) {
    return(NULL)
  }
  out <- list()
  out$qtn_table <- if (!is.null(qtn_file)) qtn_file else
    .sibling_path(file, "qtn_table", file_type)
  if (split_markers) {
    if (is.null(markers_files)) {
      out$causal    <- .sibling_path(file, "qtn_markers", file_type)
      out$noncausal <- .sibling_path(file, "noncausal_markers", file_type)
    } else {
      if (!is.character(markers_files) || length(markers_files) != 2L ||
          is.null(names(markers_files)) ||
          !setequal(names(markers_files), c("causal", "noncausal")) ||
          anyNA(markers_files) || any(!nzchar(markers_files))) {
        stop("`markers_files` must be a named character vector ",
             "c(causal = path, noncausal = path).", call. = FALSE)
      }
      out$causal    <- unname(markers_files[["causal"]])
      out$noncausal <- unname(markers_files[["noncausal"]])
    }
  }
  paths <- c(file, unlist(out, use.names = FALSE))
  key <- vapply(paths, .canonical_path, character(1), USE.NAMES = FALSE)
  if (anyDuplicated(key)) {
    dup <- duplicated(key) | duplicated(key, fromLast = TRUE)
    stop("write_phenotypes(): output paths must be distinct, but ",
         paste(sQuote(unique(paths[dup])), collapse = ", "),
         " resolve to the same file and would be written more than once.",
         call. = FALSE)
  }
  out
}

#' A comparison key for an output path that may not exist yet
#'
#' `normalizePath()` leaves a non-existent leaf alone, so a relative and an
#' absolute spelling, or a path below a directory and the same path below a
#' symlink to it, compare unequal although they name one file. The path is
#' first taken to the object it will be written to (`.write_target()`: through
#' a leaf symlink, live or dangling); then an existing leaf is resolved in
#' full and a missing one gets its longest existing ancestor normalised
#' (symlinks resolved, made absolute) with the missing tail appended. On macOS
#' the key is Unicode-normalised (APFS and HFS+ treat the NFC and NFD
#' spellings of a name as one entry) and, like Windows, case-folded. Two hard
#' links to one file have different paths and are not recognised. Used for
#' comparison only, never for writing.
#' @keywords internal
#' @noRd
.canonical_path <- function(p) {
  p <- .write_target(p)
  if (file.exists(p) && !dir.exists(p)) {
    key <- normalizePath(p, winslash = "/", mustWork = FALSE)
  } else {
    dir <- dirname(p)
    tail <- basename(p)
    while (!dir.exists(dir) && !identical(dir, dirname(dir))) {
      tail <- file.path(basename(dir), tail)
      dir <- dirname(dir)
    }
    key <- file.path(normalizePath(dir, winslash = "/", mustWork = FALSE), tail)
  }
  darwin <- identical(Sys.info()[["sysname"]], "Darwin")
  if (darwin) key <- .unicode_nfd(key)
  if (darwin || .Platform$OS.type == "windows") key <- tolower(key)
  key
}

#' Decomposed (NFD) Unicode form of a path, as macOS stores file names
#'
#' `iconv()` to `"UTF-8-MAC"` is the system's own decomposition; stringi, when
#' installed, is the fallback; otherwise the UTF-8 string is returned as is.
#' @keywords internal
#' @noRd
.unicode_nfd <- function(x) {
  out <- tryCatch(iconv(enc2utf8(x), "UTF-8", "UTF-8-MAC"),
                  error = function(e) NA_character_)
  if (length(out) == 1L && !is.na(out)) {
    return(out)
  }
  if (requireNamespace("stringi", quietly = TRUE)) {
    return(stringi::stri_trans_nfd(enc2utf8(x)))
  }
  enc2utf8(x)
}

#' The file system object a destination path denotes
#'
#' A destination that is a symbolic link -- live or dangling -- is written
#' *through*: the link is followed (a relative link target is taken relative
#' to the link's directory), the staged file is created beside the final
#' target and renamed onto it, so the link survives and points at the updated
#' (or newly created) file. A link that leads back to itself, or a chain of
#' more than 40 links, is an error, raised before anything is written.
#' Anything that is not a link is written at the path given.
#' @keywords internal
#' @noRd
.write_target <- function(p) {
  p <- path.expand(p)
  if (.Platform$OS.type == "windows") {
    return(p)
  }
  given <- p
  seen <- character(0)
  repeat {
    link <- Sys.readlink(p)
    if (is.na(link) || !nzchar(link)) {
      return(p)
    }
    if (length(seen) >= 40L) {
      stop("The output path ", sQuote(given), " goes through more than 40 ",
           "levels of symbolic links.", call. = FALSE)
    }
    key <- file.path(normalizePath(dirname(p), winslash = "/", mustWork = FALSE),
                     basename(p))
    if (key %in% seen) {
      stop("The output path ", sQuote(given), " is a symbolic link loop.",
           call. = FALSE)
    }
    seen <- c(seen, key)
    if (!grepl("^/", link)) link <- file.path(dirname(p), link)
    p <- link
  }
}

#' Check that every destination can be written before anything is written
#'
#' For each destination (the path as given, and the object it resolves to for
#' writing) the parent directory must exist and be writable, the destination
#' must not be a directory, an existing destination must be writable, and the
#' file name must fit the 255-byte component limit (a default companion name
#' is derived from `file`, so a long `file` can push it over: the message
#' points to `qtn_file =` / `markers_files =`). Errors name the path as the
#' caller gave it.
#' @keywords internal
#' @noRd
.preflight_paths <- function(paths, targets, fn) {
  for (i in seq_along(paths)) {
    p <- paths[[i]]
    t <- targets[[i]]
    if (nchar(basename(p), type = "bytes") > 255L ||
        nchar(basename(t), type = "bytes") > 255L) {
      stop(fn, "(): the file name of ", sQuote(p), " is longer than 255 ",
           "bytes, which the file system does not allow. ",
           if (identical(fn, "write_phenotypes"))
             paste0("Default companion names add a suffix to the name of ",
                    "`file`; give shorter names with `qtn_file =` and ",
                    "`markers_files =`.")
           else "Give a shorter name.",
           call. = FALSE)
    }
    parent <- dirname(t)
    if (!dir.exists(parent)) {
      stop(fn, "(): the directory of ", sQuote(p), " does not exist; create ",
           "it first.", call. = FALSE)
    }
    if (dir.exists(t)) {
      stop(fn, "(): ", sQuote(p), " is a directory, not a file path.",
           call. = FALSE)
    }
    if (file.access(parent, mode = 2L) != 0L ||
        (file.exists(t) && file.access(t, mode = 2L) != 0L)) {
      stop(fn, "(): ", sQuote(p), " is not writable.", call. = FALSE)
    }
  }
  invisible(paths)
}

#' `file.rename()` behind one name, so tests can make a rename fail
#' @keywords internal
#' @noRd
.file_rename <- function(from, to) file.rename(from, to)

#' A random token for a stage directory name (no R RNG: `tempfile()` draws
#' its own)
#' @keywords internal
#' @noRd
.random_token <- function() {
  sub("^r", "", basename(tempfile(pattern = "r", tmpdir = "")))
}

#' The extension of a destination that a writer reads its behaviour from
#'
#' `data.table::fwrite()` compresses by the trailing `.gz` / `.bz2` / `.xz` /
#' `.zip`; the extension before such a suffix is kept too (`.txt.gz`), so the
#' staged name tells the writer the same as the destination. The user's path
#' decides, not a symlink target's name.
#' @keywords internal
#' @noRd
.staging_ext <- function(path) {
  base <- basename(path)
  m <- regmatches(base, regexpr("(\\.[A-Za-z0-9]{1,10})?\\.(gz|bz2|xz|zip)$",
                                base, ignore.case = TRUE))
  if (!length(m)) {
    m <- regmatches(base, regexpr("\\.[A-Za-z0-9]{1,10}$", base))
  }
  if (length(m)) m else ""
}

#' Create a private stage directory `.<token>.stage` inside `dir`
#'
#' `dir.create()` is atomic and fails when the name exists, so the directory
#' is exclusively this export's: nothing inside it can be a file someone else
#' created or a requested destination, which is what makes the fixed names
#' used inside (`1<ext>`, `b1`, ...) safe. Up to 20 tokens are tried.
#' @keywords internal
#' @noRd
.stage_dir <- function(dir, fn, forbidden = character(0)) {
  for (i in seq_len(20L)) {
    stage <- file.path(dir, paste0(".", .random_token(), ".stage"))
    if (stage %in% forbidden) next          # a requested output has this name
    if (file.exists(stage) || !is.na(Sys.readlink(stage)) &&
        nzchar(Sys.readlink(stage))) next
    # owner-only (0700): the staged content is private until it is committed
    if (isTRUE(suppressWarnings(dir.create(stage, showWarnings = FALSE,
                                           mode = "0700")))) {
      Sys.chmod(stage, "0700", use_umask = FALSE)
      return(stage)
    }
  }
  stop(fn, "(): could not create a staging directory in ", sQuote(dir),
       " after 20 attempts.", call. = FALSE)
}

#' Remove stage directories (and everything in them), warning about leftovers
#' @keywords internal
#' @noRd
.remove_stage_dirs <- function(dirs, fn) {
  dirs <- dirs[!is.na(dirs) & dir.exists(dirs)]
  if (length(dirs)) {
    unlink(dirs, recursive = TRUE)
    left <- dirs[dir.exists(dirs)]
    if (length(left)) {
      warning(fn, "(): could not remove the staging directory(ies) ",
              paste(sQuote(left), collapse = ", "),
              "; remove them by hand.", call. = FALSE)
    }
  }
  invisible(dirs)
}

#' Write a set of files as a group, all or none
#'
#' `paths` is a named character vector of destinations. After the preflight,
#' one private stage directory is created in each destination directory (the
#' directory of the link target when the destination is a symlink, so no
#' rename crosses a file system). `writer(tmp)` receives a same-named vector
#' of files inside those directories, named `1<ext>`, `2<ext>`, ... with the
#' extension of the user's path (so the writer treats them as it would the
#' destination), and writes every file there. If the writer fails, the stage
#' directories are removed and every destination is untouched. Otherwise
#' `.commit_staged()` moves the set into place. A stage directory that cannot
#' be removed (its parent made read-only in the meantime) is named in a
#' warning.
#' @keywords internal
#' @noRd
.staged_write <- function(paths, fn, writer) {
  targets <- vapply(paths, .write_target, character(1))
  names(targets) <- names(paths)
  .preflight_paths(paths, targets, fn)
  dirs <- unique(dirname(targets))
  stage <- character(0)
  written <- FALSE
  on.exit(if (!written) .remove_stage_dirs(stage, fn), add = TRUE)
  for (d in dirs) {
    stage[[d]] <- .stage_dir(d, fn, forbidden = unname(targets))
  }
  tmp <- vapply(seq_along(paths), function(i) {
    file.path(stage[[dirname(targets[[i]])]],
              paste0(i, .staging_ext(paths[[i]])))
  }, character(1))
  names(tmp) <- names(paths)
  writer(tmp)
  written <- TRUE
  .commit_staged(tmp, targets, stage, fn)
  invisible(paths)
}

#' Move finished staged files onto their targets, all or none
#'
#' An existing target keeps its permission mode and is first moved into its
#' stage directory as `b<i>`; the staged files are then renamed onto the
#' targets. If any step fails, every placed file is replaced by its backup
#' again (an atomic rename, or unlink then rename), and the return of every
#' unlink and rename is checked: a backup that could not be put back is
#' **kept** -- its stage directory stays, with the staged temporaries removed
#' -- and the error names it together with every destination left in a mixed
#' state, so the previous content is always still on disk somewhere named in
#' the message. On success the stage directories are removed.
#' @keywords internal
#' @noRd
.commit_staged <- function(tmp, targets, stage, fn) {
  nms <- names(targets)
  existed <- file.exists(targets)
  backup <- stats::setNames(rep(NA_character_, length(nms)), nms)
  placed <- stats::setNames(logical(length(nms)), nms)
  # a replaced file keeps its permission mode (a fresh temporary has the
  # umask default, which could widen a private file)
  for (nm in nms[existed]) {
    tryCatch(Sys.chmod(tmp[[nm]], file.mode(targets[[nm]]), use_umask = FALSE),
             error = function(e) NULL, warning = function(w) NULL)
  }
  rollback <- function() {
    unrecovered <- character(0)
    mixed <- character(0)
    for (nm in rev(nms)) {
      b <- backup[[nm]]
      has_backup <- !is.na(b) && file.exists(b)
      if (has_backup) {
        # an atomic replace first; if the file system refuses to rename over
        # the placed file, remove it and rename again
        ok <- isTRUE(.file_rename(b, targets[[nm]]))
        if (!ok && placed[[nm]]) {
          removed <- unlink(targets[[nm]]) == 0L && !file.exists(targets[[nm]])
          ok <- removed && isTRUE(.file_rename(b, targets[[nm]]))
        }
        if (!ok) {
          unrecovered <- c(unrecovered, b)
          mixed <- c(mixed, targets[[nm]])
        }
      } else if (placed[[nm]]) {
        # no previous file: the destination must be absent again
        if (!(unlink(targets[[nm]]) == 0L && !file.exists(targets[[nm]]))) {
          mixed <- c(mixed, targets[[nm]])
        }
      }
    }
    list(unrecovered = unrecovered, mixed = mixed)
  }
  tryCatch({
    for (i in seq_along(nms)[existed]) {
      nm <- nms[[i]]
      b <- file.path(dirname(tmp[[nm]]), paste0("b", i))
      if (!isTRUE(.file_rename(targets[[nm]], b))) {
        stop(fn, "(): could not set aside the existing file ",
             sQuote(targets[[nm]]), " before replacing it.", call. = FALSE)
      }
      backup[[nm]] <- b
    }
    for (nm in nms) {
      if (!isTRUE(.file_rename(tmp[[nm]], targets[[nm]]))) {
        stop(fn, "(): could not move the finished file into place at ",
             sQuote(targets[[nm]]), ".", call. = FALSE)
      }
      placed[[nm]] <- TRUE
    }
  }, error = function(e) {
    rb <- rollback()
    # the staged temporaries go; a stage directory holding a backup that
    # could not be put back is kept, every other one is removed
    unlink(tmp[file.exists(tmp)])
    keep <- unique(dirname(rb$unrecovered))
    .remove_stage_dirs(setdiff(stage, keep), fn)
    msg <- conditionMessage(e)
    if (length(rb$unrecovered)) {
      msg <- paste0(
        msg, "\nThe previous content of ", paste(sQuote(rb$mixed), collapse = ", "),
        " could not be put back and is kept in ",
        paste(sQuote(rb$unrecovered), collapse = ", "),
        " (staging directory ", paste(sQuote(keep), collapse = ", "),
        "); these destinations are in a mixed state -- restore them by hand ",
        "from the backup file(s).")
    } else if (length(rb$mixed)) {
      msg <- paste0(msg, "\nThe destination(s) ",
                    paste(sQuote(rb$mixed), collapse = ", "),
                    " could not be returned to their previous state.")
    }
    stop(msg, call. = FALSE)
  })
  .remove_stage_dirs(stage, fn)
  invisible(targets)
}

# ---------------------------------------------------------------------------
# export helpers: marker metadata, causal set, marker files
# ---------------------------------------------------------------------------

#' Marker metadata of a simulation in numeric-format column order
#'
#' `snp`, `allele`, `chr`, `pos`, `cm` (and `counted` when the genotype object
#' carries it) for every marker, in map order. The foundation's `map` keeps only
#' `snp`/`chr`/`pos`; `allele` and `cm` come from the genotype object itself: a
#' numeric-format data frame has them as columns, a Population keeps them on its
#' map, and a plain matrix has neither (written as `NA`).
#' @keywords internal
#' @noRd
.marker_meta <- function(sim) {
  n <- sim$n_markers
  meta <- data.frame(
    snp = sim$map$snp, allele = rep(NA_character_, n), chr = sim$map$chr,
    pos = sim$map$pos, cm = rep(NA_real_, n), stringsAsFactors = FALSE
  )
  g <- sim$geno
  if (identical(sim$kind, "data.frame")) {
    meta$allele <- as.character(g$allele)
    meta$cm     <- as.numeric(g$cm)
  } else if (identical(sim$kind, "population")) {
    if (!is.null(g$map$allele)) meta$allele <- as.character(g$map$allele)
    meta$cm <- as.numeric(g$map$cm)
    if (!is.null(g$map$counted)) meta$counted <- as.character(g$map$counted)
  }
  meta
}

#' The causal markers of a simulation over a set of replications
#'
#' The union of the `snp` values of `qtn_table()` over `reps`, excluding the
#' gene rows of transcriptome layers (genes are not markers), mapped back to
#' marker indices and sorted into map order. Also returns the stacked table and,
#' per causal marker, the `{trait, layer, set[, rep]}` rows that use it.
#' @keywords internal
#' @noRd
.causal_markers <- function(sim, reps, tab = .qtn_table_reps(sim, reps)) {
  marker_rows <- tab[tab$layer != "transcriptome", , drop = FALSE]
  idx <- sort(unique(match(marker_rows$snp, sim$map$snp)))
  if (anyNA(idx)) {
    stop("Internal error: a QTN of the table is not in the marker map.",
         call. = FALSE)   # nocov
  }
  keep <- intersect(c("rep", "trait", "layer", "set"), names(marker_rows))
  causal_for <- split(marker_rows[, keep, drop = FALSE], marker_rows$snp)
  list(table = tab, idx = idx, causal_for = causal_for)
}

#' Write the causal and non-causal marker files of write_phenotypes()
#' @keywords internal
#' @noRd
.write_split_markers <- function(sim, cm, causal_file, noncausal_file,
                                 file_type, sep) {
  noncausal <- setdiff(seq_len(sim$n_markers), cm$idx)
  .write_marker_file(sim, cm$idx, causal_file, file_type, sep,
                     qtn_table = cm$table, causal_for = cm$causal_for)
  .write_marker_file(sim, noncausal, noncausal_file, file_type, sep)
  invisible(c(causal = causal_file, noncausal = noncausal_file))
}

#' A simulation prepared for a whole-panel export
#'
#' A data-frame or matrix backing streams through `.geno_cols()`, which reads
#' the requested columns from the user's object without copying the rest. A
#' Population stores haplotypes, and its `dosages()` rebuilds the whole
#' markers-by-individuals matrix on every call, so a chunked export through
#' `.geno_cols()` would build it once per chunk (and once more per QTN-table
#' call). For that backing the matrix is built once here, restricted to the
#' simulated individuals, and attached to this local copy of `sim` as
#' `export_dosages`; `.dosage_block()` then serves every dosage request of the
#' export (chunks and `.qtn_var()`) from it. Installed by every export that
#' computes a QTN table (`write_qtn_table()`, `write_phenotypes()` with
#' `qtn_file` or `split_markers`), since `qtn_table()` reaches the dosages once
#' per mean-effect layer and trait. The copy never leaves the export.
#' @keywords internal
#' @noRd
.export_sim <- function(sim) {
  if (identical(sim$kind, "population")) {
    D <- dosages(sim$geno)                       # markers x individuals
    storage.mode(D) <- "double"
    if (!is.null(sim$ind_idx)) D <- D[, sim$ind_idx, drop = FALSE]
    sim$export_dosages <- D
  }
  sim
}

#' Individuals-by-markers dosage block, from the export matrix when present
#' @keywords internal
#' @noRd
.dosage_block <- function(sim, idx) {
  if (is.null(sim$export_dosages)) {
    return(.geno_cols(sim, idx))
  }
  idx <- as.integer(idx)
  out <- t(sim$export_dosages[idx, , drop = FALSE])
  dimnames(out) <- list(sim$ids, sim$map$snp[idx])
  out
}

#' Write a set of markers (by map index) as a numeric-format text file or JSON
#'
#' Markers are fetched through `.dosage_block()` in chunks, so peak memory is
#' one chunk of individuals-by-markers dosages plus its text (plus, for a
#' Population, its one dosage matrix; see `.export_sim()`). `qtn_table` and
#' `causal_for` (JSON only) add the top-level `qtn_table` array and the
#' per-marker `causal_for` arrays of the causal file. JSON numbers are
#' formatted in the C numeric locale.
#' @keywords internal
#' @noRd
.write_marker_file <- function(sim, idx, file, file_type, sep,
                               qtn_table = NULL, causal_for = NULL,
                               chunk = 2000L) {
  idx <- as.integer(idx)
  meta <- .marker_meta(sim)
  chunks <- if (length(idx)) split(idx, ceiling(seq_along(idx) / chunk)) else
    list()
  if (file_type == "json") {
    con <- file(file, open = "w", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    .with_c_numeric({
      cat('{"individuals":', .json_vec(sim$ids), ',', file = con, sep = "")
      if (!is.null(qtn_table)) {
        cat('"qtn_table":',
            as.character(jsonlite::toJSON(qtn_table, dataframe = "rows",
                                          digits = I(17), na = "null",
                                          auto_unbox = TRUE)),
            ',', file = con, sep = "")
      }
      cat('"markers":[', file = con, sep = "")
      first <- TRUE
      for (ch in chunks) {
        objs <- .marker_json_objects(sim, meta, ch, causal_for)
        cat(if (first) "" else ",", paste(objs, collapse = ","),
            file = con, sep = "")
        first <- FALSE
      }
      cat(']}\n', file = con, sep = "")
    })
  } else {
    if (!length(chunks)) {
      # no marker in this set: header only, so the file is still well formed
      empty <- cbind(meta[0, , drop = FALSE],
                     as.data.frame(matrix(integer(0), 0, sim$n_ind,
                                          dimnames = list(NULL, sim$ids)),
                                   check.names = FALSE))
      data.table::fwrite(empty, file = file, sep = sep, na = "NA")
      return(invisible(file))
    }
    first <- TRUE
    for (ch in chunks) {
      G <- t(.dosage_block(sim, ch))                   # markers x individuals
      G <- .whole_to_integer(G)
      block <- cbind(meta[ch, , drop = FALSE],
                     as.data.frame(G, check.names = FALSE,
                                   stringsAsFactors = FALSE))
      rownames(block) <- NULL
      data.table::fwrite(block, file = file, sep = sep, na = "NA",
                         append = !first)
      first <- FALSE
    }
  }
  invisible(file)
}

#' Store whole-number doubles as integers (dosages are -1/0/1); leave the rest
#' @keywords internal
#' @noRd
.whole_to_integer <- function(x) {
  v <- x[is.finite(x)]
  if (all(v == floor(v)) && all(abs(v) < .Machine$integer.max)) {
    storage.mode(x) <- "integer"
  }
  x
}

#' JSON objects for a chunk of markers, one string per marker
#'
#' Vectorised over the chunk: every field is encoded column-wise and the pieces
#' are pasted, which is far faster than one `toJSON()` call per marker on a
#' 50k-marker panel. The dosage arrays follow the `individuals` order.
#' @keywords internal
#' @noRd
.marker_json_objects <- function(sim, meta, ch, causal_for = NULL) {
  m <- meta[ch, , drop = FALSE]
  fields <- vapply(names(m), function(col) {
    paste0('"', col, '":', .json_vec_elements(m[[col]]))
  }, character(length(ch)))
  if (length(ch) == 1L) fields <- matrix(fields, nrow = 1L)
  fields <- cbind(fields, paste0('"maf":', .json_vec_elements(sim$maf[ch])))
  G <- .dosage_block(sim, ch)                            # individuals x markers
  G <- .whole_to_integer(G)
  Gc <- matrix(.json_vec_elements(as.vector(G)), nrow = nrow(G))
  geno <- if (nrow(Gc) == 1L) Gc[1L, ] else
    do.call(paste, c(lapply(seq_len(nrow(Gc)), function(i) Gc[i, ]),
                     sep = ","))
  fields <- cbind(fields, paste0('"genotypes":[', geno, ']'))
  if (!is.null(causal_for)) {
    cf <- vapply(m$snp, function(s) {
      rows <- causal_for[[s]]
      if (is.null(rows)) return("[]")
      as.character(jsonlite::toJSON(rows, dataframe = "rows", digits = I(17),
                                    na = "null", auto_unbox = TRUE))
    }, character(1))
    fields <- cbind(fields, paste0('"causal_for":', cf))
  }
  paste0("{", do.call(paste, c(lapply(seq_len(ncol(fields)),
                                       function(j) fields[, j]), sep = ",")),
         "}")
}

#' Element-wise JSON encoding of a vector (strings escaped, numbers at 17
#' significant digits, NA/NaN/Inf as null)
#'
#' Plain atomic vectors are formatted here. A classed vector other than a
#' factor (`bit64::integer64`, `Date`, `POSIXct`, ...) is handed to jsonlite
#' element by element, so it is encoded exactly as `jsonlite::toJSON()` would
#' encode it through the class's own method: an `integer64` is stored as a
#' double bit pattern, and `sprintf()` on it would print a meaningless tiny
#' number instead of the integer.
#' @keywords internal
#' @noRd
.json_vec_elements <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.object(x)) {
    return(.json_classed_elements(x))
  }
  if (is.character(x)) {
    out <- .json_escape(x)
    out[is.na(x)] <- "null"
    return(out)
  }
  if (is.logical(x)) {
    out <- ifelse(x, "true", "false")
    out[is.na(x)] <- "null"
    return(out)
  }
  out <- if (is.integer(x)) as.character(x) else sprintf("%.17g", x)
  out[!is.finite(x)] <- "null"
  out
}

#' jsonlite's encoding of each element of a classed vector
#'
#' One `toJSON()` call encodes the whole vector as an array; when that array
#' holds no string (numbers, nulls, booleans) it is split on the commas, which
#' is exact and fast. Otherwise (string-valued classes such as `Date`) each
#' element is encoded on its own, so a comma inside a string cannot mislead.
#' @keywords internal
#' @noRd
.json_classed_elements <- function(x) {
  enc <- function(v) {
    as.character(jsonlite::toJSON(v, na = "null", digits = I(17)))
  }
  if (length(x) == 0L) return(character(0))
  whole <- enc(x)
  if (!grepl('"', whole, fixed = TRUE)) {
    out <- strsplit(substr(whole, 2L, nchar(whole) - 1L), ",", fixed = TRUE)[[1L]]
    if (length(out) == length(x)) return(out)
  }
  vapply(seq_along(x), function(i) {
    s <- enc(x[i])
    substr(s, 2L, nchar(s) - 1L)
  }, character(1))
}

#' A whole vector as one JSON array
#' @keywords internal
#' @noRd
.json_vec <- function(x) {
  paste0("[", paste(.json_vec_elements(x), collapse = ","), "]")
}

#' Quote and escape strings for JSON (RFC 8259: backslash, quote, control
#' characters)
#' @keywords internal
#' @noRd
.json_escape <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  x <- gsub("\b", "\\b", x, fixed = TRUE)
  x <- gsub("\f", "\\f", x, fixed = TRUE)
  ctrl <- grepl("[\x01-\x1f]", x, perl = TRUE, useBytes = TRUE)
  if (any(ctrl)) {
    x[ctrl] <- vapply(x[ctrl], function(s) {
      chars <- strsplit(s, "", fixed = TRUE)[[1L]]
      code <- utf8ToInt(s)
      bad <- code < 32L
      chars[bad] <- sprintf("\\u%04x", code[bad])
      paste(chars, collapse = "")
    }, character(1), USE.NAMES = FALSE)
  }
  paste0('"', x, '"')
}

#' Genetic values behind the simulated phenotypes
#'
#' The genetic component of each individual's phenotype: the sum of the
#' scaled additive, dominance and epistatic layers, before the residual is
#' added. This is the quantity the v1 engine wrote to `Genetic_values.txt`.
#'
#' When a phenotype includes a genome-**derived** `transcriptome()` layer, the
#' genetic-mediated part of that expression component (the share tracing to the
#' genome through gene expression) is also included, so it counts toward
#' heritability. A **real** expression source contributes nothing here: its
#' genetic content is not asserted, so it is treated as an environmental
#' predictor and excluded from the genetic value.
#'
#' Variance QTL are deliberately excluded: a `vqtl()` layer modulates the
#' spread of the residual rather than contributing a genetic value, so it has
#' no place in this matrix. By default this returns replication 1. When the
#' simulation used `vary_qtn = TRUE`, use `rep` to retrieve the architecture
#' and genetic values belonging to another replication.
#'
#' @param sim a `phenotype_sim`.
#' @param rep replication to retrieve (default 1).
#' @return An individuals-by-traits numeric matrix, with individual IDs as row
#'   names and trait names as column names.
#' @seealso [qtn_table()] for the loci and effects that produce these values,
#'   [phenotypes_long()] for the phenotypes themselves.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.5,
#'                          seed = 1) |>
#'   additive(n_qtn = 5)
#'
#' g <- genetic_values(ph)
#' head(g)
#'
#' # Realized genetic correlation between the two traits
#' round(cor(g[, 1], g[, 2]), 3)
#'
#' # Realized heritability: genetic variance over phenotypic variance
#' pheno <- phenotypes_wide(ph)
#' round(var(g[, 1]) / var(pheno$Trait_1), 3)
genetic_values <- function(sim, rep = 1L) {
  .check_sim(sim)
  .check_h2_complete(sim)
  rep <- .validate_rep(sim, rep)
  out <- .genetic_value_matrix(sim, rep)
  dimnames(out) <- list(sim$ids, paste0("Trait_", seq_len(sim$n_traits)))
  out
}

#' Genetic-mediated vs environmental split of a derived transcriptome component
#'
#' For a phenotype built with a genome-**derived** `transcriptome()` layer
#' (`simulate_phenotype(transcriptome = ...)`), the expression-mediated component
#' splits into a genetic-mediated part (the share of expression that traces to
#' the genome, which counts toward heritability and appears in
#' [genetic_values()]) and an environmental part. This returns that realized
#' split, per trait, as proportions of phenotypic variance:
#' `genetic_mediated`, `env_mediated`, and their `covariance` term
#' `2*Cov(Tx_g, Tx_e)/V_P` (finite-sample, ~0 because the generator draws the
#' genetic and non-genetic parts of expression independently). The three columns
#' sum to the realized expression-mediated share of variance.
#'
#' Returns `NULL` when the phenotype has no derived transcriptome layer: a
#' **real** expression source (`simulate_phenotype(expression = ...)`) has no
#' asserted genetic/environmental decomposition, so no split is reported.
#'
#' @param sim a `phenotype_sim`.
#' @return A data frame (`trait`, `genetic_mediated`, `env_mediated`,
#'   `covariance`), or `NULL`.
#' @seealso [genetic_values()], [transcriptome()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' tx <- simulate_transcriptome(SNP55K_maize282_maf04, n_genes = 100, seed = 1)
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.3, seed = 2,
#'                          transcriptome = tx) |>
#'   transcriptome(prop = 0.4, n_genes = 20)
#' mediation_split(ph)
mediation_split <- function(sim) {
  .check_sim(sim)
  sim$mediation
}

#' Per-QTN proportion of realized phenotypic variance for a mean-effect layer
#'
#' The realized additive/dominance layer is scaled by `k = sqrt(prop) / sd(raw)`
#' so its target variance equals `prop`. Each QTN's marginal contribution is
#' `k^2 * var(col_j)`, where `col_j` is the locus's genotypic contribution
#' (`effect * -1/0/1 dosage` for additive, `effect * heterozygote indicator` for
#' dominance, and `a * dosage + d * heterozygote` for an orthogonal layer, so its
#' dominance deviation is counted), and `var_explained` reports it as a
#' fraction of the *realized* phenotypic variance `var_p`. Dividing by the
#' realized (not the nominal, =1) V_P matters because finite-sample covariance
#' among components leaves the realized V_P slightly off 1. These are marginal
#' proportions: with LD between causal loci they do not sum exactly to `prop`,
#' because cross-locus covariances are not attributed to any single QTN. vqtl
#' layers return NA (they modulate the residual, not a genetic value).
#' @keywords internal
#' @noRd
.qtn_var <- function(sim, ly, t, idx, eff, var_p) {
  if (!ly$type %in% c("additive", "dominance")) {
    return(rep(NA_real_, length(idx)))
  }
  G <- .dosage_block(sim, idx)                       # -1/0/1 dosage
  if (isTRUE(ly$orthogonal)) {
    # Orthogonal layer: each locus contributes a * dosage + d * het, so its
    # marginal variance must include the dominance deviation, not just a.
    d_eff <- ly$d_effect[[t]]
    cols  <- sweep(G, 2L, eff, "*") + sweep((G == 0) * 1, 2L, d_eff, "*")
  } else {
    design <- if (ly$type == "dominance") (G == 0) * 1 else G
    cols   <- sweep(design, 2L, eff, "*")
  }
  raw <- rowSums(cols)
  s_raw <- .layer_sd(ly, raw, t)
  prop_t <- .expand_prop(ly$prop, sim$n_traits)[t]
  if (!is.finite(s_raw) || s_raw <= 0 || prop_t <= 0 ||
      !is.finite(var_p) || var_p <= 0) {
    return(rep(0, length(idx)))
  }
  k2 <- prop_t / s_raw^2
  vapply(seq_along(idx),
         function(j) k2 * stats::var(cols[, j]) / var_p, numeric(1))
}

#' Per-gene marginal variance share for a transcriptome layer
#'
#' The transcriptome component is `c * sum_g w_g z_g`, with `z_g` the
#' standardized expression of gene `g` (unit sample variance), `w_g` the
#' max-normalized slope, and `c = sqrt(prop) / sd(component)` the common scale
#' that realizes `prop`. Each gene's marginal contribution to the component is
#' `(c w_g)^2 Var(z_g)`; like \code{.qtn_var()} this is then divided by the
#' realized phenotypic variance `var_p`, so `var_explained` is a fraction of
#' realized phenotypic variance on the same scale as the marker rows. Like the
#' marker case these are *marginal* shares: with co-expression between causal
#' genes they do not sum exactly to `prop`, because the cross-gene covariances
#' are not attributed to any single gene.
#' @keywords internal
#' @noRd
.tx_qtn_var <- function(sim, ly, t, rep = 1L, var_p) {
  qe <- .layer_qtn_effect(ly, t, rep)
  idx <- qe$qtn; eff <- qe$effect
  if (is.null(idx) || length(idx) == 0) return(numeric(0))
  prop_t <- .expand_prop(ly$prop, sim$n_traits)[t]
  comp <- .tx_raw(ly, sim, t, rep, "total")
  s <- stats::sd(comp)
  if (!is.finite(s) || s <= 0 || prop_t <= 0 ||
      !is.finite(var_p) || var_p <= 0) return(rep(0, length(idx)))
  k <- sqrt(prop_t) / s
  w <- eff; sc <- max(abs(w))
  w <- if (is.finite(sc) && sc > 0) w / sc else rep(0, length(w))
  E <- sim$expression[idx, , drop = FALSE]
  vz <- apply(E, 1L, function(r) {                 # 1 for a standardized gene, 0 if constant
    sg <- stats::sd(r); if (is.finite(sg) && sg > 0) 1 else 0
  })
  (k * w)^2 * vz / var_p
}

#' The QTNs behind a simulation, with their effects
#'
#' One row per QTN per trait per layer, naming the marker and the effect it was
#' given. This is the v2 equivalent of the v1 `Additive_QTNs.txt` /
#' `QTN_effects_summary.txt` files, assembled from the simulation object rather
#' than written to disk.
#'
#' Epistatic QTNs come in interacting sets, so each member of a set gets its own
#' row, sharing the set's `effect` and identified by a common `set` number.
#' For every other layer type `set` is `NA`.
#'
#' Under `architecture = "ld"` the two traits have distinct causal loci that are
#' in linkage disequilibrium. `QTN_t1` and `QTN_t2` name the linked pair's causal
#' SNP for trait 1 and trait 2 (shown on both of the pair's rows so each row is
#' self-contained), and `ld_r2` is the squared correlation between them. All
#' three are `NA` for the other architectures. (For `ld_type = "indirect"` the
#' hidden cause-of-LD locus is not shown here -- it is not a QTN -- but is
#' available programmatically on the layer.)
#'
#' A `transcriptome()` layer's causal features are **genes**, not markers, so its
#' rows carry the gene identifier in the `snp` column, the per-gene slope in
#' `effect`, and the gene's marginal variance share in `var_explained`; the
#' marker-only columns `chr`, `pos` and `maf` are `NA`. Filter with
#' `subset(qtn_table(ph), layer == "transcriptome")` to isolate them.
#'
#' @param sim a `phenotype_sim`.
#' @param rep replication whose QTN architecture to report (default 1). This
#'   matters when the simulation used `vary_qtn = TRUE`.
#' @return A data frame with columns `trait`, `layer`, `set`, `snp`, `chr`,
#'   `pos`, `maf`, `effect`, `d`, `var_explained`, `QTN_t1`, `QTN_t2` and
#'   `ld_r2`, or a zero-row frame when no layers have been added. `effect` is the
#'   additive effect (`a` for an orthogonal layer); `d` is the per-locus
#'   dominance deviation of an orthogonal layer and `NA` for every other layer.
#'   For transcriptome layers `snp` holds the gene name and `chr`/`pos`/`maf` are
#'   `NA`.
#' @seealso [genetic_values()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.5, seed = 1) |>
#'   additive(prop = 0.3, n_qtn = 3) |>
#'   epistasis(prop = 0.2, n_pairs = 2)
#'
#' qtn_table(ph)
#'
#' # Which markers were used, and how large were the effects?
#' subset(qtn_table(ph), layer == "additive")
qtn_table <- function(sim, rep = 1L) {
  .check_sim(sim)
  .check_h2_complete(sim)
  rep <- .validate_rep(sim, rep)
  # Realized phenotypic variance per trait for this replication; var_explained
  # is reported as a fraction of it (SPEC realized variance-ratio convention).
  var_p <- vapply(seq_len(sim$n_traits), function(t) {
    y <- sim$pheno$value[sim$pheno$trait == paste0("Trait_", t) &
                         sim$pheno$rep == rep]
    stats::var(y)
  }, numeric(1))
  empty <- data.frame(
    trait = character(0), layer = character(0), set = integer(0),
    snp = character(0), chr = sim$map$chr[0], pos = sim$map$pos[0],
    maf = numeric(0), effect = numeric(0), d = numeric(0),
    var_explained = numeric(0),
    QTN_t1 = character(0), QTN_t2 = character(0), ld_r2 = numeric(0),
    stringsAsFactors = FALSE
  )
  if (length(sim$layers) == 0) {
    return(empty)
  }

  block <- function(trait, type, set, idx, effect, var_explained,
                    d = NA_real_,
                    QTN_t1 = NA_character_, QTN_t2 = NA_character_,
                    ld_r2 = NA_real_) {
    data.frame(
      trait         = trait,
      layer         = type,
      set           = set,
      snp           = sim$map$snp[idx],
      chr           = sim$map$chr[idx],
      pos           = sim$map$pos[idx],
      maf           = sim$maf[idx],
      effect        = effect,
      d             = d,
      var_explained = var_explained,
      QTN_t1        = QTN_t1,
      QTN_t2        = QTN_t2,
      ld_r2         = ld_r2,
      stringsAsFactors = FALSE
    )
  }

  # A transcriptome layer's causal features are genes, not markers: the feature
  # id (gene name) goes in `snp`, and the marker-only columns are NA (typed to
  # match the map so rows bind cleanly).
  gene_block <- function(trait, idx, effect, var_explained) {
    n <- length(idx)
    data.frame(
      trait         = trait,
      layer         = "transcriptome",
      set           = NA_integer_,
      snp           = rownames(sim$expression)[idx],
      chr           = rep(sim$map$chr[NA_integer_], n),
      pos           = rep(sim$map$pos[NA_integer_], n),
      maf           = rep(NA_real_, n),
      effect        = effect,
      d             = NA_real_,          # no dominance deviation for a gene predictor
      var_explained = var_explained,
      QTN_t1        = NA_character_,
      QTN_t2        = NA_character_,
      ld_r2         = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  for (ly in sim$layers) {
    for (t in seq_len(sim$n_traits)) {
      qe <- .layer_qtn_effect(ly, t, rep)
      idx <- qe$qtn
      eff <- qe$effect
      if (is.null(idx) || length(idx) == 0) {
        next
      }
      trait <- paste0("Trait_", t)
      if (identical(ly$type, "transcriptome")) {
        rows[[length(rows) + 1L]] <-
          gene_block(trait, idx, eff, .tx_qtn_var(sim, ly, t, rep, var_p[t]))
      } else if (identical(ly$type, "epistasis")) {
        # idx is an n_pairs x interaction matrix; every member of a set shares
        # the set's effect. Per-locus variance is undefined for an interaction.
        for (p in seq_len(nrow(idx))) {
          members <- idx[p, ]
          rows[[length(rows) + 1L]] <-
            block(trait, ly$type, p, members, rep(eff[p], length(members)),
                  rep(NA_real_, length(members)))
        }
      } else {
        qtn_t1 <- NA_character_
        qtn_t2 <- NA_character_
        r2 <- NA_real_
        ld <- if (!is.null(ly$qtn_reps)) attr(ly$qtn_reps[[rep]], "ld") else
          ly$ld
        if (!is.null(ld)) {
          qtn_t1 <- sim$map$snp[ld$qtn_t1]
          qtn_t2 <- sim$map$snp[ld$qtn_t2]
          r2 <- ld$r2
        }
        d_col <- if (isTRUE(ly$orthogonal)) ly$d_effect[[t]] else NA_real_
        rows[[length(rows) + 1L]] <-
          block(trait, ly$type, NA_integer_, idx, eff,
                .qtn_var(sim, ly, t, idx, eff, var_p[t]),
                d = d_col, QTN_t1 = qtn_t1, QTN_t2 = qtn_t2, ld_r2 = r2)
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
