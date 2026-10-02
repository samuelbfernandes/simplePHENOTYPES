#' Convert genomic data between formats or to a numeric dosage matrix.
#'
#' @keywords internal
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @param file Genotype input: a file path (character), an in-memory
#'   data.frame / matrix, a \code{gds.class} object, or a \code{vcfR} object.
#' @param from Character string naming the input format. When \code{NULL}
#'   (default) the format is detected automatically from the file extension or
#'   object class.
#' @param to Target format. Currently only \code{"numeric"} is implemented.
#' @param ref_allele Optional character vector (length = number of markers)
#'   specifying the reference allele for each SNP. Used with
#'   \code{method = "reference"}.
#' @param to_r Logical. Return the result as an R object? Defaults to
#'   \code{TRUE} for in-memory input, \code{FALSE} for file input.
#' @param to_file Logical. Write the result to a file? Defaults to the
#'   complement of \code{to_r}.
#' @param file_name Output file name. Auto-generated when \code{NULL}: the
#'   default is derived from the input's label (the file name, or the name of
#'   the object passed), so different inputs can share it; an existing
#'   default-named file is overwritten, with a warning.
#' @param f_name Base name used to build \code{file_name} for in-memory input.
#' @param code_as Numeric coding scheme: \code{"-101"} (major = 1, het = 0,
#'   minor = -1; default) or \code{"012"} (major = 2, het = 1, minor = 0).
#' @param hets Character vector of heterozygote codes for nucleotide tables
#'   (default: the IUPAC ambiguity codes and both orders of every digraph, the
#'   same set HapMap input uses). Matching is case-insensitive.
#' @param homo Character vector of homozygote codes for nucleotide tables
#'   (default: single letters and doubled letters). Matching is
#'   case-insensitive; a call that is in neither `hets`, `homo` nor the missing
#'   codes is set to missing and reported in a warning.
#' @param model Genetic model: \code{"Add"} (default), \code{"Dom"},
#'   \code{"Left"}, or \code{"Right"}.
#' @param impute Missing-data imputation: \code{"None"} (default),
#'   \code{"Middle"}, \code{"Minor"}, or \code{"Major"}.
#' @param method Allele orientation method: \code{"frequency"} (default) or
#'   \code{"reference"} (requires \code{ref_allele}).
#' @param verbose Logical; print progress messages when \code{TRUE}.
#' @param counted_column Logical (default \code{FALSE}). When \code{TRUE} the
#'   result, and the file written, gets a character column \code{counted}
#'   immediately after \code{cm} holding the allele coded +1 at each marker (the
#'   durable form of the \code{"counted_allele"} attribute, which a text file
#'   cannot carry). The default leaves the output unchanged.
#' @return When \code{to_r = TRUE}: a data.frame with columns
#'   \code{snp, allele, chr, pos, cm} followed by one column per sample.
#'   When \code{to_r = FALSE}: \code{invisible(NULL)}.
#' @author Samuel Fernandes
#' @seealso [as_numeric()], a shorthand for \code{to = "numeric"}.
format_conversion <- function(file,
                               from       = NULL,
                               to         = NULL,
                               ref_allele = NULL,
                               to_r       = NULL,
                               to_file    = NULL,
                               file_name  = NULL,
                               f_name     = NULL,
                               code_as    = "-101",
                               hets       = .HETS,
                               homo       = .HOMO,
                               model      = "Add",
                               impute     = "None",
                               method     = "frequency",
                               verbose    = TRUE,
                               counted_column = FALSE) {

  custom_hets <- !missing(hets)
  custom_homo <- !missing(homo)
  to <- if (is.null(to)) "numeric" else match.arg(to, "numeric")
  code_as <- match.arg(code_as, c("-101", "012"))
  model <- match.arg(model, c("Add", "Dom", "Left", "Right"))
  impute <- match.arg(impute, c("None", "Middle", "Minor", "Major"))
  method <- match.arg(method, c("frequency", "reference"))
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(counted_column) || length(counted_column) != 1L ||
      is.na(counted_column)) {
    stop("`counted_column` must be TRUE or FALSE.", call. = FALSE)
  }
  if (counted_column && identical(model, "Dom")) {
    stop("`counted_column = TRUE` records the allele counted as +1, and ",
         "model = \"Dom\" counts no allele.", call. = FALSE)
  }
  for (nm in c("to_r", "to_file")) {
    value <- get(nm)
    if (!is.null(value) &&
        (!is.logical(value) || length(value) != 1L || is.na(value))) {
      stop("`", nm, "` must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
  }

  file_class <- class(file)

  # ---- resolve to_r / to_file defaults ------------------------------------
  if (!is.null(to_r) && !is.null(to_file)) {
    if (!to_r && !to_file) {
      stop("At least one of `to_r` and `to_file` must be TRUE.",
           call. = FALSE)
    }
  } else if (is.null(to_file) && is.null(to_r)) {
    if (all(file_class == "character")) {
      to_file <- TRUE; to_r <- FALSE
    } else {
      to_file <- FALSE; to_r <- TRUE
    }
  } else if (is.null(to_file)) {
    to_file <- !to_r
  } else if (is.null(to_r)) {
    to_r <- !to_file
  }
  if (!to_file && !is.null(file_name)) {
    stop("`file_name` was supplied but `to_file` is FALSE.", call. = FALSE)
  }

  # ---- validate file exists -----------------------------------------------
  if (all(file_class == "character") && !file.exists(file)) {
    stop("file '", file, "' not found!", call. = FALSE)
  }

  # ---- auto-generate output file name -------------------------------------
  default_name <- is.null(file_name) && to_file && to == "numeric"
  if (default_name) {
    if (all(file_class == "character")) {
      # Rewrite only the extension of the file name itself (".hmp.txt", or the
      # last ".ext", after an optional .gz/.bz2); never the directory part, and
      # a name without an extension simply gains the suffix.
      dir_part <- dirname(file)
      base <- sub("(\\.hmp\\.txt|\\.[^./\\\\]*)$", "",
                  sub("\\.(gz|bz2)$", "", basename(file), ignore.case = TRUE),
                  ignore.case = TRUE)
      out_base <- paste0(base, "_numeric.txt")
      file_name <- if (identical(dir_part, ".") && !startsWith(file, "./")) {
        out_base
      } else {
        file.path(dir_part, out_base)
      }
    } else {
      # a default file name must be portable: an inline-object label such as
      # "<inline data.frame 3 x 13>" becomes "inline_data.frame_3_x_13". When
      # sanitizing changes the label (replaced characters, trimmed edges,
      # fallback or truncation to 100 characters) a short stable hash of the
      # ORIGINAL label is appended, which makes accidental clashes unlikely
      # but cannot make the name unique; labels that need no change are kept
      # exactly. An existing default-named file is overwritten with a warning.
      label <- as.character(f_name)[1L]
      safe <- gsub("[^A-Za-z0-9._-]+", "_", label)
      safe <- gsub("^_+|_+$", "", safe)
      if (is.na(safe) || !nzchar(safe)) safe <- "geno"
      changed <- is.na(label) || !identical(safe, label)
      if (nchar(safe) > 100L) {
        safe <- substr(safe, 1L, 100L)
        changed <- TRUE
      }
      if (changed) safe <- paste0(safe, "_", .label_hash(label))
      file_name <- paste0(safe, "_numeric.txt")
    }
    # A default name is derived from the label, so different inputs can share
    # it (a hash cannot make names unique, and a case-insensitive file system
    # folds case). The writer overwrites, so an existing default-named file is
    # reported rather than replaced silently; that check is made by
    # .write_numeric() just before the file is written, so a conversion that
    # fails never warns about a file it did not touch.
  }

  # ---- detect format -------------------------------------------------------
  if (is.null(from)) {
    from <- detect_format(file)
    if (from == "unknown") {
      stop(paste0("The format was not detected automatically. ",
                  'Please set "from" to one of: "hapmap", "vcf", "gds", ',
                  '"bed", "ped", "finalreport", or "table".'), call. = FALSE)
    }
  }
  from <- tolower(from)
  supported <- c("hapmap", "vcf", "vcfr", "vcfr_object", "gds",
                 "gds_object", "bed", "ped", "finalreport", "table",
                 "numeric")
  if (length(from) != 1L || is.na(from) || !from %in% supported) {
    stop("Unsupported `from` format: ", paste(from, collapse = ", "), ".",
         call. = FALSE)
  }
  if (from != "table" && (custom_hets || custom_homo)) {
    stop("`hets` and `homo` customize nucleotide-table parsing and may only ",
         "be supplied with from = \"table\".", call. = FALSE)
  }
  if (from == "table") {
    if (!is.character(hets) || !length(hets) || anyNA(hets) ||
        any(!nzchar(hets)) || !is.character(homo) || !length(homo) ||
        anyNA(homo) || any(!nzchar(homo))) {
      stop("`hets` and `homo` must be non-empty character vectors without ",
           "missing or empty codes.", call. = FALSE)
    }
    if (length(intersect(toupper(hets), toupper(homo)))) {
      stop("`hets` and `homo` must not contain overlapping genotype codes.",
           call. = FALSE)
    }
  }

  if (method == "frequency" && !is.null(ref_allele)) {
    stop("`ref_allele` is used only with method = \"reference\".",
         call. = FALSE)
  }
  if (method == "reference" && is.null(ref_allele) &&
      from %in% c("hapmap", "table")) {
    stop("method = \"reference\" requires `ref_allele`.", call. = FALSE)
  }

  # `method = "reference"` orients the coding by a user-supplied allele, which is
  # only threaded through the hapmap and table handlers. For binary/gds formats
  # the allele orientation comes from the file itself; silently ignoring
  # `ref_allele` there would mislead, so refuse rather than pretend.
  if (identical(method, "reference") &&
      !from %in% c("hapmap", "table")) {
    stop("method = \"reference\" is supported only for HapMap and nucleotide ",
         "table input; ", from, " files carry their own allele orientation. ",
         "Use method = \"frequency\" (the default).", call. = FALSE)
  }

  # ---- dispatch ------------------------------------------------------------
  # Every handler only converts (to_file = FALSE, to_r = TRUE): the optional
  # `counted` column and the file write happen once, below.
  G <- switch(
    from,
    vcf         = ,
    vcfr        = ,
    vcfr_object = handle_vcf(file, from, file_class, file_name,
                              FALSE, TRUE, to, code_as,
                              model, impute, method, verbose),

    hapmap      = handle_hapmap(file, file_class, file_name,
                                FALSE, TRUE, to, code_as, ref_allele,
                                model, impute, method, verbose),

    table       = handle_table(file, file_class, file_name,
                               FALSE, TRUE, to, code_as, ref_allele,
                               hets, homo, model, impute, method, verbose),

    gds         = ,
    gds_object  = handle_gds(file, file_name, FALSE, TRUE, to,
                              code_as, model, impute, method, verbose),

    bed         = handle_bed(file, file_name, FALSE, TRUE, to,
                             code_as, model, impute, method, verbose),

    ped         = handle_ped(file, file_name, FALSE, TRUE, to,
                             code_as, model, impute, method, verbose),

    finalreport = handle_finalreport(file, file_name, FALSE, TRUE, to,
                                     code_as, model, impute, method, verbose),

    numeric     = handle_numeric(file, file_name, FALSE, TRUE, code_as,
                                 model, impute, method, ref_allele, verbose),

    stop(paste0('Format "', from, '" is not supported. ',
                'Use one of: "hapmap", "vcf", "gds", "bed", "ped", ',
                '"finalreport", "table".'), call. = FALSE)
  )

  if (counted_column) G <- .add_counted_column(G)
  if (to_file) .write_numeric(G, file_name, default_name, verbose)

  if (verbose) message("Genotype conversion complete.")
  if (to_r) return(G)
}

#' Short stable hash of a label (8 hex digits)
#'
#' Polynomial rolling hash of the UTF-8 bytes modulo the prime 4294967291
#' (all intermediate values stay below 2^53, so it is exact in doubles and
#' identical on every platform). Used only to make accidental clashes between
#' sanitized default output names unlikely; it cannot make them unique and is
#' not a security hash.
#' @keywords internal
#' @noRd
.label_hash <- function(label) {
  label <- if (is.na(label)) "NA" else label
  bytes <- as.integer(charToRaw(enc2utf8(label)))
  h <- 5381
  for (b in bytes) h <- (h * 131 + b) %% 4294967291
  sprintf("%04x%04x", as.integer(h %/% 65536), as.integer(h %% 65536))
}
