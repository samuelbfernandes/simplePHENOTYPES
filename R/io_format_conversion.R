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
#' @param file_name Output file name. Auto-generated when \code{NULL}.
#' @param f_name Base name used to build \code{file_name} for in-memory input.
#' @param code_as Numeric coding scheme: \code{"-101"} (major = 1, het = 0,
#'   minor = -1; default) or \code{"012"} (major = 2, het = 1, minor = 0).
#' @param hets Character vector of heterozygote codes for nucleotide tables.
#' @param homo Character vector of homozygote codes for nucleotide tables.
#' @param model Genetic model: \code{"Add"} (default), \code{"Dom"},
#'   \code{"Left"}, or \code{"Right"}.
#' @param impute Missing-data imputation: \code{"None"} (default),
#'   \code{"Middle"}, \code{"Minor"}, or \code{"Major"}.
#' @param method Allele orientation method: \code{"frequency"} (default) or
#'   \code{"reference"} (requires \code{ref_allele}).
#' @param verbose Logical; print progress messages when \code{TRUE}.
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
                               hets       = c("R","Y","S","W","K","M",
                                              "AG","CT","CG","AT","GT","AC"),
                               homo       = c("A","AA","T","TT","C","CC","G","GG"),
                               model      = "Add",
                               impute     = "None",
                               method     = "frequency",
                               verbose    = TRUE) {

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
  if (is.null(file_name) && to_file && to == "numeric") {
    if (all(file_class == "character")) {
      file_name <- gsub(
        paste0(gsub(".*[.]", ".", file), "|.HMP.TXT"),
        "_numeric.txt", file, ignore.case = TRUE)
    } else {
      file_name <- paste0(f_name, "_numeric.txt")
    }
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
  G <- switch(
    from,
    vcf         = ,
    vcfr        = ,
    vcfr_object = handle_vcf(file, from, file_class, file_name,
                              to_file, to_r, to, code_as,
                              model, impute, method, verbose),

    hapmap      = handle_hapmap(file, file_class, file_name,
                                to_file, to_r, to, code_as, ref_allele,
                                model, impute, method, verbose),

    table       = handle_table(file, file_class, file_name,
                               to_file, to_r, to, code_as, ref_allele,
                               hets, homo, model, impute, method, verbose),

    gds         = ,
    gds_object  = handle_gds(file, file_name, to_file, to_r, to,
                              code_as, model, impute, method, verbose),

    bed         = handle_bed(file, file_name, to_file, to_r, to,
                             code_as, model, impute, method, verbose),

    ped         = handle_ped(file, file_name, to_file, to_r, to,
                             code_as, model, impute, method, verbose),

    finalreport = handle_finalreport(file, file_name, to_file, to_r, to,
                                     code_as, model, impute, method, verbose),

    numeric     = handle_numeric(file, file_name, to_file, to_r, code_as,
                                 model, impute, method, ref_allele, verbose),

    stop(paste0('Format "', from, '" is not supported. ',
                'Use one of: "hapmap", "vcf", "gds", "bed", "ped", ',
                '"finalreport", "table".'), call. = FALSE)
  )

  if (verbose) message("Genotype conversion complete.")
  if (to_r) return(G)
}
