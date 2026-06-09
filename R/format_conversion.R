#' Convert genomic data between formats or to a numeric dosage matrix.
#'
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn ls.gdsn
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
#' @param ... Additional arguments (currently unused).
#' @return When \code{to_r = TRUE}: a data.frame with columns
#'   \code{snp, allele, chr, pos, cm} followed by one column per sample.
#'   When \code{to_r = FALSE}: \code{invisible(NULL)}.
#' @author Samuel Fernandes
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
                               verbose    = TRUE,
                               ...) {

  file_class <- class(file)

  # ---- resolve to_r / to_file defaults ------------------------------------
  if (!is.null(to_r) && !is.null(to_file)) {
    if (!to_r && !to_file) {
      warning("to_r and to_file are both FALSE. Setting to_r = TRUE.",
              call. = FALSE, immediate. = TRUE)
      to_r <- TRUE
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

  # ---- dispatch ------------------------------------------------------------
  G <- switch(
    from,
    VCF         = ,
    vcf         = ,
    vcfR        = ,
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
    GDS         = ,
    gds_object  = handle_gds(file, file_name, to_file, to_r, to,
                              code_as, model, impute, method, verbose),

    bed         = ,
    BED         = handle_bed(file, file_name, to_file, to_r, to,
                             code_as, model, impute, method, verbose),

    ped         = ,
    PED         = handle_ped(file, file_name, to_file, to_r, to,
                             code_as, model, impute, method, verbose),

    finalreport = handle_finalreport(file, file_name, to_file, to_r, to,
                                     code_as, model, impute, method, verbose),

    stop(paste0('Format "', from, '" is not supported. ',
                'Use one of: "hapmap", "vcf", "gds", "bed", "ped", ',
                '"finalreport", "table".'), call. = FALSE)
  )

  if (to_r) return(G)
}
