# Shared format detection and HapMap character parsing utilities.
# These functions are the single source-of-truth used by both
# format_conversion() (as_numeric path) and file_loader() (create_phenotypes path).

.HMP_NAMES <- c("rs#", "alleles", "chrom", "pos", "strand",
                "assembly#", "center", "protLSID", "assayLSID",
                "panelLSID", "QCcode")

# IUPAC het and homozygote codes shared across handlers.
.HETS <- c("R","Y","S","W","K","M",
           "AG","CT","CG","AT","GT","AC",
           "GA","TC","GC","TA","TG","CA")
.HOMO <- c("A","AA","T","TT","C","CC","G","GG")
.MISS <- c("N","NN","NA","--","XX","00","+","++"," ","")

#' Detect the genotype-data format of a file path or in-memory object.
#'
#' @param file Character file path, a data.frame/matrix, or a gds.class / vcfR
#'   object.
#' @return A character string: one of "hapmap", "vcf", "gds", "bed", "ped",
#'   "finalreport", "numeric", "gds_object", "vcfr_object", or "unknown".
#' @noRd
detect_format <- function(file) {
  if (all(class(file) == "character")) {
    upper <- toupper(file[[1]])
    if (endsWith(upper, ".HMP.TXT")) return("hapmap")
    if (endsWith(upper, ".GDS"))     return("gds")
    if (endsWith(upper, ".VCF"))     return("vcf")
    if (endsWith(upper, ".BED"))     return("bed")
    if (endsWith(upper, ".PED"))     return("ped")
    if (endsWith(upper, ".TXT") || endsWith(upper, ".CSV")) {
      hdr <- tryCatch(
        data.table::fread(file[[1]], nrows = 0L, data.table = FALSE),
        error = function(e) NULL
      )
      if (is.null(hdr)) return("unknown")
      nms <- names(hdr)
      if (any(tolower(nms) == "snp name")) return("finalreport")
      if (sum(nms[seq_len(min(11L, length(nms)))] == .HMP_NAMES) > 8L)
        return("hapmap")
      if (length(nms) >= 5L &&
          all(tolower(nms[1:5]) == c("snp", "allele", "chr", "pos", "cm")))
        return("numeric")
    }
    return("unknown")
  }
  if (inherits(file, "gds.class")) return("gds_object")
  if (inherits(file, "vcfR"))      return("vcfr_object")
  if (is.data.frame(file) || is.matrix(file)) {
    return(.detect_df_format(as.data.frame(file)))
  }
  "unknown"
}

.detect_df_format <- function(df) {
  nms <- names(df)
  if (length(nms) >= 11L && sum(nms[1:11] == .HMP_NAMES) > 8L)
    return("hapmap")
  if (length(nms) >= 12L) {
    vals <- unique(as.character(df[, 12]))
    chars <- unlist(strsplit(vals, ""))
    nucleotide <- c(.HETS, .HOMO, "N", "-", "+", "0")
    if (all(chars %in% nucleotide)) return("hapmap")
  }
  if (length(nms) >= 5L &&
      all(tolower(nms[1:5]) == c("snp", "allele", "chr", "pos", "cm")))
    return("numeric")
  if (length(nms) >= 2L) {
    corner <- unlist(df[seq_len(min(3L, nrow(df))),
                        seq(max(1L, ncol(df) - 1L), ncol(df))])
    if (any(grepl("[/|]", corner))) return("vcfr_object")
  }
  "unknown"
}

# ---------------------------------------------------------------------------
# HapMap character → raw dosage (0/1/2/NA integer matrix, SNPs × samples)
#
# Returns list(raw = integer matrix, flip = logical vector).
# flip[i] = FALSE always: canonical form puts major allele as 0.
# The flip vector passed to numericalize_core() is computed separately by
# compute_flip() based on allele frequencies in the raw matrix.
# ---------------------------------------------------------------------------

#' Parse HapMap character genotype matrix to raw 0/1/2 dosage matrix.
#'
#' @param geno_mat Character matrix (SNPs × samples) from a HapMap file
#'   (columns 12: after the 11 metadata columns).
#' @return Integer matrix (SNPs × samples): 0 = hom allele-1, 1 = het,
#'   2 = hom allele-2, NA_integer_ = missing.
#' @noRd
parse_hapmap_chars_to_raw <- function(geno_mat) {
  n_snp  <- nrow(geno_mat)
  n_samp <- ncol(geno_mat)
  raw    <- matrix(NA_integer_, nrow = n_snp, ncol = n_samp)

  for (i in seq_len(n_snp)) {
    row      <- as.character(geno_mat[i, ])
    is_miss  <- row %in% .MISS | is.na(row)
    is_het   <- row %in% .HETS
    is_hom   <- !is_miss & !is_het

    hom_vals <- row[is_hom]
    if (length(hom_vals) == 0L) {
      # Only hets and/or missing — treat all present calls as het
      raw[i, is_het] <- 1L
      next
    }

    counts  <- sort(table(hom_vals), decreasing = TRUE)
    allele1 <- names(counts)[1L]  # most frequent homozygote = allele-1

    raw[i, row == allele1] <- 0L
    raw[i, is_het]         <- 1L
    raw[i, is_hom & row != allele1] <- 2L
    # missing stays NA_integer_
  }
  raw
}

#' Compute per-SNP flip flags from a raw 0/1/2 dosage matrix.
#'
#' flip[i] = TRUE when allele-2 (raw 2) is more frequent than allele-1 (raw 0).
#' This is used to ensure that numericalize_core() assigns major_val to the
#' more common homozygote.
#'
#' @param raw Integer matrix (SNPs × samples), values 0/1/2/NA.
#' @param method "frequency" (default) or "reference".
#' @param allele1 Character vector length n_snp; used only when method =
#'   "reference". When allele1 matches the desired reference, no flip.
#' @param ref Character vector length n_snp; reference alleles (method =
#'   "reference" only).
#' @noRd
compute_flip <- function(raw, method = "frequency",
                         allele1 = NULL, ref = NULL) {
  if (method == "reference") {
    if (is.null(allele1) || is.null(ref)) {
      stop("allele1 and ref must be supplied when method = 'reference'.",
           call. = FALSE)
    }
    # flip when allele1 does NOT match the reference (ref is major → allele1
    # should be treated as major, so no flip needed when allele1 == ref).
    return(allele1 != ref)
  }
  # frequency method
  count0 <- rowSums(raw == 0L, na.rm = TRUE)
  count2 <- rowSums(raw == 2L, na.rm = TRUE)
  count2 > count0
}

#' Assemble the final output data.frame from metadata and coded genotype matrix.
#'
#' @param meta  data.frame with columns snp, allele, chr, pos, cm.
#' @param sample_ids Character vector of sample names.
#' @param coded_mat Integer matrix (SNPs × samples) of coded genotype values.
#' @noRd
assemble_output <- function(meta, sample_ids, coded_mat) {
  out <- data.frame(meta, coded_mat,
                    check.names = FALSE, fix.empty.names = FALSE,
                    stringsAsFactors = FALSE)
  colnames(out) <- c("snp", "allele", "chr", "pos", "cm", sample_ids)
  out
}
