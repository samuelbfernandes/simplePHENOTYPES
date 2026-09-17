# Shared format detection and HapMap character parsing utilities.
# These functions are the single source-of-truth for numericalization, used by
# format_conversion() / as_numeric() — the one coding path, which
# create_phenotypes() reaches through genotypes().

.HMP_NAMES <- c("rs#", "alleles", "chrom", "pos", "strand",
                "assembly#", "center", "protLSID", "assayLSID",
                "panelLSID", "QCcode")

# IUPAC het and homozygote codes shared across handlers.
.HETS <- c("R","Y","S","W","K","M",
           "AG","CT","CG","AT","GT","AC",
           "GA","TC","GC","TA","TG","CA")
.HOMO <- c("A","AA","T","TT","C","CC","G","GG")
.MISS <- c("N","NN","NA","--","XX","00","+","++"," ","")

# IUPAC single-character ambiguity codes -> their two constituent allele letters.
.HET_IUPAC <- c(R = "AG", Y = "CT", S = "GC", W = "AT", K = "GT", M = "AC")

#' The two allele letters carried by a heterozygote call.
#'
#' Decodes either an IUPAC ambiguity code (`"R"` -> `c("A","G")`) or a digraph
#' (`"AG"` -> `c("A","G")`). Returns `character(0)` for anything else. Used to
#' recover the allele labels of a marker that is heterozygous in every sample, so
#' a het-only marker can still be reference-oriented.
#' @noRd
.het_to_letters <- function(code) {
  code <- toupper(as.character(code))
  if (length(code) != 1L || is.na(code)) return(character(0))
  if (nchar(code) == 1L) {
    dec <- .HET_IUPAC[[code]]
    if (is.null(dec)) return(character(0))
    strsplit(dec, "")[[1L]]
  } else if (nchar(code) == 2L) {
    strsplit(code, "")[[1L]]
  } else {
    character(0)
  }
}

#' The allele letter(s) a single genotype call resolves to.
#'
#' A heterozygote call (IUPAC ambiguity code `"R"` or digraph `"AG"`) resolves to
#' its two constituent alleles; a homozygote code (`"A"`, `"AA"`) to its own
#' letters. Used to build the set of alleles observed at a marker, so an IUPAC
#' heterozygote contributes both of its alleles rather than the ambiguity letter
#' itself (a raw `strsplit("R")` would observe `"R"`, not `"A","G"`).
#' @noRd
.call_to_letters <- function(code) {
  code <- toupper(as.character(code))
  if (length(code) != 1L || is.na(code)) return(character(0))
  het <- .het_to_letters(code)
  if (length(het)) return(het)          # het: IUPAC code or digraph
  strsplit(code, "")[[1L]]              # homozygote: "A" or "AA"
}

#' Install hint for the optional Bioconductor readers.
#'
#' SNPRelate and gdsfmt are in `Suggests` (they are Bioconductor, not CRAN), so
#' they are only needed to read GDS / VCF / PLINK BED / PLINK PED. Every read
#' site guards with `requireNamespace()` and raises this message when they are
#' absent. HapMap, the numeric object, and a plain -1/0/1 matrix never need them.
#' @noRd
.gds_needed <- function(fmt) {
  paste0("Reading ", fmt, " files needs the Bioconductor packages SNPRelate ",
         "and gdsfmt, which are not installed. Install them with:\n",
         "  if (!requireNamespace(\"BiocManager\", quietly = TRUE)) ",
         "install.packages(\"BiocManager\")\n",
         "  BiocManager::install(c(\"SNPRelate\", \"gdsfmt\"))")
}

#' Detect the genotype-data format of a file path or in-memory object.
#'
#' @param file Character file path, a data.frame/matrix, or a gds.class / vcfR
#'   object.
#' @return A character string: one of "hapmap", "vcf", "gds", "bed", "ped",
#'   "finalreport", "numeric", "gds_object", "vcfr_object", or "unknown".
#' @noRd
detect_format <- function(file) {
  if (all(class(file) == "character")) {
    # Compressed text is handled transparently by data.table::fread(); strip a
    # trailing .gz/.bz2 so the real extension is what drives detection.
    upper <- sub("\\.(GZ|BZ2)$", "", toupper(file[[1]]))
    if (endsWith(upper, ".HMP.TXT")) return("hapmap")
    if (endsWith(upper, ".GDS"))     return("gds")
    if (endsWith(upper, ".VCF"))     return("vcf")
    if (endsWith(upper, ".BED"))     return("bed")
    if (endsWith(upper, ".PED"))     return("ped")
    # Content sniff for text with a generic (or wrong) extension: VCF carries a
    # "##fileformat=VCF" preamble and Illumina FinalReport an "[Header]" block,
    # neither of which a header-name match would catch (fread would read the
    # preamble, not the real column line).
    sniff <- .sniff_text_signature(file[[1]])
    if (!is.null(sniff)) return(sniff)
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

#' Sniff the first non-empty lines of a text file for a format signature.
#'
#' Catches formats whose extension may be generic (`.txt`) but whose content is
#' unambiguous: a VCF `##fileformat=VCF...` preamble (or a `#CHROM POS ID REF ALT`
#' header line) and an Illumina FinalReport `[Header]` block. Returns the format
#' string, or `NULL` when nothing matches (so the caller falls back to
#' header-name detection). Reads at most a handful of lines and never errors.
#' @noRd
.sniff_text_signature <- function(path) {
  lines <- tryCatch(
    readLines(path, n = 40L, warn = FALSE),
    error = function(e) character(0)
  )
  if (!length(lines)) return(NULL)
  lines <- trimws(lines)
  if (any(grepl("^##fileformat=VCF", lines, ignore.case = TRUE))) return("vcf")
  if (any(grepl("^#CHROM\t|^#CHROM ", lines))) return("vcf")
  if (any(grepl("^\\[Header\\]", lines, ignore.case = TRUE))) return("finalreport")
  NULL
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

#' Parse a character genotype matrix to a raw 0/1/2 dosage matrix.
#'
#' Generalized from the HapMap parser so it also serves nucleotide tables, VCF
#' `GT` strings, and Illumina FinalReport `AB` calls: the heterozygote,
#' homozygote and missing code sets are all overridable.
#'
#' @param geno_mat Character matrix (SNPs × samples) of genotype calls.
#' @param allele1 optional allele-1 label per marker. When supplied, raw dosage
#'   0 is anchored to that allele rather than to the most frequent homozygote.
#'   Required for reference-based orientation.
#' @param hets heterozygote codes (default the IUPAC/digraph set `.HETS`).
#' @param homo optional explicit homozygote codes. When `NULL` (default) any
#'   present call that is neither het nor missing is treated as a homozygote
#'   (HapMap behaviour). When supplied, a call that is neither het, homo, nor
#'   missing is treated as missing (invalid), which is what a strict table /
#'   FinalReport parse needs.
#' @param miss missing-value codes (default `.MISS`).
#' @return Integer matrix (SNPs × samples): 0 = hom allele-1, 1 = het,
#'   2 = hom allele-2, NA_integer_ = missing. `attr(., "alleles")` is an
#'   n_snp × 2 character matrix of the (allele1, allele2) labels backing the
#'   0 / 2 codes, for building the `allele` metadata column.
#' @noRd
parse_hapmap_chars_to_raw <- function(geno_mat, allele1 = NULL,
                                      hets = .HETS, homo = NULL,
                                      miss = .MISS) {
  n_snp  <- nrow(geno_mat)
  n_samp <- ncol(geno_mat)
  raw    <- matrix(NA_integer_, nrow = n_snp, ncol = n_samp)
  alleles <- matrix(NA_character_, nrow = n_snp, ncol = 2L)

  for (i in seq_len(n_snp)) {
    row      <- as.character(geno_mat[i, ])
    is_miss  <- row %in% miss | is.na(row)
    is_het   <- row %in% hets
    if (is.null(homo)) {
      is_hom <- !is_miss & !is_het
    } else {
      is_hom <- row %in% homo
      # A present call that is neither het nor a recognised homozygote is not a
      # valid biallelic genotype: record it as missing rather than a phantom
      # third homozygote.
      is_miss <- is_miss | (!is_het & !is_hom)
    }

    hom_vals <- row[is_hom]
    if (length(hom_vals) == 0L) {
      # Only hets and/or missing: every present call is het -> raw dosage 1.
      # There is no homozygote to name the alleles from, but a reference-oriented
      # parse still needs the marker's allele letters, so recover them from the
      # heterozygote calls themselves (IUPAC code or digraph). Leaving the allele
      # metadata NA here made a valid het-only marker fail reference coding.
      raw[i, is_het] <- 1L
      het_letters <- unique(unlist(lapply(row[is_het], .het_to_letters)))
      if (length(het_letters) == 2L) {
        if (is.null(allele1)) {
          alleles[i, ] <- het_letters
        } else {
          a <- substr(as.character(allele1[[i]]), 1L, 1L)
          alleles[i, ] <- if (a %in% het_letters) {
            c(a, setdiff(het_letters, a))
          } else {
            het_letters
          }
        }
      }
      next
    }

    counts <- sort(table(hom_vals), decreasing = TRUE)

    if (length(counts) > 2L) {
      # Non-biallelic SNP: set entire row to NA, matching v1 behaviour.
      message("Non-biallelic SNP at row ", i, " set to NA.")
      next
    }

    first <- if (is.null(allele1)) {
      names(counts)[1L]
    } else {
      a <- as.character(allele1[[i]])
      candidates <- c(a, paste0(a, a))
      hit <- candidates[candidates %in% hom_vals]
      if (length(hit)) hit[[1L]] else candidates[[1L]]
    }
    second <- setdiff(names(counts), first)
    alleles[i, ] <- c(first, if (length(second)) second[[1L]] else NA_character_)

    raw[i, row == first] <- 0L
    raw[i, is_het]         <- 1L
    raw[i, is_hom & row != first] <- 2L
    # missing stays NA_integer_
  }
  attr(raw, "alleles") <- alleles
  raw
}

#' Compute per-SNP flip flags from a raw 0/1/2 dosage matrix.
#'
#' `flip[i] = TRUE` when allele-2 (raw 2) is more frequent than allele-1
#' (raw 0).
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
  method <- match.arg(method, c("frequency", "reference"))
  if (method == "reference") {
    if (is.null(ref)) {
      stop(
        'method = "reference" requires the ref_allele argument.',
        call. = FALSE
      )
    }
    if (is.null(allele1) || length(allele1) != nrow(raw) ||
        length(ref) != nrow(raw) || anyNA(allele1) || anyNA(ref)) {
      stop("`allele1` and `ref_allele` must be complete vectors with one ",
           "entry per marker.", call. = FALSE)
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
