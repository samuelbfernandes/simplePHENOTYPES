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
.MISS <- c("N","NN","NA","--","XX","00","+","++"," ","","-","0",".")

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
    # `.HET_IUPAC` is a named vector, so `[[` on a non-IUPAC letter (e.g. a
    # homozygote "A") raises "subscript out of bounds" rather than returning
    # NULL. Guard membership first so non-het single letters resolve to
    # character(0) (and .call_to_letters can fall back to splitting them).
    if (!code %in% names(.HET_IUPAC)) return(character(0))
    strsplit(.HET_IUPAC[[code]], "")[[1L]]
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
      if (.hmp_header_match(nms)) return("hapmap")
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

#' Does a header look like a HapMap header (at least 9 of the 11 standard names)?
#'
#' Case-insensitive. A HapMap table has the 11 standard metadata columns before
#' its first sample, and `handle_hapmap()` drops exactly those 11, so an object
#' with fewer than 12 columns (the 11 metadata columns plus at least one sample)
#' cannot be a usable HapMap table however its first names read: it is refused
#' here (rather than detected and then failing with a subscript error), from a
#' file path and from memory alike.
#' @noRd
.hmp_header_match <- function(nms) {
  n <- length(.HMP_NAMES)
  length(nms) > n &&
    sum(tolower(nms[seq_len(n)]) == tolower(.HMP_NAMES)) > 8L
}

.detect_df_format <- function(df) {
  nms <- names(df)
  if (.hmp_header_match(nms)) return("hapmap")
  # Header-less fallback: a nucleotide-call column in the first genotype
  # position. A numeric column (e.g. an all-zero column of a dosage matrix) is
  # never a column of nucleotide calls, and at least one real nucleotide letter
  # must be present: "0" alone (a legal missing code) is not evidence.
  if (length(nms) >= 12L && !is.numeric(df[[12L]]) && !is.logical(df[[12L]])) {
    vals <- unique(as.character(df[[12L]]))
    chars <- unlist(strsplit(vals, ""))
    nucleotide <- c(.HETS, .HOMO, "N", "-", "+", "0")
    if (all(chars %in% nucleotide) && any(chars %in% c("A", "C", "G", "T")))
      return("hapmap")
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
#' Generalized from the HapMap parser so it also serves nucleotide tables and
#' Illumina FinalReport calls: the heterozygote, homozygote and missing code
#' sets are all overridable. Calls and code sets are compared case-insensitively.
#'
#' Orientation is derived from the **observed** calls, never trusted from a
#' declared allele pair. When `allele1` (and optionally `allele2`) is supplied it
#' anchors raw 0 to that allele, but only if every allele actually observed at the
#' marker is one of the declared pair; a declaration that disagrees with the calls
#' is ignored for that marker (orientation is derived from the calls) and the
#' marker is reported in `attr(., "mismatch")` with a warning. With no usable
#' declaration, allele 1 is the more frequent homozygote (ties: alphabetically
#' first letter) or, for a heterozygote-only marker, the alphabetically first
#' letter.
#'
#' @param geno_mat Character matrix (SNPs x samples) of genotype calls.
#' @param allele1 optional declared allele-1 letter per marker (single letters
#'   only; anything else is treated as "no usable declaration").
#' @param hets heterozygote codes (default the IUPAC/digraph set `.HETS`).
#' @param homo optional explicit homozygote codes. When `NULL` (default) any
#'   present call that is neither het nor missing is treated as a homozygote
#'   (HapMap behaviour). When supplied, a call that is neither het, homo, nor
#'   missing is treated as missing (invalid) and counted in a warning.
#' @param miss missing-value codes (default `.MISS`).
#' @param allele2 optional declared allele-2 letter per marker, used with
#'   `allele1` to validate the calls against the declared pair.
#' @return Integer matrix (SNPs x samples): 0 = hom allele-1, 1 = het,
#'   2 = hom allele-2, NA_integer_ = missing. Attributes: `"alleles"`, an
#'   n_snp x 2 character matrix of the (allele1, allele2) **letters** backing the
#'   0 / 2 codes; `"nonbiallelic"` and `"mismatch"`, logical per-marker flags.
#' @noRd
parse_hapmap_chars_to_raw <- function(geno_mat, allele1 = NULL,
                                      hets = .HETS, homo = NULL,
                                      miss = .MISS, allele2 = NULL) {
  n_snp  <- nrow(geno_mat)
  n_samp <- ncol(geno_mat)
  raw    <- matrix(NA_integer_, nrow = n_snp, ncol = n_samp)
  alleles <- matrix(NA_character_, nrow = n_snp, ncol = 2L)
  nonbi    <- rep(FALSE, n_snp)
  mismatch <- rep(FALSE, n_snp)
  n_invalid <- 0L

  hets <- toupper(hets)
  miss <- toupper(miss)
  if (!is.null(homo)) homo <- toupper(homo)
  calls_up <- toupper(matrix(as.character(geno_mat), nrow = n_snp, ncol = n_samp))
  is_letter <- function(x) !is.na(x) & grepl("^[A-Z]$", x) & x != "N"
  a1 <- if (is.null(allele1)) NULL else toupper(as.character(allele1))
  a2 <- if (is.null(allele2)) NULL else toupper(as.character(allele2))

  for (i in seq_len(n_snp)) {
    row     <- calls_up[i, ]
    is_miss <- row %in% miss | is.na(row)
    is_het  <- !is_miss & row %in% hets
    if (is.null(homo)) {
      is_hom <- !is_miss & !is_het
    } else {
      is_hom <- !is_miss & !is_het & row %in% homo
      # A present call that is neither het nor a recognised homozygote is not a
      # valid biallelic genotype: record it as missing (and count it, so the
      # loss is reported) rather than as a phantom third homozygote.
      bad <- !is_miss & !is_het & !is_hom
      n_invalid <- n_invalid + sum(bad)
      is_miss <- is_miss | bad
    }
    if (all(is_miss)) next

    # Alleles observed across ALL present calls: a heterozygote contributes both
    # of its letters, so multiallelism hidden in hets (AA/AG/AT) is caught.
    u_calls <- unique(row[!is_miss])
    obs <- unique(unlist(lapply(u_calls, .call_to_letters)))
    hom_letter <- ifelse(is_hom, substr(row, 1L, 1L), NA_character_)
    hom_tab <- table(hom_letter[is_hom])
    if (length(hom_tab) > 2L || length(obs) > 2L) {
      nonbi[i] <- TRUE            # not biallelic: whole row stays missing
      next
    }

    # Orientation: anchor to the declared pair only when it is usable and
    # consistent with what is actually observed at this marker.
    pair <- NULL
    if (!is.null(a1) && is_letter(a1[[i]])) {
      if (!is.null(a2) && is_letter(a2[[i]]) && a1[[i]] != a2[[i]]) {
        if (all(obs %in% c(a1[[i]], a2[[i]]))) {
          pair <- c(a1[[i]], a2[[i]])
        } else {
          mismatch[i] <- TRUE
        }
      } else if (a1[[i]] %in% obs) {
        other <- setdiff(obs, a1[[i]])
        pair <- c(a1[[i]], if (length(other)) other[[1L]] else NA_character_)
      }
    }
    if (is.null(pair)) {
      first <- if (length(hom_tab)) {
        names(hom_tab)[order(-as.integer(hom_tab), names(hom_tab))][[1L]]
      } else {
        sort(obs)[[1L]]
      }
      other <- setdiff(obs, first)
      pair <- c(first, if (length(other)) other[[1L]] else NA_character_)
    }

    alleles[i, ] <- pair
    raw[i, is_het] <- 1L
    raw[i, is_hom & hom_letter == pair[[1L]]] <- 0L
    raw[i, is_hom & hom_letter != pair[[1L]]] <- 2L
    # missing stays NA_integer_
  }

  .rows <- function(x) {
    w <- which(x)
    paste0(paste(utils::head(w, 5L), collapse = ", "),
           if (length(w) > 5L) ", ..." else "")
  }
  if (any(nonbi)) {
    warning(sum(nonbi), " marker(s) are not biallelic (rows ", .rows(nonbi),
            ") and were set to missing.", call. = FALSE)
  }
  if (any(mismatch)) {
    warning(sum(mismatch), " marker(s) carry alleles that are not in the ",
            "declared allele pair (rows ", .rows(mismatch), "); the declared ",
            "pair was ignored for those markers and their orientation was ",
            "derived from the observed calls.", call. = FALSE)
  }
  if (n_invalid > 0L) {
    warning(n_invalid, " genotype call(s) are not among the recognised ",
            "heterozygote, homozygote or missing codes and were set to ",
            "missing.", call. = FALSE)
  }
  attr(raw, "alleles") <- alleles
  attr(raw, "nonbiallelic") <- nonbi
  attr(raw, "mismatch") <- mismatch
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
  .check_unique_ids(meta$snp, sample_ids)
  out <- data.frame(meta, coded_mat,
                    check.names = FALSE, fix.empty.names = FALSE,
                    stringsAsFactors = FALSE)
  colnames(out) <- c("snp", "allele", "chr", "pos", "cm", sample_ids)
  out
}

#' Refuse duplicated marker IDs and duplicated sample names.
#'
#' Downstream code (`qtn_table()`, `g_matrix()`, `as_population()`) indexes
#' markers and individuals by name, so a repeated name silently addresses the
#' wrong column. Missing marker IDs (`NA`) are not compared.
#' @noRd
.check_unique_ids <- function(snp, samples) {
  snp <- snp[!is.na(snp)]
  dup_snp <- unique(snp[duplicated(snp)])
  if (length(dup_snp)) {
    stop("Duplicated marker ID(s): ",
         paste(utils::head(dup_snp, 5L), collapse = ", "),
         if (length(dup_snp) > 5L) ", ..." else "",
         ". Marker IDs must be unique; rename or remove the duplicates.",
         call. = FALSE)
  }
  dup_smp <- unique(samples[duplicated(samples)])
  if (length(dup_smp)) {
    stop("Duplicated sample name(s): ",
         paste(utils::head(dup_smp, 5L), collapse = ", "),
         if (length(dup_smp) > 5L) ", ..." else "",
         ". Sample names must be unique; rename or remove the duplicates.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Fill missing VCF/PLINK marker IDs ("." / "" / NA) with `chr:pos`.
#'
#' A VCF without an ID column value writes ".", which repeats for every such
#' site; using it as a marker name would make all of them duplicates. The
#' `chr:pos` convention (as bcftools uses) keeps them unique and traceable.
#' @noRd
.fill_missing_ids <- function(ids, chr, pos) {
  ids <- as.character(ids)
  bad <- is.na(ids) | !nzchar(ids) | ids == "."
  if (any(bad)) ids[bad] <- paste0(as.character(chr)[bad], ":", as.character(pos)[bad])
  ids
}
