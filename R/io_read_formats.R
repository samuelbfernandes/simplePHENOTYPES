# Format-specific handlers for as_numeric() / format_conversion().
#
# Every handler returns the same 5-column schema:
#   data.frame(snp, allele, chr, pos, cm, <sample columns>)
# with genotype values coded per code_as / model / impute. The handlers only
# convert: format_conversion() adds the optional `counted` column and writes the
# file (.write_numeric()).
#
# The shared coding step is delegated to the Rust numericalize_core() kernel
# via .apply_coding().

# ---------------------------------------------------------------------------
# Shared coding helper
# ---------------------------------------------------------------------------

#' Apply numericalize_core() to a raw 0/1/2 matrix and wrap output.
#' @noRd
.apply_coding <- function(raw_mat, meta, sample_ids,
                          method     = "frequency",
                          ref_allele = NULL,
                          allele1    = NULL,
                          code_as    = "-101",
                          model      = "Add",
                          impute     = "None") {
  # The Rust kernel trusts its inputs (a bad shape or value can abort the R
  # process), so validate everything it relies on here, on the R side.
  if (!is.matrix(raw_mat) || !(is.numeric(raw_mat) || all(is.na(raw_mat)))) {
    stop("Internal error: the raw dosage must be a numeric matrix.",
         call. = FALSE)
  }
  bad_raw <- !is.na(raw_mat) & !(raw_mat %in% 0:2)
  if (any(bad_raw)) {
    stop("Internal error: raw dosage values must be 0, 1, 2 or NA.",
         call. = FALSE)
  }
  code_as <- match.arg(code_as, c("-101", "012"))
  model   <- match.arg(model, c("Add", "Dom", "Left", "Right"))
  impute  <- match.arg(impute, c("None", "Middle", "Minor", "Major"))
  n_snp  <- nrow(raw_mat)
  n_samp <- ncol(raw_mat)

  flip <- compute_flip(raw_mat, method = method,
                       allele1 = allele1, ref = ref_allele)
  if (length(flip) != n_snp || anyNA(flip)) {
    stop("Internal error: one non-missing orientation flag per marker is ",
         "required.", call. = FALSE)
  }

  # Rust expects row-major layout (SNP as outer dim, sample as inner).
  # R matrices are column-major, so transpose before flattening.
  coded_vec <- numericalize_core(
    raw_dosage = as.vector(t(raw_mat)),
    n_snp      = n_snp,
    n_samp     = n_samp,
    flip       = flip,
    code_as    = code_as,
    model      = model,
    impute     = impute
  )

  # Rust output is row-major; read back as n_snp x n_samp with byrow = TRUE.
  coded_mat <- matrix(coded_vec, nrow = n_snp, ncol = n_samp, byrow = TRUE)
  out <- assemble_output(meta, sample_ids, coded_mat)
  # Record which allele is the counted (+1) one. Raw 0 is allele 1, so the
  # counted allele is allele 1 unless the marker was flipped. Under the
  # dominance model ("Dom") both homozygotes share a code, so no allele is
  # counted; every other model keeps the additive homozygote coding.
  if (!identical(model, "Dom")) {
    attr(out, "counted_allele") <- .counted_allele(meta$allele, flip)
  }
  out
}

#' The allele that the numeric coding counts (+1, or 2 under `code_as = "012"`)
#'
#' The `allele` label is `"allele1/allele2"` in the raw orientation (raw 0 =
#' homozygous allele 1); a flipped marker counts allele 2. A marker whose label
#' does not name the counted allele (a single observed allele that must be
#' flipped, a missing label) is `NA`, meaning "not recorded". At a multiallelic
#' VCF site the label is `REF/ALT1,ALT2`; the coding only ever uses the first ALT.
#' @param label character vector of `allele` labels.
#' @param flip logical vector, one flag per marker (`TRUE` = allele 2 counted).
#' @noRd
.counted_allele <- function(label, flip) {
  label <- as.character(label)
  a1 <- ifelse(is.na(label), NA_character_, sub("/.*$", "", label))
  a2 <- ifelse(is.na(label) | !grepl("/", label, fixed = TRUE), NA_character_,
               sub(",.*$", "", sub("^[^/]*/", "", label)))
  a1[!is.na(a1) & !nzchar(a1)] <- NA_character_
  a2[!is.na(a2) & !nzchar(a2)] <- NA_character_
  toupper(ifelse(flip, a2, a1))
}

#' Is the sixth column of a numeric-format table the persisted `counted` record?
#'
#' `as_numeric(counted_column = TRUE)` writes the allele counted as `+1` (see
#' `.counted_allele()`) as a character column named `counted` placed immediately
#' after `cm`. A table has it when the sixth column is named `counted`
#' (case-insensitively, like the five fixed metadata names) and is not numeric:
#' a numeric column of that name is an individual called "counted" (genotype
#' columns are always numeric), so tables without the record read as before.
#' A logical column is accepted because a text round trip turns an all-`NA`
#' character column into a logical one.
#' @param df a data frame in numeric format.
#' @return `TRUE` or `FALSE`.
#' @noRd
.has_counted_col <- function(df) {
  is.data.frame(df) && ncol(df) >= 6L &&
    identical(tolower(names(df)[6L]), "counted") &&
    (is.character(df[[6L]]) || is.factor(df[[6L]]) || is.logical(df[[6L]]))
}

#' Number of leading metadata columns of a numeric-format table (5, or 6 with
#' the optional `counted` column).
#' @noRd
.n_meta <- function(df) if (.has_counted_col(df)) 6L else 5L

#' The persisted counted-allele column as a plain vector for validation
#'
#' An all-`NA` logical column (what a text round trip makes of an all-unknown
#' character column) is read as unknown; every other type is returned unchanged
#' so that `.check_counted()` can reject it.
#' @noRd
.counted_col_values <- function(x) {
  if (is.logical(x) && all(is.na(x))) rep(NA_character_, length(x)) else x
}

#' Insert the `counted` column (after `cm`) from the `"counted_allele"` record
#' @param G numeric-format data frame.
#' @return `G` with a character `counted` column as column 6 (unchanged when it
#'   already has one); the `"counted_allele"` attribute is kept.
#' @noRd
.add_counted_column <- function(G) {
  if (.has_counted_col(G)) return(G)
  cnt <- attr(G, "counted_allele", exact = TRUE)
  if (is.null(cnt)) {
    stop("`counted_column = TRUE` needs a record of the counted allele, and ",
         "this result has none (it is absent under model = \"Dom\", and for ",
         "numeric input that carries neither a `counted` column nor the ",
         "\"counted_allele\" attribute).", call. = FALSE)
  }
  cnt <- as.character(cnt)
  if (length(cnt) != nrow(G)) {
    stop("`counted_column = TRUE`: the \"counted_allele\" attribute has ",
         length(cnt), " entries for ", nrow(G), " markers.", call. = FALSE)
  }
  cols <- c(as.list(G)[1:5], list(counted = cnt), as.list(G)[-(1:5)])
  out <- data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
  attr(out, "counted_allele") <- cnt
  out
}

#' Write a numeric-format table to a text file
#'
#' The single place a converted table is written. When the file name was
#' generated (not supplied), an existing file of that name is reported with a
#' warning just before it is overwritten, so a conversion that fails earlier
#' never warns.
#' @noRd
.write_numeric <- function(G, file_name, default_name, verbose) {
  if (default_name && file.exists(file_name)) {
    warning("default output file ", file_name, " already exists and is ",
            "overwritten; pass `file_name` (or the explicit argument) to ",
            "choose another name", call. = FALSE)
  }
  suppressMessages(data.table::fwrite(
    G, file_name, row.names = FALSE, sep = "\t",
    quote = FALSE, na = NA, showProgress = FALSE, verbose = FALSE))
  if (verbose) message("Numeric file saved as '", file_name, "'.")
  invisible(file_name)
}

#' Allele label "A1/A2" from an n x 2 letter matrix; a marker whose second
#' allele was never observed is labelled by its single observed letter, and a
#' marker with no usable alleles falls back to `fallback`.
#' @noRd
.allele_label <- function(att, fallback = NA_character_) {
  fallback <- rep_len(as.character(fallback), nrow(att))
  ifelse(!is.na(att[, 1L]) & !is.na(att[, 2L]),
         paste(att[, 1L], att[, 2L], sep = "/"),
         ifelse(!is.na(att[, 1L]), att[, 1L], fallback))
}

#' Physical positions as integer when they are whole and fit in 32 bits, else
#' as double, so the `pos` type does not depend on the input format.
#' @noRd
.as_pos <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x) | (x == round(x) & abs(x) < .Machine$integer.max))) {
    as.integer(x)
  } else {
    x
  }
}

# ---------------------------------------------------------------------------
# HapMap handler (file path or in-memory data.frame)
# ---------------------------------------------------------------------------

handle_hapmap <- function(file,
                          file_class,
                          file_name,
                          to_file,
                          to_r,
                          to,
                          code_as,
                          ref_allele,
                          model,
                          impute,
                          method,
                          verbose) {
  if (all(file_class == "character")) {
    G <- try(data.table::fread(file, header = TRUE, data.table = FALSE),
             silent = TRUE)
  } else {
    G <- file
  }
  # A matrix is read by column below (`G[[j]]`), so hold it as a data frame.
  if (is.matrix(G)) G <- as.data.frame(G, stringsAsFactors = FALSE)

  if (is.data.frame(G) && ncol(G) < 12L) {
    stop("A HapMap table needs the 11 metadata columns followed by at least ",
         "one sample column; this one has ", ncol(G), " column(s).",
         call. = FALSE)
  }

  if (to == "numeric") {
    sample_ids <- colnames(G)[-(1:11)]
    geno_chars <- as.matrix(G[, -(1:11)])

    if (is.numeric(geno_chars[1L, 1L])) {
      stop("HapMap genotype columns are already numeric, so their allele ",
           "orientation cannot be inferred safely. Supply a five-metadata-",
           "column numeric-format object instead.", call. = FALSE)
    } else {
      declared <- as.character(G[[2]])
      allele1 <- toupper(sub("/.*", "", declared))
      allele2 <- toupper(sub(".*/", "", declared))
      ref_allele_vec <- if (identical(method, "reference")) {
        toupper(as.character(ref_allele))
      } else {
        NULL
      }
      if (identical(method, "reference") &&
          (length(ref_allele_vec) != nrow(G) || anyNA(ref_allele_vec))) {
        stop("`ref_allele` must give allele 1 or allele 2 for every HapMap ",
             "marker.", call. = FALSE)
      }

      # The `alleles` column is only a declaration: it is validated against the
      # observed calls, and a marker whose calls carry other alleles is oriented
      # from the calls themselves (with a warning), never from the column.
      raw <- parse_hapmap_chars_to_raw(geno_chars, allele1 = allele1,
                                       allele2 = allele2)
      att <- attr(raw, "alleles")
      eff1 <- ifelse(is.na(att[, 1L]), allele1, att[, 1L])
      eff2 <- ifelse(is.na(att[, 2L]), allele2, att[, 2L])
      if (identical(method, "reference")) {
        ok <- (ref_allele_vec == eff1) %in% TRUE | (ref_allele_vec == eff2) %in% TRUE
        if (!all(ok)) {
          stop("`ref_allele` must give allele 1 or allele 2 for every HapMap ",
               "marker (checked against the alleles observed in the calls).",
               call. = FALSE)
        }
      }
      meta <- data.frame(
        snp    = G[[1]],
        allele = .allele_label(att, fallback = declared),
        chr    = as.character(G[[3]]),
        pos    = .as_pos(G[[4]]),
        cm     = NA_real_,
        stringsAsFactors = FALSE
      )
      G_out <- .apply_coding(
        raw_mat    = raw,
        meta       = meta,
        sample_ids = sample_ids,
        method     = method,
        ref_allele = ref_allele_vec,
        allele1    = eff1,
        code_as    = code_as,
        model      = model,
        impute     = impute
      )
    }

  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# Already-numeric simplePHENOTYPES input
# ---------------------------------------------------------------------------

#' Give an already-numeric table the schema every reader emits.
#' @noRd
.normalize_numeric_schema <- function(df) {
  for (j in c(1L, 2L, 3L)) df[[j]] <- as.character(df[[j]])
  df[[4L]] <- .as_pos(df[[4L]])
  df[[5L]] <- as.numeric(df[[5L]])
  # the optional `counted` record is character (an all-NA one reads back from
  # text as logical); it is validated by handle_numeric(), never coerced to
  # a dosage type
  k <- .n_meta(df)
  if (k == 6L && is.logical(df[[6L]]) && all(is.na(df[[6L]]))) {
    df[[6L]] <- rep(NA_character_, nrow(df))
  }
  for (j in seq_along(df)[-seq_len(k)]) {
    v <- df[[j]]
    if (is.logical(v)) {
      df[[j]] <- as.integer(v)
    } else if (is.double(v) && all(is.na(v) | (v == round(v) &
                                                 abs(v) < .Machine$integer.max))) {
      # whole-number dosage codes (-1/0/1, 0/1/2, with NA) are stored integer,
      # the same type every reader emits; values are unchanged
      df[[j]] <- as.integer(v)
    }
  }
  df
}

handle_numeric <- function(file, file_name, to_file, to_r, code_as, model,
                           impute, method, ref_allele, verbose) {
  from_file <- is.character(file)
  if (from_file) {
    file <- data.table::fread(file, data.table = FALSE, showProgress = verbose)
  }
  if (!is.data.frame(file)) {
    stop("Numeric-format input must be a data frame with five metadata ",
         "columns.", call. = FALSE)
  }
  if (ncol(file) < 6L ||
      !identical(tolower(names(file)[1:5]),
                 c("snp", "allele", "chr", "pos", "cm"))) {
    stop("Numeric-format input must start with columns snp, allele, chr, pos, ",
         "and cm.", call. = FALSE)
  }
  # Restore the schema every reader emits (snp, allele, chr character; pos
  # integer when whole; cm double; sample columns integer/numeric, never
  # logical), for a file AND for an object in memory. A text round trip loses
  # types (fread reads chromosome "1" as integer and an all-NA pos or cm column
  # as logical), and a hand-built data frame may carry any of them, so
  # as_numeric(write(as_numeric(x))) reproduces as_numeric(x).
  file <- .normalize_numeric_schema(file)
  if (model != "Add" || impute != "None" || method != "frequency" ||
      !is.null(ref_allele)) {
    stop("An already-numeric input cannot be re-oriented, imputed, or changed ",
         "to another genetic model; convert from the original allele-coded ",
         "data instead.", call. = FALSE)
  }
  k <- .n_meta(file)
  if (ncol(file) <= k) {
    stop("Numeric-format input needs at least one individual column after ",
         "the metadata columns.", call. = FALSE)
  }
  if (k == 6L) {
    # the persisted counted-allele record (as_numeric(counted_column = TRUE)):
    # validated like map$counted, and mirrored into the attribute every
    # consumer reads; the column is authoritative (a stale attribute is
    # replaced, never merged)
    cc <- .check_counted(file[[6L]], file[[2L]], nrow(file),
                         "the `counted` column")
    attr(file, "counted_allele") <- cc
  }
  .check_unique_ids(file[[1L]], names(file)[-seq_len(k)])
  values <- as.matrix(file[, -seq_len(k), drop = FALSE])
  allowed <- if (code_as == "-101") c(-1, 0, 1) else c(0, 1, 2)
  if (!is.numeric(values) || any(!is.na(values) & !values %in% allowed)) {
    stop("The genotype values do not match code_as = \"", code_as, "\".",
         call. = FALSE)
  }
  if (to_r) file else invisible(NULL)
}

# ---------------------------------------------------------------------------
# Generic table handler (non-HapMap nucleotide-coded tables)
# ---------------------------------------------------------------------------

handle_table <- function(file,
                         file_class,
                         file_name,
                         to_file,
                         to_r,
                         to,
                         code_as,
                         ref_allele,
                         hets,
                         homo,
                         model,
                         impute,
                         method,
                         verbose) {
  if (all(file_class == "character")) {
    G <- try(data.table::fread(file, header = TRUE, data.table = FALSE),
             silent = TRUE)
  } else {
    G <- file
  }

  if (to == "numeric") {
    # A generic nucleotide table is markers (rows) × samples (columns); it
    # carries no metadata columns. Route it through the same Rust kernel every
    # other format uses, so genetic models, reference orientation and the
    # five-column output schema are consistent.
    G <- as.data.frame(G, stringsAsFactors = FALSE)
    sample_ids <- colnames(G)
    geno_chars <- as.matrix(G)
    snp_ids <- rownames(G)
    if (is.null(snp_ids) || all(snp_ids == as.character(seq_len(nrow(G))))) {
      snp_ids <- paste0("snp", seq_len(nrow(G)))
    }

    if (identical(method, "reference") &&
        (length(ref_allele) != nrow(G) || anyNA(ref_allele))) {
      stop("`ref_allele` must give one reference allele per marker (row) for ",
           "a nucleotide table under method = \"reference\".", call. = FALSE)
    }

    raw     <- parse_hapmap_chars_to_raw(geno_chars, hets = hets, homo = homo,
                                         miss = .MISS)
    alleles <- attr(raw, "alleles")
    # The parser reports allele *letters* ("A"), derived from every call
    # (heterozygotes included), so the metadata never carries genotype strings.
    allele1_letter <- alleles[, 1L]
    if (identical(method, "reference")) {
      # Validate the reference against the alleles actually observed at each
      # marker, decoding every non-missing call to its component alleles (so a
      # heterozygote contributes both -- the "G" in "AG", and, for an IUPAC code,
      # both alleles of "R" = A/G, not the letter "R"). A raw strsplit() would
      # observe "R" and reject a valid A/G reference; the homozygote-only check it
      # replaced missed het-only second alleles entirely.
      ref_letter <- toupper(substr(ref_allele, 1L, 1L))
      obs_ok <- vapply(seq_len(nrow(geno_chars)), function(i) {
        calls <- as.character(geno_chars[i, ])
        calls <- calls[!is.na(calls) & !toupper(calls) %in% toupper(.MISS)]
        letters_i <- unique(unlist(lapply(calls, .call_to_letters)))
        length(letters_i) == 0L || ref_letter[[i]] %in% letters_i
      }, logical(1))
      if (!all(obs_ok)) {
        bad <- which(!obs_ok)
        stop("`ref_allele` is not among the alleles observed at marker(s) ",
             paste(utils::head(bad, 5), collapse = ", "),
             if (length(bad) > 5) ", ..." else "",
             ": each reference allele must be an allele present at its marker.",
             call. = FALSE)
      }
    }
    meta <- data.frame(
      snp    = snp_ids,
      allele = .allele_label(alleles),
      chr    = NA_character_,
      pos    = NA_integer_,
      cm     = NA_real_,
      stringsAsFactors = FALSE
    )

    # compute_flip() compares allele1 (a single allele *letter*) against the
    # reference, so the reference must be a letter too: a user who writes the
    # homozygote code "AA" means allele "A". Collapse it, or a two-character
    # ref_allele would never match allele1 and would flip every marker.
    ref_letter <- if (identical(method, "reference")) {
      toupper(substr(ref_allele, 1L, 1L))
    } else {
      NULL
    }
    G_out <- .apply_coding(
      raw_mat    = raw,
      meta       = meta,
      sample_ids = sample_ids,
      method     = method,
      ref_allele = ref_letter,
      allele1    = allele1_letter,
      code_as    = code_as,
      model      = model,
      impute     = impute
    )

  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# Shared helper: open GDS file → raw 0/1/2 matrix + metadata
# ---------------------------------------------------------------------------

#' Read a SNPRelate GDS into the package's raw contract.
#'
#' SNPRelate's `snpgdsGetGeno()` returns the **number of copies of the first
#' allele** listed in `snp.allele` (2 = homozygous for the first allele). The
#' package contract (`parse_hapmap_chars_to_raw()`, `numericalize_core()`) is the
#' opposite: raw 0 = homozygous for allele 1, raw 2 = homozygous for allele 2, and
#' `allele1` is the first-listed allele. The counts are therefore reflected
#' (`2 - n`) here, at the reader, so every reader delivers the same contract and
#' a VCF gives identical dosage columns from a file path and from a data frame,
#' including on tied (MAF = 0.5) markers.
#' @noRd
.read_gds_to_raw <- function(genofile) {
  raw <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = TRUE, verbose = FALSE)
  mode(raw) <- "integer"
  raw <- 2L - raw

  sample_ids <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(genofile, "sample.id")))

  snp_node <- if ("snp.rs.id" %in% gdsfmt::ls.gdsn(genofile))
    "snp.rs.id" else "snp.id"
  snp_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, snp_node))

  alleles <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele"))
  chr     <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(genofile, "snp.chromosome")))
  pos     <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))

  snp_ids <- .fill_missing_ids(snp_ids, chr, pos)
  meta <- data.frame(snp = snp_ids, allele = alleles,
                     chr = chr, pos = pos, cm = NA_real_,
                     stringsAsFactors = FALSE)
  allele1 <- sub("/.*", "", alleles)

  list(raw = raw, meta = meta, sample_ids = sample_ids, allele1 = allele1)
}

# ---------------------------------------------------------------------------
# VCF genotype-call validation (shared by the in-memory and file-path readers)
# ---------------------------------------------------------------------------

#' Classify a matrix of VCF GT calls into the package's raw 0/1/2 contract.
#'
#' Only complete, biallelic diploid calls are recognised: `0/0` and `0|0` are raw
#' 0 (REF homozygote = allele 1), `0/1`, `0|1`, `1/0`, `1|0` raw 1, `1/1`, `1|1`
#' raw 2, and the missing codes stay `NA`. Every other non-empty call is set to
#' `NA` and reported: a call naming an allele index of 2 or more is multiallelic;
#' a call with a missing allele on one side or a single allele (haploid) is
#' partial. The reader used for a file path and the reader used for an in-memory
#' table both go through this function, so both give the same result.
#' @param gt_mat character matrix of GT strings (the `:DP:GQ...` suffix already
#'   removed).
#' @return a list: `raw` (integer matrix), `multi` and `partial` (logical
#'   matrices flagging the calls that were dropped and why).
#' @noRd
.vcf_classify_gt <- function(gt_mat) {
  nr <- nrow(gt_mat)
  hets <- c("0/1", "0|1", "1/0", "1|0")
  miss <- c("./.", ".|.", ".", "", "./", "/.")
  is_het <- matrix(gt_mat %in% hets, nrow = nr)
  is_ref <- matrix(gt_mat %in% c("0/0", "0|0"), nrow = nr)
  is_alt <- matrix(gt_mat %in% c("1/1", "1|1"), nrow = nr)
  raw <- matrix(NA_integer_, nrow = nr, ncol = ncol(gt_mat))
  raw[is_ref] <- 0L      # REF homozygote anchors raw 0 (allele 1)
  raw[is_het] <- 1L
  raw[is_alt] <- 2L
  recognized <- is_ref | is_het | is_alt | matrix(gt_mat %in% miss, nrow = nr)
  dropped <- !recognized & !is.na(gt_mat) & nzchar(gt_mat)
  multi <- dropped & matrix(grepl("[2-9]|[0-9]{2}", gt_mat), nrow = nr)
  list(raw = raw, multi = multi, partial = dropped & !multi)
}

#' Warn, with counts, about the VCF calls that `.vcf_classify_gt()` dropped.
#' @noRd
.vcf_warn_dropped <- function(n_multi, n_partial) {
  if (n_multi > 0) {
    warning(n_multi, " VCF genotype call(s) reference alleles beyond the ",
            "first ALT (multiallelic, e.g. \"0/2\" or \"2/2\"); this parser ",
            "is biallelic, so those calls were set to missing. Split ",
            "multiallelic sites (e.g. `bcftools norm -m -`) before import.",
            call. = FALSE)
  }
  if (n_partial > 0) {
    warning(n_partial, " VCF genotype call(s) are partially missing or ",
            "haploid (e.g. \"./1\", \"0/.\", \"0\"); this parser needs ",
            "complete diploid calls, so those calls were set to missing.",
            call. = FALSE)
  }
  invisible(NULL)
}

#' Scan the genotype lines of a VCF file and locate every call to be dropped.
#'
#' SNPRelate converts a VCF file without validating the calls: haploid calls are
#' kept as one allele copy (which reverses the coding of the valid diploids at the
#' same marker) and calls at multiallelic sites are kept. This reads the GT field
#' of the file in chunks (plain text, `.gz` or `.bgz`) and applies the same
#' classification as the in-memory reader, so that the two agree.
#' @param path VCF file path.
#' @param chunk number of lines read at a time.
#' @return a list: `n_snp`, `n_samp` (as seen in the text), `idx` (two-column
#'   matrix of the row and column of every call to set to missing), `n_multi` and
#'   `n_partial`.
#' @noRd
.vcf_scan_calls <- function(path, chunk = 20000L) {
  con <- gzfile(path, "r")
  on.exit(close(con), add = TRUE)
  n_snp <- 0L
  n_samp <- NA_integer_
  header_seen <- FALSE
  idx <- list()
  n_multi <- 0L
  n_partial <- 0L
  repeat {
    lines <- readLines(con, n = chunk, warn = FALSE)
    if (!length(lines)) break
    if (!header_seen) {
      h <- grep("^#CHROM", lines)
      if (!length(h)) next
      header_seen <- TRUE
      n_samp <- length(strsplit(lines[h[1L]], "\t", fixed = TRUE)[[1L]]) - 9L
      lines <- lines[-seq_len(h[1L])]
    }
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines) || n_samp < 1L) next
    fields <- strsplit(lines, "\t", fixed = TRUE)
    gt <- vapply(fields, function(f) {
      g <- rep_len(NA_character_, n_samp)
      k <- min(length(f) - 9L, n_samp)
      if (k > 0L) g[seq_len(k)] <- f[9L + seq_len(k)]
      g
    }, character(n_samp))
    gt <- matrix(sub(":.*", "", gt), nrow = length(lines), ncol = n_samp,
                 byrow = TRUE)
    cl <- .vcf_classify_gt(gt)
    bad <- which(cl$multi | cl$partial, arr.ind = TRUE)
    if (nrow(bad)) {
      bad[, 1L] <- bad[, 1L] + n_snp
      idx[[length(idx) + 1L]] <- bad
    }
    n_multi <- n_multi + sum(cl$multi)
    n_partial <- n_partial + sum(cl$partial)
    n_snp <- n_snp + length(lines)
  }
  list(n_snp = n_snp, n_samp = n_samp,
       idx = if (length(idx)) do.call(rbind, idx) else matrix(0L, 0L, 2L),
       n_multi = n_multi, n_partial = n_partial)
}

# ---------------------------------------------------------------------------
# VCF handler
# ---------------------------------------------------------------------------

handle_vcf <- function(file,
                       from,
                       file_class,
                       file_name,
                       to_file,
                       to_r,
                       to,
                       code_as,
                       model,
                       impute,
                       method,
                       verbose) {
  if (all(file_class == "character")) {
    if (!requireNamespace("SNPRelate", quietly = TRUE) ||
        !requireNamespace("gdsfmt", quietly = TRUE)) {
      stop(.gds_needed("VCF"), call. = FALSE)
    }
    # Keep the full tempfile() path: stripping the directory makes this
    # relative, which drops the intermediate GDS in the user's working
    # directory and leaves it behind.
    temp <- tempfile(fileext = ".gds")
    on.exit(unlink(temp), add = TRUE)
    SNPRelate::snpgdsVCF2GDS(
      vcf.fn = file, out.fn = temp,
      method = "copy.num.of.ref", snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    on.exit(try(SNPRelate::snpgdsClose(genofile), silent = TRUE),
            add = TRUE, after = FALSE)

    if (to == "numeric") {
      parts <- .read_gds_to_raw(genofile)
      # SNPRelate does not validate the calls, so apply the in-memory reader's
      # complete-diploid, biallelic rule to the file's own GT text: haploid,
      # partially missing and multiallelic calls become missing (with a counted
      # warning), before the orientation is chosen from the remaining calls.
      scan <- .vcf_scan_calls(file)
      if (scan$n_snp != nrow(parts$raw) || scan$n_samp != ncol(parts$raw)) {
        stop("The VCF file could not be read consistently (", scan$n_snp,
             " markers x ", scan$n_samp, " samples in the text, ",
             nrow(parts$raw), " x ", ncol(parts$raw), " after conversion); ",
             "check the file is a well-formed, tab-delimited VCF.",
             call. = FALSE)
      }
      if (nrow(scan$idx)) parts$raw[scan$idx] <- NA_integer_
      .vcf_warn_dropped(scan$n_multi, scan$n_partial)
      G_out <- .apply_coding(
        raw_mat    = parts$raw,
        meta       = parts$meta,
        sample_ids = parts$sample_ids,
        method     = method,
        allele1    = parts$allele1,
        code_as    = code_as,
        model      = model,
        impute     = impute
      )
    }
    SNPRelate::snpgdsClose(genofile)
    if (to_r) return(G_out)
    return(invisible(NULL))
  }

  # In-memory vcfR object or a plain VCF-style data.frame.
  if (inherits(file, "vcfR")) {
    fix        <- as.data.frame(file@fix, stringsAsFactors = FALSE)
    gt_cols    <- colnames(file@gt)[colnames(file@gt) != "FORMAT"]
    geno_raw   <- as.matrix(file@gt[, gt_cols, drop = FALSE])
    sample_ids <- gt_cols
    snp_ids    <- fix[["ID"]]
    chr_vec    <- fix[["CHROM"]]
    pos_vec    <- fix[["POS"]]
    ref_vec    <- fix[["REF"]]
    alt_vec    <- fix[["ALT"]]
  } else {
    G   <- as.data.frame(file, stringsAsFactors = FALSE)
    nms <- toupper(names(G))
    fmt_idx  <- match("FORMAT", nms)
    grab <- function(col) if (col %in% nms) G[[which(nms == col)[1L]]] else NULL
    snp_ids <- grab("ID"); chr_vec <- grab("#CHROM")
    if (is.null(chr_vec)) chr_vec <- grab("CHROM")
    pos_vec <- grab("POS"); ref_vec <- grab("REF"); alt_vec <- grab("ALT")
    # Sample columns are those after FORMAT when present, else every column
    # whose entries look like GT calls ("0/1", "1|1", possibly with :FORMAT).
    if (!is.na(fmt_idx)) {
      sample_ids <- names(G)[(fmt_idx + 1L):ncol(G)]
    } else {
      looks_gt <- vapply(G, function(col) {
        vals <- as.character(col)[!is.na(col) & nzchar(as.character(col))]
        # An empty column has no GT values; all(grepl(..., character(0))) is
        # vacuously TRUE, so require at least one real call before accepting it.
        length(vals) > 0L && all(grepl("^[.0-9]+[/|][.0-9]+", vals))
      }, logical(1))
      sample_ids <- names(G)[looks_gt]
    }
    if (!length(sample_ids)) {
      stop("No VCF genotype (GT) columns were found in the data frame.",
           call. = FALSE)
    }
    geno_raw <- as.matrix(G[, sample_ids, drop = FALSE])
  }

  if (to == "numeric") {
    # Keep only the GT field: strip any ":DP:GQ:..." suffix per call.
    gt_mat <- matrix(sub(":.*", "", as.character(geno_raw)),
                     nrow = nrow(geno_raw), ncol = ncol(geno_raw))
    # Only complete, biallelic diploid calls are used; anything else that is not
    # a missing code is set to missing and counted in a warning.
    cl <- .vcf_classify_gt(gt_mat)
    raw <- cl$raw
    .vcf_warn_dropped(sum(cl$multi), sum(cl$partial))

    n_snp <- nrow(raw)
    chr_out <- if (!is.null(chr_vec)) as.character(chr_vec) else rep(NA_character_, n_snp)
    pos_out <- if (!is.null(pos_vec)) .as_pos(pos_vec) else rep(NA_integer_, n_snp)
    meta <- data.frame(
      snp    = if (!is.null(snp_ids)) .fill_missing_ids(snp_ids, chr_out, pos_out) else
        paste0("snp", seq_len(n_snp)),
      allele = if (!is.null(ref_vec) && !is.null(alt_vec))
        paste(ref_vec, alt_vec, sep = "/") else NA_character_,
      chr    = chr_out,
      pos    = pos_out,
      cm     = NA_real_,
      stringsAsFactors = FALSE
    )

    G_out <- .apply_coding(
      raw_mat    = raw,
      meta       = meta,
      sample_ids = sample_ids,
      method     = "frequency",
      code_as    = code_as,
      model      = model,
      impute     = impute
    )
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# GDS handler
# ---------------------------------------------------------------------------

handle_gds <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (!requireNamespace("SNPRelate", quietly = TRUE) ||
      !requireNamespace("gdsfmt", quietly = TRUE)) {
    stop(.gds_needed("GDS"), call. = FALSE)
  }
  if (to == "numeric") {
    # `file` may be a path (open + close it here) or an already-open gds.class
    # object (use it directly; closing it is the caller's responsibility).
    if (inherits(file, "gds.class")) {
      parts <- .read_gds_to_raw(file)
    } else {
      genofile <- SNPRelate::snpgdsOpen(file)
      on.exit(SNPRelate::snpgdsClose(genofile), add = TRUE)
      parts <- .read_gds_to_raw(genofile)
    }

    G_out <- .apply_coding(
      raw_mat    = parts$raw,
      meta       = parts$meta,
      sample_ids = parts$sample_ids,
      method     = method,
      allele1    = parts$allele1,
      code_as    = code_as,
      model      = model,
      impute     = impute
    )
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# BED handler
# ---------------------------------------------------------------------------

#' Genetic distances from a PLINK .bim / .map file, when they are populated
#'
#' Column 3 of both formats is the genetic distance in centiMorgans. PLINK
#' writes 0 there when it is unknown, which is the common case, so an all-zero
#' column is reported as missing rather than as a map where every marker sits
#' at position 0.
#' @keywords internal
#' @noRd
.plink_cm <- function(path, n_markers) {
  if (!file.exists(path)) {
    return(rep(NA_real_, n_markers))
  }
  cm <- tryCatch({
    tab <- data.table::fread(path, header = FALSE, select = 3L,
                             data.table = FALSE, showProgress = FALSE)
    as.numeric(tab[[1L]])
  }, error = function(e) rep(NA_real_, n_markers))
  if (length(cm) != n_markers || all(!is.finite(cm)) ||
      all(cm[is.finite(cm)] == 0)) {
    return(rep(NA_real_, n_markers))
  }
  cm
}

handle_bed <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (!requireNamespace("SNPRelate", quietly = TRUE) ||
      !requireNamespace("gdsfmt", quietly = TRUE)) {
    stop(.gds_needed("PLINK BED"), call. = FALSE)
  }
  if (to == "numeric") {
    temp <- tempfile(fileext = ".gds")
    on.exit(unlink(temp), add = TRUE)
    base <- sub("\\.bed$", "", file, ignore.case = TRUE)
    SNPRelate::snpgdsBED2GDS(
      bed.fn = file, fam.fn = paste0(base, ".fam"),
      bim.fn = paste0(base, ".bim"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    parts    <- .read_gds_to_raw(genofile)
    SNPRelate::snpgdsClose(genofile)

    # BED: SNPRelate lists the .bim alleles as A1/A2 and counts copies of A1;
    # .read_gds_to_raw() has already reflected the counts, so allele 1 = A1 (the
    # first-listed allele), exactly as for every other reader.
    parts$meta$cm <- .plink_cm(paste0(base, ".bim"), nrow(parts$meta))

    G_out <- .apply_coding(
      raw_mat    = parts$raw,
      meta       = parts$meta,
      sample_ids = parts$sample_ids,
      method     = method,
      allele1    = parts$allele1,
      code_as    = code_as,
      model      = model,
      impute     = impute
    )
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# PED handler
# ---------------------------------------------------------------------------

handle_ped <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (!requireNamespace("SNPRelate", quietly = TRUE) ||
      !requireNamespace("gdsfmt", quietly = TRUE)) {
    stop(.gds_needed("PLINK PED"), call. = FALSE)
  }
  if (to == "numeric") {
    temp <- tempfile(fileext = ".gds")
    on.exit(unlink(temp), add = TRUE)
    base <- sub("\\.ped$", "", file, ignore.case = TRUE)
    SNPRelate::snpgdsPED2GDS(
      ped.fn = file, map.fn = paste0(base, ".map"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    parts    <- .read_gds_to_raw(genofile)
    SNPRelate::snpgdsClose(genofile)
    parts$meta$cm <- .plink_cm(paste0(base, ".map"), nrow(parts$meta))

    G_out <- .apply_coding(
      raw_mat    = parts$raw,
      meta       = parts$meta,
      sample_ids = parts$sample_ids,
      method     = method,
      allele1    = parts$allele1,
      code_as    = code_as,
      model      = model,
      impute     = impute
    )
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# Illumina FinalReport handler
#
# FinalReport is a long-format file (one row per sample × SNP) produced by
# Illumina GenomeStudio — the standard delivery format for livestock
# genotyping chips (BovineSNP50, PorcineSNP50, OvineSNP50, EquineSNP50).
#
# Expected structure after the [Data] section header:
#   SNP Name, Sample ID, Allele1-<conv>, Allele2-<conv>, ...
# where <conv> is Top, Forward, Design, or AB.
# Chr and Position columns are used when present.
# ---------------------------------------------------------------------------

handle_finalreport <- function(file, file_name, to_file, to_r, to,
                               code_as, model, impute, method, verbose) {
  if (to != "numeric") {
    if (to_r) return(invisible(NULL))
    return(invisible(NULL))
  }

  # Locate [Data] section
  con   <- file(file, "r")
  lines <- readLines(con, n = 2000L, warn = FALSE)
  close(con)
  data_line  <- which(grepl("^\\[Data\\]", lines, ignore.case = TRUE))
  skip_n     <- if (length(data_line) > 0L) data_line[[1L]] else 0L

  long <- data.table::fread(file, skip = skip_n, header = TRUE,
                             data.table = TRUE, sep = "\t")

  # Normalise column names for robust matching
  orig_names  <- names(long)
  names(long) <- tolower(gsub("[ -]", "_", orig_names))
  nms         <- names(long)

  # First column matching `pattern`, or NA when there is none (so the required-
  # column diagnostic below is reachable rather than a subscript error).
  find_col <- function(pattern) {
    i <- grep(pattern, nms)
    if (length(i) > 0L) nms[[i[[1L]]]] else NA_character_
  }
  snp_col    <- find_col("^snp_name$|^snp$")
  sample_col <- find_col("^sample_id$|^sample$")
  a1_col     <- find_col("allele1")
  a2_col     <- find_col("allele2")
  chr_col    <- find_col("^chr$|^chromosome$")
  pos_col    <- find_col("^position$|^pos$")

  if (any(is.na(c(snp_col, sample_col, a1_col, a2_col)))) {
    stop("FinalReport: cannot find required columns (a SNP name, a sample ID, ",
         "and Allele1/Allele2 columns). Found: ",
         paste(orig_names, collapse = ", "), call. = FALSE)
  }

  # An allele field that is NA (a literal "NA" in the file) or blank on either
  # side is a no-call: concatenating would otherwise give the string "NANA",
  # which the parser would take for a called homozygote, or turn "A" + blank
  # into a homozygote.
  al1 <- as.character(long[[a1_col]])
  al2 <- as.character(long[[a2_col]])
  no_call <- is.na(al1) | is.na(al2) | !nzchar(trimws(al1)) | !nzchar(trimws(al2))
  long[, geno := ifelse(no_call, NA_character_, paste0(al1, al2))]

  wide <- data.table::dcast(
    long,
    formula     = paste(snp_col, "~", sample_col),
    value.var   = "geno",
    fun.aggregate = function(x) x[[1L]]
  )

  snp_ids    <- wide[[snp_col]]
  sample_ids <- setdiff(names(wide), snp_col)
  geno_chars <- as.matrix(wide[, sample_ids, with = FALSE])

  chr_vec <- pos_vec <- NULL
  if (!is.na(chr_col) && !is.na(pos_col)) {
    cp  <- unique(long[, c(snp_col, chr_col, pos_col), with = FALSE])
    idx <- match(snp_ids, cp[[snp_col]])
    chr_vec <- as.character(cp[[chr_col]][idx])
    pos_vec <- as.integer(cp[[pos_col]][idx])
  }

  # FinalReport allele conventions include Illumina AB, whose heterozygote is
  # "AB"/"BA"; add those to the IUPAC/digraph set so AB calls are not mistaken
  # for a third homozygote and dropped as non-biallelic. Orientation is derived
  # from the calls (a FinalReport declares no allele pair).
  raw   <- parse_hapmap_chars_to_raw(geno_chars,
                                     hets = c(.HETS, "AB", "BA"))
  meta <- data.frame(
    snp    = snp_ids,
    allele = .allele_label(attr(raw, "alleles")),
    chr    = if (!is.null(chr_vec)) chr_vec else NA_character_,
    pos    = if (!is.null(pos_vec)) pos_vec else NA_integer_,
    cm     = NA_real_,
    stringsAsFactors = FALSE
  )
  G_out <- .apply_coding(
    raw_mat    = raw,
    meta       = meta,
    sample_ids = sample_ids,
    method     = method,
    code_as    = code_as,
    model      = model,
    impute     = impute
  )

  if (to_r) return(G_out)
}
