# Format-specific handlers for as_numeric() / format_conversion().
#
# Every handler returns the same 5-column schema:
#   data.frame(snp, allele, chr, pos, cm, <sample columns>)
# with genotype values coded per code_as / model / impute.
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
  n_snp  <- nrow(raw_mat)
  n_samp <- ncol(raw_mat)

  flip <- compute_flip(raw_mat, method = method,
                       allele1 = allele1, ref = ref_allele)

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
  assemble_output(meta, sample_ids, coded_mat)
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

  if (to == "numeric") {
    meta <- data.frame(
      snp    = G[[1]],
      allele = G[[2]],
      chr    = G[[3]],
      pos    = G[[4]],
      cm     = NA_real_,
      stringsAsFactors = FALSE
    )
    sample_ids <- colnames(G)[-(1:11)]
    geno_chars <- as.matrix(G[, -(1:11)])

    if (is.numeric(geno_chars[1L, 1L])) {
      # Already numeric — just wrap
      G_out <- assemble_output(meta, sample_ids, geno_chars)
    } else {
      if (is.null(ref_allele)) {
        ref_allele_vec <- gsub("/.", "", G[[2]])
        rlang::inform(
          paste0('"ref_allele" was not provided. ',
                 'The first allele in the HapMap file will be used.'),
          .frequency = "once",
          .frequency_id = "hapmap_ref_allele_default"
        )
      } else {
        ref_allele_vec <- ref_allele
      }
      allele1 <- gsub("/.*", "", G[[2]])

      raw   <- parse_hapmap_chars_to_raw(geno_chars)
      G_out <- .apply_coding(
        raw_mat    = raw,
        meta       = meta,
        sample_ids = sample_ids,
        method     = method,
        ref_allele = if (method == "reference") ref_allele_vec else NULL,
        allele1    = allele1,
        code_as    = code_as,
        model      = model,
        impute     = impute
      )
    }

    if (to_file) {
      data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                         quote = FALSE, na = NA, showProgress = FALSE)
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
  }
  if (to_r) return(G_out)
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
    G_out <- table_to_numeric(
      G,
      code_as    = code_as,
      ref_allele = ref_allele,
      hets       = hets,
      homo       = homo,
      model      = model,
      impute     = impute,
      method     = method,
      verbose    = verbose
    )
    if (to_file) {
      suppressMessages(data.table::fwrite(
        G_out, file_name, row.names = FALSE, sep = "\t",
        quote = FALSE, na = NA, showProgress = FALSE, verbose = FALSE))
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# Shared helper: open GDS file → raw 0/1/2 matrix + metadata
# ---------------------------------------------------------------------------

.read_gds_to_raw <- function(genofile) {
  raw <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = TRUE, verbose = FALSE)
  mode(raw) <- "integer"

  sample_ids <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(genofile, "sample.id")))

  snp_node <- if ("snp.rs.id" %in% gdsfmt::ls.gdsn(genofile))
    "snp.rs.id" else "snp.id"
  snp_ids <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, snp_node))

  alleles <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele"))
  chr     <- as.character(gdsfmt::read.gdsn(
    gdsfmt::index.gdsn(genofile, "snp.chromosome")))
  pos     <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position"))

  meta <- data.frame(snp = snp_ids, allele = alleles,
                     chr = chr, pos = pos, cm = NA_real_,
                     stringsAsFactors = FALSE)
  allele1 <- sub("/.*", "", alleles)

  list(raw = raw, meta = meta, sample_ids = sample_ids, allele1 = allele1)
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
    temp <- paste0(gsub(".*/", "", tempfile()), ".gds")
    SNPRelate::snpgdsVCF2GDS(
      vcf.fn = file, out.fn = temp,
      method = "copy.num.of.ref", snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)

    if (to == "numeric") {
      parts <- .read_gds_to_raw(genofile)
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
      if (to_file) {
        data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                           quote = FALSE, na = NA, showProgress = FALSE)
        cat("\nNumeric file saved as '", file_name, "'\n")
      }
    }
    SNPRelate::snpgdsClose(genofile)
    if (to_r) return(G_out)
    return(invisible(NULL))
  }

  # In-memory vcfR / VCF data.frame
  if (from == "vcfR") {
    G          <- data.frame(file@gt[, colnames(file@gt) != "FORMAT"],
                             stringsAsFactors = FALSE)
    ref_allele <- file@fix[, "REF"]
  } else {
    G          <- file[, -1:-which(colnames(file) == "FORMAT")]
    ref_allele <- file$REF
  }
  if (!all(grepl("[/]|[|]", G[, 1]))) G <- G[, -1]
  if (!any(class(G) %in% "data.frame")) G <- data.frame(G, stringsAsFactors = TRUE)

  if (to == "numeric") {
    G_out <- table_to_numeric(
      G,
      code_as    = code_as,
      hets       = c("0/1", "0|1", "1/0", "1|0"),
      homo       = c("0/0", "0|0", "1/1", "1|1"),
      ref_allele = ref_allele,
      model      = model,
      impute     = impute,
      method     = method,
      verbose    = verbose
    )
    if (to_file) {
      suppressMessages(data.table::fwrite(
        G_out, file_name, row.names = FALSE, sep = "\t",
        quote = FALSE, na = NA, showProgress = FALSE, verbose = FALSE))
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# GDS handler
# ---------------------------------------------------------------------------

handle_gds <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (to == "numeric") {
    genofile <- SNPRelate::snpgdsOpen(file)
    parts    <- .read_gds_to_raw(genofile)
    SNPRelate::snpgdsClose(genofile)

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
    if (to_file) {
      data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                         quote = FALSE, na = NA, showProgress = FALSE)
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# BED handler
# ---------------------------------------------------------------------------

handle_bed <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (to == "numeric") {
    temp <- tempfile(fileext = ".gds")
    base <- sub("\\.bed$", "", file, ignore.case = TRUE)
    SNPRelate::snpgdsBED2GDS(
      bed.fn = file, fam.fn = paste0(base, ".fam"),
      bim.fn = paste0(base, ".bim"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    parts    <- .read_gds_to_raw(genofile)
    SNPRelate::snpgdsClose(genofile)

    # BED: SNPRelate stores A1/A2 where A2 is typically the major allele.
    # Swap the allele string so allele1 in meta = A2 (more common).
    allele_parts <- strsplit(parts$meta$allele, "/")
    parts$meta$allele <- vapply(allele_parts,
      function(x) if (length(x) == 2L) paste(x[2L], x[1L], sep = "/") else x[[1L]],
      character(1L))
    parts$allele1 <- sub("/.*", "", parts$meta$allele)

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
    if (to_file) {
      data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                         quote = FALSE, na = NA, showProgress = FALSE)
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
  }
  if (to_r) return(G_out)
}

# ---------------------------------------------------------------------------
# PED handler
# ---------------------------------------------------------------------------

handle_ped <- function(file, file_name, to_file, to_r, to,
                       code_as, model, impute, method, verbose) {
  if (to == "numeric") {
    temp <- tempfile(fileext = ".gds")
    base <- sub("\\.ped$", "", file, ignore.case = TRUE)
    SNPRelate::snpgdsPED2GDS(
      ped.fn = file, map.fn = paste0(base, ".map"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    parts    <- .read_gds_to_raw(genofile)
    SNPRelate::snpgdsClose(genofile)

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
    if (to_file) {
      data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                         quote = FALSE, na = NA, showProgress = FALSE)
      cat("\nNumeric file saved as '", file_name, "'\n")
    }
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
  lines <- readLines(con, n = 200L, warn = FALSE)
  close(con)
  data_line  <- which(grepl("^\\[Data\\]", lines, ignore.case = TRUE))
  skip_n     <- if (length(data_line) > 0L) data_line[[1L]] else 0L

  long <- data.table::fread(file, skip = skip_n, header = TRUE,
                             data.table = TRUE, sep = "\t")

  # Normalise column names for robust matching
  orig_names  <- names(long)
  names(long) <- tolower(gsub("[ -]", "_", orig_names))
  nms         <- names(long)

  snp_col    <- nms[grep("^snp_name$|^snp$",     nms)[[1L]]]
  sample_col <- nms[grep("^sample_id$|^sample$", nms)[[1L]]]
  a1_col     <- nms[grep("allele1", nms)[[1L]]]
  a2_col     <- nms[grep("allele2", nms)[[1L]]]

  chr_idx <- grep("^chr$|^chromosome$", nms)
  pos_idx <- grep("^position$|^pos$",   nms)
  chr_col <- if (length(chr_idx) > 0L) nms[[chr_idx[[1L]]]] else NA_character_
  pos_col <- if (length(pos_idx) > 0L) nms[[pos_idx[[1L]]]] else NA_character_

  if (any(is.na(c(snp_col, sample_col, a1_col, a2_col)))) {
    stop("FinalReport: cannot find required columns. Found: ",
         paste(orig_names, collapse = ", "), call. = FALSE)
  }

  long[, geno := paste0(get(a1_col), get(a2_col))]

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

  meta <- data.frame(
    snp    = snp_ids,
    allele = NA_character_,
    chr    = if (!is.null(chr_vec)) chr_vec else NA_character_,
    pos    = if (!is.null(pos_vec)) pos_vec else NA_integer_,
    cm     = NA_real_,
    stringsAsFactors = FALSE
  )

  raw   <- parse_hapmap_chars_to_raw(geno_chars)
  G_out <- .apply_coding(
    raw_mat    = raw,
    meta       = meta,
    sample_ids = sample_ids,
    method     = method,
    code_as    = code_as,
    model      = model,
    impute     = impute
  )

  if (to_file) {
    data.table::fwrite(G_out, file_name, row.names = FALSE, sep = "\t",
                       quote = FALSE, na = NA, showProgress = FALSE)
    cat("\nNumeric file saved as '", file_name, "'\n")
  }
  if (to_r) return(G_out)
}
