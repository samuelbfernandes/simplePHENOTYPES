#' Load marker data for create_phenotypes().
#'
#' Returns list(GT, GD, GI, input_format, out_name, temp) where GD is the
#' numeric genotype matrix (samples × SNPs) used by the simulation engine.
#' The return contract is identical to v1.3.0 — this function is internal.
#'
#' @keywords internal
#' @param geno_obj   In-memory data.frame (HapMap or numeric).
#' @param geno_file  Path to a single genotype file.
#' @param geno_path  Directory path containing multiple HapMap/VCF/BED/PED/GDS files.
#' @param nrows      Maximum rows to read (Inf = all).
#' @param na_string  String representing missing values.
#' @param prefix     File-name prefix filter for geno_path.
#' @param SNP_impute Imputation method: "Middle", "Minor", "Major", "None".
#' @param SNP_effect Genetic model: "Add", "Dom", "Left", "Right".
#' @param verbose    Print progress messages.
#' @param chr_prefix Chromosome prefix string for VCF (passed to SNPRelate).
#' @return list(GT, GD, GI, input_format, out_name, temp).
#' @author Samuel Fernandes
file_loader <- function(geno_obj   = NULL,
                        geno_file  = NULL,
                        geno_path  = NULL,
                        nrows      = Inf,
                        na_string  = "NA",
                        prefix     = NULL,
                        SNP_impute = "Middle",
                        SNP_effect = "Add",
                        verbose    = TRUE,
                        chr_prefix = "chr") {

  temp         <- tempfile(pattern = "", fileext = ".gds")
  input_format <- NULL
  out_name     <- NULL

  # ---------------------------------------------------------------------------
  # Helper: HapMap character matrix → numeric GD (samples × SNPs)
  # Uses parse_hapmap_chars_to_raw() + numericalize_core() Rust kernel.
  # ---------------------------------------------------------------------------
  .hapmap_to_GD <- function(G) {
    geno_chars <- as.matrix(G[, -(1:11)])
    if (is.numeric(geno_chars[1L, 1L])) {
      # Pre-coded numeric HapMap
      return(t(as.matrix(G[, -(1:11)])))
    }
    raw   <- parse_hapmap_chars_to_raw(geno_chars)   # SNPs × samples
    flip  <- compute_flip(raw)
    n_snp  <- nrow(raw)
    n_samp <- ncol(raw)
    coded <- numericalize_core(
      raw_dosage = as.vector(t(raw)),   # row-major: SNP outer, sample inner
      n_snp      = n_snp,
      n_samp     = n_samp,
      flip       = flip,
      code_as    = "-101",
      model      = SNP_effect,
      impute     = SNP_impute
    )
    # Rust output is row-major. Reading as (n_samp x n_snp) column-major gives
    # samples as rows and SNPs as columns — the GD convention.
    matrix(coded, nrow = n_samp, ncol = n_snp)
  }

  # ---------------------------------------------------------------------------
  # in-memory geno_obj
  # ---------------------------------------------------------------------------
  if (!is.null(geno_obj)) {
    fmt <- detect_format(geno_obj)

    if (fmt == "hapmap") {
      input_format <- "hapmap"
      if (verbose) message("File (geno_obj) loaded from memory.")
      if (verbose) message("Converting HapMap format to numerical under model of ", SNP_impute)
      GT <- as.matrix(colnames(geno_obj)[-(1:11)])
      colnames(GT) <- "taxa"
      GI <- geno_obj[, 1:4]
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      GD <- .hapmap_to_GD(geno_obj)
    } else if (fmt == "numeric") {
      input_format <- "numeric"
      if (verbose) message("File (geno_obj) loaded from memory.")
      GT <- as.matrix(colnames(geno_obj)[-(1:5)])
      colnames(GT) <- "taxa"
      GI <- geno_obj[, 1:4]
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      GD <- as.matrix(geno_obj[, -(1:5)])
      # Normalise 0/1/2 → -1/0/1 if needed
      probe <- unique(as.vector(GD[1:min(2, nrow(GD)), ]))
      if (!any(probe == -1, na.rm = TRUE) && any(probe == 2, na.rm = TRUE)) {
        GD <- GD - 1L
      }
      if (any(is.na(GD))) {
        GD[is.na(GD)] <- switch(SNP_impute,
                                 Middle = 0L, Minor = -1L, Major = 1L, 0L)
      }
    } else {
      stop("Format of 'geno_obj' was not recognised. ",
           "Provide a numeric or HapMap data.frame.", call. = FALSE)
    }
    return(list(GT = GT, GD = GD, GI = GI,
                input_format = input_format, out_name = out_name, temp = temp))
  }

  # ---------------------------------------------------------------------------
  # Single file: geno_file
  # ---------------------------------------------------------------------------
  if (is.null(geno_path)) {
    if (!file.exists(geno_file)) {
      stop("File ", geno_file, " not found.", call. = FALSE)
    }
    out_name     <- tools::file_path_sans_ext(basename(geno_file))
    out_name     <- gsub("\\.hmp$", "", out_name, ignore.case = TRUE)
    input_format <- detect_format(geno_file)

    if (input_format == "hapmap") {
      if (verbose) message("Performing numericalization")
      G <- data.table::fread(geno_file, header = TRUE,
                             nrows = nrows, na.strings = na_string,
                             data.table = FALSE)
      GT <- as.matrix(colnames(G)[-(1:11)])
      colnames(GT) <- "taxa"
      GI <- G[, 1:4]
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      GD <- .hapmap_to_GD(G)

    } else if (input_format == "numeric") {
      G <- data.table::fread(geno_file, header = TRUE,
                             nrows = nrows, na.strings = na_string,
                             data.table = FALSE)
      GT <- as.matrix(colnames(G)[-(1:5)])
      colnames(GT) <- "taxa"
      GI <- G[, 1:4]
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      GD <- as.matrix(G[, -(1:5)])
      probe <- unique(as.vector(GD[1:min(2, nrow(GD)), ]))
      if (!any(probe == -1, na.rm = TRUE) && any(probe == 2, na.rm = TRUE)) {
        GD <- GD - 1L
      }
      if (any(is.na(GD))) {
        GD[is.na(GD)] <- switch(SNP_impute,
                                 Middle = 0L, Minor = -1L, Major = 1L, 0L)
      }
      return(list(GT = GT, GD = GD, GI = GI,
                  input_format = input_format, out_name = out_name, temp = temp))

    } else if (input_format %in% c("vcf", "VCF")) {
      SNPRelate::snpgdsVCF2GDS(
        vcf.fn = geno_file, out.fn = temp,
        method = "biallelic.only", snpfirstdim = FALSE,
        verbose = FALSE, ignore.chr.prefix = chr_prefix)
      genofile <- SNPRelate::snpgdsOpen(temp)
      GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                     verbose = FALSE) - 1L
      GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
      colnames(GT) <- "taxa"
      GI <- data.frame(
        SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
        allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
        Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
        Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
        stringsAsFactors = FALSE)
      SNPRelate::snpgdsClose(genofile)

    } else if (input_format == "gds") {
      genofile <- SNPRelate::snpgdsOpen(geno_file)
      GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                     verbose = FALSE) - 1L
      GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
      colnames(GT) <- "taxa"
      snp_node <- if ("snp.rs.id" %in% gdsfmt::ls.gdsn(genofile))
        "snp.rs.id" else "snp.id"
      GI <- data.frame(
        SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, snp_node)),
        allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
        Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
        Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
        stringsAsFactors = FALSE)
      SNPRelate::snpgdsClose(genofile)

    } else if (input_format == "bed") {
      SNPRelate::snpgdsBED2GDS(
        bed.fn = geno_file,
        fam.fn = paste0(sub("\\.bed$", "", geno_file, ignore.case = TRUE), ".fam"),
        bim.fn = paste0(sub("\\.bed$", "", geno_file, ignore.case = TRUE), ".bim"),
        out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
      genofile <- SNPRelate::snpgdsOpen(temp)
      GD <- (SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                      verbose = FALSE) - 1L) * -1L
      GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
      colnames(GT) <- "taxa"
      raw_allele <- unlist(strsplit(
        gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")), "/"))
      l <- length(raw_allele)
      GI <- data.frame(
        SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
        allele     = paste0(raw_allele[seq(2, l, 2)], "/", raw_allele[seq(1, l, 2)]),
        Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
        Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
        stringsAsFactors = FALSE)
      SNPRelate::snpgdsClose(genofile)

    } else if (input_format == "ped") {
      SNPRelate::snpgdsPED2GDS(
        ped.fn = geno_file,
        map.fn = paste0(sub("\\.ped$", "", geno_file, ignore.case = TRUE), ".map"),
        out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
      genofile <- SNPRelate::snpgdsOpen(temp)
      GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                     verbose = FALSE) - 1L
      GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
      colnames(GT) <- "taxa"
      GI <- data.frame(
        SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
        allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
        Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
        Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
        stringsAsFactors = FALSE)
      SNPRelate::snpgdsClose(genofile)

    } else {
      stop("Format of '", geno_file, "' was not recognised. ",
           "Provide one of: HapMap, VCF, GDS, BED, or PED.", call. = FALSE)
    }

    # non-HapMap imputation
    if (input_format != "hapmap" && input_format != "numeric" && any(is.na(GD))) {
      GD[is.na(GD)] <- switch(SNP_impute,
                               Middle = 0L, Minor = -1L, Major = 1L, 0L)
    }
    return(list(GT = GT, GD = GD, GI = GI,
                input_format = input_format, out_name = out_name, temp = temp))
  }

  # ---------------------------------------------------------------------------
  # Directory: geno_path
  # ---------------------------------------------------------------------------
  files <- if (is.null(prefix)) {
    file.path(geno_path, dir(geno_path))
  } else {
    file.path(geno_path, dir(geno_path)[grepl(prefix, dir(geno_path))])
  }
  files <- sort(files)

  # Determine common prefix for out_name
  nn     <- strsplit(files[c(1L, length(files))], "")
  l_nn   <- lengths(nn)
  if (length(unique(l_nn)) > 1L) {
    pad <- max(l_nn) - min(l_nn)
    if (l_nn[1L] < l_nn[2L]) nn[[1L]] <- c(nn[[1L]], rep(".", pad)) else
      nn[[2L]] <- c(nn[[2L]], rep(".", pad))
  }
  do_nn  <- do.call("==", nn)
  nn_com <- match(FALSE, if (all(do_nn)) FALSE else do_nn) - 1L
  if (nn_com <= 0L) {
    out_name <- "out_geno"
  } else {
    out_name <- substr(files[1L], 1L, nn_com)
    if (length(files) > 2L) {
      nn2 <- strsplit(out_name, "")[[1L]]
      if (length(nn2) != length(nn[[1L]]))
        nn2 <- c(nn2, rep(".", length(nn[[1L]]) - length(nn2)))
      nn_com <- match(FALSE, nn[[1L]] == nn2) - 1L
      out_name <- substr(out_name, 1L, nn_com)
    }
    out_name <- basename(out_name)
  }

  input_format <- detect_format(files[1L])

  if (input_format == "hapmap") {
    if (verbose) { message("Reading HapMap files:"); message(files, sep = "\n") }
    G_list <- lapply(files, function(f) {
      data.table::fread(f, header = TRUE, nrows = nrows,
                        na.strings = na_string, data.table = FALSE)
    })
    G  <- do.call(rbind, G_list)
    GT <- as.matrix(colnames(G)[-(1:11)])
    colnames(GT) <- "taxa"
    GI <- G[, 1:4]; colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
    GD <- .hapmap_to_GD(G)

  } else if (input_format %in% c("vcf", "VCF")) {
    if (verbose) { message("Reading VCF files:"); message(files, sep = "\n") }
    SNPRelate::snpgdsVCF2GDS(vcf.fn = files, out.fn = temp,
                              method = "biallelic.only", snpfirstdim = FALSE,
                              verbose = FALSE, ignore.chr.prefix = chr_prefix)
    genofile <- SNPRelate::snpgdsOpen(temp)
    GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE, verbose = FALSE) - 1L
    GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
    colnames(GT) <- "taxa"
    GI <- data.frame(
      SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
      allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
      Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
      Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
      stringsAsFactors = FALSE)
    SNPRelate::snpgdsClose(genofile)

  } else if (input_format == "bed") {
    if (verbose) { message("Reading BED files:"); message(files, sep = "\n") }
    SNPRelate::snpgdsBED2GDS(
      bed.fn = files,
      fam.fn = paste0(sub("\\.bed$", "", files, ignore.case = TRUE), ".fam"),
      bim.fn = paste0(sub("\\.bed$", "", files, ignore.case = TRUE), ".bim"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    GD <- (SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE, verbose = FALSE) - 1L) * -1L
    GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
    colnames(GT) <- "taxa"
    raw_allele <- unlist(strsplit(
      gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")), "/"))
    l <- length(raw_allele)
    GI <- data.frame(
      SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
      allele     = paste0(raw_allele[seq(2, l, 2)], "/", raw_allele[seq(1, l, 2)]),
      Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
      Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
      stringsAsFactors = FALSE)
    SNPRelate::snpgdsClose(genofile)

  } else if (input_format == "ped") {
    if (verbose) { message("Reading PED files:"); message(files, sep = "\n") }
    SNPRelate::snpgdsPED2GDS(
      ped.fn = files,
      map.fn = paste0(sub("\\.ped$", "", files, ignore.case = TRUE), ".map"),
      out.gdsfn = temp, snpfirstdim = FALSE, verbose = FALSE)
    genofile <- SNPRelate::snpgdsOpen(temp)
    GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE, verbose = FALSE) - 1L
    GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
    colnames(GT) <- "taxa"
    GI <- data.frame(
      SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
      allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
      Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
      Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
      stringsAsFactors = FALSE)
    SNPRelate::snpgdsClose(genofile)

  } else if (input_format == "gds") {
    if (verbose) { message("Reading GDS files:"); message(files, sep = "\n") }
    genofile <- SNPRelate::snpgdsOpen(files)
    GD <- SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE, verbose = FALSE) - 1L
    GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
    colnames(GT) <- "taxa"
    snp_node <- if ("snp.rs.id" %in% gdsfmt::ls.gdsn(genofile)) "snp.rs.id" else "snp.id"
    GI <- data.frame(
      SNP        = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, snp_node)),
      allele     = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
      Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
      Position   = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
      stringsAsFactors = FALSE)
    SNPRelate::snpgdsClose(genofile)

  } else {
    stop("Format found in '", geno_path, "' was not recognised.", call. = FALSE)
  }

  if (input_format != "hapmap" && any(is.na(GD))) {
    GD[is.na(GD)] <- switch(SNP_impute, Middle = 0L, Minor = -1L, Major = 1L, 0L)
  }

  list(GT = GT, GD = GD, GI = GI,
       input_format = input_format, out_name = out_name, temp = temp)
}
