#' Load and numericalize marker data for create_phenotypes().
#'
#' Reads HapMap / VCF / GDS / PLINK BED-PED / FinalReport input and returns it
#' as a numeric-format data frame (`snp, allele, chr, pos, cm` + one column per
#' sample). Numericalization is delegated to [as_numeric()] /
#' [format_conversion()], the single coding implementation in the package, so
#' create_phenotypes() and the modern grammar code genotypes identically.
#' @keywords internal
#' @param geno_obj In-memory data.frame (e.g. HapMap).
#' @param geno_file Path to a single genotype file.
#' @param geno_path Directory containing multiple genotype files.
#' @param nrows Maximum rows to read from HapMap text (Inf = all).
#' @param na_string String representing missing values in HapMap text.
#' @param prefix File-name prefix filter for geno_path.
#' @param maf_cutoff Drop markers with minor-allele frequency below this.
#' @param SNP_effect Genetic model: "Add", "Dom", "Left", "Right".
#' @param SNP_impute Imputation: "Middle", "Minor", "Major", "None".
#' @param verbose Print progress messages.
#' @param chr_prefix Chromosome prefix string (accepted for backward
#'   compatibility; VCF chromosome handling is managed by SNPRelate).
#' @return list(geno_obj, input_format, out_name, temp).
#' @author Samuel Fernandes and Alexander Lipka.
#'
genotypes <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           nrows = Inf,
           na_string = "NA",
           prefix = NULL,
           maf_cutoff = NULL,
           SNP_effect = "Add",
           SNP_impute = "Middle",
           verbose = TRUE,
           chr_prefix = "chr") {
    #---------------------------------------------------------------------------
    # Numericalize a non-numeric input through the shared as_numeric() pipeline.
    numify <- function(input, from = NULL) {
      format_conversion(file = input, from = from, to = "numeric",
                        to_r = TRUE, to_file = FALSE,
                        model = SNP_effect, impute = SNP_impute,
                        verbose = verbose)
    }
    # Already-numeric input needs normalization (0/1/2 -> -1/0/1) and NA
    # imputation, NOT numericalization -- and the as_numeric numeric handler
    # deliberately refuses imputation/model changes -- so handle it directly,
    # matching create_phenotypes()'s own in-memory numeric path.
    numeric_df <- function(G) {
      G <- as.data.frame(G, check.names = FALSE, stringsAsFactors = FALSE)
      meta <- G[, 1:5, drop = FALSE]
      names(meta) <- c("snp", "allele", "chr", "pos", "cm")
      vals <- as.matrix(G[, -(1:5), drop = FALSE])
      probe <- unique(as.vector(vals))
      if (!any(probe == -1, na.rm = TRUE) && any(probe == 2, na.rm = TRUE)) {
        vals <- vals - 1L
      }
      if (any(is.na(vals))) {
        vals[is.na(vals)] <- switch(SNP_impute, Middle = 0L, Minor = -1L,
                                    Major = 1L, 0L)
      }
      cbind(meta, as.data.frame(vals, check.names = FALSE))
    }
    # HapMap text is read here so nrows / na_string still apply, then coded in
    # memory; binary formats are read straight from disk by SNPRelate.
    read_one <- function(f, fmt) {
      if (fmt == "numeric") {
        numeric_df(data.table::fread(f, header = TRUE, nrows = nrows,
                                     na.strings = na_string, data.table = FALSE))
      } else if (fmt == "hapmap") {
        G <- data.table::fread(f, header = TRUE, nrows = nrows,
                               na.strings = na_string, data.table = FALSE)
        numify(G, from = "hapmap")
      } else {
        numify(f, from = fmt)
      }
    }

    out_name <- NULL
    if (!is.null(geno_obj)) {
      input_format <- detect_format(geno_obj)
      df <- if (input_format == "numeric") numeric_df(geno_obj) else
        numify(geno_obj)
    } else if (!is.null(geno_file)) {
      input_format <- detect_format(geno_file)
      df <- read_one(geno_file, input_format)
      out_name <- gsub("\\.hmp$", "",
                       tools::file_path_sans_ext(basename(geno_file)),
                       ignore.case = TRUE)
    } else if (!is.null(geno_path)) {
      files <- if (is.null(prefix)) {
        file.path(geno_path, dir(geno_path))
      } else {
        file.path(geno_path, dir(geno_path)[grepl(prefix, dir(geno_path))])
      }
      files <- sort(files)
      input_format <- detect_format(files[[1L]])
      df <- do.call(rbind, lapply(files, read_one, fmt = input_format))
      out_name <- "out_geno"
    } else {
      stop("genotypes(): supply one of geno_obj, geno_file, or geno_path.",
           call. = FALSE)
    }

    # maf_cutoff: keep markers at or above the threshold. MAF is computed from
    # the -1/0/1 dosages: allele frequency p = mean((g + 1) / 2), maf = min(p, 1 - p).
    if (!is.null(maf_cutoff)) {
      vals <- as.matrix(df[, -(1:5), drop = FALSE])
      p    <- rowMeans((vals + 1) / 2, na.rm = TRUE)
      maf  <- pmin(p, 1 - p)
      df   <- df[maf >= maf_cutoff, , drop = FALSE]
    }

    list(
      geno_obj     = df,
      input_format = input_format,
      out_name     = out_name,
      temp         = tempfile(pattern = "", fileext = ".gds")
    )
  }
