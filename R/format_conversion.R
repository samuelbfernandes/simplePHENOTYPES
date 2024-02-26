#' architectures.
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen  snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom data.table fwrite fread
#' @param file file
#' @param from = NULL,
#' @param to = NULL,
#' @param ref_allele = NULL,
#' @param to_r = NULL,
#' @param to_file = NULL,
#' @param file_name = NULL,
#' @param f_name = NULL,
#' @param code_as = "-101",
#' @param model = "Add",
#' @param impute = "None",
#' @param method = "frequency",
#' @param ... ...
#' @return Phenotypes for ntraits traits
#' @author Samuel Fernandes
#' Last update: Apr 20, 2021
#'
format_conversion <-
  function(file,
           from = NULL,
           to = NULL,
           ref_allele = NULL,
           to_r = NULL,
           to_file = NULL,
           file_name = NULL,
           f_name = NULL,
           code_as = "-101",
           hets = c("R", "Y", "S", "W", "K", "M", "AG", "CT", "CG", "AT", "GT", "AC"),
           homo = c("A", "AA", "T", "TT", "C", "CC", "G", "GG"),
           model = "Add",
           impute = "None",
           method = "frequency",
           verbose = TRUE,
           ...) {
    file_class <- class(file)
    if (!is.null(to_r) & !is.null(to_file)) {
      if (!to_r & !to_file) {
        warning(
          "to_r and to_file are both FALSE! Setting to_r = TRUE",
          call. = F,
          immediate. = T
        )
        to_r <- TRUE
      }
    } else if (is.null(to_file) & is.null(to_r)) {
      if (all(file_class == "character")) {
        to_file <- TRUE
        to_r <- FALSE
      } else {
        to_file <- FALSE
        to_r <- TRUE
      }
    } else if (is.null(to_file)) {
      if (to_r) {
        to_file <- FALSE
      } else{
        to_file <- TRUE
      }
    } else if (is.null(to_r)) {
      if (to_file) {
        to_r <- FALSE
      } else{
        to_r <- TRUE
      }
    }
    if (all(file_class == "character")) {
      file_upper <- toupper(file)
      if (!file.exists(file)) {
        stop(paste0("file \'", file, "\' not found!"), call. = F)
      }
      if (is.null(file_name) & to_file) {
        if (to == "numeric") {
          file_name <-  gsub(paste0(gsub(".*[.]", ".", file), "|.HMP.TXT"), "_numeric.txt", file, ignore.case = T)
        }
        
      }
    } else {
      if (is.null(file_name) & to_file) {
        if (to == "numeric") {
          file_name <- paste0(f_name, "_numeric.txt")
        }
      }
    }
    if (is.null(from)) {
      if (all(file_class == "character")) {
        if (endsWith(file_upper, ".HMP.TXT")) {
          from <- "hapmap"  
        } else if (endsWith(file_upper, ".GDS")) {
          from <- "gds"
        } else if (endsWith(file_upper, ".VCF")) {
          from <- "VCF"
        } else if (endsWith(file_upper, ".BED")) {
          from <- "bed"
        } else if (endsWith(file_upper, ".PED")) {
          from <- "ped"
        }
      }  else if (any(file_class == "gds.class")) {
        from <- "GDS"
      }  else if (any(file_class == "vcfR")) {
        from <- "vcfR"
      } else if (any(file_class == "data.frame") |
                 any(file_class == "matrix")) {
        hap_names <- c("rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center", "protLSID", "assayLSID", "panelLSID", "QCcode")
        nucleotide <- c(hets, homo, "+","++","-", "--", "0", "00")
        if_hap <- try((sum(colnames(file)[1:11] == hap_names) > 8 |
                         all(unlist(ifelse(is.character(file[, 12]), strsplit(unique(file[, 12]), ""), unique(file[, 12]))) %in% nucleotide) &
                         !all(unlist(ifelse(is.character(file[, 1]), strsplit(unique(file[, 1]), ""), unique(file[, 1]))) %in% nucleotide)),
                      silent = TRUE)
        if (class(if_hap) == "try-error"){
          stop("The original data format was not detected. Please provide it with the argument \"from\".",
               call. = F)
        } else if (if_hap) {
          from <- "hapmap" 
        } else if (all(unlist(ifelse(is.character(file[, 1]), strsplit(unique(file[, 1]), ""),unique(file[, 1]))) %in% nucleotide)) {
          from <- "table" 
        } else if (all(grepl("[/]|[|]", tail(file)[1:2, (ncol(file)-1):ncol(file)]))) {
          from <- "VCF"
        }
      }
    }
    if (from == "VCF" | from == "vcfR") {
      G <- handle_vcf(file,
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
                      verbose)
    } else if (from == "hapmap") {
      G <- handle_hapmap(file,
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
                         verbose)
    } else if (from == "table") {
      G <- handle_table(file,
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
                        verbose)
    } else {
      stop(
        "The format was not recognized automatically. Please set \"from\" to one of: \"hapmap\", \"VCF\", \"GDS\", \"BED\", \"PED\", or \"table\"",
        call. = F
      )
    }
    if (to_r) {
      return(G)
    }
  }
