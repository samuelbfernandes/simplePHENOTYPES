#' Function used by create_phenotypes to load marker dataset.
#' @param geno_obj = NULL,
#' @param geno_file = NULL,
#' @param geno_path = NULL,
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param prefix = NULL,
#' @param SNP_impute = "Middle",
#' @keywords internal
#' @param SNP_effect = "Add",
#' @param verbose = TRUE
#' @param chr_prefix = "chrm"
#' @return genotype data sampled
#' @author  Samuel Fernandes
#' Last update: Apr 20, 2020
#'------------------------------------------------------------------------------
file_loader <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           nrows = Inf,
           na_string = "NA",
           prefix = NULL,
           SNP_impute = "Middle",
           SNP_effect = "Add",
           verbose = TRUE,
           chr_prefix = "chr") {
    #'--------------------------------------------------------------------------
    hap_names <- c(
      "rs#",
      "alleles",
      "chrom",
      "pos",
      "strand",
      "assembly#",
      "center",
      "protLSID",
      "assayLSID",
      "panelLSID",
      "QCcode"
    )
    nucleotide <-
      c("A", "C", "T", "G", "R", "Y", "S", "W", "K", "M", "N", "+", "-", "0")
    temp <- tempfile(pattern = "", fileext = ".gds")
    input_format <- NULL
    out_name <- NULL
    if (!is.null(geno_obj)) {
      if (sum(colnames(geno_obj)[1:11] == hap_names) > 8 |
          all(unlist(ifelse(
            is.character(geno_obj[, 12]),
            strsplit(unique(geno_obj[, 12]), ""),
            unique(geno_obj[, 12])
          )) %in% nucleotide)) {
        input_format <- "hapmap"
      } else if (any(grepl("[/]|[|]", tail(geno_obj)[1, - (1:5)]))) {
        stop("Please read VCF files as \'geno_file\' or \'geno_path\'.",
             call. = F)
      } else {
        stop(
          "File format provied by \'geno_obj\' was not recognized! Please provied one of: numeric or HapMap.",
          call. = F
        )
      }
      if (verbose)
        message("File (geno_obj) loaded from memory.")
      if (input_format == "hapmap") {
        GT <- as.matrix(colnames(geno_obj)[- (1:11)])
        GI <- geno_obj[, c(1, 2, 3, 4)]
        if (verbose)
          message(paste0(
            "Converting HapMap format to numerical under model of ",
            SNP_impute
          ))
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        GD <- NULL
        bit <- nchar(as.character(geno_obj[2, 12]))
        if (verbose)
          message("Performing numericalization")
        GD <- apply(geno_obj[, - (1:11)], 1, function(one)
          numericalization(
            one,
            bit = bit,
            effect = SNP_effect,
            impute = SNP_impute
          ))
        if (is.null(GD)) {
          GT <- NULL
          GI <- NULL
        }
        return(list(
          GT = GT,
          GD = GD,
          GI = GI,
          input_format = input_format,
          out_name = out_name
        ))
      }
    } else if (is.null(geno_path)) {
      if (!geno_file %in% dir() &
          (!sub(".*/", "", geno_file) %in% list.files(dirname(geno_file)))) {
        stop(paste("File ", geno_file, " not found."), call. = F)
      }
      out_name <- gsub(".vcf|.hmp.txt|.bed|.txt|.gds|.ped",
                       "",
                       sub(".*/", "", geno_file))
      if (grepl(".VCF$", toupper(geno_file))) {
        input_format <- "VCF"
      } else if (grepl(".TXT$", toupper(geno_file))) {
        data_type <- try(data.table::fread(
          file = geno_file,
          header = T,
          skip = 0,
          nrows = 1,
          na.strings = na_string,
          data.table = F
        ),
        silent = TRUE)
        if (sum(hap_names %in% colnames(data_type)) > 8 |
            all(unlist(ifelse(
              is.character(data_type[, 12]),
              strsplit(unique(data_type[, 12]), ""),
              "numeric"
            )) %in% nucleotide)) {
          input_format <- "hapmap"
        } else if (all(colnames(data_type)[1:5] == c("snp", "allele", "chr", "pos", "cm"))) {
          G <-
            try(data.table::fread(
              file = geno_file,
              header = TRUE,
              skip = 0,
              nrows = nrows,
              na.strings = na_string,
              data.table = F
            ),
            silent = TRUE)
          GT <- as.matrix(colnames(G)[- (1:5)])
          GI <- G[, c(1, 2, 3, 4)]
          colnames(GT) <- "taxa"
          colnames(GI) <-
            c("SNP", "allele", "Chromosome", "Position")
          GD <- G[, - c(1:5)]
          dose <- 0
          counter <- 1
          while (all(dose != 2) & all(dose != -1)) {
            dose <- unique(GD[, counter])
            counter <- counter + 1
          }
          if (all(dose != -1) | any(dose == 2)) {
            GD <- GD - 1
          }
          isna <- is.na(GD)
          if (any(isna)) {
            if (SNP_impute == "Middle") {
              GD[isna] <- 0
            } else
              if (SNP_impute == "Minor") {
                GD[isna] <- -1
              } else
                if (SNP_impute == "Major") {
                  GD[isna] <- 1
                }
          }
          return(list(
            GT = GT,
            GD = GD,
            GI = GI,
            input_format = input_format,
            out_name = out_name
          ))
        } else {
          stop(
            "File format provied by \'geno_file\' was not recognized! Please provied one of: Numeric, VCF, HapMap, gds, or plink bed or ped files.",
            call. = F
          )
        }
      }  else if (grepl(".GDS$", toupper(geno_file))) {
        input_format <- "gds"
      } else if (grepl(".BED$", toupper(geno_file))) {
        input_format <- "bed"
      } else if (grepl(".PED$", toupper(geno_file))) {
        input_format <- "ped"
      } else {
        data_type <- try(data.table::fread(
          file = geno_file,
          header = F,
          skip = 0,
          nrows = 1,
          na.strings = na_string,
          data.table = F
        ),
        silent = TRUE)
        if (any(grepl("VCF", data_type))) {
          input_format <- "VCF"
        } else if (sum(hap_names %in% data_type) > 8 |
                   all(unlist(ifelse(
                     is.character(data_type[, 12]),
                     strsplit(unique(data_type[, 12]), ""),
                     "numeric"
                   )) %in% nucleotide)) {
          input_format <- "hapmap"
        } else {
          stop(
            "File format provied by \'geno_file\' was not recognized! Please provied one of: Numeric, VCF, HapMap, gds, or plink bed or ped files.",
            call. = F
          )
        }
      }
      if (input_format == "hapmap") {
        G <-
          try(data.table::fread(
            file = geno_file,
            header = TRUE,
            skip = 0,
            nrows = nrows,
            na.strings = na_string,
            data.table = F
          ),
          silent = TRUE)
        GT <- as.matrix(colnames(G)[- (1:11)])
        GI <- G[, c(1, 2, 3, 4)]
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        GD <- NULL
        bit <- nchar(as.character(G[2, 12]))
        if (verbose)
          message("Performing numericalization")
        GD <- apply(G[, - (1:11)], 1, function(one)
          numericalization(
            one,
            bit = bit,
            effect = SNP_effect,
            impute = SNP_impute
          ))
      } else if (input_format == "VCF") {
        SNPRelate::snpgdsVCF2GDS(
          vcf.fn = geno_file,
          out.fn = temp,
          method = "biallelic.only",
          snpfirstdim = FALSE,
          verbose = FALSE,
          ignore.chr.prefix = chr_prefix
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                   verbose = FALSE) - 1
        GT <-  as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "bed") {
        SNPRelate::snpgdsBED2GDS(
          bed.fn = geno_file,
          fam.fn = paste0(gsub(".bed", "", geno_file), ".fam"),
          bim.fn = paste0(gsub(".bed", "", geno_file), ".bim"),
          out.gdsfn = temp,
          snpfirstdim = FALSE,
          verbose = FALSE
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          (SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                    verbose = FALSE) - 1) * -1
        GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        allele <- unlist(strsplit(gdsfmt::read.gdsn(
          gdsfmt::index.gdsn(genofile, "snp.allele")
        ), "/"))
        l <- length(allele)
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
          allele = paste0(allele[seq(2, l, 2)], "/", allele[seq(1, l, 2)]),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "ped") {
        SNPRelate::snpgdsPED2GDS(
          ped.fn = geno_file,
          map.fn = paste0(gsub(".ped", "", geno_file), ".map"),
          out.gdsfn = temp,
          snpfirstdim = FALSE,
          verbose = FALSE
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                   verbose = FALSE) - 1
        GT <-  as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "gds") {
        genofile <- SNPRelate::snpgdsOpen(geno_file)
        GD <-
          SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                   verbose = FALSE) - 1
        GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      }
    } else {
      if (is.null(prefix)) {
        files <- paste0(geno_path, "/", dir(geno_path))
      } else{
        files <- paste0(geno_path, "/",
                        dir(geno_path)[grepl(prefix,
                                             dir(geno_path))])
      }
      files <- sort(files)
      nn <- strsplit(files[c(1, length(files))], "")
      l_nn <- lengths(nn)
      if (length(unique(l_nn)) > 1) {
        if (l_nn[1] - l_nn[2] > 0) {
          nn[[2]] <- c(nn[[2]], rep(".", l_nn[1] - l_nn[2])) 
        } else {
          nn[[1]] <- c(nn[[1]], rep(".", l_nn[2] - l_nn[1]))
        }
      }
      do_nn <- do.call("==", nn)
      nn_com <- match(FALSE, if(all(do_nn)) { FALSE } else {do_nn}) - 1
      if (nn_com == 0) {
        out_name <- "out_geno"
      } else {
        out_name <- substr(files[1], 1, nn_com)
        if (length(files) > 2) {
          nn[[2]] <- unlist(strsplit(out_name, ""))
          if (length(nn[[2]]) != length(nn[[1]])) {
              nn[[2]] <- c(nn[[2]], rep(".", length(nn[[1]]) - length(nn[[2]]))) 
          }
          nn_com <- match(FALSE, do.call("==", nn)) - 1
          out_name <- substr(out_name, 1, nn_com)
        }
      }
      if (grepl(".VCF$", toupper(files[1]))) {
        input_format <- "VCF"
      } else if (grepl(".TXT$", toupper(files[1]))) {
        data_type <- try(data.table::fread(
          file = files[1],
          header = T,
          skip = 0,
          nrows = 1,
          na.strings = na_string,
          data.table = F
        ),
        silent = TRUE)
        if (sum(hap_names %in% colnames(data_type)) > 8 |
            all(unlist(ifelse(
              is.character(data_type[, 12]),
              strsplit(unique(data_type[, 12]), ""),
              "numeric"
            )) %in% nucleotide)) {
          input_format <- "hapmap"
        } else if (all(colnames(data_type)[1:5] == c("snp", "allele", "chr", "pos", "cm"))) {
          G <-
            try(data.table::fread(
              file = files[1],
              header = TRUE,
              skip = 0,
              nrows = nrows,
              na.strings = na_string,
              data.table = F
            ),
            silent = TRUE)
          GT <- as.matrix(colnames(G)[- (1:5)])
          GI <- G[, c(1, 2, 3, 4)]
          colnames(GT) <- "taxa"
          colnames(GI) <-
            c("SNP", "allele", "Chromosome", "Position")
          GD <- G[, - c(1:5)]
          return(list(
            GT = GT,
            GD = GD,
            GI = GI,
            input_format = input_format,
            out_name = out_name
          ))
        } else {
          stop(
            "File format provied by \'geno_file\' was not recognized! Please provied one of: Numeric, VCF, HapMap, gds, or plink bed or ped files.",
            call. = F
          )
        }
      } else if (grepl(".GDS$", toupper(files[1]))) {
        input_format <- "gds"
      } else if (grepl(".BED$", toupper(files[1]))) {
        input_format <- "bed"
      } else if (grepl(".PED$", toupper(files[1]))) {
        input_format <- "ped"
      } else {
        data_type <- try(data.table::fread(
          file = files[1],
          header = F,
          skip = 0,
          nrows = 1,
          na.strings = na_string,
          data.table = F
        ),
        silent = TRUE)
        if (any(grepl("VCF", data_type))) {
          input_format <- "VCF"
        } else if (sum(hap_names %in% data_type) > 8 |
                   all(unlist(ifelse(
                     is.character(data_type[, 12]),
                     strsplit(unique(data_type[, 12]), ""),
                     "numeric"
                   )) %in% nucleotide)) {
          input_format <- "hapmap"
        } else {
          stop(
            "File format found in \'geno_path\' was not recognized! Please provied one of: Numeric, VCF, HapMap, gds, or plink bed or ped files.",
            call. = F
          )
        }
      }
      if (input_format == "hapmap") {
        if (verbose)
          message("Reading the following HapMap files: \n")
        if (verbose)
          message(files, sep = "\n")
        G <- vector("list", length(files))
        count <- 1
        for (i in files) {
          G[[count]] <-
            try(data.table::fread(
              file = i,
              header = TRUE,
              skip = 0,
              nrows = nrows,
              na.strings = na_string,
              data.table = F
            ),
            silent = TRUE)
          count <- count + 1
        }
        if (all(lengths(G) != 0)) {
          G <- do.call(rbind, G)
        } else {
          stop("Problems reading multiple HapMap files!
               Check your input data set.")
        }
        GT <- as.matrix(colnames(G)[- (1:11)])
        GI <- G[, c(1, 2, 3, 4)]
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        GD <- NULL
        bit <- nchar(as.character(G[2, 12]))
        if (verbose)
          message("Performing numericalization")
        GD <- apply(G[, - (1:11)], 1, function(one)
          numericalization(
            one,
            bit = bit,
            effect = SNP_effect,
            impute = SNP_impute
          ))
      } else if (input_format == "VCF") {
        if (verbose)
          message("Reading the following VCF files: \n")
        if (verbose)
          message(files, sep = "\n")
        SNPRelate::snpgdsVCF2GDS(
          vcf.fn = files,
          out.fn = temp,
          method = "biallelic.only",
          snpfirstdim = FALSE,
          verbose = FALSE,
          ignore.chr.prefix = chr_prefix
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                   verbose = FALSE) - 1
        GT <- as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "bed") {
        SNPRelate::snpgdsBED2GDS(
          bed.fn = files,
          fam.fn = paste0(gsub(".bed", "", files), ".fam"),
          bim.fn = paste0(gsub(".bed", "", files), ".bim"),
          out.gdsfn = temp,
          snpfirstdim = FALSE,
          verbose = FALSE
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          (SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                    verbose = FALSE) - 1) * -1
        GT <-  as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        allele <- unlist(strsplit(gdsfmt::read.gdsn(
          gdsfmt::index.gdsn(genofile, "snp.allele")
        ), "/"))
        l <- length(allele)
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id")),
          allele = paste0(allele[seq(2, l, 2)], "/", allele[seq(1, l, 2)]),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "ped") {
        SNPRelate::snpgdsPED2GDS(
          ped.fn = files,
          map.fn = paste0(gsub(".ped", "", files), ".map"),
          out.gdsfn = temp,
          snpfirstdim = FALSE,
          verbose = FALSE
        )
        genofile <- SNPRelate::snpgdsOpen(temp)
        GD <-
          SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                   verbose = FALSE) - 1
        GT <-  as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      } else if (input_format == "gds") {
        genofile <- SNPRelate::snpgdsOpen(files)
        GD <-  SNPRelate::snpgdsGetGeno(genofile, snpfirstdim = FALSE,
                                        verbose = FALSE) - 1
        GT <-  as.matrix(gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id")))
        GI <- data.frame(
          SNP = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.rs.id")),
          allele = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.allele")),
          Chromosome = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.chromosome")),
          Position = gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.position")),
          stringsAsFactors = FALSE
        )
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
      }
    }
    if (is.null(GD)) {
      GT <- NULL
      GI <- NULL
    }
    if (input_format != "hapmap") {
      if (any(is.na(GD))) {
        if (SNP_impute == "Middle") {
          GD[is.na(GD)] <- 0
        } else
          if (SNP_impute == "Minor") {
            GD[is.na(GD)] <- -1
          } else
            if (SNP_impute == "Major") {
              GD[is.na(GD)] <- 1
            }
      }
    }
    return(list(
      GT = GT,
      GD = GD,
      GI = GI,
      input_format = input_format,
      out_name = out_name,
      temp = temp
    ))
  }
