#' Function used by create_phenotypes to load marker dataset.
#' @param geno_obj = NULL,
#' @param geno_file = NULL,
#' @param geno_path = NULL,
#' @param input_format = "hapmap",
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param prefix = NULL,
#' @param SNP_impute = "Middle",
#' @keywords internal
#' @param SNP_effect = "Add",
#' @param major_allele_zero = FALSE
#' @param verbose = TRUE
#' @return genotype data sampled
#' @author  Samuel Fernandes
#' Last update: Nov 05, 2019
#'------------------------------------------------------------------------------
file_loader <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           input_format = "hapmap",
           nrows = Inf,
           na_string = "NA",
           prefix = NULL,
           SNP_impute = "Middle",
           SNP_effect = "Add",
           major_allele_zero = FALSE,
           verbose = TRUE) {
    #'--------------------------------------------------------------------------
    if (is.null(geno_obj) &&
        is.null(geno_file) &&
        is.null(geno_path)){
      stop("Please provide one of: \'geno_obj\', \'geno_file\' or \'geno_path\'", call. = F)
    }
    if (!is.null(geno_obj)){
      if (verbose) cat("File loaded from memory. \n")
      if (input_format == "hapmap") {
        GT <- as.matrix(colnames(geno_obj)[- (1:11)])
        GI <- geno_obj[, c(1, 2, 3, 4)]
        if (verbose) print(paste0(
          "Converting HapMap format to numerical under model of ",
          SNP_impute
        ))
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        GD <- NULL
        bit <- nchar(as.character(geno_obj[2, 12]))
        if (verbose) print("Performing numericalization")
        GD <- apply(geno_obj[, - (1:11)], 1, function(one)
          GAPIT_numericalization(
            one,
            bit = bit,
            effect = SNP_effect,
            impute = SNP_impute,
            major_allele_zero = major_allele_zero
          ))
        if (is.null(GD)) {
          GT <- NULL
          GI <- NULL
        }
      }
      return(list(GT = GT, GD = GD, GI = GI))
    }else{
      if (input_format == "hapmap") {
        if (is.null(geno_path)) {
          if (!geno_file %in% dir()) {
            stop(paste("File ",geno_file," not found."), call. = F)
          }
          G <-
            try(data.table::fread(
              file = geno_file,
              head = TRUE,
              skip = 0,
              nrows = nrows,
              na.strings = na_string,
              data.table = F
            ),
            silent = TRUE)
        } else {
          if (is.null(prefix)) {
            files <- paste0(geno_path, "/", dir(geno_path))
          } else{
            files <- paste0(geno_path, "/",
                            dir(geno_path)[grepl(prefix,
                                                      dir(geno_path))])
          }
          if (verbose) cat("Reading the following HapMap files: \n")
          if (verbose) cat( files, sep = "\n")
          G <- vector("list", length(files))
          count <- 1
          for (i in files) {
            G[[count]] <-
              try(data.table::fread(
                file = i,
                head = TRUE,
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
        }
      }
      GT <- as.matrix(colnames(G)[- (1:11)])
      GI <- G[, c(1, 2, 3, 4)]
      colnames(GT) <- "taxa"
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      GD <- NULL
      bit <- nchar(as.character(G[2, 12]))
      GD <- apply(G[, - (1:11)], 1, function(one)
        GAPIT_numericalization(
          one,
          bit = bit,
          effect = SNP_effect,
          impute = SNP_impute,
          major_allele_zero = major_allele_zero
        ))
      if (is.null(GD)) {
        GT <- NULL
        GI <- NULL
      }
      return(list(GT = GT, GD = GD, GI = GI))
    }
  }
