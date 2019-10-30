#' ...
#' @export
#' @param genotypes_object = NULL,
#' @param genotypes_file = NULL,
#' @param genotypes_path = NULL,
#' @param input_format = "hapmap",
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param prefix = NULL,
#' @param SNP_impute = "Middle",
#' @param SNP_effect = "Add",
#' @param major_allele_zero = FALSE
#' @return genotype data sampled
#' @author  Samuel Fernandes
#' Last update: Sep 19, 2019
#'------------------------------------------------------------------------------
file_loader <-
  function(genotypes_object = NULL,
           genotypes_file = NULL,
           genotypes_path = NULL,
           input_format = "hapmap",
           nrows = Inf,
           na_string = "NA",
           prefix = NULL,
           SNP_impute = "Middle",
           SNP_effect = "Add",
           major_allele_zero = FALSE) {
    #'--------------------------------------------------------------------------
    if (is.null(genotypes_object) &&
        is.null(genotypes_file) &&
        is.null(genotypes_path)){
      stop("Please provide one of: \'genotypes_object\',
           \'genotypes_file\' or \'genotypes_path\'")
    }
    if (!is.null(genotypes_object)){
      cat("File loaded from memory. \n")
      if (input_format == "hapmap") {
        GT <- as.matrix(colnames(genotypes_object)[- (1:11)])
        GI <- genotypes_object[, c(1, 2, 3, 4)]
        print(paste0(
          "Converting HapMap format to numerical under model of ",
          SNP_impute
        ))
        # Set column names
        colnames(GT) <- "taxa"
        colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
        # Initial GD
        GD <- NULL
        # to determine number of bits of genotypes_object
        bit <- nchar(as.character(genotypes_object[2, 12]))
        print("Performing numericalization")
        GD <- apply(genotypes_object[, - (1:11)], 1, function(one)
          GAPIT_numericalization(
            one,
            bit = bit,
            effect = SNP_effect,
            impute = SNP_impute,
            major_allele_zero = major_allele_zero
          ))
        # set GT and GI to NULL in case of null GD
        if (is.null(GD)) {
          GT <- NULL
          GI <- NULL
        }
      }
      return(list(GT = GT, GD = GD, GI = GI))
    }else{
      if (input_format == "hapmap") {
        if (is.null(genotypes_path)) {
          G <-
            try(data.table::fread(
              file = genotypes_file,
              head = TRUE,
              skip = 0,
              nrows = nrows,
              na.strings = na_string,
              data.table = F
            ),
            silent = TRUE)
        } else {
          if (is.null(prefix)) {
            files <- paste0(genotypes_path,"/", dir(genotypes_path))
          } else{
            files <- paste0(genotypes_path,"/", 
                            dir(genotypes_path)[grepl(prefix,
                                                      dir(genotypes_path))])
          }
          cat("Reading the following HapMap files: \n")
          cat( files, sep = "\n")
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
      print(paste0(
        "Converting HapMap format to numerical under model of ",
        SNP_impute
      ))
      GT <- as.matrix(colnames(G)[- (1:11)])
      GI <- G[, c(1, 2, 3, 4)]
      # Set column names
      colnames(GT) <- "taxa"
      colnames(GI) <- c("SNP", "allele", "Chromosome", "Position")
      # Initial GD
      GD <- NULL
      # to determine number of bits of genotypes_object
      bit <- nchar(as.character(G[2, 12]))
      print("Performing numericalization")
      GD <- apply(G[, - (1:11)], 1, function(one)
        GAPIT_numericalization(
          one,
          bit = bit,
          effect = SNP_effect,
          impute = SNP_impute,
          major_allele_zero = major_allele_zero
        ))
      # set GT and GI to NULL in case of null GD
      if (is.null(GD)) {
        GT <- NULL
        GI <- NULL
      }
      return(list(GT = GT, GD = GD, GI = GI))
    }
  }
