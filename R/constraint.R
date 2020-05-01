#' Select SNPs to be assigned as QTNs: To be included LD, FST
#' @keywords internal
#' @param genotypes a numericalized genotype object (geno_obj).
#' @param maf_above Threshold for the minimum value of minor allele frequency.
#' @param maf_below Threshold for the maximum value of minor allele frequency.
#' @param hets Option of including (\'include\') and removing (\'remove\') only
#' heterozygotes.
#' @param verbose = verbose
#' @return Return a filtered dataset to be used when selecting QTNs.
#' Last update: Apr 20, 2020
#'
#'----------------------------- constrain ---------------------------------
constraint <-
  function(genotypes = NULL,
           maf_above = NULL,
           maf_below = NULL,
           hets = NULL,
           verbose = verbose) {
    GD <- genotypes[, - (1:5)]
    list_h <- NULL
    list_maf <- NULL
    if (!is.null(hets)) {
      if (hets != "include" & hets != "remove") {
        stop("hets option must be either \'include\' or \'remove\'.",
             call. = F)
      } else if (hets == "include") {
        if (verbose)
          message("* Filtering variants without heterozygotes.")
        list_h <- apply(GD, 1, function(x)
          any(unique(x) == 0))
      } else {
        if (verbose)
          message("* Filtering heterozygote variants.")
        list_h <- !apply(GD, 1, function(x)
          any(unique(x) == 0))
      }
    }
    if (!is.null(maf_above) |
        !is.null(maf_below)) {
      if (verbose)
        message("* Filtering variants based on MAF.")
      ns <- ncol(GD)
      maf_calc <- apply(GD, 1, function(x) {
        sumx <- ((sum(x) + ns) / ns * 0.5)
        min(sumx,  (1 - sumx))
      })
      if (!is.null(maf_above) &
          !is.null(maf_below)) {
        list_maf <- (maf_calc > maf_above & maf_calc < maf_below)
      } else if (!is.null(maf_above)) {
        list_maf <- (maf_calc > maf_above)
      } else if (!is.null(maf_below)) {
        list_maf <- (maf_calc < maf_below)
      }
    }
    if (!is.null(list_h) & !is.null(list_maf)) {
      selected_snps <- which(list_h & list_maf)
    } else if (!is.null(list_h)) {
      selected_snps <- which(list_h)
    } else if (!is.null(list_maf)) {
      selected_snps <- which(list_maf)
    }
    return(selected_snps)
  }
