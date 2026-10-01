#' Select SNPs to be assigned as QTNs: To be included LD, FST
#' @keywords internal
#' @param genotypes a numericalized genotype object (geno_obj).
#' @param maf_above Threshold for the minimum value of minor allele frequency.
#' The comparison is strict: a marker is kept only if its MAF is *greater than*
#' `maf_above`.
#' @param maf_below Threshold for the maximum value of minor allele frequency.
#' The comparison is strict: a marker is kept only if its MAF is *smaller than*
#' `maf_below`.
#' @param hets Option of including (\'include\') and removing (\'remove\') only
#' heterozygotes.
#' @param verbose = verbose
#' @return Row indices of the markers that satisfy the constraints (all rows
#' when no constraint is given). The frozen selection code applies them to the
#' randomly drawn QTNs only: in the LD architectures the linked partner markers
#' are found by walking along the marker order and are *not* filtered.
#' Last update: Apr 20, 2020
#'
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
      ns <- ncol(GD)
      maf_calc <- apply(GD, 1, function(x) {
        sumx <- ((sum(x) + ns) / ns * 0.5)
        min(sumx,  (1 - sumx))
      })
      if (!is.null(maf_above) &
          !is.null(maf_below)) {
        if (verbose)
          message(paste("* Removing variants with MAF above", maf_below, "or below", maf_above,"!"))
        list_maf <- (maf_calc > maf_above & maf_calc < maf_below)
      } else if (!is.null(maf_above)) {
        if (verbose)
          message(paste("* Removing variants with MAF below", maf_above,"!"))
        list_maf <- (maf_calc > maf_above)
      } else if (!is.null(maf_below)) {
        if (verbose)
          message(paste("* Removing variants with MAF above", maf_below,"!"))
        list_maf <- (maf_calc < maf_below)
      }
    }
    if (!is.null(list_h) & !is.null(list_maf)) {
      selected_snps <- which(list_h & list_maf)
    } else if (!is.null(list_h)) {
      selected_snps <- which(list_h)
    } else if (!is.null(list_maf)) {
      selected_snps <- which(list_maf)
    } else {
      # no constraint requested: every marker is eligible
      selected_snps <- seq_len(nrow(genotypes))
    }
    return(selected_snps)
  }
