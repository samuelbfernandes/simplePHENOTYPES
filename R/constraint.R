#' Select SNPs to be assigned as QTNs: To be included LD, FST
#' @keywords internal
#' @param genotypes a numericalized genotype object (geno_obj).
#' @param maf_above Threshold for the minimum value of minor allele frequency.
#' @param maf_below Threshold for the maximum value of minor allele frequency.
#' @param hets Option of including (\'include\') and removing (\'remove\') only heterozygotes.
#' @return Return a filtered dataset to be used when selecting QTNs.
#' @author Samuel Fernandes
#' Last update: Nov 25, 2019
#'
#'----------------------------- constrain ---------------------------------
constraint <-
  function(genotypes = NULL,
           maf_above = NULL,
           maf_below = NULL,
           hets = NULL) {
    GD <- genotypes[,-(1:5)]
    list_h <- NULL
    list_maf <- NULL
    if (!is.null(hets)) {
      if (hets != "include" & hets != "remove") {
        stop("hets option must be either \'include\' or \'remove\'.",
             call. = F)
      } else if (hets == "include") {
        h <- apply(GD, 1, table)
        list_h <-  unlist(lapply(h, function(x) {
          any(names(x) == "0")
        }))
      } else {
        h <- apply(GD, 1, table)
        list_h <-  !unlist(lapply(h, function(x) {
          any(names(x) == "0")
        }))
      }
    }
    if (!is.null(maf_above) |
        !is.null(maf_below)) {
      GD <- GD + 1
      ns <- ncol(GD)
      ss <- apply(GD, 1, sum)
      maf_matrix <- rbind((0.5 * ss / ns), (1 - (0.5 * ss / ns)))
      maf_calc <- apply(maf_matrix, 2, min)
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
      selected_snps <-which(list_maf)
    }
    return(selected_snps)
  }
