#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param maf_above = NULL
#' @param maf_below = NULL
#' @return ...
#' @author Samuel Fernandes
#' Last update: Sep 23, 2019
#'
#'----------------------------- constrain ---------------------------------
constrain <-
  function(genotypes = NULL,
           maf_above = NULL,
           maf_below = NULL) {
    GD <- genotypes[, - (1:5)]
    ns <- ncol(GD)
    ss <- apply(GD, 1, sum)
    maf_matrix <- rbind( (0.5 * ss / ns), (1 - (0.5 * ss / ns)))
    maf_calc <- apply(maf_matrix, 2, min)
    if (!is.null(maf_above) &&
        !is.null(maf_below)) {
      selected_snps <- which(maf_calc > maf_above && maf_calc < maf_below)
    } else if (!is.null(maf_above)) {
      selected_snps <- which(maf_calc > maf_above)
    } else if (!is.null(maf_below)){
      selected_snps <- which(maf_calc < maf_below)
    }
    return(selected_snps)
  }
