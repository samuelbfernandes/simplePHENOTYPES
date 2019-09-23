#' Generate a numeric (dosaje) HapMap file
#' @export
#' @param genotypes_object = NULL,
#' @param genotypes_file = NULL,
#' @param input_format = "hapmap",
#' @param skip = 0,
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param shared_name = NULL,
#' @param genotypes_path = NULL,
#' @param maf_cutoff = NULL,
#' @param SNP_effect = 'Add',
#' @param SNP_impute = 'Middle',
#' @param major_allele_zero = FALSE,
#' @return A numeric HapMap
#' @author Alex lipka and Samuel Fernandes
#' Last update: Sep 19, 2019
#'
#'------------------------------------------------------------------------------
genotypes <-
  function(genotypes_object = NULL,
           genotypes_file = NULL,
           genotypes_path = NULL,
           input_format = "hapmap",
           skip = 0,
           nrows = Inf,
           na_string = "NA",
           shared_name = NULL,
           maf_cutoff = NULL,
           SNP_effect = "Add",
           SNP_impute = "Middle",
           major_allele_zero = FALSE) {
    #---------------------------------------------------------------------------
    hmp <- file_loader(
      genotypes_object = genotypes_object,
      genotypes_file = genotypes_file,
      genotypes_path = genotypes_path,
      input_format = input_format,
      skip = skip,
      nrows = nrows,
      na_string = na_string,
      shared_name = shared_name
    )
if (!is.null(maf_cutoff)) {
  hm <- list(GT = hmp$GT,
             GD = hmp$GD,
             GI = hmp$GI)
  # Obtain the mafs of all SNPs
  #-------------------------------------------------------------------------
  # Total number of lines
  ns <- nrow(hm$GD)
  # Sum of the allele scores for each SNP
  ss <- apply(hm$GD, 2, sum)
  # Combine two situations: one where the allele coded as '2' is major;
  # one where '0' is coded as major.
  maf_matrix <- rbind( (0.5 * ss / ns), (1 - (0.5 * ss / ns)))
  # Copy the minor allele frequencies for all SNPs
  maf <- apply(maf_matrix, 2, min)
  # Find out which SNPs have MAF < maf_cutoff
  snps_below_maf <- which(maf < maf_cutoff)
  # Remove these SNPs from hm$GD
  hm_GD_without_snps_below_maf <- hm$GD[, -snps_below_maf]
  genotypes_object <-
    data.frame(hm$GI,
               rep(NA, nrow(hm$GI)),
               t(hm_GD_without_snps_below_maf))
  colnames(genotypes_object) <-
    c("snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
} else {
  genotypes_object <-
    data.frame(hmp$GI,
               rep(NA, nrow(hmp$GI)),
               t(hmp$GD))
  colnames(genotypes_object) <-
    c("snp", "allele", "chr", "pos", "cm", t(as.character(hmp$GT)))
}
return(genotypes_object)
}
