#' Generate a numeric (dosage) HapMap file
#' @keywords internal
#' @param geno_obj = NULL,
#' @param geno_file = NULL,
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param prefix = NULL,
#' @param geno_path = NULL,
#' @param maf_cutoff = NULL,
#' @param SNP_effect = 'Add',
#' @param SNP_impute = 'Middle',
#' @param verbose = verbose
#' @param chr_prefix = "chrm"
#' @return A numeric HapMap
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Apr 20, 2020
#'
#'------------------------------------------------------------------------------
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
    hmp <- file_loader(
      geno_obj = geno_obj,
      geno_file = geno_file,
      geno_path = geno_path,
      nrows = nrows,
      na_string = na_string,
      prefix = prefix,
      verbose = verbose,
      SNP_impute = SNP_impute,
      chr_prefix = chr_prefix
    )
    if (!is.null(maf_cutoff)) {
      hm <- list(GT = hmp$GT,
                 GD = hmp$GD,
                 GI = hmp$GI)
      ns <- nrow(hm$GD)
      ss <- apply(hm$GD, 2, sum)
      maf_matrix <- rbind((0.5 * ss / ns), (1 - (0.5 * ss / ns)))
      maf <- apply(maf_matrix, 2, min)
      snps_below_maf <- which(maf < maf_cutoff)
      hm_GI_filtered <- hm$GI[-snps_below_maf, ]
      geno_obj <-
        data.frame(
          hm$GI[-snps_below_maf, ],
          rep(NA, nrow(hm_GI_filtered)),
          t(hm$GD[, -snps_below_maf]),
          check.names = FALSE,
          fix.empty.names = FALSE
        )
      colnames(geno_obj) <-
        c("snp", "allele", "chr", "pos", "cm", t(as.character(hm$GT)))
    } else {
      if (nrow(hmp$GD) == nrow(hmp$GI)) {
        geno_obj <-
          data.frame(
            hmp$GI,
            rep(NA, nrow(hmp$GI)),
            hmp$GD,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
      } else{
        geno_obj <-
          data.frame(
            hmp$GI,
            rep(NA, nrow(hmp$GI)),
            t(hmp$GD),
            check.names = FALSE,
            fix.empty.names = FALSE
          )
      }
      colnames(geno_obj) <-
        c("snp", "allele", "chr", "pos", "cm", t(as.character(hmp$GT)))
    }
    return(list(
      geno_obj = geno_obj,
      input_format =  hmp$input_format,
      out_name =  hmp$out_name,
      temp = hmp$temp
    ))
  }
