#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param additive_QTN_number = NULL,
#' @param ld = NULL,
#' @param gdsfile NULL
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_linkage ---------------------------------
QTN_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           additive_QTN_number = NULL,
           ld = NULL,
           gdsfile = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL)) {
    #---------------------------------------------------------------------------
    # Randomly select (without replacement)
    # k additive QTN, and assign an effect size
    if (!is.null(seed)) {
      set.seed(seed)
    }
    if (any(lengths(constrains)>0)) { 
      index <- constrain(genotypes = genotypes, 
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      # First SNP at column 6
      index <- 6:nrow(genotypes)
    }
    vector_of_add_QTN <-
      sample(index, additive_QTN_number, replace = FALSE)
    genofile <- SNPRelate::snpgdsOpen(gdsfile)
    inf <- c()
    sup <- c()
    x <- 1
    for (k in 1:length(ld)){
      for (j in vector_of_add_QTN) {
      ldsup <- 1
      i <- j + 1
      while (ldsup >= ld[k]) {
        snp1 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(j, 1),
            count = c(1, -1)
          )
        snp2 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(i, 1),
            count = c(1, -1)
          )
        ldsup <-
          abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
        if (is.nan(ldsup)) {
          SNPRelate::snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        i <- i + 1
      }
      sup[x] <- i
      ldinf <- 1
      i2 <- j - 1
      while (ldinf >= ld[k]) {
        snp3 <-
          gdsfmt::read.gdsn(
            gdsfmt::index.gdsn(genofile, "genotype"),
            start = c(i2, 1),
            count = c(1, -1)
          )
        ldinf <-
          abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = "composite"))
        if (is.nan(ldinf)) {
          SNPRelate::snpgdsClose(genofile)
          stop("Monomorphic SNPs are not accepted", call. = F)
        }
        i2 <- i2 - 1
      }
      inf[x] <- i2
      x <- x + 1
      }
    }
    # close the genotype file
    SNPRelate::snpgdsClose(genofile)
    add_QTN_genotypic_information_sup <- genotypes[sup, ]
    add_QTN_genotypic_information_inf <- genotypes[inf, ]
    # Create an output file that gives the chromosome, bp,
    # and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set_seed not assigned"
      ),
      paste0("seed_number_for_", additive_QTN_number,
             "Add_QTN", ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      add_QTN_genotypic_information_sup,
      paste0(
        "Genotypic_information_for_",
        additive_QTN_number,
        "LD_",
        paste(ld, collapse = "_"),
        "_SUP_Additive_QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    data.table::fwrite(
      add_QTN_genotypic_information_inf,
      paste0(
        "Genotypic_information_for_",
        additive_QTN_number,
        "LD_",
        paste(ld, collapse = "_"),
        "_INF_Additive_QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    # Create a 'base line' trait, which is basically just the additive effects;
    # this is what we would see if the heritability of the simulated
    # trait were 1
    additive_effect_trait_object_sup <-
      t(add_QTN_genotypic_information_sup[, -c(1:5)])
    additive_effect_trait_object_inf <-
      t(add_QTN_genotypic_information_inf[, -c(1:5)])
    colnames(additive_effect_trait_object_sup) <-
      paste0(
        "Chr_",
        unlist(add_QTN_genotypic_information_sup[, 3]),
        "_",
        unlist(add_QTN_genotypic_information_sup[, 4])
      )
    colnames(additive_effect_trait_object_inf) <-
      paste0(
        "Chr_",
        unlist(add_QTN_genotypic_information_inf[, 3]),
        "_",
        unlist(add_QTN_genotypic_information_inf[, 4])
      )
    return(
      list(
        additive_effect_trait_object_sup = additive_effect_trait_object_sup,
        additive_effect_trait_object_inf = additive_effect_trait_object_inf
      )
    )
  }