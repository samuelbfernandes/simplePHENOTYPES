#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param overlap = NULL,
#' @param overlap_e = NULL,
#' @param specific_QTN_number = NULL,
#' @param specific_e_QTN_number = NULL,
#' @param ntraits = NULL
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_partially_pleiotropic ----------------------
QTN_partially_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           overlap = NULL,
           overlap_e = NULL,
           specific_QTN_number = NULL,
           specific_e_QTN_number = NULL,
           ntraits = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL)) {
    #---------------------------------------------------------------------------
    # Randomly select (without replacement) k additive QTN, and
    # assign an effect size using data.table::fwrite from data.table package
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
    vector_of_pleiotropic_add_QTN <-
      sample(index, overlap, replace = FALSE)
    add_pleiotropic_QTN_genotypic_info <-
      genotypes[vector_of_pleiotropic_add_QTN, ]
    data.table::fwrite(
      add_pleiotropic_QTN_genotypic_info,
      paste0(
        "Genotypic_information_for_",
        overlap,
        "_pleiotropic_Additive_QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    snps <-
      setdiff(index, vector_of_pleiotropic_add_QTN)
    vector_of_specific_add_QTN <- list()
    add_specific_QTN_genotypic_info <- list()
    ss <- c()
    for (i in 1:ntraits) {
      if (!is.null(seed)) {
        ss[i] <- seed + i
        set.seed(seed + i)
      }
      vector_of_specific_add_QTN[[i]] <-
        sample(snps, specific_QTN_number[i], replace = FALSE)
      snps <- setdiff(snps, vector_of_specific_add_QTN[[i]])
      add_specific_QTN_genotypic_info[[i]] <-
        genotypes[vector_of_specific_add_QTN[[i]], ]
      data.table::fwrite(
        add_specific_QTN_genotypic_info[[i]],
        paste0(
          "Trait",
          i,
          "_Additive_Genotypic_info_for_",
          specific_QTN_number[i],
          "_specific_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    }
    # Create an output file that gives the chromosome, bp, and
    # additive effect of the effect sizes, as well as the seed
    write.table(
      c(seed, ss),
      paste0(
        "seed_number_for_",
        paste0(specific_QTN_number + overlap, collapse = "_"),
        "Add_QTN",
        ".txt"
      ),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    if (!is.null(overlap_e) | !is.null(specific_e_QTN_number)) {
      # Randomly select (without replacement) 2*k epistatic QTN,
      # and assign an effect size
      if (!is.null(seed)) {
        set.seed(seed + seed)
      }
      vector_of_pleiotropic_epi_QTN <-
        sample(index, (2 * overlap_e), replace = FALSE)
      epi_pleiotropic_QTN_genotypic_info <-
        genotypes[vector_of_pleiotropic_epi_QTN, ]
      # Create an output file that gives the chromosome, bp, and
      # additive effect of the effect sizes, as well as the seed
      data.table::fwrite(
        epi_pleiotropic_QTN_genotypic_info,
        paste0(
          "Genotypic_information_for_",
          overlap_e,
          "_pleiotropic_Epistatic_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      snps_e <-
        setdiff(index, vector_of_pleiotropic_epi_QTN)
      vector_of_specific_epi_QTN <- list()
      epi_specific_QTN_genotypic_info <- list()
      sse <- c()
      for (i in 1:ntraits) {
        if (!is.null(seed)) {
          sse[i] <- (seed + i) + seed
          set.seed( (seed + i) + seed)
        }
        vector_of_specific_epi_QTN[[i]] <-
          sample(snps_e, (2 * specific_e_QTN_number[i]), replace = FALSE)
        snps_e <- setdiff(snps_e, vector_of_specific_epi_QTN[[i]])
        epi_specific_QTN_genotypic_info[[i]] <-
          genotypes[vector_of_specific_epi_QTN[[i]], ]
        data.table::fwrite(
          epi_specific_QTN_genotypic_info[[i]],
          paste0(
            "Trait",
            i,
            "_Epistatic_Genotypic_info_for_",
            specific_e_QTN_number[i],
            "_specific_QTN",
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      write.table(
        c(seed + seed, sse),
        paste0(
          "seed_number_for_",
          paste0(specific_e_QTN_number + overlap_e, collapse = "_"),
          "pleiotropic_Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      # Create a 'base line' trait, which is basically just the additive effects
      # this is what we would see if the heritability of the simulated
      # trait were 1
      additive_effect_trait_object <-
        as.data.frame(add_pleiotropic_QTN_genotypic_info)
      epistatic_effect_trait_object <-
        as.data.frame(epi_pleiotropic_QTN_genotypic_info)
      for (i in 1:ntraits) {
        additive_effect_trait_object <-
          rbind(additive_effect_trait_object,
                as.data.frame(add_specific_QTN_genotypic_info[[i]]))
        epistatic_effect_trait_object <-
          rbind(epistatic_effect_trait_object,
                as.data.frame(epi_specific_QTN_genotypic_info[[i]]))
      }
      rownames(additive_effect_trait_object) <-
        paste0("Chr_",
               additive_effect_trait_object$chr,
               "_",
               additive_effect_trait_object$pos)
      rownames(epistatic_effect_trait_object) <-
        paste0("Chr_",
               epistatic_effect_trait_object$chr,
               "_",
               epistatic_effect_trait_object$pos)
      additive_effect_trait_object <-
        list(t(additive_effect_trait_object[, -c(1:5)]))
      epistatic_effect_trait_object <-
        list(t(epistatic_effect_trait_object[, -c(1:5)]))
      return(
        list(
          additive_effect_trait_object = additive_effect_trait_object,
          epistatic_effect_trait_object = epistatic_effect_trait_object
        )
      )
    } else {
      # Create a 'base line' trait, which is basically just the additive effects
      # this is what we would see if the heritability of the simulated
      # trait were 1
      additive_effect_trait_object <-
        as.data.frame(add_pleiotropic_QTN_genotypic_info)
      for (i in 1:ntraits) {
        additive_effect_trait_object <-
          rbind(additive_effect_trait_object,
          as.data.frame(add_specific_QTN_genotypic_info[[i]]))
      }
      rownames(additive_effect_trait_object) <-
        paste0("Chr_",
              additive_effect_trait_object$chr,
              "_",
              additive_effect_trait_object$pos)
      additive_effect_trait_object <-
        list(t(additive_effect_trait_object[, -c(1:5)]))
      return(list(additive_effect_trait_object = additive_effect_trait_object))
    }
  }