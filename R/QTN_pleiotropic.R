#' Select SNPs to be assigned as QTNs.
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param additive_QTN_number = NULL,
#' @param epistatic_QTN_number = NULL
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#------------------------------  QTN_pleiotropic -------------------------------
QTN_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           additive_QTN_number = NULL,
           epistatic_QTN_number = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL)) {
    #---------------------------------------------------------------------------
    # Randomly select (without replacement) k additive QTN,
    # and assign an effect size
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
    add_QTN_genotypic_information <- genotypes[vector_of_add_QTN, ]
    # Create an output file that gives the chromosome, bp, and
    # additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set.seed not assigned"
      ),
      paste0("seed_number_for_", additive_QTN_number,
             "Add_QTN", ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      add_QTN_genotypic_information,
      paste0(
        "Genotypic_information_for_",
        additive_QTN_number,
        "_Additive_QTN",
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
    additive_effect_trait_object <-
      t(add_QTN_genotypic_information[, -c(1:5)])
    colnames(additive_effect_trait_object) <-
      paste0(
        "Chr_",
        unlist(add_QTN_genotypic_information[, 3]),
        "_",
        unlist(add_QTN_genotypic_information[, 4])
      )
    if (!is.null(epistatic_QTN_number)){
      # Randomly select (without replacement) 2*k epistatic QTN, and
      # assign an effect size
      if (!is.null(seed)) {
        set.seed(seed * seed)
      }
      vector_of_epi_QTN <-
        sample(index, (2 * epistatic_QTN_number), replace = FALSE)
      epi_QTN_genotypic_information <-
        genotypes[vector_of_epi_QTN, ]
      # Create an output file that gives the chromosome, bp, and
      # additive effect of the effect sizes, as well as the seed
      write.table(
        ifelse(
          !is.null(seed),
          paste0("Here_is_the_seed_number: ", seed * seed),
          "set.seed not assigned"
        ),
        paste0(
          "seed_number_for_",
          epistatic_QTN_number,
          "Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        epi_QTN_genotypic_information,
        paste0(
          "Genotypic_information_for_",
          epistatic_QTN_number,
          "_Epistatic_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      epistatic_effect_trait_object <-
        t(epi_QTN_genotypic_information[, -c(1:5)])
      colnames(epistatic_effect_trait_object) <-
        paste0(
          "Chr_",
          unlist(epi_QTN_genotypic_information[, 3]),
          "_",
          unlist(epi_QTN_genotypic_information[, 4])
        )
      return(
        list(
          additive_effect_trait_object = list(additive_effect_trait_object),
          epistatic_effect_trait_object = list(epistatic_effect_trait_object)
        )
      )
    } else {
      return(list(
        additive_effect_trait_object = list(additive_effect_trait_object)
      ))
    }
  }