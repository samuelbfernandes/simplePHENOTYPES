#' Select SNPs to be assigned as QTNs.
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param Additive.QTN.number = NULL,
#' @param Epistatic.QTN.number = NULL
#' @export
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#------------------------------  QTN_pleiotropic ---------------------------------
`QTN_pleiotropic` <-
  function(genotypes = NULL,
           seed = NULL,
           Additive.QTN.number = NULL,
           Epistatic.QTN.number = NULL) {

    #---------------------------------------------------------------------------
    #Randomly select (without replacement) k additive QTN, and assign an effect size
    if (!is.null(seed)) {
      set.seed(seed)
    }

    vector.of.add.QTN <-
      sample(1:nrow(genotypes), Additive.QTN.number, replace = FALSE)
    Add.QTN.genotypic.information <- genotypes[vector.of.add.QTN, ]

    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed),
        "set.seed not assigned"
      ),
      paste0("seed.number.for.",
             Additive.QTN.number,
             "Add.QTN",
             ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      Add.QTN.genotypic.information,
      paste0(
        "Genotypic.information.for.",
        Additive.QTN.number,
        ".Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )

    #Create a "base line" trait, which is basically just the additive effects;
    #this is what we would see if the heritability of the simulated trait were 1
    additive.effect.trait.object <-
      t(Add.QTN.genotypic.information[, -c(1:5)])

    colnames(additive.effect.trait.object) <-
      paste0("Chr_",
             Add.QTN.genotypic.information[, 3],
             "_",
             Add.QTN.genotypic.information[, 4])

    if ( !is.null(Epistatic.QTN.number) ){
    #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
    if (!is.null(seed)) {
      set.seed(seed * seed)
    }

    vector.of.epi.QTN <-
      sample(1:nrow(genotypes), (2 * Epistatic.QTN.number), replace = FALSE)
    Epi.QTN.genotypic.information <- genotypes[vector.of.epi.QTN, ]

    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      ifelse(
        !is.null(seed),
        paste0("Here_is_the_seed_number: ", seed * seed),
        "set.seed not assigned"
      ),
      paste0("seed.number.for.",
             Epistatic.QTN.number,
             "Epi.QTN",
             ".txt"),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    data.table::fwrite(
      Epi.QTN.genotypic.information,
      paste0(
        "Genotypic.information.for.",
        Epistatic.QTN.number,
        ".Epistatic.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )

    epistatic.effect.trait.object <-
      t(Epi.QTN.genotypic.information[, -c(1:5)])

    colnames(epistatic.effect.trait.object) <-
      paste0("Chr_",
             Epi.QTN.genotypic.information[, 3],
             "_",
             Epi.QTN.genotypic.information[, 4])
    return(
      list(
        additive.effect.trait.object = list(additive.effect.trait.object),
        epistatic.effect.trait.object = list(epistatic.effect.trait.object)
      )
    )
    } else {
      return(
        list(
          additive.effect.trait.object = list(additive.effect.trait.object)
        )
      )
    }
  }

