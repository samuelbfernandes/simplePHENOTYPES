#' Select SNPs to be assigned as QTNs
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param overlap = NULL,
#' @param overlapE = NULL,
#' @param specific.QTN.number = NULL,
#' @param specific.E.QTN.number = NULL,
#' @param ntraits = NULL
#' @export
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_partially_pleiotropic ---------------------------------
`QTN_partially_pleiotropic` <-
  function(genotypes = NULL,
           seed = NULL,
           overlap = NULL,
           overlapE = NULL,
           specific.QTN.number = NULL,
           specific.E.QTN.number = NULL,
           ntraits = NULL) {
    #---------------------------------------------------------------------------
    #Randomly select (without replacement) k additive QTN, and assign an effect size
    # using data.table::fwrite from data.table package
    if (!is.null(seed)) {
      set.seed(seed)
    }

    vector.of.pleiotropic.add.QTN <-
      sample(1:nrow(genotypes), overlap, replace = FALSE)

    Add.pleiotropic.QTN.genotypic.info <-
      genotypes[vector.of.pleiotropic.add.QTN,]

    data.table::fwrite(
      Add.pleiotropic.QTN.genotypic.info,
      paste0(
        "Genotypic.information.for.",
        overlap,
        ".pleiotropic.Additive.QTN",
        ".txt"
      ),
      row.names = FALSE,
      sep = "\t",
      quote = FALSE,
      na = NA
    )

    snps <-
      setdiff(1:nrow(genotypes), vector.of.pleiotropic.add.QTN)
    vector.of.specific.add.QTN <- list()
    Add.specific.QTN.genotypic.info <- list()
    ss <- c()
    for (i in 1:ntraits) {
      if (!is.null(seed)) {
        ss[i] <- seed + i
        set.seed(seed + i)
      }
      vector.of.specific.add.QTN[[i]] <-
        sample(snps, specific.QTN.number[i], replace = FALSE)

      snps <- setdiff(snps, vector.of.specific.add.QTN[[i]])

      Add.specific.QTN.genotypic.info[[i]] <-
        genotypes[vector.of.specific.add.QTN[[i]], ]

      data.table::fwrite(
        Add.specific.QTN.genotypic.info[[i]],
        paste0(
          "Trait",
          i,
          ".Additive.Genotypic.info.for.",
          specific.QTN.number[i],
          ".specific.QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    }

    #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed
    write.table(
      c(seed, ss),
      paste0(
        "seed.number.for.",
        paste0(specific.QTN.number + overlap, collapse = "_"),
        "Add.QTN",
        ".txt"
      ),
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      quote = FALSE
    )

    if (!is.null(overlapE) | !is.null(specific.E.QTN.number)) {
      #Randomly select (without replacement) 2*k epistatic QTN, and assign an effect size
      if (!is.null(seed)) {
        set.seed(seed * seed)
      }

      vector.of.pleiotropic.epi.QTN <-
        sample(1:nrow(genotypes), (2 * overlapE), replace = FALSE)
      Epi.pleiotropic.QTN.genotypic.info <-
        genotypes[vector.of.pleiotropic.epi.QTN,]

      #Create an output file that gives the chromosome, bp, and additive effect of the effect sizes, as well as the seed

      data.table::fwrite(
        Epi.pleiotropic.QTN.genotypic.info,
        paste0(
          "Genotypic.information.for.",
          overlapE,
          ".pleiotropic.Epistatic.QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )

      snpse <-
        setdiff(1:nrow(genotypes), vector.of.pleiotropic.epi.QTN)
      vector.of.specific.epi.QTN <- list()
      Epi.specific.QTN.genotypic.info <- list()
      sse <- c()
      for (i in 1:ntraits) {
        if (!is.null(seed)) {
          sse[i] <- (seed + i) * seed
          set.seed((seed + i) * seed)
        }
        vector.of.specific.epi.QTN[[i]] <-
          sample(snps, (2 * specific.E.QTN.number[i]) , replace = FALSE)

        snps <- setdiff(snps, vector.of.specific.epi.QTN[[i]])

        Epi.specific.QTN.genotypic.info[[i]] <-
          genotypes[vector.of.specific.epi.QTN[[i]], ]

        data.table::fwrite(
          Epi.specific.QTN.genotypic.info[[i]],
          paste0(
            "Trait",
            i,
            ".Epistatic.Genotypic.info.for.",
            specific.E.QTN.number[i],
            ".specific.QTN",
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }

      write.table(
        c(seed * seed, sse),
        paste0(
          "seed.number.for.",
          paste0(specific.E.QTN.number + overlapE, collapse = "_"),
          "pleiotropic.Epi.QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )

      #Create a "base line" trait, which is basically just the additive effects;
      #this is what we would see if the heritability of the simulated trait were 1
      additive.effect.trait.object <- list()
      epistatic.effect.trait.object <- list()
      for (i in 1:ntraits) {
        additive.effect.trait.object[[i]] <-
          t(
            rbind(
              Add.pleiotropic.QTN.genotypic.info,
              Add.specific.QTN.genotypic.info[[i]]
            )[,-c(1:5)]
          )
        epistatic.effect.trait.object[[i]] <-
          t(
            rbind(
              Epi.pleiotropic.QTN.genotypic.info,
              Epi.specific.QTN.genotypic.info[[i]]
            )[,-c(1:5)]
          )

        colnames(additive.effect.trait.object[[i]]) <-
          c(
            paste0(
              "Chr_",
              unlist(Add.pleiotropic.QTN.genotypic.info[, 3]),
              "_",
              unlist(Add.pleiotropic.QTN.genotypic.info[, 4])
            ),
            paste0(
              "Chr_",
              unlist(Add.specific.QTN.genotypic.info[[i]][, 3]),
              "_",
              unlist(Add.specific.QTN.genotypic.info[[i]][, 4])
            )
          )

        colnames(epistatic.effect.trait.object[[i]]) <-
          c(
            paste0(
              "Chr_",
              unlist(Epi.pleiotropic.QTN.genotypic.info[, 3]),
              "_",
              unlist(Epi.pleiotropic.QTN.genotypic.info[, 4])
            ),
            paste0(
              "Chr_",
              unlist(Epi.specific.QTN.genotypic.info[[i]][, 3]),
              "_",
              unlist(Epi.specific.QTN.genotypic.info[[i]][, 4])
            )
          )
      }

      return(
        list(
          additive.effect.trait.object = additive.effect.trait.object,
          epistatic.effect.trait.object = epistatic.effect.trait.object
        )
      )
    } else{
      #Create a "base line" trait, which is basically just the additive effects;
      #this is what we would see if the heritability of the simulated trait were 1
      additive.effect.trait.object <- list()

      for (i in 1:ntraits) {
        additive.effect.trait.object[[i]] <-
          t(
            rbind(
              Add.pleiotropic.QTN.genotypic.info,
              Add.specific.QTN.genotypic.info[[i]]
            )[,-c(1:5)]
          )

        colnames(additive.effect.trait.object[[i]]) <-
          c(
            paste0(
              "Chr_",
              unlist(Add.pleiotropic.QTN.genotypic.info[, 3]),
              "_",
              unlist(Add.pleiotropic.QTN.genotypic.info[, 4])
            ),
            paste0(
              "Chr_",
              unlist(Add.specific.QTN.genotypic.info[[i]][, 3]),
              "_",
              unlist(Add.specific.QTN.genotypic.info[[i]][, 4])
            )
          )
      }

      return(list(additive.effect.trait.object = additive.effect.trait.object))
    }
  }
