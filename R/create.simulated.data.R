#' Main function for simulation of ntraits phenotypes based on a SNP file
#' @export
#' @importFrom data.table fwrite fread
#' @import utils
#' @import stats
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes jjj
#' @param file.G = NULL,
#' @param file.Ext.G = NULL,
#' @param file.fragment = Inf,
#' @param file.from = 1,
#' @param file.to = 1,
#' @param maf_cutoff = NULL,
#' @param SNP.effect = "Add",
#' @param SNP.impute = "Middle",
#' @param Create.indicator = FALSE,
#' @param Major.allele.zero = FALSE,
#' @param Additive.QTN.number = NULL,
#' @param Epistatic.QTN.number = NULL,
#' @param additive.effect = NULL,
#' @param epistatic.effect = NULL,
#' @param big.additive.QTN.effect = NULL,
#' @param model = c("pleiotropic", "partially", "LD"),
#' @param overlap = NULL,
#' @param overlapE = NULL,
#' @param specific.QTN.number = NULL,
#' @param specific.E.QTN.number = NULL,
#' @param rep = NULL,
#' @param ntraits = 1,
#' @param h2 = NULL,
#' @param h2_MT = NULL,
#' @param set.cor = TRUE,
#' @param correlation = NULL,
#' @param seed = NULL,
#' @param home.dir = getwd(),
#' @param output.dir = NULL,
#' @param format = "multi-file",
#' @param out.geno = FALSE
#' @return Numeric hapmap, selected QTNs, phenotypes for ntraits traits
#' @author Alex Lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'-------------------------------------------------------------------------------
create.simulated.data <-
  function(genotypes = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
           #'default for numericalization#
           file.fragment = Inf,
           file.from = 1,
           file.to = 1,
           maf_cutoff = NULL,
           SNP.effect = "Add",
           SNP.impute = "Middle",
           Create.indicator = FALSE,
           Major.allele.zero = FALSE,
           Additive.QTN.number = NULL,
           Epistatic.QTN.number = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           model = c("pleiotropic", "partially", "LD"),
           overlap = NULL,
           overlapE = NULL,
           specific.QTN.number = NULL,
           specific.E.QTN.number = NULL,
           rep = NULL,
           ntraits = 1,
           h2 = NULL,
           h2_MT = NULL,
           set.cor = TRUE,
           correlation = NULL,
           seed = NULL,
           home.dir = getwd(),
           output.dir = NULL,
           format = "multi-file",
           out.geno = FALSE) {
    #' -----------------------------------------------------------------------------
    try({
      .onAttach <- function(libname, pkgname) {
        packageStartupMessage("Welcome to my package")
      }
      setwd(home.dir)
      on.exit(setwd(home.dir), add = TRUE)
      if (is.null(genotypes)) {
        genotypes <- Genotypes(
          file.path = paste0(home.dir, "/"),
          maf_cutoff = maf_cutoff,
          seed = 123,
          file.G = file.G,
          file.Ext.G = file.Ext.G,
          SNP.effect = SNP.effect,
          SNP.impute = SNP.impute,
          file.fragment = file.fragment,
          Create.indicator = Create.indicator,
          Major.allele.zero = Major.allele.zero,
          file.from = file.from,
          file.to = file.to
        )
      }

      #' Create a working directory for the output results:
      if (!is.null(output.dir)) {
        dir.create(paste0(home.dir, "/", output.dir))
      }

      #' Set the working directory
      setwd(paste0(home.dir, "/", output.dir, "/"))
      sink("Log_Sim.txt", type = "output")

      if (ntraits > 1) {
        mm <- ifelse(
          model == "pleiotropic",
          "pleiotropic",
          ifelse(
            model == "partially",
            "partially pleiotropic",
            ifelse(model == "LD", "Linkage Desiquilibrium",
                   {
                     stop("model used is not valid!", call. = F)
                     sink()
                   })
          )
        )

        cat(paste("Simulation of a", mm, "genetic model \n"))
      }
      if (is.null(file.G)) {
        cat(paste("genotype object:", deparse(substitute(genotypes))))
      } else{
        cat(paste("genotype file:", file.G , "*", file.Ext.G))
      }

      if (ntraits > 1) {
        if (model == "pleiotropic") {
          cat(
            paste(
              "\nNumber of additive QTNs:",
              Additive.QTN.number,
              "\nNumber of epistatic QTNs:",
              Epistatic.QTN.number,
              "\n"
            )
          )
        }

        if (model == "partially") {
          if (is.null(overlapE) | is.null(specific.E.QTN.number)) {
            cat(
              paste(
                "\nNumber of pleiotropic additive QTNs:",
                overlap,
                "\nNumber of trait specific additive QTNs:",
                paste(specific.QTN.number, collapse = ", "),
                "\n \n"
              )
            )
          } else{
            cat(
              paste(
                "\nNumber of pleiotropic additive QTNs:",
                overlap,
                "\nNumber of trait specific additive QTNs:",
                paste(specific.QTN.number, collapse = ", "),
                "\nNumber of pleiotropic epistatic QTNs:",
                overlapE,
                "\nNumber of trait specific epistatic QTNs:",
                paste(specific.E.QTN.number, collapse = ", "),
                "\n \n"
              )
            )
          }
        }
      }
      if (out.geno) {
        data.table::fwrite(
          genotypes,
          paste0(ifelse(is.null(output.dir), "", "../"), file.G, "NUM.txt"),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        print(paste0("Numeric Genotypes saved at:",  home.dir))
      }
      if (ntraits > 1) {
        if (set.cor &
            !all(ntraits == dim(correlation))) {
          stop("ntraits and Correlation matrix do not match!", call. = F)
          sink()
        }

        if (ntraits != (length(h2_MT) + 1)) {
          if (length(h2_MT) > (ntraits - 1)) {
            cat(
              paste(
                "Length of h2_MT > number of correlated traits! using the first",
                ntraits - 1,
                "(",
                paste(h2_MT[1:(ntraits - 1)], collapse = ", "),
                ")",
                "values\n"
              )
            )
            h2_MT <- h2_MT[1:(ntraits - 1)]
          } else{
            cat(
              paste(
                "Heritability not assigned for all correlated traits! \nSetting all traits with the same heritability (",
                h2[1],
                ")\n"
              )
            )
            h2_MT <- rep(h2[1], (ntraits - 1))
          }
        }

        if (model == "partially" &
            !is.null(specific.QTN.number) &
            ntraits != length(specific.QTN.number)) {
          if (length(specific.QTN.number) > ntraits) {
            cat(
              paste(
                "Length of additive.effect > ntraits! using the first",
                ntraits,
                "(",
                paste(specific.QTN.number[1:ntraits], collapse = ", "),
                ")",
                "values\n"
              )
            )
            specific.QTN.number <- specific.QTN.number[1:ntraits]
          } else{
            stop("Please set up [additive] specific.QTN.number",
                 call. = F)
            sink()
          }
        }

        if (model == "partially" &
            !is.null(specific.E.QTN.number) &
            ntraits != length(specific.E.QTN.number)) {
          if (length(specific.E.QTN.number) > ntraits) {
            cat(
              paste(
                "Length of specific.E.QTN.number > ntraits! using the first",
                ntraits,
                "(",
                paste(specific.E.QTN.number[1:ntraits], collapse = ", "),
                ")",
                "values\n"
              )
            )
            specific.E.QTN.number <- specific.E.QTN.number[1:ntraits]
          } else{
            stop("Please set up specific.E.QTN.number", call. = F)
            sink()
          }
        }

        if (ntraits != length(additive.effect)) {
          if (length(additive.effect) > ntraits) {
            cat(
              paste(
                "Length of additive.effect > ntraits! using the first",
                ntraits,
                "(",
                paste(additive.effect[1:ntraits], collapse = ", "),
                ")",
                "values\n"
              )
            )
            additive.effect <- additive.effect[1:ntraits]
          } else{
            cat(
              paste(
                "Length additive.effect < ntraits! \nSetting all traits with the additive.effect = ",
                additive.effect[1],
                "\n"
              )
            )
            additive.effect <- rep(additive.effect[1], ntraits)
          }
        }

        if (!is.null(epistatic.effect) &
            ntraits != length(epistatic.effect)) {
          if (length(epistatic.effect) > ntraits) {
            cat(
              paste(
                "Length of epistatic.effect > ntraits! using the first",
                ntraits,
                "(",
                paste(epistatic.effect[1:ntraits], collapse = ", "),
                ")",
                "values\n"
              )
            )
            epistatic.effect <- epistatic.effect[1:ntraits]
          } else{
            cat(
              paste(
                "Length epistatic.effect < ntraits! \nSetting all traits with the epistatic.effect = ",
                epistatic.effect[1],
                "\n"
              )
            )
            epistatic.effect <- rep(epistatic.effect[1], ntraits)
          }
        }

        if (ntraits != length(big.additive.QTN.effect)) {
          if (length(big.additive.QTN.effect) > ntraits) {
            cat(
              paste(
                "Length of big.additive.QTN.effect > ntraits! using the first",
                ntraits,
                "(",
                paste(big.additive.QTN.effect[1:ntraits], collapse = ", "),
                ")",
                "values\n"
              )
            )
            big.additive.QTN.effect <-
              big.additive.QTN.effect[1:ntraits]
          } else{
            cat(
              paste(
                "Length big.additive.QTN.effect < ntraits! \nSetting all traits with the big.additive.QTN.effect = ",
                big.additive.QTN.effect[1],
                "\n"
              )
            )
            big.additive.QTN.effect <-
              rep(big.additive.QTN.effect[1], ntraits)
          }
        }
      }

      if (ntraits == 1 | !any(model != "pleiotropic")) {
        QTN <- QTN_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          Additive.QTN.number = Additive.QTN.number,
          Epistatic.QTN.number = Epistatic.QTN.number
        )
      }

      if (ntraits > 1 & !any(model != "partially")) {
        QTN <- QTN_partially_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          overlap = overlap,
          overlapE = overlapE,
          specific.QTN.number = specific.QTN.number,
          specific.E.QTN.number = specific.E.QTN.number,
          ntraits = ntraits
        )
      }

      if (ntraits == 1) {
        Genetic_value <- Base_line_single_trait(
          seed = seed,
          additive.object = QTN$additive.effect.trait.object,
          epistatic.object = QTN$epistatic.effect.trait.object,
          additive.effect = additive.effect,
          epistatic.effect = epistatic.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
      } else {
        Genetic_value <- Base_line_multi_traits(
          seed = seed,
          set.cor = set.cor,
          ntraits = ntraits,
          correlation = correlation,
          additive.object = QTN$additive.effect.trait.object,
          epistatic.object = QTN$epistatic.effect.trait.object,
          additive.effect = additive.effect,
          epistatic.effect = epistatic.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
        cat("\nDIAGNOSTICS: \n")
        if (set.cor) {
          cat("Populational Correlation \n")
          colnames(correlation) <-
            c("Target", c(paste0("Trait_", 2:ntraits)))
          rownames(correlation) <-
            c("Target", c(paste0("Trait_", 2:ntraits)))
          print(correlation)
        }
        cat("\nSample Correlation \n")
        sample.cor <- Genetic_value$sample.cor
        colnames(sample.cor) <-
          c("Target", c(paste0("Trait_", 2:ntraits)))
        rownames(sample.cor) <-
          c("Target", c(paste0("Trait_", 2:ntraits)))
        print(sample.cor)
      }

      Phenotypes(
        seed = seed,
        base.line.trait = Genetic_value$base.line,
        h2 = h2,
        rep = rep,
        ntraits = ntraits,
        h2_MT =  h2_MT,
        format = format
      )

      closeAllConnections()
      setwd(home.dir)
      file.show(paste0(output.dir, "/Log_Sim.txt"))
    }, silent = TRUE)

  }#'end "create.simluated.data()"
