#' Simulation of ntraits phenotypes based on a SNP file
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes A HapMap format dataset from which simulated data will
#' be generated.
#' @param file.G Optional file name (without the ".hmp.txt" extension).
#' @param file.Ext.G Optional file extension.
#' @param file.fragment Optional argument for loading only a subset of rows
#' of the genotipic dataset.
#' @param file.from Option for files saved by chromosome
#' (minimum chromosome number).
#' @param file.to Option for files saved by chromosome
#' (maximum chromosome number).
#' @param maf_cutoff Optional filter for minor allele frequency.
#' @param SNP.effect Following GAPIT implementation. Default 'Add'.
#' @param SNP.impute Following GAPIT implementation. Default 'Middle'.
#' @param Create.indicator Following GAPIT implementation. Default FALSE.
#' @param Major.allele.zero Following GAPIT implementation. Default FALSE.
#' @param Additive.QTN.number Number of additive quantitative trait nucleotide
#' to be simulated.
#' @param Epistatic.QTN.number Number of epistatic (additive x additive)
#' quantitative trait nucleotide to be simulated.
#' @param additive.effect Additive effect size to be simulated. Follows a
#' geometric series.
#' @param epistatic.effect Epistatic (additive x additive) effect size to be
#' simulated. Follows a geometric series.
#' @param big.additive.QTN.effect Additive effect size for one possible major
#' effect quantitative trait nucleotide.
#' @param model Genetic architecture to be simulated. Possible options are:
#' 'pleiotropic', for traits being controled by all the same QTNs;
#' 'partially', for traits being controled by pleiotropic and trait specific QTNs;
#' 'LD', for traits being exclusively controled QTNs in linkage disequilibrium.
#' @param overlap Number of pleiotropic additive QTNs if model = 'partially'.
#' @param overlapE Number of pleiotropic epistatic QTNs if model = 'partially'.
#' @param specific.QTN.number Number of trait specific additive QTNs if
#' model = 'partially'.
#' @param specific.E.QTN.number Number of trait specific epistatic QTNs if
#' model = 'partially'.
#' @param ld Linkage disequilibrium between selected marker two adjacent markers
#' to be used as QTN. Default is ld = 05.
#' @param rep Number of experiments to be simulated.
#' @param ntraits Number of traits to be simulated under pleitropic and
#' partially pleiotropic models. The default for linkage disequilibrium models is two.
#' @param h2 Heritability of target trait (if multiple traits are simulated).
#' If a vector, the simulation will "loop" over it and generate one file for each
#' combination of genetic settings.
#' @param h2_MT Heritability of correlated traits (in case it should be different).
#' than target trait. It should be length ntraits-1.
#' @param set.cor Should correlation among traits be controled for? Default TRUE.
#' @param correlation Trait correlation matrix to be simulated.
#' To be used if set.cor = TRUE.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 10e5) will be used.
#' @param home.dir Home directory. Default is current working directory.
#' @param output.dir Name to be used to create folder and save output files.
#' @param format Three options for saving outputs: 'multi-file'
#' (default for multiple traits), saves one simulation setting in a separate file;
#' 'long', appends each experiment (rep) to the last one (by row); 'wide', saves
#' experiments by column (default for single trait).
#' @param out.geno Saves numericalized genotype either as "hapmap", "plink" or "gds",
#' @param gdsfile gds file (in case there is one already created) to be used
#' with option model = "LD". Default is NULL.
#' @return Numericalized marker dataset, selected QTNs, phenotypes for 'ntraits' traits.
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Aug 2nd, 2019
#' @examples
#' # Simulate 50 replications of a single phenotype.
#'
#' create.simulated.data(
#'   genotypes = SNP55K_maize282,
#'   Additive.QTN.number = 3,
#'   additive.effect = c(0.1, 0.2),
#'   big.additive.QTN.effect = 0.9,
#'   rep = 50,
#'   h2 = 0.7
#' )
#'

create.simulated.data <-
  function(genotypes = NULL,
           file.G = NULL,
           file.Ext.G = NULL,
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
           model = "pleiotropic",
           overlap = NULL,
           overlapE = NULL,
           specific.QTN.number = NULL,
           specific.E.QTN.number = NULL,
           ld = 0.5,
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
           out.geno = "none",
           gdsfile = NULL) {
    # -------------------------------------------------------------------------
    .onAttach <- function(libname, pkgname) {
      packageStartupMessage("Welcome to the simplePHENOTYPES package :)")
    }
    .onAttach()
    setwd(home.dir)
    on.exit(setwd(home.dir), add = TRUE)

    if (is.null(genotypes) ||
        (class(unlist(genotypes[, 12])) != "numeric" &&
         class(unlist(genotypes[, 12])) != "integer")) {
      genotypes <-
        Genotypes(
          geno = genotypes,
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

    if (is.null(seed))
      seed <- round(runif(1, 0, 10e5))
    
    # Createing and setting a working directory for the output results
    if (!is.null(output.dir)) {
      tempdir <- paste0(home.dir, "/", output.dir)
      if (dir.exists(tempdir)) {
        j <- 1
        while (dir.exists(tempdir)) {
          tempdir <- paste0(home.dir, "/", output.dir, "(", j, ")")
          j <- j + 1
        }
        message("Directory alredy exists! Creating ", tempdir)
        path_out <- tempdir
        dir.create(path_out)
        setwd(path_out)
      } else {
        path_out <- tempdir
        dir.create(path_out)
        setwd(path_out)
      }
    } else {
      path_out <- home.dir
    }
    
    zz <- file("Log_Sim.txt", "w")
    sink(zz, type = "output")

    if (ntraits > 1 | model == "LD") {
      mm <-
        ifelse(
          model == "pleiotropic",
          "pleiotropic",
          ifelse(
            model == "partially",
            "partially pleiotropic",
            ifelse(model == "LD", "Linkage Desiquilibrium",
                   {
                     stop("model used is not valid!", call. = F)
                     sink()
                     close(zz)
                   })
          )
        )
      cat(paste("Simulation of a", mm, "genetic model \n"))
    }
    
    if (model == "LD" ||
        out.geno == "plink" ||
        out.geno == "gds") {
      
      if (model == "LD") {
        ntraits <- 2
        if (length(ld) > Additive.QTN.number) {
          message("Length of ld object > Additive.QTN.number. Using first ", Additive.QTN.number)
          ld <- ld[1:Additive.QTN.number]
        }
        if (length(ld) > 1 & length(ld) < Additive.QTN.number) {
          message("Length of ld object < Additive.QTN.number. Using ld[1]=", ld[1])
          ld <- ld[1]
        }
      }
      
      genotypes <- genotypes[!duplicated(genotypes$Snp), ]
      if (is.null(gdsfile))
        gdsfile <- "geno"
      SNPRelate::snpgdsCreateGeno(
        paste0(gdsfile, ".gds"),
        genmat = as.matrix(genotypes[, -c(1:5)]),
        sample.id = colnames(genotypes)[-c(1:5)],
        snp.id = genotypes$Snp,
        snp.chromosome = genotypes$chr,
        snp.position = genotypes$pos,
        snp.allele = genotypes$allele,
        snpfirstdim = TRUE
      )
      gdsfile <- paste0(home.dir, "/", gdsfile, ".gds")
      cat(paste0("GDS files saved at:", home.dir, "\n"))
    }
    
    if (!is.null(genotypes)) {
      cat("genotype object inputed form R \n")
    } else {
      cat(paste("genotype file:", file.G, "*", file.Ext.G))
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
              paste(specific.QTN.number,
                    collapse = ", "),
              "\n \n"
            )
          )
        } else {
          cat(
            paste(
              "\nNumber of pleiotropic additive QTNs:",
              overlap,
              "\nNumber of trait specific additive QTNs:",
              paste(specific.QTN.number,
                    collapse = ", "),
              "\nNumber of pleiotropic epistatic QTNs:",
              overlapE,
              "\nNumber of trait specific epistatic QTNs:",
              paste(specific.E.QTN.number,
                    collapse = ", "),
              "\n \n"
            )
          )
        }
      }
    }

    cat("Number of traits:", ntraits, "\n")
      
    if (!is.null(out.geno)) {
      if (out.geno == "numeric") {
        data.table::fwrite(
          genotypes,
          paste0(
            ifelse(is.null(output.dir), "", "../"),
            file.G,
            "Numericalized_Genotypes.txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        print(paste0("Numeric Genotypes saved at:", home.dir))
      } else if (out.geno == "plink") {
        # Open the GDS file
        genofile <- SNPRelate::snpgdsOpen(gdsfile)
        
        snpset <-
          SNPRelate::snpgdsSelectSNP(genofile, missing.rate = 0.95)
        SNPRelate::snpgdsGDS2BED(genofile, bed.fn = "geno", snp.id = snpset)
        SNPRelate::snpgdsClose(genofile)
        
        cat(paste0("Plink files saved at:", home.dir, "\n"))
      }
    }
    if (ntraits > 1) {
      if (set.cor & !all(ntraits == dim(correlation))) {
        stop("ntraits and Correlation matrix do not match!", call. = F)
        sink()
        close(zz)
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
        } else {
          cat(
            paste(
              "Heritability not assigned for all correlated traits!
              \nSetting all traits with the same heritability (",
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
          specific.QTN.number <-
            specific.QTN.number[1:ntraits]
        } else {
          stop("Please set up [additive] specific.QTN.number", call. = F)
          sink()
          close(zz)
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
              paste(specific.E.QTN.number[1:ntraits],
                    collapse = ", "),
              ")",
              "values\n"
            )
          )
          specific.E.QTN.number <-
            specific.E.QTN.number[1:ntraits]
        } else {
          stop("Please set up specific.E.QTN.number", call. = F)
          sink()
          close(zz)
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
        } else {
          cat(
            paste(
              "Length additive.effect < ntraits!
              \nSetting all traits with the additive.effect = ",
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
        } else {
          cat(
            paste(
              "Length epistatic.effect < ntraits!
              \nSetting all traits with the epistatic.effect = ",
              epistatic.effect[1],
              "\n"
            )
          )
          epistatic.effect <-
            rep(epistatic.effect[1], ntraits)
        }
      }
      
      if (ntraits != length(big.additive.QTN.effect)) {
        if (length(big.additive.QTN.effect) > ntraits) {
          cat(
            paste(
              "Length of big.additive.QTN.effect > ntraits! using the first",
              ntraits,
              "(",
              paste(big.additive.QTN.effect[1:ntraits],
                    collapse = ", "),
              ")",
              "values\n"
            )
          )
          big.additive.QTN.effect <-
            big.additive.QTN.effect[1:ntraits]
        } else {
          cat(
            paste(
              "Length big.additive.QTN.effect < ntraits!
              \nSetting all traits with the big.additive.QTN.effect = ",
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
      QTN <-
        QTN_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          Additive.QTN.number = Additive.QTN.number,
          Epistatic.QTN.number = Epistatic.QTN.number
        )
    }
    
    if (ntraits > 1 & !any(model != "partially")) {
      QTN <-
        QTN_partially_pleiotropic(
          genotypes = genotypes,
          seed = seed,
          overlap = overlap,
          overlapE = overlapE,
          specific.QTN.number = specific.QTN.number,
          specific.E.QTN.number = specific.E.QTN.number,
          ntraits = ntraits
        )
    }
    
    if (ntraits > 1 & !any(model != "LD")) {
      QTN <-
        QTN_linkage(
          genotypes = genotypes,
          seed = seed,
          Additive.QTN.number = Additive.QTN.number,
          ld = ld,
          gdsfile = gdsfile
        )
    }
    
    if (ntraits == 1) {
      Genetic_value <-
        Base_line_single_trait(
          seed = seed,
          additive.object = QTN$additive.effect.trait.object,
          epistatic.object = QTN$epistatic.effect.trait.object,
          additive.effect = additive.effect,
          epistatic.effect = epistatic.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
    } else if (!any(model != "LD")) {
      cat("LD between SNP and QTN:", ld)
      Genetic_value_sup <-
        Base_line_single_trait(
          seed = seed,
          additive.object = QTN$additive.effect.trait.object.sup,
          additive.effect = additive.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
      
      Genetic_value_inf <-
        Base_line_single_trait(
          seed = seed,
          additive.object = QTN$additive.effect.trait.object.inf,
          additive.effect = additive.effect,
          big.additive.QTN.effect = big.additive.QTN.effect
        )
      
      Genetic_value <-
        list(
          base.line = cbind(
            Genetic_value_sup$base.line,
            Genetic_value_inf$base.line
          ),
          VA = c(Genetic_value_sup$VA, Genetic_value_inf$VA)
        )
      
      cat("\nDIAGNOSTICS: \n")
      
      cat("\nGenetic Correlation \n")
      print(cor(Genetic_value$base.line)[1, 2])
    } else {
      Genetic_value <-
        Base_line_multi_traits(
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
      h2_MT = h2_MT,
      format = format
    )
    cat(paste("\n\nResults are saved at:", getwd()))
    sink()
    close(zz)
    file.show(paste0(path_out, "/Log_Sim.txt"))
    
  }  # end 'create.simluated.data()'
