#' Simulation of ntraits phenotypes based on a SNP file
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes_object A HapMap format dataset from which simulated data will
#' be generated.
#' @param genotypes_file = NULL,
#' @param input_format = "hapmap",
#' @param skip = 0,
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param shared_name = NULL,
#' @param genotypes_path = NULL,
#' @param maf_cutoff Optional filter for minor allele frequency.
#' @param SNP_effect Following GAPIT implementation. Default 'Add'.
#' @param SNP_impute Following GAPIT implementation. Default 'Middle'.
#' @param major_allele_zero Following GAPIT implementation. Default FALSE.
#' @param additive_QTN_number Number of additive quantitative trait nucleotide
#' to be simulated.
#' @param epistatic_QTN_number Number of epistatic (additive x additive)
#' quantitative trait nucleotide to be simulated.
#' @param additive_effect Additive effect size to be simulated. Follows a
#' geometric series.
#' @param epistatic_effect Epistatic (additive x additive) effect size to be
#' simulated. Follows a geometric series.
#' @param big_additive_QTN_effect Additive effect size for one possible major
#' effect quantitative trait nucleotide.
#' @param model Genetic architecture to be simulated. Possible options are:
#' 'pleiotropic', for traits being controled by all the same QTNs;
#' 'partially', for traits being controled by pleiotropic and trait specific QTNs;
#' 'LD', for traits being exclusively controled QTNs in linkage disequilibrium.
#' @param overlap Number of pleiotropic additive QTNs if model = 'partially'.
#' @param overlap_e Number of pleiotropic epistatic QTNs if model = 'partially'.
#' @param specific_QTN_number Number of trait specific additive QTNs if
#' model = 'partially'.
#' @param specific_e_QTN_number Number of trait specific epistatic QTNs if
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
#' @param set_cor Should correlation among traits be controled for? Default TRUE.
#' @param correlation Trait correlation matrix to be simulated.
#' To be used if set_cor = TRUE.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 10e5) will be used.
#' @param home_dir Home directory. Default is current working directory.
#' @param output_dir Name to be used to create folder and save output files.
#' @param format Three options for saving outputs: 'multi-file'
#' (default for multiple traits), saves one simulation setting in a separate file;
#' 'long', appends each experiment (rep) to the last one (by row); 'wide', saves
#' experiments by column (default for single trait).
#' @param out_geno Saves numericalized genotype either as "hapmap", "plink" or "gds",
#' @param gdsfile gds file (in case there is one already created) to be used
#' with option model = "LD". Default is NULL.
#' @return Numericalized marker dataset, selected QTNs, phenotypes for 'ntraits' traits.
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Aug 2nd, 2019
#' @examples
#' # Simulate 50 replications of a single phenotype.
#'
#' create_simulated_data(
#'   genotypes_object = SNP55K_maize282,
#'   additive_QTN_number = 3,
#'   additive_effect = c(0.1, 0.2),
#'   big_additive_QTN_effect = 0.9,
#'   rep = 50,
#'   h2 = 0.7
#' )
#'
create_simulated_data <-
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
           major_allele_zero = FALSE,
           additive_QTN_number = NULL,
           epistatic_QTN_number = NULL,
           additive_effect = NULL,
           epistatic_effect = NULL,
           big_additive_QTN_effect = NULL,
           model = "pleiotropic",
           overlap = NULL,
           overlap_e = NULL,
           specific_QTN_number = NULL,
           specific_e_QTN_number = NULL,
           ld = 0.5,
           rep = NULL,
           ntraits = 1,
           h2 = NULL,
           h2_MT = NULL,
           set_cor = TRUE,
           correlation = NULL,
           seed = NULL,
           home_dir = getwd(),
           output_dir = NULL,
           format = "multi-file",
           out_geno = "none",
           gdsfile = NULL) {
    # -------------------------------------------------------------------------
    .onAttach <- function(libname, pkgname) {
      packageStartupMessage("Thank you for using the simplePHENOTYPES package!")
    }
    .onAttach()
    setwd(home_dir)
    on.exit(setwd(home_dir), add = TRUE)
    if (!is.null(genotypes_path) ||
        !is.null(genotypes_file) ||
      (class(unlist(genotypes_object[, 12])) != "numeric" &&
         class(unlist(genotypes_object[, 12])) != "integer")) {
      genotypes_object <-
        genotypes(genotypes_object = genotypes_object,
        genotypes_path = genotypes_path,
        genotypes_file = genotypes_file,
        input_format = input_format,
        skip = skip,
        nrows = nrows,
        na_string = na_string,
        shared_name = shared_name,
        maf_cutoff = maf_cutoff,
        SNP_impute = SNP_impute,
        major_allele_zero = major_allele_zero)
    }
    if (is.null(seed))
      seed <- round(runif(1, 0, 10e5))
    # Createing and setting a working directory for the output results
    if (!is.null(output_dir)) {
      tempdir <- paste0(home_dir, "/", output_dir)
      if (dir.exists(tempdir)) {
        j <- 1
        while (dir.exists(tempdir)) {
          tempdir <- paste0(home_dir, "/", output_dir, "(", j, ")")
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
      path_out <- home_dir
    }
    zz <- file("Log_Sim.txt", "w")
    sink(zz, type = "output")
    if (ntraits > 1 | model == "LD") {
      mm <-
        ifelse(
          model == "pleiotropic",
          "Pleiotropic",
          ifelse(
            model == "partially",
            "Partially Pleiotropic",
            ifelse(model == "LD", "Linkage Desiquilibrium", {
                     stop("model used is not valid!", call. = F)
                     sink()
                     close(zz)
                   })
          )
        )
      cat(paste("Simulation of a", mm, "genetic model \n"))
    } else {
      cat("Simulation of a Single Trait genetic model \n")
      }
    if (model == "LD" ) ntraits <- 2
    cat("Number of traits:", ntraits, "\n")
      if (model == "pleiotropic" || model == "LD" || ntraits == 1 ) {
        cat(
          paste(
            "\nNumber of additive QTNs:",
            additive_QTN_number,
            "\nNumber of epistatic QTNs:",
            ifelse(is.null(epistatic_QTN_number), 0, epistatic_QTN_number ),
            "\n"
          )
        )
      }
      if (model == "partially") {
        if (is.null(overlap_e) | is.null(specific_e_QTN_number)) {
          cat(
            paste(
              "\nNumber of pleiotropic additive QTNs:",
              overlap,
              "\nNumber of trait specific additive QTNs:",
              paste(specific_QTN_number,
                    collapse = ", "),
              "\n"
            )
          )
        } else {
          cat(
            paste(
              "\nNumber of pleiotropic additive QTNs:",
              overlap,
              "\nNumber of trait specific additive QTNs:",
              paste(specific_QTN_number,
                    collapse = ", "),
              "\nNumber of pleiotropic epistatic QTNs:",
              overlap_e,
              "\nNumber of trait specific epistatic QTNs:",
              paste(specific_e_QTN_number,
                    collapse = ", "),
              "\n"
            )
          )
        }
      }
    if (model == "LD" ||
        out_geno == "plink" ||
        out_geno == "gds"){
      if ( out_geno == "plink" || out_geno == "gds")
        cat("Creating GDS files...\n")
      if (model == "LD") {
        if (length(ld) > additive_QTN_number) {
          message("Length of ld object > additive_QTN_number. Using first ",
                  additive_QTN_number, "values")
          ld <- ld[1:additive_QTN_number]
        }
        if (length(ld) > 1 & length(ld) < additive_QTN_number) {
          message("Length of ld object < additive_QTN_number.
                  Using ld[1]=", ld[1])
          ld <- ld[1]
        }
      }
      genotypes_object <- genotypes_object[!duplicated(genotypes_object$snp), ]
      if (any(grepl("geno.gds", dir(home_dir)))){
        message("A file named geno.gds is already present in this folder")
        gdsfile <- "geno2"
        }
      if (is.null(gdsfile))  gdsfile <- "geno"
      if (!is.numeric(genotypes_object$chr)) {
        genotypes_object$chr <- as.numeric(
          gsub("\\D+", "", genotypes_object$chr)
        )
        }
      SNPRelate::snpgdsCreateGeno(
        paste0(home_dir, "/", gdsfile, ".gds"),
        genmat = as.matrix(genotypes_object[, -c(1:5)]),
        sample.id = colnames(genotypes_object)[-c(1:5)],
        snp.id = genotypes_object$snp,
        snp.chromosome = genotypes_object$chr,
        snp.position = genotypes_object$pos,
        snp.allele = genotypes_object$allele,
        snpfirstdim = TRUE
      )
      gdsfmt::showfile.gds(closeall = TRUE)
      gdsfile <- paste0(home_dir, "/", gdsfile, ".gds")
      cat(paste0("GDS files saved at:", home_dir, "\n"))
    }
    if (!is.null(out_geno)) {
      if (out_geno == "numeric") {
        data.table::fwrite(
          genotypes_object,
          paste0(
            ifelse(is.null(output_dir), "", "../"),
            "Numeric_SNP_File"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        print(paste0("Numeric Genotypes saved at:", home_dir))
      } else if (out_geno == "plink") {
        # Open the GDS file
        genofile <- SNPRelate::snpgdsOpen(gdsfile)
        snpset <-
          SNPRelate::snpgdsSelectSNP(genofile, missing.rate = 0.95)
        SNPRelate::snpgdsGDS2BED(genofile, bed.fn = "geno", snp.id = snpset)
        SNPRelate::snpgdsClose(genofile)
        cat(paste0("Plink files saved at:", home_dir, "\n"))
      }
    }
    if (ntraits > 1) {
      if (set_cor & !all(ntraits == dim(correlation))) {
        stop("ntraits do not match with dim of correlation matrix!", call. = F)
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
          !is.null(specific_QTN_number) &
          ntraits != length(specific_QTN_number)) {
        if (length(specific_QTN_number) > ntraits) {
          cat(
            paste(
              "Length of additive_effect > ntraits! using the first",
              ntraits,
              "(",
              paste(specific_QTN_number[1:ntraits], collapse = ", "),
              ")",
              "values\n"
            )
          )
          specific_QTN_number <-
            specific_QTN_number[1:ntraits]
        } else {
          stop("Please set up [additive] specific_QTN_number", call. = F)
          sink()
          close(zz)
        }
      }
      if (model == "partially" &
          !is.null(specific_e_QTN_number) &
          ntraits != length(specific_e_QTN_number)) {
        if (length(specific_e_QTN_number) > ntraits) {
          cat(
            paste(
              "Length of specific_e_QTN_number > ntraits! using the first",
              ntraits,
              "(",
              paste(specific_e_QTN_number[1:ntraits],
                    collapse = ", "),
              ")",
              "values\n"
            )
          )
          specific_e_QTN_number <-
            specific_e_QTN_number[1:ntraits]
        } else {
          stop("Please set up specific_e_QTN_number", call. = F)
          sink()
          close(zz)
        }
      }
      if (ntraits != length(additive_effect)) {
        if (length(additive_effect) > ntraits) {
          cat(
            paste(
              "Length of additive_effect > ntraits! using the first",
              ntraits,
              "(",
              paste(additive_effect[1:ntraits], collapse = ", "),
              ")",
              "values\n"
            )
          )
          additive_effect <- additive_effect[1:ntraits]
        } else {
          cat(
            paste(
              "Length additive_effect < ntraits!
              \nSetting all traits with the additive_effect = ",
              additive_effect[1],
              "\n"
            )
          )
          additive_effect <- rep(additive_effect[1], ntraits)
        }
      }
      if (!is.null(epistatic_effect) &
          ntraits != length(epistatic_effect)) {
        if (length(epistatic_effect) > ntraits) {
          cat(
            paste(
              "Length of epistatic_effect > ntraits! using the first",
              ntraits,
              "(",
              paste(epistatic_effect[1:ntraits], collapse = ", "),
              ")",
              "values\n"
            )
          )
          epistatic_effect <- epistatic_effect[1:ntraits]
        } else {
          cat(
            paste(
              "Length epistatic_effect < ntraits!
              \nSetting all traits with the epistatic_effect = ",
              epistatic_effect[1],
              "\n"
            )
          )
          epistatic_effect <-
            rep(epistatic_effect[1], ntraits)
        }
      }
      if (ntraits != length(big_additive_QTN_effect)) {
        if (length(big_additive_QTN_effect) > ntraits) {
          cat(
            paste(
              "Length of big_additive_QTN_effect > ntraits! using the first",
              ntraits,
              "(",
              paste(big_additive_QTN_effect[1:ntraits],
                    collapse = ", "),
              ")",
              "values\n"
            )
          )
          big_additive_QTN_effect <-
            big_additive_QTN_effect[1:ntraits]
        } else {
          cat(
            paste(
              "Length big_additive_QTN_effect < ntraits!
              \nSetting all traits with the big_additive_QTN_effect = ",
              big_additive_QTN_effect[1],
              "\n"
            )
          )
          big_additive_QTN_effect <-
            rep(big_additive_QTN_effect[1], ntraits)
        }
      }
    }
    cat("\nDIAGNOSTICS: \n\n")
    if (ntraits == 1 | !any(model != "pleiotropic")) {
      QTN <-
        QTN_pleiotropic(
          genotypes = genotypes_object,
          seed = seed,
          additive_QTN_number = additive_QTN_number,
          epistatic_QTN_number = epistatic_QTN_number
        )
    }
    if (ntraits > 1 & !any(model != "partially")) {
      QTN <-
        QTN_partially_pleiotropic(
          genotypes = genotypes_object,
          seed = seed,
          overlap = overlap,
          overlap_e = overlap_e,
          specific_QTN_number = specific_QTN_number,
          specific_e_QTN_number = specific_e_QTN_number,
          ntraits = ntraits
        )
    }
    if (ntraits > 1 & !any(model != "LD")) {
      QTN <-
        QTN_linkage(
          genotypes = genotypes_object,
          seed = seed,
          additive_QTN_number = additive_QTN_number,
          ld = ld,
          gdsfile = gdsfile
        )
    }
    if (ntraits == 1) {
      genetic_value <-
        base_line_single_trait(
          seed = seed,
          additive_object = QTN$additive_effect_trait_object,
          epistatic_object = QTN$epistatic_effect_trait_object,
          additive_effect = additive_effect,
          epistatic_effect = epistatic_effect,
          big_additive_QTN_effect = big_additive_QTN_effect
        )
    } else if (!any(model != "LD")) {
      cat("LD between SNP and QTN:", ld, "\n")
      genetic_value_sup <-
        base_line_single_trait(
          seed = seed,
          additive_object = QTN$additive_effect_trait_object_sup,
          additive_effect = additive_effect,
          big_additive_QTN_effect = big_additive_QTN_effect
        )
      genetic_value_inf <-
        base_line_single_trait(
          seed = seed,
          additive_object = QTN$additive_effect_trait_object_inf,
          additive_effect = additive_effect,
          big_additive_QTN_effect = big_additive_QTN_effect
        )
      genetic_value <-
        list(
          base.line = cbind(
            genetic_value_sup$base_line,
            genetic_value_inf$base_line
          ),
          VA = c(genetic_value_sup$VA, genetic_value_inf$VA)
        )
      cat("\nGenetic Correlation \n")
      print(cor(genetic_value$base_line)[1, 2])
    } else {
      genetic_value <-
        base_line_multi_traits(
          seed = seed,
          set_cor = set_cor,
          ntraits = ntraits,
          correlation = correlation,
          additive_object = QTN$additive_effect_trait_object,
          epistatic_object = QTN$epistatic_effect_trait_object,
          additive_effect = additive_effect,
          epistatic_effect = epistatic_effect,
          big_additive_QTN_effect = big_additive_QTN_effect
        )
      if (set_cor) {
        cat("Populational Correlation \n")
        colnames(correlation) <-
          c("Target", c(paste0("Trait_", 2:ntraits)))
        rownames(correlation) <-
          c("Target", c(paste0("Trait_", 2:ntraits)))
        print(correlation)
      }
      cat("\nSample Correlation \n")
      sample_cor <- genetic_value$sample_cor
      colnames(sample_cor) <-
        c("Target", c(paste0("Trait_", 2:ntraits)))
      rownames(sample_cor) <-
        c("Target", c(paste0("Trait_", 2:ntraits)))
      print(sample_cor)
    }
    phenotypes(
      seed = seed,
      base_line_trait = genetic_value$base_line,
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
    if (out_geno != "gds") {
      unlink(gdsfile, force = TRUE)
      }
  }