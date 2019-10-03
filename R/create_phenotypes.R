#' Simulation of ntraits phenotypes based on a SNP file
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes_object Marker data set loaded as an R object 
#' (currently accepting either HapMap or numericalized files). 
#' Only one of `genotypes_object`, `genotypes_file` or `genotypes_path` should be used.
#' @param genotypes_file Name of a marker data set to be read from file.
#' @param genotypes_path Path to a folder with multiple marker data set files 
#' (e.g. separated by chromosome).
#' @param input_format Currently only "hapmap" is implemented.
#' @param nrows Please see data.table::fread for details.
#' @param na_string Sets missing data as "NA".
#' @param shared_name If files are saved by chromosome and in a folder with 
#' other files, a portion of the name may be used to match files 
#' (e.g. shared_name = "Chr" for files Chr1.hmp.txt, ..., Chr10.hmp.txt).
#' @param maf_cutoff Optional filter for minor allele frequency 
#' (The dataset will be filtered. Not to be confunded with the constrain option
#' which will only filter possible QTNs).
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
#' @param rep_by Simulate same QTNs but with different residuals ('residuals'). 
#' Simulate each replication with different QTNs (Defalut 'QTN').
#' @param export_gt if genotypes of selected QTNs should be exported.
#' @param ntraits Number of traits to be simulated under pleitropic and
#' partially pleiotropic models. The default for linkage disequilibrium models is two.
#' @param h2 Heritability of trait 1 (possible the target trait) if multiple traits are simulated.
#' If a vector, the simulation will "loop" over it and generate one file for each
#' combination of genetic settings.
#' @param correlation Trait correlation matrix to be simulated. 
#' Should have nrow == ncol == ntraits.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 10e5) will be used.
#' @param home_dir Home directory. Default is current working directory.
#' @param output_dir Name to be used to create folder and save output files.
#' @param to_r Option for saving simulated results to R in addition to saving to file. 
#' If TRUE, results need to be assinged to an R object (see vign).
#' @param output_format Four options for saving outputs: 'multi-file', 
#' saves one simulation setting in a separate file;
#' 'long' (default for multiple traits), appends each experiment (rep) to the last one (by row); 'wide', saves
#' experiments by column (default for single trait) and 'gemma', saves .fam files 
#' to be used by gemma with plink bed files (renaming .fam file might be necessary).
#' @param out_geno Saves numericalized genotype either as "numeric", "plink" or "gds".
#' Default is NULL.
#' @param gdsfile gds file (in case there is one already created) to be used
#' with option model = "LD". Default is NULL.
#' @param constrains Set constrains for QTN selection. 
#' Currently only minor allelic frequency is implemented. 
#' Eiter one or both of the following options may be non-null: 'list(maf_above = NULL, maf_below = NULL)'.
#' @return Numericalized marker dataset, selected QTNs, phenotypes for 'ntraits' traits.
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Aug 2nd, 2019
#' @examples
#' # Simulate 50 replications of a single phenotype.
#'
#' pheno <- create_phenotypes(
#'   genotypes_object = SNP55K_maize282_maf04,
#'   additive_QTN_number = 3,
#'   additive_effect = 0.1,
#'   big_additive_QTN_effect = 0.9,
#'   rep = 10,
#'   h2 = 0.7,
#'   to_r = TRUE
#' )
#'
create_phenotypes <-
  function(genotypes_object = NULL,
           genotypes_file = NULL,
           genotypes_path = NULL,
           input_format = "hapmap",
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
           rep_by = "QTN",
           export_gt = FALSE,
           ntraits = 1,
           h2 = NULL,
           correlation = NULL,
           seed = NULL,
           home_dir = getwd(),
           output_dir = NULL,
           to_r = FALSE,
           output_format = "long",
           out_geno = NULL,
           gdsfile = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL)) {
    # -------------------------------------------------------------------------
    .onAttach <- function(libname, simplePHENOTYPES) {
      packageStartupMessage("Thank you for using the simplePHENOTYPES package!")
    }
    .onAttach()
    if (is.null(out_geno)) out_geno <- "none"
    if (model == "LD" ) {
      ntraits <- 2
      if (length(h2) == 1) h2 <- cbind(h2, h2)
      if (length(h2) > 2) h2 <- matrix(h2, ncol = 1 )
      if (is.matrix(h2)){
        if (ncol(h2)>2) 
          stop("In the LD model the parameter \'h2\' should two columns, one for each trait!", call. = F)
        if (ncol(h2) ==1 ) h2 <- cbind(h2, h2)
      }
      if (length(big_additive_QTN_effect)!= 2 |
          length(additive_effect)!= 2) {
        stop("In the LD model, parameters \'additive_effect\' and \'big_additive_QTN_effect\' should have length == 2", call. = F)
      }
    }
    if (is.vector(h2)) {
      if(ntraits > 1){
        h2 <- matrix(h2, nrow = 1)
      } else {
        h2 <- matrix(h2, ncol = 1)
      }
    }
    if (rep_by != "QTN" &
        rep_by != "experiment")
      stop("Parameter \'rep_by\' should be either \'QTN\' or \'experiment\'!", call. = F)
    if (ntraits > 1) {
      if (ntraits != ncol(h2)) {
        stop("Parameter \'h2\' should have length == \'ntraits\'", call. = F)
      }
      mm <-
        ifelse(
          model == "pleiotropic",
          "Pleiotropic",
          ifelse(
            model == "partially",
            "Partially Pleiotropic",
            ifelse(model == "LD", "Linkage Disequilibrium", {
              stop("The genetic model used is not valid! 
                   Please choose one of: \'pleiotropic\', \'partially\' or \'LD\' ", call. = F)
            })
          )
        )
      if (!is.null(correlation) & !all(ntraits == dim(correlation))) {
        stop("ntraits do not match with dim of correlation matrix!", call. = F)
      }
      if (model == "partially"){
        if  (!is.null(specific_QTN_number) |
             !is.null(overlap)) {
          if (any(is.null(specific_QTN_number),
                  is.null(overlap)))
            stop("For Partially Pleiotropic additive model, both \'overlap\' and \'specific_QTN_number\' must be provided!", call. = F)
          if (ntraits != length(specific_QTN_number)) 
            stop("Parameter [additive] \'specific_QTN_number\' should have length == \'ntraits\'", call. = F)
        }
        if (!is.null(specific_e_QTN_number) |
            !is.null(overlap_e)) {
          if (any(is.null(specific_e_QTN_number),
                  is.null(overlap_e))) 
            stop("For Partially Pleiotropic epistatic model, both \'overlap_e\' and \'specific_e_QTN_number\' must be provided!", call. = F)
          if (ntraits != length(specific_e_QTN_number)) {
            stop("Parameter \'specific_e_QTN_number\' should have length == \'ntraits\'", call. = F)
          }
        }
      }
      if (!is.null(epistatic_effect) &
          ntraits != length(epistatic_effect)) {
        stop("Parameter \'epistatic_effect\' should have length == \'ntraits\'", call. = F)
      }
      if (ntraits != length(additive_effect) |
          ntraits != length(big_additive_QTN_effect) ) {
        stop("Parameters \'additive_effect\' and \'big_additive_QTN_effect\' should have length == \'ntraits\'", call. = F)
      }
      if (model == "partially" |
          model == "pleiotropic"){
        if (is.null(correlation) &
            var(additive_effect) == 0) {
          stop("In order to simulate correlated traits, please provide either a 
               correlation matrix or different genetic effects.", call. = F)
        }
      }
    }
    setwd(home_dir)
    on.exit(setwd(home_dir), add = TRUE)
    if (!is.null(genotypes_path) |
        !is.null(genotypes_file) |
        (class(unlist(genotypes_object[, 12])) != "numeric" &
         class(unlist(genotypes_object[, 12])) != "integer")) {
      genotypes_object <-
        genotypes(genotypes_object = genotypes_object,
                  genotypes_path = genotypes_path,
                  genotypes_file = genotypes_file,
                  input_format = input_format,
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
    if (ntraits > 1) {
      if (!is.null(epistatic_QTN_number)) {
        cat("Simulation of a", mm, ", Additive and Epistatic Genetic Model \n") 
      } else{
        cat("Simulation of a", mm, "Additive Genetic Model \n") 
      }
    } else {
      if (!is.null(epistatic_QTN_number)) {
        cat("Simulation of a Single Trait, Additive and Epistatic Genetic Model \n") 
      } else {
        cat("Simulation of a Single Trait Additive Genetic Model \n") 
      }
    }
    cat("\nSIMULATION PARAMETERS: \n\n")
    cat("Number of traits:", ntraits)
    if (model == "pleiotropic" | model == "LD" | ntraits == 1 ) {
      cat("\nNumber of additive QTNs:",
          additive_QTN_number,
          "\nNumber of epistatic QTNs:",
          ifelse(is.null(epistatic_QTN_number), 0, epistatic_QTN_number )
      )
    }
    if (model == "partially") {
      if (is.null(overlap_e) | is.null(specific_e_QTN_number)) {
        cat("\nNumber of pleiotropic additive QTNs:",
            overlap,
            "\nNumber of trait specific additive QTNs:",
            paste(specific_QTN_number,
                  collapse = ", ")
        )
      } else {
        cat("\nNumber of pleiotropic additive QTNs:",
            overlap,
            "\nNumber of trait specific additive QTNs:",
            paste(specific_QTN_number,
                  collapse = ", "),
            "\nNumber of pleiotropic epistatic QTNs:",
            overlap_e,
            "\nNumber of trait specific epistatic QTNs:",
            paste(specific_e_QTN_number,
                  collapse = ", ")
        )
      }
    }
    if(ntraits == 1) {
      cat("\nPopulation Heritability:", h2)
      cat("\nAdditive genetic effect:", additive_effect)
      cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
      if (!is.null(epistatic_QTN_number)) {
        cat("\nEpistatic genetic effect:", epistatic_effect)
      }
      cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
    } else {
      if(nrow(h2) > 1){
        cat("\nPopulation Heritability:")
        colnames(h2) <- paste0("Trait_", 1:ntraits)
        print(h2)
      } else {
        cat("\nPopulation Heritability", paste0("(Traits 1:", ntraits,", respectively)") , h2)
      }
      cat("\nAdditive genetic effects:", additive_effect)
      cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
      if (!is.null(epistatic_QTN_number)) {
        cat("\nEpistatic genetic effects:", epistatic_effect)
      }
      cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
    }
    if (model == "LD" |
        out_geno == "plink" |
        out_geno == "gds" |
        output_format == "gemma"){
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
      dup <- duplicated(genotypes_object$snp)
      if(any(dup)){
        message("Removing ", sum(dup), " markers for being duplicated!")
        genotypes_object <- genotypes_object[!dup, ]      
      }
      if (any(grepl("geno.gds", dir(home_dir)))){
        message("A file named geno.gds is already present in this folder, creating \'geno2\'.")
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
      if (out_geno == "gds") cat("GDS files saved at:", home_dir, "\n")
      if (out_geno == "plink" |
          output_format == "gemma") {
        genofile <- SNPRelate::snpgdsOpen(gdsfile)
        snpset <-
          SNPRelate::snpgdsSelectSNP(genofile, remove.monosnp=F, verbose=F)
        SNPRelate::snpgdsGDS2BED(genofile, bed.fn = "geno", snp.id = snpset, verbose=F)
        SNPRelate::snpgdsClose(genofile)
        cat("\nPlink bed files saved at:", home_dir, "\n")
      }
    }
    if (output_format == "gemma") {
      fam <- data.table::fread("geno.fam", drop = 6)
      fam$V1 <- fam$V2
    }
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
    }
    if (ntraits > 1) {
      if (ntraits != length(big_additive_QTN_effect)) {
        if (length(big_additive_QTN_effect) > ntraits) {
          cat("\nLength of big_additive_QTN_effect > ntraits! using the first",
              ntraits,
              "(",
              paste(big_additive_QTN_effect[1:ntraits],
                    collapse = ", "),
              ")",
              "values\n"
          )
          big_additive_QTN_effect <-
            big_additive_QTN_effect[1:ntraits]
        } else {
          cat("\nLength big_additive_QTN_effect < ntraits!",
              "Setting all traits with big_additive_QTN_effect = ",
              big_additive_QTN_effect[1],
              "\n"
          )
          big_additive_QTN_effect <-
            rep(big_additive_QTN_effect[1], ntraits)
        }
      }
    }
    cat("\n\nDIAGNOSTICS:\n\n")
    if (ntraits == 1 | !any(model != "pleiotropic")) {
      QTN <-
        QTN_pleiotropic(
          genotypes = genotypes_object,
          seed = seed,
          additive_QTN_number = additive_QTN_number,
          epistatic_QTN_number = epistatic_QTN_number,
          constrains = constrains,
          rep = rep,
          rep_by = rep_by,
          export_gt = export_gt
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
          ntraits = ntraits,
          constrains = constrains,
          rep = rep,
          rep_by = rep_by,
          export_gt = export_gt
        )
    }
    if (ntraits > 1 & !any(model != "LD")) {
      QTN <-
        QTN_linkage(
          genotypes = genotypes_object,
          seed = seed,
          additive_QTN_number = additive_QTN_number,
          ld = ld,
          gdsfile = gdsfile,
          constrains = constrains,
          rep = rep,
          rep_by = rep_by,
          export_gt = export_gt
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
          big_additive_QTN_effect = big_additive_QTN_effect,
          rep = rep,
          rep_by = rep_by
        )
    } else {
      genetic_value <-
        base_line_multi_traits(
          seed = seed,
          ntraits = ntraits,
          correlation = correlation,
          additive_object = QTN$additive_effect_trait_object,
          epistatic_object = QTN$epistatic_effect_trait_object,
          additive_effect = additive_effect,
          epistatic_effect = epistatic_effect,
          big_additive_QTN_effect = big_additive_QTN_effect,
          rep = rep,
          rep_by = rep_by,
          model = model
        )
      #tirar a media das correlacoes quando rep_by = "QTN"
      
      if (!is.null(correlation)) {
        cat("Population Correlation \n")
        colnames(correlation) <- paste0("Trait_", 1:ntraits)
        rownames(correlation) <- paste0("Trait_", 1:ntraits)
        print(correlation)
      }
      if (rep_by == "QTN") {
        sample_cor <- matrix(0, ntraits, ntraits)
        for(v in 1:rep) {
          sample_cor <- (sample_cor + genetic_value[[v]]$sample_cor)
        }
        sample_cor <- sample_cor/rep
      } else {
        sample_cor <- genetic_value[[1]]$sample_cor
      }
      cat("\nSample Correlation \n")
      colnames(sample_cor) <- paste0("Trait_", 1:ntraits)
      rownames(sample_cor) <- paste0("Trait_", 1:ntraits)
      print(sample_cor)
    }
    results <- phenotypes(
      seed = seed,
      base_line_trait = genetic_value,
      h2 = h2,
      rep = rep,
      ntraits = ntraits,
      output_format = output_format,
      fam = fam,
      to_r = to_r,
      rep_by = rep_by
    )
    cat("\n\nResults are saved at:", getwd())
    sink()
    close(zz)
    file.show(paste0(path_out, "/Log_Sim.txt"))
    if (out_geno != "gds") {
      unlink(gdsfile, force = TRUE)
    }
    if(to_r) return(results)
  }