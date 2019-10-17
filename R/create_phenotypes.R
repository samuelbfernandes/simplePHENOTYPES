#' Simulation of ntraits phenotypes based on a SNP file
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes_object Marker dataset loaded as an R object.
#' Currently either HapMap or numericalized files (e.g. `data("SNP55K_maize282_maf04")`) are accepted. 
#' Only one of `genotypes_object`, `genotypes_file` or `genotypes_path` should be provided.
#' @param genotypes_file Name of a marker data set to be read from file.
#' @param genotypes_path Path to a folder containing the marker dataset file/files 
#' (e.g. separated by chromosome).
#' @param nrows Option for loading only part of a dataset. Please see data.table::fread for details.
#' @param na_string Sets missing data as "NA".
#' @param shared_name If `genotypes_path` points to a folder with files other than 
#' the marker dataset, a part of the dataset name may be used to select the desired files
#' (e.g. shared_name = "Chr" would read files Chr1.hmp.txt, ..., Chr10.hmp.txt but not HapMap.hmp.txt).
#' @param maf_cutoff Optional filter for minor allele frequency 
#' (The dataset will be filtered. Not to be confunded with the constrain option
#' which will only filter possible QTNs).
#' @param SNP_effect Following GAPIT implementation. Default 'Add'.
#' @param SNP_impute Following GAPIT implementation. Default 'Middle'.
#' @param major_allele_zero Following GAPIT implementation. Default FALSE.
#' @param additive_QTN_number Number of additive quantitative trait nucleotide
#' to be simulated.
#' @param dominance_QTN_number Number of dominance quantitative trait nucleotide
#' to be simulated.
#' @param epistatic_QTN_number Number of epistatic (additive x additive)
#' quantitative trait nucleotide to be simulated.
#' @param additive_effect Additive effect size to be simulated. If `ntraits` > 1, 
#' additive_effect should have length equals to `ntraits`. For multiple QTNs a
#' geometric series is simulated, i.e. if the additive_effect = 0.5, 
#' the effect size of the first QTNs is 0.2, and the effect size of the second is 0.5^2 and the 
#' effect of the n^th QTN is 0.5^n.
#' @param epistatic_effect Epistatic (additive x additive) effect size to be
#' simulated. Follows a geometric series similar to the additive_effect. If `ntraits` > 1, 
#' epistatic_effect should have length equals to `ntraits`.
#' @param big_additive_QTN_effect Additive effect size for one possible major
#' effect QTN. If `ntraits` > 1, big_additive_QTN_effect should have length equals `ntraits`.
#' @param architecture Genetic architecture to be simulated. Should be provided if `ntraits` > 1. Possible options are:
#' 'pleiotropic', for traits being controled by the same QTNs;
#' 'partially', for traits being controled by pleiotropic and trait specific QTNs;
#' 'LD', for traits being exclusively controled QTNs in linkage disequilibrium (controled by parameter `ld`). 
#' Currently the only option for `architecture = "LD"` is `ntraits = 2`.
#' @param overlap Number of pleiotropic additive QTNs if architecture = 'partially'.
#' @param overlap_e Number of pleiotropic epistatic QTNs if architecture = 'partially'.
#' @param specific_QTN_number Number of trait specific additive QTNs if
#' architecture = 'partially'. It should have length equals to `ntraits`.
#' @param specific_e_QTN_number Number of trait specific epistatic QTNs if
#' architecture = 'partially'. It should have length equals to `ntraits`.
#' @param ld Linkage disequilibrium between selected marker two adjacent markers
#' to be used as QTN. Default is ld = 05.
#' @param rep Number of experiments to be simulated.
#' @param vary_QTN TRUE or FALSE ... If rep_by = "QTN" (default), at each replication it selects a differet set of SNPs to be the QTNs. 
#' If rep_by = "experiment", the same set of QTNs are used to generate phenotypes on the "rep" replications.
#' @param export_gt If TRUE genotypes of selected QTNs will be saved at file. If FALSE (default), only the QTN information will be saved.
#' @param ntraits Number of traits to be simulated under pleitropic,
#' partially and LD models. If not assigned, a single trait will be simulated. 
#' The only option for the LD architecture is two.
#' @param h2 Heritability of all traits being simulated. If could be either a vector with length equals to `ntraits`, 
#' or a matrix with ncol equals to ntraits. If the later is used, the simulation will loop over the number of rows and will
#' generate a result row. If a single trait is being simulated and h2 is a vector, one simulation of each value heritability will be
#' conducted.
#' @param correlation Trait correlation matrix to be simulated. 
#' Should have nrow == ncol == ntraits.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 10e5) will be used.
#' @param home_dir Home directory. Default is current working directory.
#' @param output_dir Name to be used to create a folder and save output files.
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
#' with option architecture = "LD". Default is NULL.
#' @param constrains Set constrains for QTN selection. 
#' Currently only minor allelic frequency is implemented. 
#' @param model Options: "A" (additive), "D" (dominance), "E" (epistatic) 
#' as well as any combination of those models such as "AE", "DE" or "ADE".
#' @param sim_method Either "geometric" or "user".
#' @param dominance_effect ...
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
           nrows = Inf,
           na_string = "NA",
           shared_name = NULL,
           maf_cutoff = NULL,
           SNP_effect = "Add",
           SNP_impute = "Middle",
           major_allele_zero = FALSE,
           model = NULL,
           additive_QTN_number = NULL,
           dominance_QTN_number = NULL,
           epistatic_QTN_number = NULL,
           sim_method = "geometric",
           additive_effect = NULL,
           dominance_effect = NULL,
           epistatic_effect = NULL,
           same_add_dom_QTN = TRUE,
           average_degree_of_dom = NULL,
           big_additive_QTN_effect = NULL,
           architecture = "pleiotropic",
           overlap = NULL,
           overlap_e = NULL,
           specific_QTN_number = NULL,
           specific_e_QTN_number = NULL,
           ld = 0.5,
           rep = NULL,
           vary_QTN = TRUE,
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
    
    if(vary_QTN) {rep_by <-  "QTN"} else {rep_by <- "experiment"}
    
    if (!is.null(model)){
      if (!(grepl("A", model) |
          grepl("D", model) |
          grepl("E", model))) {
        
      }
    }else{
      stop("Please assign a \'model\' parameter. Options:\'A\', \'D\', \'E\' or combinations such as \'ADE\'.", call. = F)
    }
    if (grepl("A", model)) { add <- TRUE } else {add <- FALSE  }
    if (grepl("D", model)) { dom <- TRUE } else {dom <- FALSE  }
    if (grepl("E", model)) { epi <- TRUE } else {epi <- FALSE  }
    if (grepl("E", model) &
        is.null(epistatic_QTN_number) &
        is.null(overlap_e) &
        is.null(specific_e_QTN_number))
      stop("In order to simulate an epistatic model, please provide either \'epistatic_QTN_number\' or \'overlap_e\' and \'specific_e_QTN_number\'.", call. = F)
    if (is.null(out_geno)) out_geno <- "none"
    if (architecture == "LD" ) ntraits <- 2
    if (sim_method != "geometric") {
      if (ntraits > 1){
        if (add){
          if (is.matrix(additive_effect)) {
            if (is.null(big_additive_QTN_effect)) {
              if (nrow(additive_effect) != additive_QTN_number |
                  ncol(additive_effect) != ntraits){
                stop("Please provide an \'additive_effect\' matrix of nrow = \'additive_QTN_number\' and ncol = \'ntraits\'.", call. = F)
              }
            } else {
              if (nrow(additive_effect) != (additive_QTN_number-1) |
                  ncol(additive_effect) != ntraits){
                stop("Please provide an \'additive_effect\' matrix of nrow = (\'additive_QTN_number\'- 1) and ncol = \'ntraits\'.", call. = F)
              }
            }
          } else {
            stop("Please provide an \'additive_effect\' matrix of nrow = \'additive_QTN_number\' and ncol = \'ntraits\'.", call. = F)
          }
        }
        if (dom){
          if (is.matrix(dominance_effect)) {
              if (nrow(dominance_effect) != dominance_QTN_number |
                  ncol(dominance_effect) != ntraits){
                stop("Please provide an \'dominance_effect\' matrix of nrow = \'dominance_QTN_number\' and ncol = \'ntraits\'.", call. = F)
              }
          } else {
            stop("Please provide an \'dominance_effect\' matrix of nrow = \'dominance_QTN_number\' and ncol = \'ntraits\'.", call. = F)
          }
        }
        if (epi){
          if (is.matrix(epistatic_effect)) {
            if (nrow(epistatic_effect) != epistatic_QTN_number |
                ncol(epistatic_effect) != ntraits){
              stop("Please provide an \'epistatic_effect\' matrix of nrow = \'epistatic_QTN_number\' and ncol = \'ntraits\'.", call. = F)
            }
          } else {
            stop("Please provide an \'epistatic_effect\' matrix of nrow = \'epistatic_QTN_number\' and ncol = \'ntraits\'.", call. = F)
          }
        }
      } else {
        if (add){
          if (is.vector(additive_effect) |
              is.matrix(additive_effect)) {
            additive_effect <- matrix(additive_effect, ncol = 1)
          }
          if(nrow(additive_effect) > additive_QTN_number)
            stop("Please provide an \'additive_effect\' object of nrow/length = \'additive_QTN_number\'.", call. = F)
        }
        if (dom){
          if (is.vector(dominance_effect) |
              is.matrix(dominance_effect)) {
            dominance_effect <- matrix(dominance_effect, ncol = 1)
          }
          if(nrow(dominance_effect) > dominance_QTN_number)
            stop("Please provide a \'dominance_effect\' object of nrow/length  = \'dominance_QTN_number\'.", call. = F)
        }
        if (epi){
          if (is.vector(epistatic_effect) |
              is.matrix(epistatic_effect)) {
            epistatic_effect <- matrix(epistatic_effect, ncol = 1)
          }
          if(nrow(epistatic_effect) > epistatic_QTN_number)
            stop("Please provide an \'epistatic_effect\' object of nrow/length = \'epistatic_QTN_number\'.", call. = F)
        }
      }
    } else {
      if (add){
        if (is.vector(additive_effect) |
            is.matrix(additive_effect)) {
          additive_effect <- matrix(additive_effect, nrow = 1)
        }
        if(ncol(additive_effect) != ntraits)
          stop("sim_method = \'geometric\' requires an \'additive_effect\' object of ncol/length = \'ntraits\'.", call. = F)
      }
      if (dom){
        if (is.vector(dominance_effect) |
            is.matrix(dominance_effect)) {
          dominance_effect <- matrix(dominance_effect, nrow = 1)
        }
        if(ncol(dominance_effect) != ntraits)
          stop("sim_method = \'geometric\' requires a \'dominance_effect\' object of ncol/length  = \'ntraits\'.", call. = F)
      }
      if (epi){
        if (is.vector(epistatic_effect) |
            is.matrix(epistatic_effect)) {
          epistatic_effect <- matrix(epistatic_effect, nrow = 1)
        }
        if(ncol(epistatic_effect) > ntraits)
          stop("sim_method = \'geometric\' requires an \'epistatic_effect\' object of ncol/length = \'ntraits\'.", call. = F)
      }
    }
    if (architecture == "LD" ) {
      if (length(h2) == 1) h2 <- cbind(h2, h2)
      if (length(h2) > 2) h2 <- matrix(h2, ncol = 1 )
      if (is.matrix(h2)){
        if (ncol(h2)>2) 
          stop("In the LD architecture the parameter \'h2\' should two columns, one for each trait!", call. = F)
        if (ncol(h2) ==1 ) h2 <- cbind(h2, h2)
      }
      if (length(big_additive_QTN_effect)!= 2 |
          ncol(additive_effect)!= 2) {
        stop("In the LD architecture, parameters \'additive_effect\' and \'big_additive_QTN_effect\' should have length/ncol = 2", call. = F)
      }
    }
    if (is.vector(h2)) {
      if(ntraits > 1){
        h2 <- matrix(h2, nrow = 1)
      } else {
        h2 <- matrix(h2, ncol = 1)
      }
    }
    colnames(h2) <- paste0("Trait_", 1:ntraits)
    if (rep_by != "QTN" & rep_by != "experiment")
      stop("Parameter \'rep_by\' should be either \'QTN\' or \'experiment\'!", call. = F)
    if (sim_method != "geometric" & sim_method != "user")
      stop("Parameter \'sim_method\' should be either \'geometric\' or \'user\'!", call. = F)
    if (sum(c(!is.null(genotypes_object),
              !is.null(genotypes_file),
              !is.null(genotypes_path))) > 1)
      stop("Only one of `genotypes_object`, `genotypes_file` or `genotypes_path` should be provided!", call. = F)
    if (ntraits > 1) {
      if (ntraits != ncol(h2)) {
        stop("Parameter \'h2\' should have length/ncol = \'ntraits\'", call. = F)
      }
      mm <-
        ifelse(
          architecture == "pleiotropic",
          "Pleiotropic",
          ifelse(
            architecture == "partially",
            "Partially Pleiotropic",
            ifelse(architecture == "LD", "Linkage Disequilibrium", {
              stop("The genetic architecture used is not valid! 
                   Please choose one of: \'pleiotropic\', \'partially\' or \'LD\' ", call. = F)
            })
          )
        )
      if (!is.null(correlation) & !all(ntraits == dim(correlation))) {
        stop("ntraits do not match with dim of correlation matrix!", call. = F)
      }
      if (architecture == "partially"){
        if  (!is.null(specific_QTN_number) |
             !is.null(overlap)) {
          if (any(is.null(specific_QTN_number),
                  is.null(overlap)))
            stop("For Partially Pleiotropic additive architecture, both \'overlap\' and \'specific_QTN_number\' must be provided!", call. = F)
          if (ntraits != length(specific_QTN_number)) 
            stop("Parameter [additive] \'specific_QTN_number\' should have length = \'ntraits\'", call. = F)
        }
        if (!is.null(specific_e_QTN_number) |
            !is.null(overlap_e)) {
          if (any(is.null(specific_e_QTN_number),
                  is.null(overlap_e))) 
            stop("For Partially Pleiotropic epistatic architecture, both \'overlap_e\' and \'specific_e_QTN_number\' must be provided!", call. = F)
          if (ntraits != length(specific_e_QTN_number)) {
            stop("Parameter \'specific_e_QTN_number\' should have length = \'ntraits\'", call. = F)
          }
        }
      }
      if (!is.null(epistatic_effect) &
          ntraits != length(epistatic_effect)) {
        stop("Parameter \'epistatic_effect\' should have length = \'ntraits\'", call. = F)
      }
      # if (ntraits != length(additive_effect) |
      #     ntraits != length(big_additive_QTN_effect) ) {
      #   stop("Parameters \'additive_effect\' and \'big_additive_QTN_effect\' should have length = \'ntraits\'", call. = F)
      # }
      if (ntraits > 1 & 
          (architecture == "partially" |
          architecture == "pleiotropic")){
        if (is.null(correlation)) {
          if (add == TRUE & dom == TRUE & epi == TRUE) {
            if (var(additive_effect[1,]) == 0 &
                var(dominance_effect[1,]) == 0 &
                var(epistatic_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'additive_effect\', \'dominance_effect\', \'epistatic_effect\').", call. = F)
          } else if (add == TRUE & epi == TRUE){
            if (var(additive_effect[1,]) == 0 &
                var(epistatic_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'additive_effect\', \'epistatic_effect\').", call. = F)
          } else if (add == TRUE & dom == TRUE){
            if (var(additive_effect[1,]) == 0 &
                var(dominance_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'additive_effect\', \'dominance_effect\').", call. = F)
          } else if (dom == TRUE & epi == TRUE){
            if (var(dominance_effect[1,]) == 0 &
                var(epistatic_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'dominance_effect\', \'epistatic_effect\').", call. = F)
          } else if (add == TRUE){
            if (var(additive_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'additive_effect\').", call. = F)
          } else if (dom == TRUE){
            if (var(dominance_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'dominance_effect\').", call. = F)
          } else if (epi == TRUE){
            if (var(epistatic_effect[1,]) == 0)
              stop("In order to simulate correlated traits, please provide either a 
                   correlation matrix or different genetic effects (\'epistatic_effect\').", call. = F)
          }
        }
      }
    }
    #TODO VCF and other formats.
    input_format <- "hapmap"
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
      if (add & !dom & !epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Additive Effects\n")
      } else if (!add & dom & !epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Dominance Effects\n") 
      } else if (!add & !dom & epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Epistatic Effects \n") 
      } else if (add & dom & !epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Additive and Dominance Effects \n") 
      } else if (add & !dom & epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Additive and Epistatic Effects \n")
      } else if (!add & dom & epi) {
        cat("Simulation of a", mm, "Genetic Architecture with Dominance and Epistatic Effects \n")
      } else {
        cat("Simulation of a", mm, "Genetic Architecture with Additive, Dominance and Epistatic Effects \n")
      }
    } else {
      if (add & !dom & !epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Additive Effects\n")
      } else if (!add & dom & !epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Dominance Effects\n") 
      } else if (!add & !dom & epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Epistatic Effects \n") 
      } else if (add & dom & !epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Additive and Dominance Effects \n") 
      } else if (add & !dom & epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Additive and Epistatic Effects \n")
      } else if (!add & dom & epi) {
        cat("Simulation of a Single Trait Genetic Architecture with Dominance and Epistatic Effects \n")
      } else {
        cat("Simulation of a Single Trait Genetic Architecture with Additive, Dominance and Epistatic Effects \n")
      }
    }
    cat("\nSIMULATION PARAMETERS: \n\n")
    cat("Number of traits:", ntraits)
    if (architecture == "pleiotropic" | architecture == "LD" | ntraits == 1 ) {
      if (add) cat("\nNumber of additive QTNs:", additive_QTN_number)
      if (dom) cat("\nNumber of dominance QTNs:", dominance_QTN_number)
      if (epi) cat("\nNumber of epistatic QTNs:", epistatic_QTN_number)
    }
    if (architecture == "partially") {
      if (add) {
        cat("\nNumber of pleiotropic additive QTNs:", overlap)
        cat("\nNumber of trait specific additive QTNs:",
            paste(specific_QTN_number,
                  collapse = ", ")
        )
      }
      if (dom) {
        cat("\nNumber of pleiotropic dominant QTNs:", overlap)
        cat("\nNumber of trait specific dominant QTNs:",
            paste(specific_QTN_number,
                  collapse = ", "))
      }
      if (epi) {
        cat("\nNumber of pleiotropic epistatic QTNs:", overlap_e)
        cat("\nNumber of trait specific epistatic QTNs:",
            paste(specific_e_QTN_number,
                  collapse = ", "))
      }
    }
    cat("\nReplication by:", rep_by)
    if(ntraits == 1) {
      cat("\nPopulation Heritability:", h2)
      if (add) {
        colnames(additive_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nAdditive genetic effect:\n")
        print(additive_effect)
        cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
      }
      if (dom) {
        colnames(dominance_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nDominance genetic effect:\n")
        print(dominance_effect)
      }
      if (epi) {
        colnames(epistatic_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nEpistatic genetic effect:\n")
        print(epistatic_effect)
      }
      cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
    } else {
      if(nrow(h2) > 1){
        cat("\nPopulation Heritability:")
        #colnames(h2) <- paste0("Trait_", 1:ntraits)
        print(h2)
      } else {
        cat("\nPopulation Heritability:\n")
        print(h2)
      }
      if (add) {
        colnames(additive_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nAdditive genetic effects:\n")
        print(additive_effect)
        cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
      }
      if (dom) {
        colnames(dominance_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nDominance genetic effects:\n")
        print(dominance_effect)
      }
      if (epi) {
        colnames(epistatic_effect) <- paste0("Trait_", 1:ntraits)
        cat("\nEpistatic genetic effects:\n")
        print(epistatic_effect)
      }
      cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
    }
    if (architecture == "LD" |
        out_geno == "plink" |
        out_geno == "gds" |
        output_format == "gemma"){
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
    cat("\n\nDIAGNOSTICS:\n\n")
    if (ntraits == 1 | !any(architecture != "pleiotropic")) {
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
    if (ntraits > 1 & !any(architecture != "partially")) {
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
    if (ntraits > 1 & !any(architecture != "LD")) {
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
      if (dom) {
        if (same_add_dom_QTN) {
          genetic_value <-
            base_line_single_trait(
              seed = seed,
              additive_object = QTN$additive_effect_trait_object,
              dominance_object = QTN$additive_effect_trait_object,
              epistatic_object = QTN$epistatic_effect_trait_object,
              additive_effect = additive_effect,
              dominance_effect = dominance_effect,
              epistatic_effect = epistatic_effect,
              big_additive_QTN_effect = big_additive_QTN_effect,
              rep = rep,
              rep_by = rep_by,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method
            )
        }else{
          genetic_value <-
            base_line_single_trait(
              seed = seed,
              additive_object = QTN$additive_effect_trait_object,
              dominance_object = QTN$dominance_effect_trait_object,
              epistatic_object = QTN$epistatic_effect_trait_object,
              additive_effect = additive_effect,
              dominance_effect = dominance_effect,
              epistatic_effect = epistatic_effect,
              big_additive_QTN_effect = big_additive_QTN_effect,
              rep = rep,
              rep_by = rep_by,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method
            )
        }
      } else {
        genetic_value <-
          base_line_single_trait(
            seed = seed,
            additive_object = QTN$additive_effect_trait_object,
            epistatic_object = QTN$epistatic_effect_trait_object,
            additive_effect = additive_effect,
            epistatic_effect = epistatic_effect,
            big_additive_QTN_effect = big_additive_QTN_effect,
            rep = rep,
            rep_by = rep_by,
            add = add,
            dom = dom,
            epi = epi,
            sim_method = sim_method
          )
      }
    } else {
      if(dom){
        #including same_add_dom_QTN
        if(same_add_dom_QTN){
          genetic_value <-
            base_line_multi_traits(
              seed = seed,
              ntraits = ntraits,
              correlation = correlation,
              additive_object = QTN$additive_effect_trait_object,
              dominance_object = QTN$additive_effect_trait_object,
              epistatic_object = QTN$epistatic_effect_trait_object,
              additive_effect = additive_effect,
              dominance_effect = dominance_effect,
              epistatic_effect = epistatic_effect,
              big_additive_QTN_effect = big_additive_QTN_effect,
              rep = rep,
              rep_by = rep_by,
              architecture = architecture,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method
            )
        }else{
          genetic_value <-
            base_line_multi_traits(
              seed = seed,
              ntraits = ntraits,
              correlation = correlation,
              additive_object = QTN$additive_effect_trait_object,
              dominance_object = QTN$dominance_effect_trait_object,
              epistatic_object = QTN$epistatic_effect_trait_object,
              additive_effect = additive_effect,
              dominance_effect = dominance_effect,
              epistatic_effect = epistatic_effect,
              big_additive_QTN_effect = big_additive_QTN_effect,
              rep = rep,
              rep_by = rep_by,
              architecture = architecture,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method
            )
        }
      }else{
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
            architecture = architecture,
            add = add,
            dom = dom,
            epi = epi,
            sim_method = sim_method
          )
      }
      
      

      
      
      
      
      
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