#' Simulation of ntraits phenotypes based on a SNP file
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param genotypes_object Marker dataset loaded as an R object.
#' Currently either HapMap or numericalized files 
#' (e.g. `data("SNP55K_maize282_maf04")`) are accepted. Only one of 
#' `genotypes_object`, `genotypes_file` or `genotypes_path` should be provided.
#' @param genotypes_file Name of a marker data set to be read from file.
#' @param genotypes_path Path to a folder containing the marker dataset 
#' file/files (e.g. separated by chromosome).
#' @param rep Number of experiments to be simulated.
#' @param ntraits Number of traits to be simulated under pleitropic,
#' partially and LD architectures (see below). If not assigned, 
#' a single trait will be simulated. Currently, the only option for the 
#' LD architecture is `ntraits = 2`.
#' @param h2 Heritability of all traits being simulated. 
#' It could be either a vector with length equals to `ntraits`, 
#' or a matrix with ncol equals to ntraits. If the later is used, the simulation
#' will loop over the number of rows and will generate a result for each row. 
#' If a single trait is being simulated and h2 is a vector, 
#' one simulation of each heritability value will be conducted.
#' @param model The genetic model to be assumed. The options are:
#' "A" (additive), "D" (dominance), "E" (epistatic) 
#' as well as any combination of those models such as "AE", "DE" or "ADE".
#' @param additive_QTN_number Number of additive quantitative trait nucleotide
#' to be simulated.
#' @param dominance_QTN_number Number of dominance quantitative trait nucleotide
#' to be simulated.
#' @param epistatic_QTN_number Number of epistatic (Currently, only additive x 
#' additive epistasis are simulated) quantitative trait nucleotide to be simulated.
#' @param additive_effect Additive effect size to be simulated. It may be either 
#' a vector (assuming `ntraits` = 1 or one allelic effect per trait) or a list 
#' of length = `ntraits`, i.e., if `ntraits` > 1, one set of additive effects 
#' should be provided for each trait. In that case, each component should be a 
#' vector of either lenght one, if `sim_method = "geometric"` (see below), or 
#' length equal to the number of additive QTNs being simulated. 
#' @param same_add_dom_QTN A bolean for having the same quantitative trait
#' nucleotide having additive and dominance effects.
#' @param dominance_effect Similar to the `additive_effect`, it could be either 
#' a vector or a list. Optional if `same_add_dom_QTN = TRUE`. 
#' @param degree_of_dominance If the same set of quantitative trait nucleotide 
#' are being used for simulating additive and dominance effects, the dominance 
#' allelic effect could be a proportion of the additive allelic effect.
#' In other words, `degree_of_dominance` equals to 0.5, 1, 1.5 will simulate,
#' partial dominance, complete dominance and overdominance, respectivelly. 
#' @param epistatic_effect Epistatic (additive x additive) effect size to be
#' simulated. Similar to the `additive_effect`, it could be either a vector or 
#' a list.
#' @param architecture Genetic architecture to be simulated. Should be provided 
#' if `ntraits` > 1. Possible options are: 'pleiotropic', for traits being 
#' controled by the same QTNs; 'partially', for traits being controled by 
#' pleiotropic and trait specific QTNs; 'LD', for traits being exclusively 
#' controled QTNs in linkage disequilibrium (controled by parameter `ld`). 
#' Currently the only option for `architecture = "LD"` is `ntraits = 2`.
#' @param pleitropic_a Number of pleiotropic additive QTNs to be 
#' used if `architecture = "partially"`.
#' @param pleitropic_d Number of pleiotropic dominance QTNs to be 
#' used if `architecture = "partially"`.
#' @param pleitropic_e Number of pleiotropic epistatic QTNs to be 
#' used if `architecture = "partially"`.
#' @param trait_specific_a_QTN_number Number of trait specific additive QTNs if
#'`architecture = "partially"`.. It should have length equals to `ntraits`.
#' @param trait_specific_d_QTN_number Number of trait specific dominance QTNs if
#' `architecture = "partially"`. It should have length equals to `ntraits`.
#' @param trait_specific_e_QTN_number Number of trait specific epistatic QTNs if
#' `architecture = "partially"`. It should have length equals to `ntraits`.
#' @param ld Linkage disequilibrium between selected marker two adjacent markers
#' to be used as QTN. Default is `ld = 0.5`.
#' @param sim_method Provide the method of simulating allelic effects. 
#' The options available are "geometric" and "custom". For multiple QTNs, 
#' a geometric series may be simulated, i.e. if the additive_effect = 0.5, 
#' the effect size of the first QTNs is 0.2, and the effect size of the second 
#' is 0.5^2 and the effect of the n^th QTN is 0.5^n.
#' @param vary_QTN A bolean to determine if the same set of quantitative trait 
#' nucleotide (QTN) should be used to generate genetic effects for each 
#' experiment (`vary_QTN = FALSE`) or  if a different set of QTNs should be 
#' used for each experiment (`vary_QTN = TRUE`).
#' @param big_additive_QTN_effect Additive effect size for one possible major
#' effect quantitative trait nucleotide. If `ntraits` > 1, 
#' big_additive_QTN_effect should have length equals `ntraits`.
#' If `additive_QTN_number` > 1, the fist QTN will have the large effect.
#' @param correlation Option to simulate traits with a pre-defined correlation. 
#' It should be a square matrix with number of rows = `ntraits`.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 10e5) will be used. Notice that at each sampling step, 
#' a different seed generated based on the `seed` parameter is used. For example, 
#' when simulating the 10th replication of trait 1, the seed to be used is 
#' `round( (seed * 10 * 10) * 1)`. On the other hand, for simulating the 21st 
#' replication of trait 2, the seed to be used will be 
#' `round( (seed * 21 * 21) * 2)`.The actual seed used in every simulation is 
#' exported along with simulated phenotypes.
#' @param export_gt If TRUE genotypes of selected QTNs will be saved at file.
#' If FALSE (default), only the QTN information will be saved.
#' @param home_dir Home directory. Default is current working directory.
#' @param output_dir Name to be used to create a folder and save output files.
#' @param to_r Option for saving simulated results into R in addition to saving 
#' it to file. If TRUE, results need to be assinged to an R object (see vignette).
#' @param output_format Four options for saving outputs: 'multi-file', 
#' saves one simulation setting in a separate file; 'long' (default for multiple
#' traits), appends each experiment (rep) to the last one (by row); 'wide', saves
#' experiments by column (default for single trait) and 'gemma', saves .fam files 
#' to be used by gemma with plink bed files (renaming .fam file might be necessary).
#' @param out_geno Saves numericalized genotype either as "numeric", "plink" or
#' "gds". Default is NULL.
#' @param gdsfile Points to a gds file (in case there is one already created) to
#' be used with option architecture = "LD". Default is NULL.
#' @param constrains Set constrains for QTN selection. Currently, only minor 
#' allelic frequency is implemented. Eiter one or both of the following options 
#' may be non-null: 'list(maf_above = NULL, maf_below = NULL)'.
#' @param prefix If `genotypes_path` points to a folder with files other than the
#' marker dataset, a part of the dataset name may be used to select the desired 
#' files (e.g. prefix = "Chr" would read files Chr1.hmp.txt, ..., Chr10.hmp.txt 
#' but not HapMap.hmp.txt).
#' @param maf_cutoff Optional filter for minor allele frequency 
#' (The dataset will be filtered. Not to be confunded with the constrain option
#' which will only filter possible QTNs).
#' @param nrows Option for loading only part of a dataset. Please see 
#' data.table::fread for details.
#' @param na_string Sets missing data as "NA".
#' @param SNP_effect Parameter used for numericalization. Following GAPIT 
#' implementation, the default is 'Add'.
#' @param SNP_impute Parameter used for numericalization. Following GAPIT 
#' implementation, the default is 'Middle'.
#' @param major_allele_zero Parameter used for numericalization. Following 
#' GAPIT implementation, the default is FALSE.
#' @return Numericalized marker dataset, selected QTNs, phenotypes for 'ntraits'
#'  traits, log file.
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Oct 29th, 2019
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
#'   to_r = TRUE,
#'   model = "A"
#' )
#'
#' # For more examples, please run the following:
#' # vignette("simplePHENOTYPES")
#'
create_phenotypes <-
  function(genotypes_object = NULL,
           genotypes_file = NULL,
           genotypes_path = NULL,
           nrows = Inf,
           na_string = "NA",
           prefix = NULL,
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
           degree_of_dominance = NULL,
           big_additive_QTN_effect = NULL,
           architecture = "pleiotropic",
           pleitropic_a = NULL,
           pleitropic_d = NULL,
           pleitropic_e = NULL,
           trait_specific_a_QTN_number = NULL,
           trait_specific_d_QTN_number = NULL,
           trait_specific_e_QTN_number = NULL,
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
    #TODO:
    # list all set.seed values
    x <- try({
      .onAttach <- function(libname, simplePHENOTYPES) {
        packageStartupMessage("Thank you for using the simplePHENOTYPES package!")
      }
      .onAttach()
      if (grepl("A", model)) { add <- TRUE } else {add <- FALSE  }
      if (grepl("D", model)) { dom <- TRUE } else {dom <- FALSE  }
      if (grepl("E", model)) { epi <- TRUE } else {epi <- FALSE  }
      if (architecture == "LD" ) ntraits <- 2
      if (is.null(out_geno)) out_geno <- "none"
      if (vary_QTN) {rep_by <-  "QTN"} else {rep_by <- "experiment"}
      if (is.vector(h2)) {
        if (ntraits > 1){
          h2 <- matrix(h2, nrow = 1)
        } else {
          h2 <- matrix(h2, ncol = 1)
        }
      }
      colnames(h2) <- paste0("Trait_", 1:ntraits)
      if (add) {
        if (is.vector(additive_effect)) {
          additive_effect <- as.list(additive_effect)
        } else if (!is.list(additive_effect)) {
          stop("\'additive_effect\' should be either a vector or a list of length = ntraits.", call. = F)
        }
      }
      if (epi) {if (is.vector(epistatic_effect)) {
        epistatic_effect <- as.list(epistatic_effect)
      } else if (!is.list(epistatic_effect)) {
        stop("\'epistatic_effect\' should be either a vector or a list of length = ntraits.", call. = F)
      }
      }
      if (dom & same_add_dom_QTN & add){
        if (!is.null(degree_of_dominance) & is.null(dominance_effect)) {
          dominance_effect <- lapply(additive_effect, function(x) x *
                                       degree_of_dominance)
        }
        if (architecture == "partially") {
          if (is.null(pleitropic_d) & is.null(trait_specific_d_QTN_number)) {
            trait_specific_d_QTN_number <- trait_specific_a_QTN_number
            pleitropic_d <- pleitropic_a
          }
        } else if (is.null(dominance_QTN_number)){
          dominance_QTN_number <- additive_QTN_number
        }
      }
      if (dom) {if (is.null(dominance_effect)) {
        stop("Please set \'same_add_dom_QTN\'=TRUE and a \'degree_of_dominance\' between -2 and 2, or provide \'dominance_effect\'.", call. = F)
      }
      }
      if (sum(c(!is.null(genotypes_object),
                !is.null(genotypes_file),
                !is.null(genotypes_path))) != 1) {
        stop("Please provide one of `genotypes_object`, `genotypes_file` or `genotypes_path`.", call. = F)
      }
      if (!is.null(model)){
        if (!(grepl("A", model) | grepl("D", model) | grepl("E", model))) {
          stop("Please assign a \'model\'. Options:\'A\', \'D\', \'E\' or combinations such as \'ADE\'.", call. = F)
        }
      }else{
        stop("Please assign a \'model\'. Options:\'A\', \'D\', \'E\' or combinations such as \'ADE\'.", call. = F)
      }
      if (sim_method != "geometric" & sim_method != "custom") {
        stop("Parameter \'sim_method\' should be either \'geometric\' or \'custom\'!", call. = F)
      }
      if (sim_method == "geometric") {
        if (add) {
          temp_add <- additive_effect
          if (architecture == "partially"){
            a_qtns <- (trait_specific_a_QTN_number + pleitropic_a)
          } else {
            a_qtns <-  rep(additive_QTN_number, ntraits)
          }
          additive_effect <- vector("list", ntraits)
          if (!is.null(big_additive_QTN_effect)) {
            a_qtns <- a_qtns-1
            for(i in 1:ntraits){
              additive_effect[[i]] <- c(big_additive_QTN_effect[i],
                                        rep(temp_add[[i]], a_qtns[i]) ^
                                          (1:a_qtns[i]))
            }
          }else {
            for(i in 1:ntraits){
              additive_effect[[i]] <- 
                rep(temp_add[[i]], a_qtns[i]) ^
                (1:a_qtns[i])
            }
          }
        }
        if (dom) {
          if (architecture == "partially"){ 
            d_qtns <- (trait_specific_d_QTN_number + pleitropic_d) 
          } else {
            d_qtns <-  rep(dominance_QTN_number, ntraits)
          }
          temp_dom <- dominance_effect
          dominance_effect <- vector("list", ntraits)
          for(i in 1:ntraits){
            dominance_effect[[i]] <- 
              rep(temp_dom[[i]], d_qtns[i]) ^
              (1:d_qtns[i])
          }
        }
        if (epi) {
          if (architecture == "partially"){
            e_qtns <- (trait_specific_e_QTN_number + pleitropic_e) 
          } else {
            e_qtns <-  rep(epistatic_QTN_number, ntraits)
          }
          temp_epi <- epistatic_effect
          epistatic_effect <- vector("list", ntraits)
          for(i in 1:ntraits){
            epistatic_effect[[i]] <- 
              rep(temp_epi[[i]], e_qtns[i]) ^
              (1:e_qtns[i])
          }
        }
      } else {
        if (add){
          if (!is.null(big_additive_QTN_effect)) {
            for(i in 1:ntraits){
              additive_effect[[i]] <- 
                c(big_additive_QTN_effect[i], 
                  additive_effect[[i]])
            }
            if (architecture != "partially") {
              if(any(lengths(additive_effect) != additive_QTN_number)){
                stop("When simulating big effect QTNs, \'additive_effect\' must be of length = \'additive_QTN_number\'-1.", call. = F)
              } else if (any(lengths(additive_effect) !=
                             (trait_specific_a_QTN_number + pleitropic_a))){
                stop("When simulating big effect QTNs, \'additive_effect\' must be of length = (\'trait_specific_a_QTN_number\' + \'pleitropic_a\') -1.", call. = F)
              }
            }
          } else if(any(lengths(additive_effect) != additive_QTN_number)){
            stop("Please provide an \'additive_effect\' object of length = \'additive_QTN_number\'.", call. = F)
          }
        }
        if (dom){
          if(any(lengths(dominance_effect) !=  dominance_QTN_number))
            stop("Please provide a \'dominance_effect\' object of length  = \'dominance_QTN_number\'.", call. = F)
        }
        if (epi){
          if(any(lengths(epistatic_effect) !=  epistatic_QTN_number))
            stop("Please provide an \'epistatic_effect\' object of length = \'epistatic_QTN_number\'.", call. = F)
        }
      }
      if (ntraits > 1) {
        if (!is.null(correlation) & !all(ntraits == dim(correlation))) {
          stop("ntraits do not match with dim of correlation matrix!", call. = F)
        }
        if (architecture == "LD") {
          if (ncol(h2) > 2) {
            stop("For the LD architecture, \'h2\' should have two columns, one for each trait!", call. = F) 
          }
          if (ncol(h2) ==1 ) {h2 <- cbind(h2, h2)}
        }
        if (ntraits != ncol(h2)) {
          stop("Parameter \'h2\' should have length/ncol = \'ntraits\'.", call. = F)
        }
        if (architecture == "partially"){
          
          if (add) {
            if (is.null(trait_specific_a_QTN_number) | is.null(pleitropic_a)) {
              stop("For Partially Pleiotropic additive architecture, both \'pleitropic_a\' and \'trait_specific_a_QTN_number\' must be provided!", call. = F)
            }
            if (ntraits != length(trait_specific_a_QTN_number)) {
              stop("Parameter \'trait_specific_a_QTN_number\' should have length = \'ntraits\'.", call. = F) 
            }
          }
          if (dom) {
            if (same_add_dom_QTN) { 
              if (!add) {
                stop("\'same_add_dom_QTN\' is only valid if model includes \'A\'.", call. = F)   
              }
            } else if (is.null(pleitropic_d) | is.null(trait_specific_d_QTN_number)) {
              stop("Please provide \'pleitropic_d\' and \'trait_specific_d_QTN_number\'.", call. = F)  
            }
          }
          if (epi) {
            if (is.null(trait_specific_e_QTN_number) | is.null(pleitropic_e)) {
              stop("For Partially Pleiotropic epistatic architecture, both \'pleitropic_e\' and \'trait_specific_e_QTN_number\' must be provided!", call. = F)
            }
            if (ntraits != length(trait_specific_e_QTN_number)) {
              stop("Parameter \'trait_specific_e_QTN_number\' should have length = \'ntraits\'", call. = F)
            }
          }
        } else {
          if (add) {
            if (is.null(additive_QTN_number)) {
              stop("Please provide a valid \'additive_QTN_number\'.", call. = F)
            }
          }
          if (dom) {
            if (is.null(dominance_QTN_number)) {
              stop("Please provide a valid \'dominance_QTN_number\'.", call. = F)
            }
            if (same_add_dom_QTN) {
              if (!add) {
                stop("\'same_add_dom_QTN\' is only valid if \'additive_QTN_number\' is not NULL.", call. = F)  
              }
            }
          }
          if (epi) {
            if (is.null(epistatic_QTN_number)) {
              stop("Please provide a valid \'epistatic_QTN_number\'.", call. = F)
            }
          }
        }
        
        if (architecture == "partially" | architecture == "pleiotropic"){
          if (is.null(correlation)) {
            if (!is.null(additive_effect)) {
              a_var <- all(lapply(additive_effect, var) == 0)
            }
            if (!is.null(dominance_effect)) {
              d_var <- all(lapply(dominance_effect, var) == 0)
            }
            if (!is.null(epistatic_effect)) {
              e_var <- all(lapply(epistatic_effect, var) == 0)
            }
            w <- c()
            if (exists("a_var")) w <- a_var
            if (exists("d_var")) w <- c(w, d_var)
            if (exists("e_var")) w <- c(w, e_var) 
            if (all(w)) warning("If \'correlation = NULL\', allelic effect must be different to generate different correlated traits!")
          }
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
      }
      
      
      #TODO VCF and other formats.
      input_format <- "hapmap"
      setwd(home_dir)
      on.exit(setwd(home_dir), add = TRUE)
      if (!is.null(genotypes_path) | !is.null(genotypes_file) |
          (class(unlist(genotypes_object[, 12])) != "numeric" &
           class(unlist(genotypes_object[, 12])) != "integer")) {
        genotypes_object <-
          genotypes(genotypes_object = genotypes_object,
                    genotypes_path = genotypes_path,
                    genotypes_file = genotypes_file,
                    input_format = input_format,
                    nrows = nrows,
                    na_string = na_string,
                    prefix = prefix,
                    maf_cutoff = maf_cutoff,
                    SNP_impute = SNP_impute,
                    major_allele_zero = major_allele_zero)
      }
      if (is.null(seed)){ seed <- round(runif(1, 0, 10e5)) }
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
      zz <- file("Log_Sim.txt", open = "wt")
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
      if (architecture == "pleiotropic" |
          architecture == "LD" |
          ntraits == 1 ) {
        if (add) cat("\nNumber of additive QTNs:", additive_QTN_number)
        if (dom) {
          if(same_add_dom_QTN) {
            cat("\nNumber of dominance QTNs: Same QTNs used for the additive model!")
            if (!is.null(degree_of_dominance))
              cat("\nDegree of dominance:", degree_of_dominance)
          } else {
            cat("\nNumber of dominance QTNs:", dominance_QTN_number)
          }
        }
        if (epi) cat("\nNumber of epistatic QTNs:", epistatic_QTN_number)
      }
      if (architecture == "partially") {
        if (add) {
          cat("\nNumber of pleiotropic additive QTNs:", pleitropic_a)
          cat("\nNumber of trait specific additive QTNs:",
              paste(trait_specific_a_QTN_number,
                    collapse = ", ")
          )
        }
        if (dom) {
          if(same_add_dom_QTN & add) { 
            cat("\nNumber of pleiotropic and trait specific dominant QTNs: Same as for the additive model!")
          } else {
            cat("\nNumber of pleiotropic dominant QTNs:", pleitropic_d)
            cat("\nNumber of trait specific dominant QTNs:",
                paste(trait_specific_d_QTN_number,
                      collapse = ", "))
          }
        }
        if (epi) {
          cat("\nNumber of pleiotropic epistatic QTNs:", pleitropic_e)
          cat("\nNumber of trait specific epistatic QTNs:",
              paste(trait_specific_e_QTN_number,
                    collapse = ", "))
        }
      }
      cat("\nReplicating set of QTNs every simulation (vary_QTN): ", vary_QTN)
      if(ntraits == 1) {
        cat("\nPopulation Heritability:", h2)
        if (add) {
          names(additive_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nAdditive genetic effect:\n")
          print(additive_effect)
          #cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
        }
        if (dom) {
          names(dominance_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nDominance genetic effect:\n")
          print(dominance_effect)
        }
        if (epi) {
          names(epistatic_effect) <- paste0("Trait_", 1:ntraits)
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
          names(additive_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nAdditive genetic effects:\n")
          print(additive_effect)
          #cat("\nBig effect Additive QTN:", big_additive_QTN_effect)
        }
        if (dom) {
          names(dominance_effect) <- paste0("Trait_", 1:ntraits)
          if (same_add_dom_QTN & !is.null(degree_of_dominance)) {
            cat("\nDominance genetic effects (degree_of_dominance = ",degree_of_dominance,"):\n")
          } else {
            cat("\nDominance genetic effects:\n")
          }
          print(dominance_effect) 
        }
        if (epi) {
          names(epistatic_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nEpistatic genetic effects:\n")
          print(epistatic_effect)
        }
        cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
      }
      if (architecture == "LD" | out_geno == "plink" |
          out_geno == "gds" | output_format == "gemma"){
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
          SNPRelate::snpgdsGDS2BED(genofile,
                                   bed.fn = "geno",
                                   snp.id = snpset,
                                   verbose=F)
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
            same_add_dom_QTN = same_add_dom_QTN,
            additive_QTN_number = additive_QTN_number,
            dominance_QTN_number = dominance_QTN_number,
            epistatic_QTN_number = epistatic_QTN_number,
            constrains = constrains,
            rep = rep,
            rep_by = rep_by,
            export_gt = export_gt,
            add = add,
            dom = dom,
            epi = epi
          )
      }
      if (ntraits > 1 & !any(architecture != "partially")) {
        QTN <-
          QTN_partially_pleiotropic(
            genotypes = genotypes_object,
            seed = seed,
            pleitropic_a = pleitropic_a,
            pleitropic_e = pleitropic_e,
            trait_specific_a_QTN_number = trait_specific_a_QTN_number,
            trait_specific_e_QTN_number = trait_specific_e_QTN_number,
            ntraits = ntraits,
            constrains = constrains,
            rep = rep,
            rep_by = rep_by,
            export_gt = export_gt,
            same_add_dom_QTN = same_add_dom_QTN,
            add = add,
            dom = dom,
            epi = epi
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
            export_gt = export_gt,
            same_add_dom_QTN = same_add_dom_QTN,
            add = add,
            dom = dom)
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
          if(same_add_dom_QTN){
            #include pleitropic for partially and QTN number
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
      
    }, 
    silent = FALSE)
    if (class(x)== "try-error") {
      sink()
      close(zz)
    }
  }