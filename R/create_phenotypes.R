#' Simulation of single/multiple traits under different models and genetic architectures.
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom lqmm is.positive.definite make.positive.definite
#' @param geno_obj Marker dataset loaded as an R object.
#' Currently either HapMap or numericalized files 
#' (code as aa = 0, Aa = 1 and AA = 2, e.g. `data("SNP55K_maize282_maf04")`) are accepted. Only one of 
#' `geno_obj`, `geno_file` or `geno_path` should be provided.
#' @param geno_file Name of a marker data set to be read from file.
#' @param geno_path Path to a folder containing the marker dataset 
#' file/files (e.g. separated by chromosome).
#' @param rep Number of experiments to be simulated.
#' @param ntraits Number of traits to be simulated under pleiotropic,
#' partially and LD architectures (see below). If not assigned, 
#' a single trait will be simulated. Currently, the only option for the 
#' LD architecture is `ntraits = 2`.
#' @param h2 Heritability of all traits being simulated. 
#' It could be either a vector with length equals to `ntraits`, 
#' or a matrix with ncol equals to ntraits. If the later is used, the simulation
#' will loop over the number of rows and will generate a result for each row. 
#' If a single trait is being simulated and h2 is a vector, 
#' one simulation of each heritability value will be conducted. Either none or all
#' traits are expected to have a `h2 = 0`.
#' @param model The genetic model to be assumed. The options are:
#' "A" (additive), "D" (dominance), "E" (epistatic) 
#' as well as any combination of those models such as "AE", "DE" or "ADE".
#' @param add_QTN_num Number of additive quantitative trait nucleotide
#' to be simulated.
#' @param dom_QTN_num Number of dominance quantitative trait nucleotide
#' to be simulated.
#' @param epi_QTN_num Number of epistatic (Currently, only additive x 
#' additive epistasis are simulated) quantitative trait nucleotide to be simulated.
#' @param add_effect Additive effect size to be simulated. It may be either 
#' a vector (assuming `ntraits` = 1 or one allelic effect per trait) or a list 
#' of length = `ntraits`, i.e., if `ntraits` > 1, one set of additive effects 
#' should be provided for each trait. In that case, each component should be a 
#' vector of either length one, if `sim_method = "geometric"` (see below), or 
#' length equal to the number of additive QTNs being simulated. 
#' @param same_add_dom_QTN A boolean for having the same quantitative trait
#' nucleotide having additive and dominance effects.
#' @param dom_effect Similar to the `add_effect`, it could be either 
#' a vector or a list. Optional if `same_add_dom_QTN = TRUE`. 
#' @param degree_of_dom If the same set of quantitative trait nucleotide 
#' are being used for simulating additive and dominance effects, the dominance 
#' allelic effect could be a proportion of the additive allelic effect.
#' In other words, `degree_of_dom` equals to 0.5, 1, 1.5 will simulate,
#' partial dominance, complete dominance and overdominance, respectively. 
#' @param epi_effect Epistatic (additive x additive) effect size to be
#' simulated. Similar to the `add_effect`, it could be either a vector or 
#' a list.
#' @param architecture Genetic architecture to be simulated. Should be provided 
#' if `ntraits` > 1. Possible options are: 'pleiotropic', for traits being 
#' controlled by the same QTNs; 'partially', for traits being controlled by 
#' pleiotropic and trait specific QTNs; 'LD', for traits being exclusively 
#' controlled QTNs in linkage disequilibrium (controlled by parameter `ld`). 
#' Currently the only option for `architecture = "LD"` is `ntraits = 2`.
#' @param pleio_a Number of pleiotropic additive QTNs to be 
#' used if `architecture = "partially"`.
#' @param pleio_d Number of pleiotropic dominance QTNs to be 
#' used if `architecture = "partially"`.
#' @param pleio_e Number of pleiotropic epistatic QTNs to be 
#' used if `architecture = "partially"`.
#' @param trait_spec_a_QTN_num Number of trait specific additive QTNs if
#'`architecture = "partially"`. It should have length equals to `ntraits`.
#' @param trait_spec_d_QTN_num Number of trait specific dominance QTNs if
#' `architecture = "partially"`. It should have length equals to `ntraits`.
#' @param trait_spec_e_QTN_num Number of trait specific epistatic QTNs if
#' `architecture = "partially"`. It should have length equals to `ntraits`.
#' @param ld Linkage disequilibrium between selected marker two adjacent markers
#' to be used as QTN. Default is `ld = 0.5`.
#' @param sim_method Provide the method of simulating allelic effects. 
#' The options available are "geometric" and "custom". For multiple QTNs, 
#' a geometric series may be simulated, i.e. if the add_effect = 0.5, 
#' the effect size of the first QTNs is 0.2, and the effect size of the second 
#' is 0.5^2 and the effect of the n^th QTN is 0.5^n.
#' @param vary_QTN A boolean to determine if the same set of quantitative trait 
#' nucleotide (QTN) should be used to generate genetic effects for each 
#' experiment (`vary_QTN = FALSE`) or  if a different set of QTNs should be 
#' used for each experiment (`vary_QTN = TRUE`).
#' @param big_add_QTN_effect Additive effect size for one possible major
#' effect quantitative trait nucleotide. If `ntraits` > 1, 
#' big_add_QTN_effect should have length equals `ntraits`.
#' If `add_QTN_num` > 1, the fist QTN will have the large effect.
#' @param cor Option to simulate traits with a pre-defined cor. 
#' It should be a square matrix with number of rows = `ntraits`.
#' @param seed Value to be used by set.seed. If NULL (default),
#'  runif(1, 0, 1000000) will be used. Notice that at each sampling step, 
#' a different seed generated based on the `seed` parameter is used. For example, 
#' if one uses `seed = 123`, when simulating the 10th replication of trait 1, 
#' the seed to be used is `round( (123 * 10 * 10) * 1)`. On the other hand, 
#' for simulating the 21st replication of trait 2, the seed to be used will be 
#' `round( (123 * 21 * 21) * 2)`.The actual seed used in every simulation is 
#' exported along with simulated phenotypes.
#' @param export_gt If TRUE genotypes of selected QTNs will be saved at file.
#' If FALSE (default), only the QTN information will be saved.
#' @param home_dir Directory where files will be saved. It might be home_dir = getwd().
#' @param output_dir Name to be used to create a folder and save output files.
#' @param to_r Option for saving simulated results into R in addition to saving 
#' it to file. If TRUE, results need to be assigned to an R object (see vignette).
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
#' allelic frequency is implemented. Either one or both of the following options 
#' may be non-null: 'list(maf_above = NULL, maf_below = NULL)'.
#' @param prefix If `geno_path` points to a folder with files other than the
#' marker dataset, a part of the dataset name may be used to select the desired 
#' files (e.g. prefix = "Chr" would read files Chr1.hmp.txt, ..., Chr10.hmp.txt 
#' but not HapMap.hmp.txt).
#' @param maf_cutoff Optional filter for minor allele frequency 
#' (The dataset will be filtered. Not to be confounded with the constrain option
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
#' @param quiet Whether or not the log file should be opened once the simulation is done.
#' @param verbose if FALSE, suppress prints.
#' @return Numericalized marker dataset, selected QTNs, phenotypes for 'ntraits'
#'  traits, log file.
#' @references Rice, B., Lipka, A. E. (2019). Evaluation of RR-BLUP genomic selection models that incorporate peak genome-wide association study signals in maize and sorghum. Plant Genome 12, 1–14.\doi{10.3835/plantgenome2018.07.0052} \cr
#' 
#' Alexander E. Lipka, Feng Tian, Qishan Wang, Jason Peiffer, Meng Li, Peter J. Bradbury, Michael A. Gore, Edward S. Buckler, Zhiwu Zhang, GAPIT: genome association and prediction integrated tool, Bioinformatics, Volume 28, Issue 18, 15 September 2012, Pages 2397–2399, \doi{10.1093/bioinformatics/bts444}
#' @author Samuel B Fernandes and Alexander E Lipka
#' Last update: Nov 14, 2019
#' @examples
#' # Simulate 50 replications of a single phenotype.
#'\dontrun{
#' pheno <- create_phenotypes(
#'   geno_obj = SNP55K_maize282_maf04,
#'   add_QTN_num = 3,
#'   add_effect = 0.1,
#'   big_add_QTN_effect = 0.9,
#'   rep = 10,
#'   h2 = 0.7,
#'   to_r = TRUE,
#'   model = "A",
#'   home_dir = tempdir()
#' )
#'}
#' # For more examples, please run the following:
#' # vignette("simplePHENOTYPES")
#'
create_phenotypes <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           rep = NULL,
           ntraits = 1,
           h2 = NULL,
           model = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           epi_QTN_num = NULL,
           add_effect = NULL,
           same_add_dom_QTN = FALSE,
           dom_effect = NULL,
           degree_of_dom = 1,
           epi_effect = NULL,
           architecture = "pleiotropic",
           pleio_a = NULL,
           pleio_d = NULL,
           pleio_e = NULL,
           trait_spec_a_QTN_num = NULL,
           trait_spec_d_QTN_num = NULL,
           trait_spec_e_QTN_num = NULL,
           ld = 0.5,
           sim_method = "geometric",
           vary_QTN = FALSE,
           big_add_QTN_effect = NULL,
           cor = NULL,
           seed = NULL,
           export_gt = FALSE,
           home_dir = NULL,
           output_dir = NULL,
           to_r = FALSE,
           output_format = "long",
           out_geno = NULL,
           gdsfile = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL),
           prefix = NULL,
           maf_cutoff = NULL,
           nrows = Inf,
           na_string = "NA",
           SNP_effect = "Add",
           SNP_impute = "Middle",
           major_allele_zero = FALSE,
           quiet = FALSE,
           verbose = TRUE) {
    # -------------------------------------------------------------------------
    x <- try({
      packageStartupMessage("Thank you for using the simplePHENOTYPES package!")
      sunk <- FALSE
      if (is.null(home_dir)) {
        stop("Please provide a path to output results (It may be getwd())!.", call. = F)
      }
      if (grepl("A", model)) {
        add <- TRUE
      } else {
        add <- FALSE
      }
      if (grepl("D", model)) {
        dom <- TRUE
      } else {
        dom <- FALSE
      }
      if (grepl("E", model)) {
        epi <- TRUE
      } else {
        epi <- FALSE
      }
      if (architecture == "LD") {
        ntraits <- 2
      }
      if (is.null(out_geno)) {
        out_geno <- "none"
      }
      if (vary_QTN) {
        rep_by <-  "QTN"
      } else {
        rep_by <- "experiment"
      }
      if (is.vector(h2)) {
        if (ntraits > 1){
          h2 <- matrix(h2, nrow = 1)
          if (ntraits != ncol(h2)) {
            stop("Parameter \'h2\' should have length/ncol = \'ntraits\'.", call. = F)
          }
        } else {
          h2 <- matrix(h2, ncol = 1)
        }
      }
      colnames(h2) <- paste0("Trait_", 1:ntraits)
      if (add) {
        if (is.vector(add_effect)) {
          add_effect <- as.list(add_effect)
        } else if (!is.list(add_effect)) {
          stop("\'add_effect\' should be either a vector or a list of length = ntraits.", call. = F)
        }
      }
      if (epi) {
        if (is.vector(epi_effect)) {
          epi_effect <- as.list(epi_effect)
        } else if (!is.list(epi_effect)) {
          stop("\'epi_effect\' should be either a vector or a list of length = ntraits.", call. = F)
        }
      }
      if (dom) {
        if(same_add_dom_QTN & add){
          if (is.null(dom_effect)) {
            dom_effect <- lapply(add_effect, function(x) x * degree_of_dom)
          } else {
            if (!is.list(dom_effect)) {
              dom_effect <- as.list(dom_effect)
            }
          }
          if (architecture == "partially") {
            if (is.null(pleio_d) & is.null(trait_spec_d_QTN_num)) {
              trait_spec_d_QTN_num <- trait_spec_a_QTN_num
              pleio_d <- pleio_a
            } else {
              if (pleio_d != pleio_a |
                  any(trait_spec_d_QTN_num != trait_spec_a_QTN_num)) {
                stop("If \'same_add_dom_QTN = TRUE\', please use \'pleio_d\' = \'pleio_a\' and \'trait_spec_d_QTN_num\' = \'trait_spec_a_QTN_num\'!", call. = F)
              }
            }
          } else if (is.null(dom_QTN_num)) {
            if (add_QTN_num != dom_QTN_num ) {
              stop("If \'same_add_dom_QTN = TRUE\', please use \'dom_QTN_num\' = \'add_QTN_num\'!", call. = F)
            }
            dom_QTN_num <- add_QTN_num
          }
        } else {
          if (is.null(dom_effect)) {
            stop("Please provide a set of \'dom_effect\'.", call. = F)
          } else {
            if (!is.list(dom_effect)) {
              dom_effect <- as.list(dom_effect)
            }
          }
        }
      }
      if (dom) {
        if (is.null(dom_effect)) {
          stop("Please set \'same_add_dom_QTN\'=TRUE and a \'degree_of_dom\' between -2 and 2, or provide \'dom_effect\'.", call. = F)
        }
      }
      if (sum(c(!is.null(geno_obj),
                !is.null(geno_file),
                !is.null(geno_path))) != 1) {
        stop("Please provide (only) one of `geno_obj`, `geno_file` or `geno_path`.", call. = F)
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
      if (ntraits > 1) {
        if (add) {
          if (!is.null(big_add_QTN_effect)){
            if(length(big_add_QTN_effect) != ntraits) {
              stop("Parameter \'big_add_QTN_effect\' should be of length ntraits", call. = F)
            }
          }
          if (length(add_effect) != ntraits){
            stop("Parameter \'add_effect\' should be of length ntraits", call. = F)
          }
        }
        if (dom & length(dom_effect) != ntraits){
          stop("Parameter \'dom_effect\' should be of length ntraits", call. = F)
        }
        if (epi & length(epi_effect) != ntraits){
          stop("Parameter \'epi_effect\' should be of length ntraits", call. = F)
        }
      }
      if (sim_method == "geometric") {
        if (add) {
          temp_add <- add_effect
          if (architecture == "partially"){
            if (length(trait_spec_a_QTN_num) != ntraits){
              stop("Parameter \'trait_spec_a_QTN_num\' should be of length ntraits", call. = F)
            }
            a_qtns <- (trait_spec_a_QTN_num + pleio_a)
            if (all(a_qtns == 0) ) {
              a_qtns <- rep(0, ntraits)
            }
          } else {
            a_qtns <-  rep(add_QTN_num, ntraits)
          }
          add_effect <- vector("list", ntraits)
          if (!is.null(big_add_QTN_effect)) {
            a_qtns <- ifelse(a_qtns == 0, 0,  a_qtns - 1)
            for (i in 1:ntraits) {
              if (a_qtns[i] == 0) {
                add_effect[[i]] <- 0
              } else {
                add_effect[[i]] <- c(big_add_QTN_effect[i],
                                     rep(temp_add[[i]], a_qtns[i]) ^
                                       (1:a_qtns[i]))
              }
            }
          } else {
            for (i in 1:ntraits) {
              if (a_qtns[i] == 0) {
                add_effect[[i]] <- 0
              } else {
              add_effect[[i]] <-
                rep(temp_add[[i]], a_qtns[i]) ^
                (1:a_qtns[i])
              }
            }
          }
        }
        if (dom) {
          if (architecture == "partially") {
            if (length(trait_spec_d_QTN_num) != ntraits){
              stop("Parameter \'trait_spec_d_QTN_num\' should be of length ntraits", call. = F)
            }
            d_qtns <- (trait_spec_d_QTN_num + pleio_d)
            if (all(d_qtns == 0 )) {
              d_qtns <- rep(0, ntraits)
            }
          } else {
            d_qtns <- rep(dom_QTN_num, ntraits)
          }
          temp_dom <- dom_effect
          dom_effect <- vector("list", ntraits)
          for (i in 1:ntraits){
            if (d_qtns[i] == 0) {
              dom_effect[[i]] <- 0
            } else {
            dom_effect[[i]] <-
              rep(temp_dom[[i]], d_qtns[i]) ^
              (1:d_qtns[i])
            }
          }
        }
        if (epi) {
          if (architecture == "partially") {
            if (length(trait_spec_e_QTN_num) != ntraits){
              stop("Parameter \'trait_spec_e_QTN_num\' should be of length ntraits", call. = F)
            }
            e_qtns <- (trait_spec_e_QTN_num + pleio_e)
            if (all(e_qtns == 0 )) {
              e_qtns <- rep(0, ntraits)
            }
          } else {
            e_qtns <-  rep(epi_QTN_num, ntraits)
          }
          temp_epi <- epi_effect
          epi_effect <- vector("list", ntraits)
          for (i in 1:ntraits) {
            if (e_qtns[i] == 0) {
              epi_effect[[i]] <- 0
            } else {
            epi_effect[[i]] <-
              rep(temp_epi[[i]], e_qtns[i]) ^
              (1:e_qtns[i])
            }
          }
        }
      } else {
        if (add){
          if (!is.null(big_add_QTN_effect)) {
            for (i in 1:ntraits){
              add_effect[[i]] <-
                c(big_add_QTN_effect[i],
                  add_effect[[i]])
            }
            if (architecture != "partially") {
              if (any(lengths(add_effect) != add_QTN_num)) {
                stop("When simulating big effect QTNs, \'add_effect\' must be of length = \'add_QTN_num\'-1.", call. = F)
              } else if (any(lengths(add_effect) !=
                             (trait_spec_a_QTN_num + pleio_a))){
                stop("When simulating big effect QTNs, \'add_effect\' must be of length = (\'trait_spec_a_QTN_num\' + \'pleio_a\') -1.", call. = F)
              }
            }
          } else if (any(lengths(add_effect) != add_QTN_num)) {
            stop("Please provide an \'add_effect\' object of length = \'add_QTN_num\'.", call. = F)
          }
        }
        if (dom) {
          if (any(lengths(dom_effect) !=  dom_QTN_num))
            stop("Please provide a \'dom_effect\' object of length  = \'dom_QTN_num\'.", call. = F)
        }
        if (epi) {
          if (any(lengths(epi_effect) !=  epi_QTN_num))
            stop("Please provide an \'epi_effect\' object of length = \'epi_QTN_num\'.", call. = F)
        }
      }
      if (ntraits > 1) {
        if (!is.null(cor) & !all(ntraits == dim(cor))) {
          stop("ntraits do not match with dim of cor matrix!", call. = F)
        }
        if (architecture == "LD") {
          if (ncol(h2) > 2) {
            stop("For the LD architecture, \'h2\' should have two columns, one for each trait!", call. = F)
          }
          if (ncol(h2) == 1 ) {
            h2 <- cbind(h2, h2)
          }
        }
        if (architecture == "partially"){
          if (add) {
            if (is.null(trait_spec_a_QTN_num) | is.null(pleio_a)) {
              stop("For Partially Pleiotropic additive architecture, both \'pleio_a\' and \'trait_spec_a_QTN_num\' must be provided!", call. = F)
            }
            if (ntraits != length(add_effect)) {
              stop("Parameter \'add_effect\' should have length = \'ntraits\'.", call. = F)
            }
          }
          if (dom) {
            if (!add & same_add_dom_QTN) {
              stop("\'same_add_dom_QTN\' is only valid if model includes \'A\'.", call. = F)
            } else if (is.null(pleio_d) | is.null(trait_spec_d_QTN_num)) {
              stop("Please provide \'pleio_d\' and \'trait_spec_d_QTN_num\' or set \'same_add_dom_QTN = TRUE\'.", call. = F)
            }
          }
          if (epi) {
            if (is.null(trait_spec_e_QTN_num) | is.null(pleio_e)) {
              stop("For Partially Pleiotropic epistatic architecture, both \'pleio_e\' and \'trait_spec_e_QTN_num\' must be provided!", call. = F)
            }
            if (ntraits != length(epi_effect)) {
              stop("Parameter \'epi_effect\' should have length = \'ntraits\'", call. = F)
            }
          }
        } else {
          if (add) {
            if (is.null(add_QTN_num)) {
              stop("Please provide a valid \'add_QTN_num\'.", call. = F)
            }
          }
          if (dom) {
            if (is.null(dom_QTN_num)) {
              stop("Please provide a valid \'dom_QTN_num\'.", call. = F)
            }
          }
          if (epi) {
            if (is.null(epi_QTN_num)) {
              stop("Please provide a valid \'epi_QTN_num\'.", call. = F)
            }
          }
        }
        if (architecture == "pleiotropic") {
          if (is.null(cor)) {
            if (add) {
              if (!is.null(add_effect) &
                add_QTN_num > 0) {
                af <- matrix(unlist(add_effect), add_QTN_num, ntraits)
                a_var <- all(apply(af, 1, var) == 0)
              }
            }
           if (dom) {
             if (!is.null(dom_effect) &
                dom_QTN_num > 0) {
               df <- matrix(unlist(dom_effect), dom_QTN_num, ntraits)
               d_var <- all(apply(df, 1, var) == 0)
            }
          }
          if (epi) {
            if (!is.null(epi_effect) &
                epi_QTN_num > 0) {
              ef <- matrix(unlist(epi_effect), epi_QTN_num, ntraits)
              e_var <- all(apply(ef, 1, var) == 0)
            }
          }
            w <- c()
            if (exists("a_var")) w <- a_var
            if (exists("d_var")) w <- c(w, d_var)
            if (exists("e_var")) w <- c(w, e_var)
            if (all(w)) warning("If \'cor = NULL\', allelic effect must be different to generate different correlated traits!",call. = F, immediate. = T)
          }
        }
        if (architecture == "partially") {
          if (is.null(cor)) {
            if (add) {
              if (!is.null(add_effect) &
                  all(pleio_a > 1)) {
                af <- matrix(NA, pleio_a, ntraits)
                for (i in 1:ntraits) {
                  af[,i] <- add_effect[[i]][1:pleio_a]
                }
                a_var <- all(apply(af, 2, var) == 0)
              }
            }
            if (dom) {
              if (!is.null(dom_effect) &
                  all(pleio_d > 1)) {
                df <- matrix(NA, pleio_d, ntraits)
                for (i in 1:ntraits) {
                  df[,i] <- dom_effect[[i]][1:pleio_d]
                }
                d_var <- all(apply(df, 2, var) == 0)
              }
            }
            if (epi) {
              if (!is.null(epi_effect) &
                  all(pleio_e > 1)) {
                ef <- matrix(NA, pleio_e, ntraits)
                for (i in 1:ntraits) {
                  ef[,i] <- epi_effect[[i]][1:pleio_e]
                }
                e_var <- all(apply(ef, 2, var) == 0)
              }
            }
            w <- c()
            if (exists("a_var")) w <- a_var
            if (exists("d_var")) w <- c(w, d_var)
            if (exists("e_var")) w <- c(w, e_var)
            if (all(w)) warning("If \'cor = NULL\', allelic effects must be different to generate different correlated traits!",call. = F, immediate. = T)
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
      input_format <- "hapmap"
      setwd(home_dir)
      on.exit(setwd(home_dir), add = TRUE)
      if (!is.null(geno_path) | !is.null(geno_file) |
          (class(unlist(geno_obj[, 12])) != "numeric" &
           class(unlist(geno_obj[, 12])) != "integer")) {
        geno_obj <-
          genotypes(geno_obj = geno_obj,
                    geno_path = geno_path,
                    geno_file = geno_file,
                    input_format = input_format,
                    nrows = nrows,
                    na_string = na_string,
                    prefix = prefix,
                    maf_cutoff = maf_cutoff,
                    SNP_impute = SNP_impute,
                    major_allele_zero = major_allele_zero,
                    verbose = verbose)
      }
      if (is.null(seed)) {
        seed <- as.integer(runif(1, 0, 1000000))
      }
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
      sunk <- TRUE
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
          ntraits == 1) {
        if (add) cat("\nNumber of additive QTNs:", add_QTN_num)
        if (dom) {
          if (add & same_add_dom_QTN) {
            cat("\nNumber of dominance QTNs: Same QTNs used for the additive model!")
            if (!is.null(degree_of_dom)){
              cat("\nDegree of dominance:", degree_of_dom)
            if (degree_of_dom < -2 | degree_of_dom > + 2){
              cat("Note: suggested values should range between -2 and 2.
            Please see appropriate literature for more information.")
            }
            }
          } else {
            cat("\nNumber of dominance QTNs:", dom_QTN_num)
          }
        }
        if (epi) cat("\nNumber of epistatic QTNs:", epi_QTN_num)
      }
      if (architecture == "partially") {
        if (add) {
          cat("\nNumber of pleiotropic additive QTNs:", pleio_a)
          cat("\nNumber of trait specific additive QTNs:",
              paste(trait_spec_a_QTN_num,
                    collapse = ", ")
          )
        }
        if (dom) {
          if (same_add_dom_QTN & add) {
            cat("\nNumber of pleiotropic and trait specific dominant QTNs: Same as for the additive model!")
          } else {
            cat("\nNumber of pleiotropic dominant QTNs:", pleio_d)
            cat("\nNumber of trait specific dominant QTNs:",
                paste(trait_spec_d_QTN_num,
                      collapse = ", ")
            )
          }
        }
        if (epi) {
          cat("\nNumber of pleiotropic epistatic QTNs:", pleio_e)
          cat("\nNumber of trait specific epistatic QTNs:",
              paste(trait_spec_e_QTN_num,
                    collapse = ", "))
        }
      }
      cat("\nReplicating set of QTNs every simulation (vary_QTN): ", vary_QTN)
      if (ntraits == 1) {
        cat("\nPopulation Heritability:", h2)
        if (add) {
          names(add_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nAdditive genetic effect:\n")
          print(add_effect)
        }
        if (dom) {
          names(dom_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nDominance genetic effect:\n")
          print(dom_effect)
        }
        if (epi) {
          names(epi_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nEpistatic genetic effect:\n")
          print(epi_effect)
        }
        cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
      } else {
        if (nrow(h2) > 1){
          cat("\nPopulation Heritability:")
          print(h2)
        } else {
          cat("\nPopulation Heritability:\n")
          print(h2)
        }
        if (add) {
          names(add_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nAdditive genetic effects:\n")
          print(add_effect)
        }
        if (dom) {
          names(dom_effect) <- paste0("Trait_", 1:ntraits)
          if (same_add_dom_QTN & !is.null(degree_of_dom)) {
            cat("\nDominance genetic effects (degree_of_dom = ", degree_of_dom, "):\n")
            if (degree_of_dom < -2 | degree_of_dom > + 2){
              cat("Note: suggested values should range between -2 and 2.
            Please see appropriate literature for more information.")
            }
          } else {
            cat("\nDominance genetic effects:\n")
          }
          print(dom_effect)
        }
        if (epi) {
          names(epi_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nEpistatic genetic effects:\n")
          print(epi_effect)
        }
        cat(paste0("\nOutput file format: \'", output_format, "\'\n"))
      }
      if (architecture == "LD" | out_geno == "plink" |
          out_geno == "gds" | output_format == "gemma") {
        dup <- duplicated(geno_obj$snp)
        if (any(dup)) {
          message("Removing ", sum(dup), " markers for being duplicated!")
          geno_obj <- geno_obj[!dup, ]
        }
        if (any(grepl("geno.gds", dir(home_dir)))) {
          message("A file named geno.gds is already present in this folder, creating \'geno2\'.")
          gdsfile <- "geno2"
        }
        if (is.null(gdsfile))  gdsfile <- "geno"
        if (!is.numeric(geno_obj$chr)) {
          geno_obj$chr <- as.numeric(
            gsub("\\D+", "", geno_obj$chr)
          )
        }
        al_na <- is.na(geno_obj$allele)
        if (any(al_na)) {
          stop("Allele information must be provided to create GDS file.",
                  call. = F, immediate. = T)
        }
        SNPRelate::snpgdsCreateGeno(
          paste0(home_dir, "/", gdsfile, ".gds"),
          genmat = t(geno_obj[, -c(1:5)]),
          sample.id = colnames(geno_obj)[-c(1:5)],
          snp.id = as.character(geno_obj$snp),
          snp.chromosome = geno_obj$chr,
          snp.position = geno_obj$pos,
          snp.allele = as.character(geno_obj$allele),
          snpfirstdim = FALSE
        )
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
        gdsfile <- paste0(home_dir, "/", gdsfile, ".gds")
        if (out_geno == "gds") cat("GDS files saved at:", home_dir, "\n")
        if (out_geno == "plink" |
            output_format == "gemma") {
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          snpset <-
            SNPRelate::snpgdsSelectSNP(genofile, remove.monosnp = F, verbose = F)
          SNPRelate::snpgdsGDS2BED(genofile,
                                   bed.fn = "geno",
                                   snp.id = snpset,
                                   verbose = F, 
                                   snpfirstdim = F)
          gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
          cat("\nPlink bed files saved at:", home_dir, "\n")
        }
      }
      if (output_format == "gemma") {
        fam <- data.table::fread("geno.fam", drop = 6)
        fam$V1 <- fam$V2
      }
      if (out_geno == "numeric") {
        data.table::fwrite(
          geno_obj,
          paste0(
            ifelse(is.null(output_dir), "", "../"),
            "Numeric_SNP_File.txt"
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
            genotypes = geno_obj,
            seed = seed,
            same_add_dom_QTN = same_add_dom_QTN,
            add_QTN_num = add_QTN_num,
            dom_QTN_num = dom_QTN_num,
            epi_QTN_num = epi_QTN_num,
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
            genotypes = geno_obj,
            seed = seed,
            pleio_a = pleio_a,
            pleio_d = pleio_d,
            pleio_e = pleio_e,
            trait_spec_a_QTN_num = trait_spec_a_QTN_num,
            trait_spec_d_QTN_num = trait_spec_d_QTN_num,
            trait_spec_e_QTN_num = trait_spec_e_QTN_num,
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
            genotypes = geno_obj,
            seed = seed,
            add_QTN_num = add_QTN_num,
            dom_QTN_num = dom_QTN_num,
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
      hets <- NULL
      if (dom & !is.null(dom_effect)) {
        h_num <- if(!is.null(dom_QTN_num)) {
          dom_QTN_num > 0 & any(unlist(dom_effect) > 0)
        } else if (!is.null(pleio_d) &
                   !is.null(trait_spec_d_QTN_num)) {
          all(pleio_d > 0 & trait_spec_d_QTN_num > 0)
        }
        if (h_num & any(unlist(dom_effect) > 0) ) {
          if (same_add_dom_QTN & add){
            if (rep_by == "QTN" |
                architecture == "partially" |
                architecture == "LD") {
              if ( class(QTN$add_ef_trait_obj[[1]]) == "matrix") {
                hets <- lapply(
                  QTN$add_ef_trait_obj,
                  function (x) {
                    f <- apply(x, 2, function (b) {
                      b == 1
                    })
                    hrow <- sum(apply(f, 1, sum)) > 0
                    return(hrow)
                  })
              } else {
                hets <- lapply(
                  QTN$add_ef_trait_obj,
                  function (x) {
                    lapply(x, function (x2) {
                      f <- apply(x2, 2, function (b) {
                        b == 1
                      })
                      hrow <- sum(apply(f, 1, sum)) > 0
                      return(hrow)
                    }) })
              }
            } else {
              hets <- lapply(
                QTN$add_ef_trait_obj,
                function(x) {
                  f <- apply(x, 2, function(b) {
                    b == 1
                  })
                  hrow <- sum(apply(f, 1, sum)) > 0 
                  #TODO count number of het QTNs
                  #hcol <- apply(f, 2, sum)
                  #sum(hcol > 0)
                  return(hrow)
                })
            }
          } else {
            if (rep_by == "QTN"|
                architecture == "partially"|
                architecture == "LD") {
              if ( class(QTN$dom_ef_trait_obj[[1]]) == "matrix") {
                hets <- lapply(
                  QTN$dom_ef_trait_obj,
                  function (x) {
                    f <- apply(x, 2, function (b) {
                      b == 1
                    })
                    hrow <- sum(apply(f, 1, sum)) > 0
                    return(hrow)
                  })
              } else {
                hets <- lapply(
                  QTN$dom_ef_trait_obj,
                  function (x) {
                    lapply(x, function (x2) {
                      f <- apply(x2, 2, function (b) {
                        b == 1
                      })
                      hrow <- sum(apply(f, 1, sum)) > 0
                      return(hrow)
                    }) })
              }
            } else {
              hets <- lapply(
                QTN$dom_ef_trait_obj,
                function(x) {
                  f <- apply(x, 2, function(b) {
                    b == 1
                  })
                  hrow <- sum(apply(f, 1, sum)) > 0 
                  return(hrow)
                })
            }
          }
          if (any(!unlist(hets))) {
            if (!add & !epi) {
              stop("All individuals are homozygote for the selected QTNs. Dominance effect will be zero! Consider using a different seed number to select new QTNs.", call. = F)
            } else {
              warning("Most individuals are homozygote for the selected QTNs. Dominance effect will be zero! Consider using a different seed number to select new QTNs.",
                      call. = F, immediate. = T)
            }
          }
        }
      }
      if (ntraits == 1) {
        if (dom) {
          if (add & same_add_dom_QTN) {
            genetic_value <-
              base_line_single_trait(
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$add_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                add_effect = add_effect,
                dom_effect = dom_effect,
                epi_effect = epi_effect,
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
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$dom_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                add_effect = add_effect,
                dom_effect = dom_effect,
                epi_effect = epi_effect,
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
              add_obj = QTN$add_ef_trait_obj,
              epi_obj = QTN$epi_ef_trait_obj,
              add_effect = add_effect,
              epi_effect = epi_effect,
              rep = rep,
              rep_by = rep_by,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method
            )
        }
      } else {
        if (dom){
          if (add & same_add_dom_QTN) {
            genetic_value <-
              base_line_multi_traits(
                ntraits = ntraits,
                cor = cor,
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$add_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                add_effect = add_effect,
                dom_effect = dom_effect,
                epi_effect = epi_effect,
                rep = rep,
                rep_by = rep_by,
                architecture = architecture,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method,
                verbose = verbose
              )
          }else{
            genetic_value <-
              base_line_multi_traits(
                ntraits = ntraits,
                cor = cor,
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$dom_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                add_effect = add_effect,
                dom_effect = dom_effect,
                epi_effect = epi_effect,
                rep = rep,
                rep_by = rep_by,
                architecture = architecture,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method,
                verbose = verbose
              )
          }
        }else{
          genetic_value <-
            base_line_multi_traits(
              ntraits = ntraits,
              cor = cor,
              add_obj = QTN$add_ef_trait_obj,
              epi_obj = QTN$epi_ef_trait_obj,
              add_effect = add_effect,
              epi_effect = epi_effect,
              rep = rep,
              rep_by = rep_by,
              architecture = architecture,
              add = add,
              dom = dom,
              epi = epi,
              sim_method = sim_method,
              verbose = verbose
            )
        }
        if (!is.null(cor)) {
          cat("Population cor \n")
          colnames(cor) <- paste0("Trait_", 1:ntraits)
          rownames(cor) <- paste0("Trait_", 1:ntraits)
          print(cor)
        }
        if (rep_by == "QTN") {
          sample_cor <- matrix(0, ntraits, ntraits)
          for (v in 1:rep) {
            sample_cor <- (sample_cor + genetic_value[[v]]$sample_cor)
          }
          sample_cor <- sample_cor / rep
        } else {
          sample_cor <- genetic_value[[1]]$sample_cor
        }
        if (!is.null(sample_cor) & length(sample_cor) > 0){
        cat("\nSample cor \n")
        colnames(sample_cor) <- paste0("Trait_", 1:ntraits)
        rownames(sample_cor) <- paste0("Trait_", 1:ntraits)
        print(sample_cor)
      }
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
        rep_by = rep_by,
        hets = unlist(hets),
        verbose = verbose
      )
      cat("\n\nResults are saved at:", home_dir)
      sink()
      close(zz)
      if (!quiet) {file.show(paste0(path_out, "/Log_Sim.txt"))}
      if (out_geno != "gds") {
        unlink(gdsfile, force = TRUE)
      }
      if (to_r) return(results)
    },
    silent = FALSE)
    if (class(x) == "try-error" & sunk) {
      sink()
      close(zz)
    }
  }