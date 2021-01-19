#' Simulation of single/multiple traits under different models and genetic
#' architectures.
#' @export
#' @import utils
#' @import stats
#' @importFrom data.table fwrite fread
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpair snpgdsClose snpgdsCreateGeno
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @param geno_obj Marker data set loaded as an R object.
#' Currently either HapMap or numericalized files
#' (code as aa = -1, Aa = 0 and AA = 1, e.g. `data("SNP55K_maize282_maf04")`)
#' are accepted. These and other file formats (VCF, GDS, and Plink Bed/Ped
#' files) may be read from file with `geno_file` or `geno_path`. Only one of
#' `geno_obj`, `geno_file` or `geno_path` should be provided.
#' @param geno_file Name of a marker data set to be read from a file. If in a 
#' different folder, the whole path should be provided. Formats accepted are
#' Numeric, HapMap, VCF, GDS, and Plink Bed/Ped files. Notice that the major
#' allele will always be 1 and the minor allele -1. Thus, when using Plink Bed
#' files, the dosage information will be converted to the opposite value.
#' @param geno_path Path to a folder containing the marker data set
#' file/files (e.g., separated by chromosome). Formats accepted are:
#' Numeric, HapMap, VCF, GDS, and Plink Bed/Ped files
#' @param QTN_list A list of specific markers to be used as QTNs. If one wants to specify the QTNs instead of selecting them randomly, at least one of the following elements should be provided: `QTN_list$add`, `QTN_list$dom`, and/or `QTN_list$epi`. The element `$add`, `$dom`, and `$epi` are lists containing a vector of markers for each of the traits to be simulated. For example, to simulate 2 traits controlled by 1 pleiotropic and 2 trait-specific additive QTNs, the user would create a list of marker names `marker_list <- list(add = list(trait1 = c("marker1", "marker2", "marker3"), trait2 = c("marker1", "marker4", "marker5")))` and set `QTN_list = marker_list`. On the other hand, to simulate a single trait controlled by 1 additive and 2 dominance QTNs, the marker list would be `marker_list <- list(add = list("marker1"), dom = list(c("marker2", "marker3")))`. Notice that these vectors with maker names is used in the order they appear. For instance, in the list `marker_list <- list(add = list(trait9 = c("marker1"), trait4 = c("marker5")))`, the vector names itself ("trait9" and "trait4") are ignored and "trait9" will be the vector of markers used to simulate the first trait and "trait4" will be the vector of markers used to simulate the second trait. Also, when using `QTN_list`, many parameters used for selecting QTNs will be ignored (e.g., `constraints`).
#' @param prefix If `geno_path` points to a folder with files other than the
#' marker data set, a part of the data set name may be used to select the desired
#' files (e.g., prefix = "Chr" would read files Chr1.hmp.txt, ..., Chr10.hmp.txt
#' but not HapMap.hmp.txt).
#' @param rep The number of experiments (replicates of a trait with the same
#' genetic architecture) to be simulated.
#' @param ntraits The number of multi-trait phenotypes to simulate under
#' pleiotropic, partially [pleiotropic], and LD (spurious pleiotropy)
#' architectures (see `architecture`). If not assigned, a single trait will be
#' simulated. Currently, the only option for the LD architecture is
#' `ntraits = 2`.
#' @param h2 The heritability for each traits being simulated.
#' It could be either a vector with length equals to `ntraits`,
#' or a matrix with ncol equals to `ntraits`. If the later is used, the simulation
#' will loop over the number of rows and will generate a result for each row.
#' If a single trait is being simulated and h2 is a vector,
#' one simulation of each heritability value will be conducted. Either none or
#' all traits are expected to have `h2 = 0`.
#' @param mean A vector with the mean (intercept) value for each of the simulated traits. If omitted, the simulated traits will be centered to zero. 
#' @param model The genetic model to be assumed. The options are
#' "A" (additive), "D" (dominance), "E" (epistatic)
#' as well as any combination of those models such as "AE", "DE" or "ADE".
#' @param architecture The genetic architecture to be simulated. Should be provided
#' if `ntraits` > 1. Possible options are: 'pleiotropic' (default), for traits being
#' controlled by the same QTNs; 'partially', for traits being controlled by
#' pleiotropic and trait-specific QTNs; 'LD', for traits being exclusively
#' controlled by different QTNs in "direct" or "indirect" (See `type_of_ld`, `ld_min`, and
#' `ld_max` below) linkage disequilibrium. Currently the
#' only option for `architecture = "LD"` is `ntraits = 2`.
#' @param add_QTN_num The number of additive quantitative trait nucleotides
#' (QTNs) to be simulated.
#' @param dom_QTN_num The number of dominance QTNs to be simulated.
#' @param epi_QTN_num The number of epistatic (Currently, only additive x
#' additive epistasis are simulated) QTNs to be simulated.
#' @param epi_type to be implemented...
#' @param epi_interaction Number of markers that compose an epistatic QTN. 
#' If `epi_interaction = 2` (default), a 2-way interaction (marker1 x marker2) will be 
#' used to simulate epistatic QTNs. If `epi_interaction = 3` a 3-way interaction (marker1 x marker2 x marker3) will be used instead.
#' @param pleio_a The number of pleiotropic additive QTNs to be
#' used if `architecture = "partially"`. When `sim_method = custom` (see below),
#' the first effects will be assigned to the pleiotropic QTNs and the last to
#' the trait-specific ones. For instance, in a scenario where ntraits = 2,
#' pleio_a = 2, trait_spec_a_QTN_num = 1, and add_effect = list(
#' trait1 = c(0.1, 0.2, 0.3), trait2 = c(0.4, 0.5, 0.6)), the trait-specific
#' QTNs for trait 1 and trait 2 will be 0.3 and 0.6, respectively. The first
#' two allelic effects will be assigned to the pleiotropic QTNs.
#' @param pleio_d The number of pleiotropic dominance QTNs to be
#' used if `architecture = "partially"` (See pleio_a for details).
#' @param pleio_e The number of pleiotropic epistatic QTNs to be
#' used if `architecture = "partially"` (See pleio_a for details).
#' @param trait_spec_a_QTN_num The number of trait-specific additive QTNs if
#'`architecture = "partially"`. It should be a vector of length equals
#' to `ntraits`.
#' @param trait_spec_d_QTN_num The number of trait-specific dominance QTNs if
#' `architecture = "partially"`. It should be a vector of length equals
#' to `ntraits`.
#' @param trait_spec_e_QTN_num The number of trait-specific epistatic QTNs if
#' `architecture = "partially"`. It should be a vector of length equals
#' to `ntraits`.
#' @param add_effect Additive effect size to be simulated. It may be either
#' a vector (assuming `ntraits` = 1 or one allelic effect per trait to create a
#' geometric series [`sim_method = "geometric"`]) or a list
#' of length = `ntraits`, i.e., if `ntraits` > 1, a list with one vector of
#' additive effects should be provided for each trait. Unless
#' `big_add_QTN_effect` is provided, the length of each vector
#' should be equal to the number of additive QTNs being simulated.
#' @param big_add_QTN_effect Additive effect size for one possible major
#' effect quantitative trait nucleotide. If `ntraits` > 1,
#' big_add_QTN_effect should have length equals `ntraits`.
#' If `add_QTN_num` > 1, this large effect will be assigned to the fist QTN.
#' @param same_add_dom_QTN A boolean for selecting markers to be both additive
#' and dominance QTNs.
#' @param dom_effect Similar to the `add_effect`, it could be either
#' a vector or a list. Optional if `same_add_dom_QTN = TRUE`.
#' @param degree_of_dom If the same set of QTNs are being used for simulating
#' additive and dominance effects, the dominance allelic effect could be a
#' proportion of the additive allelic effect.
#' In other words, `degree_of_dom` equals to 0.5, 1, 1.5 will simulate,
#' partial dominance, complete dominance and overdominance, respectively.
#' @param epi_effect Epistatic (additive x additive) effect size to be
#' simulated. Similar to the `add_effect`, it could be either a vector or
#' a list.
#' @param type_of_ld Type of LD used to simulate spurious pleiotropy. If
#' "indirect" (default), an intermediate marker is selected from which two
#' adjacent markers (one upstream and another downstream) will be chosen based
#' on its LD with the intermediate marker to be the QTNs. Optionally,
#' in the "direct" method, one marker is selected to be a QTN for trait 1, and
#' a second marker is selected based on its LD with the first selected marker to
#' be the QTN for trait 2.
#' @param ld_min Minimum Linkage disequilibrium for selecting QTNs when
#' `architecture = LD`. The default is `ld_min = 0.2` (markers should have a minimum LD of
#' 0.2 to be used as QTNs).
#' @param ld_max Maximum Linkage disequilibrium for selecting QTNs when
#' `architecture = LD`. The default is `ld_max = 0.8` (markers should have an LD of
#' at maximum 0.8 to be used as QTNs).
#' @param ld_method Four methods can be used to calculate linkage disequilibrium values: "composite" for LD composite measure (Default), "r" for R coefficient (by EM algorithm assuming HWE, it could be negative), "dprime" for D', and "corr" for correlation coefficient (see snpgdsLDpair from package SNPRelate).
#' @param sim_method Provide the method of simulating allelic effects.
#' The options available are "geometric" and "custom". For multiple QTNs,
#' a geometric series may be simulated, i.e., if add_effect = 0.5,
#' the effect size of the first QTNs will be 0.5, the effect size of the second
#' QTN will be 0.5^2, and the effect of the n^th QTN will be 0.5^n.
#' @param vary_QTN A boolean that determines if the same set of quantitative trait
#' nucleotide (QTN) should be used to generate genetic effects for each
#' experiment (`vary_QTN = FALSE`) or  if a different set of QTNs should be
#' used for each replication (`vary_QTN = TRUE`).
#' @param cor Option to simulate traits with a predefined genetic correlation.
#' It should be a correlation matrix with a number of rows = `ntraits`.
#' Default = NULL. Notice that when opting for controlling the correlation, the
#' genetic effects are transformed using Cholesky decomposition. In this case,
#' the correlation of genetic effects for different traits will be as provided, 
#' but due to the transformation, the actual allelic effects of correlated
#' traits may be different than the input allelic effect.
#' @param cor_res Option to simulate traits with a predefined residual
#' correlation. It should be a correlation matrix with number of
#' rows = `ntraits`. If NULL, an identity matrix (independent residuals)
#' will be used.
#' @param QTN_variance Whether or not the percentage of the phenotypic variance
#' explained by each QTN (QTN variance / phenotypic variance) should be
#' exported. The default is FALSE. Notice that this is calculated prior to any
#' transformation, such as the whitening/coloring transformation used to assign
#' user-specified correlation to the genetic effect. In may not reflect the
#' actual variance explained when the data is transformed.
#' @param seed Value to be used by set.seed. If NULL (default),
#' runif(1, 0, 1000000) will be used. Notice that at each sampling step,
#' a different seed generated based on the `seed` parameter used.
#' For example, if one uses `seed = 123`, when simulating the 10th replication
#' of trait 1, the seed to be used is `round( (123 * 10 * 10) * 1)`. On the
#' other hand, for simulating the 21st replication of trait 2, the seed to be
#' used will be `round( (123 * 21 * 21) * 2)`. The master seed (unique value required to reproduce  results) is saved at the top of the log file. Unless verbose = FALSE the actual seed used in every  simulation is exported along with simulated phenotypes.
#' @param home_dir Directory where files should be saved. It may be
#' home_dir = getwd().
#' @param output_dir Name to be used to create a folder inside `home_dir` and
#' save output files.
#' @param export_gt If TRUE genotypes of selected QTNs will be saved at file.
#' If FALSE (default), only the QTN information will be saved.
#' @param output_format Four options are available for saving simulated
#' phenotypes: 'multi-file', saves each simulation in a separate file;
#' 'long' (default for multiple traits), appends each experiment (rep) to the
#' last one (by row); 'wide', saves experiments by column (default for single
#' trait) and 'gemma', saves .fam files to be used by gemma or other software
#' that uses plink bed files. (renaming .fam file with the same name of the bim
#' and bed files is necessary).
#' @param to_r Option for outputting the simulated results as an R data.frame in
#' addition to saving it to file. If TRUE, results need to be assigned to an
#' R object (see vignette).
#' @param out_geno Optionally saves the numericalized genotype either as "numeric" (see
#' vignettes for an example data), "BED" or "gds". The default is NULL.
#' @param chr_prefix If input file format is VCF and out_geno = "BED", and a prefix
#' is used in the chromosomes names, chr_prefix may be used to avoid issues in
#' converting to bed files (e.g., chr_prefix = "chr" in "chr01").
#' @param remove_QTN Whether or not a copy of the genotypic file should be saved
#' without the simulated QTNs. The default is FALSE. If `vary_QTN = TRUE`, the
#' question "Are you sure that you want to save one genotypic file/rep
#' (remove_QTN = TRUE and vary_QTN = TRUE) [type yes or no] ?" will pop up to
#' avoid saving multiple large files unintentionally
#' @param warning_file_saver Skips the interactive question and saves all files
#' when `remove_QTN = TRUE` and `vary_QTN = TRUE`.
#' @param constraints Set constraints for QTN selection. Currently, the options
#' are maf_above (the minimum value of minor allele frequency, a double between
#' 0 - 0.5), maf_below (the maximum value of minor allele frequency, a double
#' between 0 - 0.5),and hets ('include' and 'remove'). All of these options
#' are NULL by default ('list(maf_above = NULL, maf_below = NULL, hets = NULL )'
#' ). For instance, if the parameters used are 
#' `constraints = list(maf_above = 0.3, maf_below = 0.44, hets = "include")`,
#' only heterozygote markers with minor allele frequency between 0.3 and 0.44 will
#' be selected to be QTNs. The option "remove" would only select homozygote
#' markers to be QTNs.
#' @param maf_cutoff Option for filtering the data set based on minor allele
#' frequency (Not to be confounded with the constraints option which will only
#' filter possible QTNs). It may be useful when outputting the genotypic data
#' set.
#' @param nrows Option for loading only part of a data set. Used when marker
#' data is in numeric or HapMap format. Please see data.table::fread for details.
#' @param na_string Tell create_phenotypes what character represents missing
#' data (default is "NA"). Used when the input marker data is numeric or
#' HapMap.
#' @param SNP_effect Parameter used for numericalization. The options are: Add
#' (AA = 1, Aa = 0, aa = -1),  Dom (AA = -1, Aa = 0, aa = -1), Left (AA = 1,
#' Aa = -1, aa = -1), Right (AA = 1, Aa = 1, aa = -1). The default option is Add.
#' @param SNP_impute Naive imputation for HapMap numericalization. The options
#' are: Major (NA <- 1), Middle (NA <- 0), and Minor (NA <- -1).
#' @param quiet Whether or not the log file should pop up into R once the
#' simulation is done.
#' @param verbose If FALSE, suppress all prints and suppress individual seed numbers from being saved to file. The master seed (unique value required to reproduce results) is saved at the top of the log file.
#' @param RNGversion Parameter to set the random number generator. Different R versions may be selected, the default value is `3.5.1`. 
#' @return Single or multi-trait phenotypes in one of many formats.
#' Numericalized marker data set with or without the selected QTNs.
#' Diagnostic files (log, QTN information, summary of LD between QTNs,
#' proportion of phenotypic variance explained by each QTN).
#' @references Fernandes, S.B., and Lipka, A.E., 2020 simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic
#' SIMulation of Pleiotropic, Linked and Epistatic PHENOTYPES. BMC Bioinformatics 21(1):491,
#' \doi{https://doi.org/10.1186/s12859-020-03804-y} \cr
#' @author Samuel B Fernandes and Alexander E Lipka
#' Last update: Jan 19, 2021
#' @examples
#' # Simulate 50 replications of a single phenotype.
#' data("SNP55K_maize282_maf04")
#' pheno <- 
#'   create_phenotypes(
#'     geno_obj = SNP55K_maize282_maf04,
#'     add_QTN_num = 3,
#'     add_effect = 0.2,
#'     big_add_QTN_effect = 0.9,
#'     rep = 10,
#'     h2 = 0.7,
#'     model = "A",
#'     to_r = TRUE,
#'     home_dir = tempdir(),
#'     quiet = T
#'     )
#' # For more examples, please run the following:
#' # vignette("simplePHENOTYPES")
#'
create_phenotypes <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           QTN_list = list(add = list(NULL), dom = list(NULL), epi = list(NULL)),
           prefix = NULL,
           rep = NULL,
           ntraits = 1,
           h2 = NULL,
           mean = NULL,
           model = NULL,
           architecture = "pleiotropic",
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           epi_QTN_num = NULL,
           epi_type = NULL,
           epi_interaction = 2,
           pleio_a = NULL,
           pleio_d = NULL,
           pleio_e = NULL,
           trait_spec_a_QTN_num = NULL,
           trait_spec_d_QTN_num = NULL,
           trait_spec_e_QTN_num = NULL,
           add_effect = NULL,
           big_add_QTN_effect = NULL,
           same_add_dom_QTN = FALSE,
           dom_effect = NULL,
           degree_of_dom = 1,
           epi_effect = NULL,
           type_of_ld = "indirect",
           ld_min = 0.2,
           ld_max = 0.8,
           ld_method = "composite",
           sim_method = "geometric",
           vary_QTN = FALSE,
           cor = NULL,
           cor_res = NULL,
           QTN_variance = FALSE,
           seed = NULL,
           home_dir = NULL,
           output_dir = NULL,
           export_gt = FALSE,
           output_format = "long",
           to_r = FALSE,
           out_geno = NULL,
           chr_prefix = "chr",
           remove_QTN = FALSE,
           warning_file_saver = TRUE,
           constraints = list(maf_above = NULL,
                              maf_below = NULL,
                              hets = NULL),
           maf_cutoff = NULL,
           nrows = Inf,
           na_string = "NA",
           SNP_effect = "Add",
           SNP_impute = "Middle",
           quiet = FALSE,
           verbose = TRUE,
           RNGversion = '3.5.1'
           ) {
    # -------------------------------------------------------------------------
    home_exit <- getwd()
    tryCatch({
      sunk <- FALSE
      gdsfile <- NULL
      suppressWarnings(RNGversion(RNGversion))
      if (is.null(home_dir)) {
        stop("Please provide a path to output results (It may be getwd())!.",
             call. = F)
      } else {
        files_in_dir <- dir(home_dir, full.names = T)
      }
      if (!is.null(model)) {
        if (any(!toupper(unlist(strsplit(model, ""))) %in% c("A", "D", "E")) |
            nchar(model) > 3) {
          stop(
            "Please assign a \'model\'. Options:\'A\', \'D\', \'E\' or combinations such as \'ADE\'.",
            call. = F
          )
        }
      } else{
        stop(
          "Please assign a \'model\'. Options:\'A\', \'D\', \'E\' or combinations such as \'ADE\'.",
          call. = F
        )
      }
      mm <-
        ifelse(
          architecture == "pleiotropic",
          "Pleiotropic",
          ifelse(
            architecture == "partially",
            "Partially Pleiotropic",
            ifelse(architecture == "LD", "Linkage Disequilibrium", {
              stop(
                "The genetic architecture used is not valid! Please choose one of: \'pleiotropic\', \'partially\' or \'LD\' ",
                call. = F
              )
            })
          )
        )
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
      if (!add) {
        add_QTN_num <- NULL
        add_effect <- NULL
        pleio_a <- NULL
        trait_spec_a_QTN_num <- NULL
        big_add_QTN_effect <- NULL
        same_add_dom_QTN <- FALSE
      }
      if (!dom) {
        dom_QTN_num <- NULL
        dom_effect <- NULL
        pleio_d <- NULL
        trait_spec_d_QTN_num <- NULL
      }
      if (!epi) {
        epi_QTN_num <- NULL
        epi_effect <- NULL
        pleio_e <- NULL
        trait_spec_e_QTN_num <- NULL
      }
      if (architecture == "LD") {
        ntraits <- 2
        if (type_of_ld != "indirect" & type_of_ld != "direct") {
          stop(
            "Parameter \'type_of_ld\' should be either \'direct\' or \'indirect\'.",
            call. = F
          )
        }
      }
      if (is.null(out_geno)) {
        out_geno <- "none"
      }
      if (out_geno != "none" & out_geno != "numeric" & out_geno != "BED" & out_geno !=  "gds") {
        stop(
          "Parameter \'out_geno\' should be either \'numeric\', \'BED\' or \'gds\'.",
          call. = F
        )
      }
      if (vary_QTN) {
        rep_by <-  "QTN"
      } else {
        rep_by <- "experiment"
      }
      if (is.null(mean)) {
        mean <- rep(0, ntraits)
      } else {
        if (length(mean) != ntraits) {
          stop("Parameter \'mean\' should have length = \'ntraits\'.",
               call. = F)
        }
      }
      if (is.vector(h2)) {
        if (ntraits > 1) {
          h2 <- matrix(h2, nrow = 1)
          if (ntraits != ncol(h2)) {
            stop("Parameter \'h2\' should have length/ncol = \'ntraits\'.",
                 call. = F)
          }
        } else {
          h2 <- matrix(h2, ncol = 1)
        }
      }
      colnames(h2) <- paste0("Trait_", 1:ntraits)
      if (!is.null(unlist(QTN_list))) {
        if (rep_by == "QTN") {
          stop("The option for using user inputted QTNs is only valid if \'vary_QTN = FALSE\'.",
               call. = F)
        }
        if (same_add_dom_QTN & !is.null(unlist(QTN_list$dom))) {
          stop("If \'same_add_dom_QTN = TRUE\', \'QTN_list$dom\' should not be provided. Instead, the QTNs in \'QTN_list$add\' will be used.",
               call. = F)
        }
        if (add) {
          qadd <-  suppressWarnings(do.call(cbind, QTN_list$add))
          if (ntraits != ncol(qadd) & ncol(qadd) != 1) {
            stop("The length of \'QTN_list$add\' should either be 1 (if architecture == \'pleiotropic\') or \'ntraits\'",
                 call. = F)
          }
          if (architecture == "pleiotropic") {
            if (any(apply(qadd, 2, duplicated))) {
              stop("If architecture == \'pleiotropic\', \'QTN_list$add\' should have length 1 or all of its elements should have the same length. For traits controlled by different QTN numbers, please use architecture == \'partially\'.",
                   call. = F)
            } else if (ncol(qadd) >1){
              qaddu <- apply(qadd, 1, unique)
              if (is.matrix(qaddu)) {
                if (any(nrow(qaddu) > 1 )) {
                  stop("If architecture == \'pleiotropic\' and \'length(QTN_list$add)\' = \'ntraits\' all elements of \'QTN_list$add\' should be identical.",
                       call. = F)
                }
              } else {
                if (any(lengths(qaddu) > 1 )) {
                  stop("If architecture == \'pleiotropic\' and \'length(QTN_list$add)\' = \'ntraits\' all elements of \'QTN_list$add\' should be identical.",
                       call. = F)
                } else {
                  QTN_list$add <- list(QTN_list$add[[1]])
                }
              }
            }
          }
        }
        if (dom) {
          qdom <-  suppressWarnings(do.call(cbind, QTN_list$dom))
          if (ntraits != ncol(qdom) & ncol(qdom) != 1) {
            stop("The length of \'QTN_list$dom\' should either be 1 (if architecture == \'pleiotropic\') or \'ntraits\'",
                 call. = F)
          }
          if (architecture == "pleiotropic") {
            if (any(apply(qdom, 2, duplicated))) {
              stop("If architecture == \'pleiotropic\', \'QTN_list$dom\' should have length 1 or all of its elements should have the same length. For traits controlled by different QTN numbers, please use architecture == \'partially\'.",
                   call. = F)
            } else if (ncol(qdom) >1){
              qdomu <- apply(qdom, 1, unique)
              if (is.matrix(qdomu)) {
                if (any(nrow(qdomu) > 1 )) {
                  stop("If architecture == \'pleiotropic\' and \'length(QTN_list$dom)\' = \'ntraits\' all elements of \'QTN_list$dom\' should be identical.",
                       call. = F)
                }
              } else {
                if (any(lengths(qdomu) > 1 )) {
                  stop("If architecture == \'pleiotropic\' and \'length(QTN_list$dom)\' = \'ntraits\' all elements of \'QTN_list$dom\' should be identical.",
                       call. = F)
                } else {
                  QTN_list$dom <- list(QTN_list$dom[[1]])
                }
              }
            }
          }
        }
        if (epi) {
        }
      }
      if (add) {
        if (is.vector(add_effect)) {
          add_effect <- as.list(add_effect)
        } else if (!is.list(add_effect)) {
          stop("\'add_effect\' should be either a vector or a list of length = ntraits.",
               call. = F)
        }
      }
      if (epi) {
        if (is.vector(epi_effect)) {
          epi_effect <- as.list(epi_effect)
        } else if (!is.list(epi_effect)) {
          stop("\'epi_effect\' should be either a vector or a list of length = ntraits.",
               call. = F)
        }
      }
      if (dom) {
        if (same_add_dom_QTN & add) {
          if (is.null(dom_effect)) {
            dom_effect <- lapply(add_effect, function(x)
              x * degree_of_dom)
          } else {
            if (!is.list(dom_effect)) {
              dom_effect <- as.list(dom_effect)
            }
          }
          if (architecture == "partially") {
            if (is.null(pleio_d) & is.null(trait_spec_d_QTN_num)) {
              trait_spec_d_QTN_num <- trait_spec_a_QTN_num
              pleio_d <- pleio_a
            } else if (is.null(pleio_d) | is.null(trait_spec_d_QTN_num)){
              stop(
                "If \'same_add_dom_QTN = TRUE\', please either use \'pleio_d\' = \'pleio_a\' and \'trait_spec_d_QTN_num\' = \'trait_spec_a_QTN_num\' or remove \'pleio_d\' and \'trait_spec_d_QTN_num\'",
                call. = F
              )
            } else {
              if (pleio_d != pleio_a |
                  any(trait_spec_d_QTN_num != trait_spec_a_QTN_num)) {
                stop(
                  "If \'same_add_dom_QTN = TRUE\', please either use \'pleio_d\' = \'pleio_a\' and \'trait_spec_d_QTN_num\' = \'trait_spec_a_QTN_num\' or remove \'pleio_d\' and \'trait_spec_d_QTN_num\'",
                  call. = F
                )
              }
            }
          } else if (is.null(dom_QTN_num)) {
            if (add_QTN_num != dom_QTN_num) {
              stop(
                "If \'same_add_dom_QTN = TRUE\', please use \'dom_QTN_num\' = \'add_QTN_num\'!",
                call. = F
              )
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
          stop(
            "Please set \'same_add_dom_QTN\'=TRUE and a \'degree_of_dom\' between -2 and 2, or provide \'dom_effect\'.",
            call. = F
          )
        }
      }
      if (sum(c(
        !is.null(geno_obj), !is.null(geno_file), !is.null(geno_path)
      )) != 1) {
        stop("Please provide (only) one of `geno_obj`, `geno_file` or `geno_path`.",
             call. = F)
      }
      if (sim_method != "geometric" & sim_method != "custom") {
        stop("Parameter \'sim_method\' should be either \'geometric\' or \'custom\'!",
             call. = F)
      }
      if (ntraits > 1) {
        if (add) {
          if (!is.null(big_add_QTN_effect)) {
            if (length(big_add_QTN_effect) != ntraits) {
              stop("Parameter \'big_add_QTN_effect\' should be of length ntraits",
                   call. = F)
            }
          }
          if (length(add_effect) != ntraits) {
            stop("Parameter \'add_effect\' should be of length ntraits",
                 call. = F)
          }
        }
        if (dom & length(dom_effect) != ntraits) {
          stop("Parameter \'dom_effect\' should be of length ntraits",
               call. = F)
        }
        if (epi & length(epi_effect) != ntraits) {
          stop("Parameter \'epi_effect\' should be of length ntraits",
               call. = F)
        }
      }
      if (sim_method == "geometric") {
        if (add) {
          temp_add <- add_effect
          if (architecture == "partially") {
            if (length(trait_spec_a_QTN_num) != ntraits) {
              stop("Parameter \'trait_spec_a_QTN_num\' should be of length ntraits",
                   call. = F)
            }
            a_qtns <- (trait_spec_a_QTN_num + pleio_a)
            if (all(a_qtns == 0)) {
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
            if (length(trait_spec_d_QTN_num) != ntraits) {
              stop("Parameter \'trait_spec_d_QTN_num\' should be of length ntraits",
                   call. = F)
            }
            d_qtns <- (trait_spec_d_QTN_num + pleio_d)
            if (all(d_qtns == 0)) {
              d_qtns <- rep(0, ntraits)
            }
          } else {
            d_qtns <- rep(dom_QTN_num, ntraits)
          }
          temp_dom <- dom_effect
          dom_effect <- vector("list", ntraits)
          for (i in 1:ntraits) {
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
            if (length(trait_spec_e_QTN_num) != ntraits) {
              stop("Parameter \'trait_spec_e_QTN_num\' should be of length ntraits",
                   call. = F)
            }
            e_qtns <- (trait_spec_e_QTN_num + pleio_e)
            if (all(e_qtns == 0)) {
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
        if (add) {
          if (!is.null(big_add_QTN_effect)) {
            for (i in 1:ntraits) {
              add_effect[[i]] <-
                c(big_add_QTN_effect[i],
                  add_effect[[i]])
            }
            if (architecture != "partially") {
              if (any(lengths(add_effect) != add_QTN_num)) {
                stop(
                  "When simulating big effect QTNs, \'add_effect\' must be of length = \'add_QTN_num\'-1.",
                  call. = F
                )
              } else if (any(lengths(add_effect) !=
                             (trait_spec_a_QTN_num + pleio_a))) {
                stop(
                  "When simulating big effect QTNs, \'add_effect\' must be of length = (\'trait_spec_a_QTN_num\' + \'pleio_a\') -1.",
                  call. = F
                )
              }
            }
          } else if (any(lengths(add_effect) != add_QTN_num)) {
            stop(
              "Please provide an \'add_effect\' object of length = \'add_QTN_num\'.",
              call. = F
            )
          }
        }
        if (dom) {
          if (any(lengths(dom_effect) !=  dom_QTN_num))
            stop("Please provide a \'dom_effect\' object of length  = \'dom_QTN_num\'.",
                 call. = F)
        }
        if (epi) {
          if (any(lengths(epi_effect) !=  epi_QTN_num))
            stop("Please provide an \'epi_effect\' object of length = \'epi_QTN_num\'.",
                 call. = F)
        }
      }
      if (ntraits > 1) {
        if (!is.null(cor)) {
          if (!all(ntraits == dim(cor))) {
            stop("The dimension of the \'cor\' matrix must be equal to \'ntraits\'!",
                 call. = F)
          }
          if (unique(diag(cor)) != 1 | any(abs(cor) > 1)) {
            stop("\'cor\' must be a correlation matrix!", call. = F)
          }
        }
        if (!is.null(cor_res)) {
          if (!all(ntraits == dim(cor_res))) {
            stop("The dimension of the \'cor_res\' matrix must be equal to ntraits!",
                 call. = F)
          }
          if (unique(diag(cor_res)) != 1 | any(abs(cor_res) > 1)) {
            stop("\'cor_res\' must be a correlation matrix!", call. = F)
          }
        }
        if (architecture == "LD") {
          if (ncol(h2) > 2) {
            stop(
              "For the LD architecture, \'h2\' should have two columns, one for each trait!",
              call. = F
            )
          }
          if (ncol(h2) == 1) {
            h2 <- cbind(h2, h2)
          }
        }
        if (architecture == "partially") {
          if (add) {
            if (is.null(trait_spec_a_QTN_num) | is.null(pleio_a)) {
              stop(
                "For Partially Pleiotropic additive architecture, both \'pleio_a\' and \'trait_spec_a_QTN_num\' must be provided!",
                call. = F
              )
            }
            if (ntraits != length(add_effect)) {
              stop("Parameter \'add_effect\' should have length = \'ntraits\'.",
                   call. = F)
            }
          }
          if (dom) {
            if (!add & same_add_dom_QTN) {
              stop("\'same_add_dom_QTN\' is only valid if model includes \'A\'.",
                   call. = F)
            } else if (is.null(pleio_d) |
                       is.null(trait_spec_d_QTN_num)) {
              stop(
                "Please provide \'pleio_d\' and \'trait_spec_d_QTN_num\' or set \'same_add_dom_QTN = TRUE\'.",
                call. = F
              )
            }
          }
          if (epi) {
            if (is.null(trait_spec_e_QTN_num) | is.null(pleio_e)) {
              stop(
                "For Partially Pleiotropic epistatic architecture, both \'pleio_e\' and \'trait_spec_e_QTN_num\' must be provided!",
                call. = F
              )
            }
            if (ntraits != length(epi_effect)) {
              stop("Parameter \'epi_effect\' should have length = \'ntraits\'",
                   call. = F)
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
            if (exists("a_var"))
              w <- a_var
            if (exists("d_var"))
              w <- c(w, d_var)
            if (exists("e_var"))
              w <- c(w, e_var)
            if (all(w))
              warning(
                "If \'cor = NULL\', allelic effect must be different among traits to generate different correlated traits!",
                call. = F,
                immediate. = T
              )
          }
        }
      }
      setwd(home_dir)
      on.exit({
        setwd(home_exit)
        RNGversion(getRversion())
        if (sunk) {
          sink()
          close(zz)
        }
        gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
        gc()
      }, add = TRUE)
      if (is.null(geno_obj) &
          is.null(geno_file) &
          is.null(geno_path)) {
        stop("Please provide one of: \'geno_obj\', \'geno_file\' or \'geno_path\'",
             call. = F)
      }
      out_name <- NULL
      input_format <- NULL
      if (!is.null(geno_obj)) {
        out_name <- deparse(substitute(geno_obj))
        if (any(class(geno_obj) != "data.frame")) {
          geno_obj <- as.data.frame(geno_obj)
        }
        if (is.numeric(unlist(geno_obj[, 6:7]))) {
          if (any(colnames(geno_obj)[1:5] != c("snp", "allele", "chr",  "pos",  "cm"))) {
            stop(
              "If a numeric format is provided, the first 5 columns of \'geno_obj\' should have the following names:\n       c(\"snp\", \"allele\", \"chr\",  \"pos\",  \"cm\").\n       Please see data(SNP55K_maize282_maf04) for an example. ",
              call. = F
            )
          } else {
            nonnumeric <- FALSE
          }
        } else {
          nonnumeric <- TRUE
        }
      } else {
        nonnumeric <- TRUE
      }
      path_out <- NULL
      if (!is.null(geno_path) | !is.null(geno_file) | nonnumeric) {
        geno_obj <-
          genotypes(
            geno_obj = geno_obj,
            geno_path = geno_path,
            geno_file = geno_file,
            nrows = nrows,
            na_string = na_string,
            prefix = prefix,
            maf_cutoff = maf_cutoff,
            SNP_impute = SNP_impute,
            verbose = verbose,
            chr_prefix = chr_prefix
          )
        input_format <- geno_obj$input_format
        temp <- geno_obj$temp
        if (is.null(out_name))
          out_name <- geno_obj$out_name
        geno_obj <-  geno_obj$geno_obj
      } else {
        temp <- tempfile(pattern = "", fileext = ".gds")
        if (verbose)
          message(paste0("File ", "\'",deparse(substitute(geno_obj)),"\'", " loaded from memory."))
        dose <- 0
        counter <- 6
        while (all(dose != 2) & all(dose != -1)) {
          dose <- unique(geno_obj[, counter])
          counter <- counter + 1
        }
        if (all(dose != -1) | any(dose == 2)) {
          geno_obj[, -c(1:5)] <- geno_obj[, -c(1:5)] - 1
        }
        isna <- is.na(geno_obj[, -c(1:5)])
        if (any(isna)) {
          if (SNP_impute == "Middle") {
            geno_obj[, -c(1:5)][isna] <- 0
          } else
            if (SNP_impute == "Minor") {
              geno_obj[, -c(1:5)][isna] <- -1
            } else
              if (SNP_impute == "Major") {
                geno_obj[, -c(1:5)][isna] <- 1
              }
        }
      }
      if (is.null(input_format))
        input_format <- "numeric"
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
          message("Directory name provided by \'output_dir\' alredy exists! \nCreating: ",
                  tempdir)
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
      sunk <- TRUE
      if (ntraits > 1) {
        if (add & !dom & !epi) {
          cat("Simulation of a",
              mm,
              "Genetic Architecture with Additive Effects\n")
        } else if (!add & dom & !epi) {
          cat("Simulation of a",
              mm,
              "Genetic Architecture with Dominance Effects\n")
        } else if (!add & !dom & epi) {
          cat("Simulation of a",
              mm,
              "Genetic Architecture with Epistatic Effects \n")
        } else if (add & dom & !epi) {
          cat(
            "Simulation of a",
            mm,
            "Genetic Architecture with Additive and Dominance Effects \n"
          )
        } else if (add & !dom & epi) {
          cat(
            "Simulation of a",
            mm,
            "Genetic Architecture with Additive and Epistatic Effects \n"
          )
        } else if (!add & dom & epi) {
          cat(
            "Simulation of a",
            mm,
            "Genetic Architecture with Dominance and Epistatic Effects \n"
          )
        } else {
          cat(
            "Simulation of a",
            mm,
            "Genetic Architecture with Additive, Dominance and Epistatic Effects \n"
          )
        }
      } else {
        if (add & !dom & !epi) {
          cat("Simulation of a Single Trait Genetic Architecture with Additive Effects\n")
        } else if (!add & dom & !epi) {
          cat("Simulation of a Single Trait Genetic Architecture with Dominance Effects\n")
        } else if (!add & !dom & epi) {
          cat("Simulation of a Single Trait Genetic Architecture with Epistatic Effects \n")
        } else if (add & dom & !epi) {
          cat(
            "Simulation of a Single Trait Genetic Architecture with Additive and Dominance Effects \n"
          )
        } else if (add & !dom & epi) {
          cat(
            "Simulation of a Single Trait Genetic Architecture with Additive and Epistatic Effects \n"
          )
        } else if (!add & dom & epi) {
          cat(
            "Simulation of a Single Trait Genetic Architecture with Dominance and Epistatic Effects \n"
          )
        } else {
          cat(
            "Simulation of a Single Trait Genetic Architecture with Additive, Dominance and Epistatic Effects \n"
          )
        }
      }
      cat("\nSIMULATION PARAMETERS: \n\n")
      cat("\nMaster Seed:", seed, "\n\n")
      cat("Number of traits:", ntraits)
      if (architecture == "pleiotropic" |
          architecture == "LD" |
          ntraits == 1) {
        if (add)
          cat("\nNumber of additive QTNs:", add_QTN_num)
        if (dom) {
          if (add & same_add_dom_QTN) {
            cat("\nNumber of dominance QTNs: Same QTNs used for the additive model!")
            if (!is.null(degree_of_dom)) {
              cat("\nDegree of dominance:", degree_of_dom)
              if (degree_of_dom < -2 | degree_of_dom > +2) {
                cat(
                  "Note: suggested values should range between -2 and 2.
            Please see appropriate literature for more information."
                )
              }
            }
          } else {
            cat("\nNumber of dominance QTNs:", dom_QTN_num)
          }
        }
        if (epi)
          cat("\nNumber of epistatic QTNs:", epi_QTN_num)
      }
      if (architecture == "partially") {
        if (add) {
          cat("\nNumber of pleiotropic additive QTNs:", pleio_a)
          cat(
            "\nNumber of trait specific additive QTNs:",
            paste(trait_spec_a_QTN_num,
                  collapse = ", ")
          )
        }
        if (dom) {
          if (same_add_dom_QTN & add) {
            cat(
              "\nNumber of pleiotropic and trait specific dominant QTNs: Same as for the additive model!"
            )
          } else {
            cat("\nNumber of pleiotropic dominant QTNs:", pleio_d)
            cat(
              "\nNumber of trait specific dominant QTNs:",
              paste(trait_spec_d_QTN_num,
                    collapse = ", ")
            )
          }
        }
        if (epi) {
          cat("\nNumber of pleiotropic epistatic QTNs:", pleio_e)
          cat(
            "\nNumber of trait specific epistatic QTNs:",
            paste(trait_spec_e_QTN_num,
                  collapse = ", ")
          )
        }
      }
      cat("\nReplicating set of QTNs at each simulation (vary_QTN): ",
          vary_QTN)
      if (ntraits == 1) {
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
        if (add) {
          names(add_effect) <- paste0("Trait_", 1:ntraits)
          cat("\nAdditive genetic effects:\n")
          print(add_effect)
        }
        if (dom) {
          names(dom_effect) <- paste0("Trait_", 1:ntraits)
          if (same_add_dom_QTN & !is.null(degree_of_dom)) {
            cat("\nDominance genetic effects (degree_of_dom = ",
                degree_of_dom,
                "):\n")
            if (degree_of_dom < -2 | degree_of_dom > +2) {
              cat(
                "Note: suggested values should range between -2 and 2.
            Please see appropriate literature for more information."
              )
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
      if (architecture == "LD" | out_geno == "BED" |
          out_geno == "gds" |
          (output_format == "gemma" & remove_QTN == TRUE)) {
        if (!exists("temp")) {
          temp <- tempfile(pattern = "", fileext = ".gds")
          } else if (is.null(temp)) {
            temp <- tempfile(pattern = "", fileext = ".gds")
            }
        if (input_format == "hapmap" |
            input_format == "numeric") {
          dup <- duplicated(geno_obj$snp)
          if (any(dup)) {
            message("Removing ", sum(dup), " markers for being duplicated!")
            geno_obj <- geno_obj[!dup, ]
          }
          if (!is.numeric(geno_obj$chr)) {
            geno_obj$chr <- as.numeric(gsub("\\D+", "", geno_obj$chr))
          }
          al_na <- is.na(geno_obj$allele)
          if (any(al_na)) {
            stop(
              "Allele information must be provided to create GDS file.",
              call. = F
            )
          }
          SNPRelate::snpgdsCreateGeno(
            temp,
            genmat = t(geno_obj[, -c(1:5)]) + 1,
            sample.id = colnames(geno_obj)[-c(1:5)],
            snp.id = as.character(geno_obj$snp),
            snp.chromosome = geno_obj$chr,
            snp.position = geno_obj$pos,
            snp.allele = as.character(geno_obj$allele),
            snpfirstdim = FALSE
          )
          gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
        }
        gdsfile <- temp
        if ((out_geno == "BED" |
             output_format == "gemma") & remove_QTN == FALSE) {
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          snpset <-
            SNPRelate::snpgdsSelectSNP(genofile,
                                       remove.monosnp = F,
                                       verbose = F,
                                       autosome.only = F)
          try_bed <- try(
            SNPRelate::snpgdsGDS2BED(
            genofile,
            bed.fn = out_name,
            snp.id = snpset,
            verbose = F,
            snpfirstdim = F
          ), silent = TRUE)
          if (class(try_bed) == "try-error") {
            stop(
              "Conversion to Bed files failed, probably because of chromosome names. Try using \'chr_prefix\' to remove the prefix and have names as numbers.",
              call. = F
            )
          }
          gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
        }
      }
      if (output_format == "gemma") {
        fam <- data.frame(
          colnames(geno_obj)[- (1:5)],
          colnames(geno_obj)[- (1:5)],
          0,
          0,
          0,
          check.names = FALSE,
          fix.empty.names = FALSE
        )
        colnames(fam) <- paste0("V", 1:5)
      }
      if (out_geno == "numeric" & remove_QTN == FALSE) {
        data.table::fwrite(
          geno_obj,
          paste0(out_name,
                 "_numeric.txt"),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA,
          showProgress = FALSE
        )
      }
      cat("\n\nDIAGNOSTICS:\n\n")
      if (is.null(unlist(QTN_list))){
        if (ntraits == 1 | !any(architecture != "pleiotropic")) {
          QTN <-
            qtn_pleiotropic(
              genotypes = geno_obj,
              seed = seed,
              same_add_dom_QTN = same_add_dom_QTN,
              add_QTN_num = add_QTN_num,
              dom_QTN_num = dom_QTN_num,
              epi_QTN_num = epi_QTN_num,
              epi_type = epi_type,
              epi_interaction = epi_interaction,
              constraints = constraints,
              rep = rep,
              rep_by = rep_by,
              export_gt = export_gt,
              add = add,
              dom = dom,
              epi = epi,
              verbose = verbose
            )
        }
        if (ntraits > 1 & !any(architecture != "partially")) {
          QTN <-
            qtn_partially_pleiotropic(
              genotypes = geno_obj,
              seed = seed,
              pleio_a = pleio_a,
              pleio_d = pleio_d,
              pleio_e = pleio_e,
              trait_spec_a_QTN_num = trait_spec_a_QTN_num,
              trait_spec_d_QTN_num = trait_spec_d_QTN_num,
              trait_spec_e_QTN_num = trait_spec_e_QTN_num,
              ntraits = ntraits,
              constraints = constraints,
              rep = rep,
              rep_by = rep_by,
              export_gt = export_gt,
              same_add_dom_QTN = same_add_dom_QTN,
              add = add,
              dom = dom,
              epi = epi,
              verbose = verbose
            )
        }
        if (ntraits > 1 & !any(architecture != "LD")) {
          QTN <-
            qtn_linkage(
              genotypes = geno_obj,
              seed = seed,
              add_QTN_num = add_QTN_num,
              dom_QTN_num = dom_QTN_num,
              ld_min = ld_min,
              ld_max =  ld_max,
              ld_method = ld_method,
              gdsfile = gdsfile,
              constraints = constraints,
              rep = rep,
              rep_by = rep_by,
              export_gt = export_gt,
              same_add_dom_QTN = same_add_dom_QTN,
              add = add,
              dom = dom,
              type_of_ld = type_of_ld,
              verbose = verbose
            )
        }
      } else {
        QTN <- qtn_from_user(
          genotypes = geno_obj,
          QTN_list = QTN_list,
          export_gt = export_gt,
          architecture = architecture,
          same_add_dom_QTN = same_add_dom_QTN,
          add = add,
          dom = dom,
          epi = epi,
          ntraits = ntraits,
          type_of_ld = type_of_ld,
          ld_method = ld_method,
          gdsfile = gdsfile,
          ld_min = ld_min,
          ld_max = ld_max,
          verbose = verbose
        )
      }
      if (remove_QTN) {
        if (add) {
          if (dom & same_add_dom_QTN) {
            selected_add_QTN <-
              data.table::fread("Additive_and_Dominance_Selected_QTNs.txt",
                                data.table = F)
            if (architecture == "LD") {
              selected_add_QTN <- selected_add_QTN[selected_add_QTN$type != "cause_of_LD",]
            }
          } else {
            selected_add_QTN <-
              data.table::fread("Additive_Selected_QTNs.txt", data.table = F)
            if (architecture == "LD") {
              selected_add_QTN <- selected_add_QTN[selected_add_QTN$type != "cause_of_LD",]
            }
          }
        }
        if (dom & !same_add_dom_QTN) {
          selected_dom_QTN <-
            data.table::fread("Dominance_Selected_QTNs.txt", data.table = F)
          if (architecture == "LD") {
            selected_dom_QTN <-
              selected_dom_QTN[selected_dom_QTN$type != "cause_of_LD", ]
          }
        }
        if (epi) {
          selected_epi_QTN <-
            data.table::fread("Epistatic_Selected_QTNs.txt", data.table = F)
          if (architecture == "LD") {
            selected_epi_QTN <-
              selected_epi_QTN[selected_epi_QTN$type != "cause_of_LD", ]
          }
        }
        sel_a <- NULL
        sel_d <- NULL
        sel_e <- NULL
        if (rep_by == "QTN") {
          if (warning_file_saver) {
            yes_no <- "NO"
            yes_no <-
              readline(prompt = "Are you sure that you want to save one genotypic file/rep (remove_QTN = TRUE and vary_QTN = TRUE) [type yes or no] ?\n")
            if (toupper(yes_no) != "YES" &
                toupper(yes_no) != "NO") {
              yes_no <- readline(prompt = "Please answer yes or no: \n")
            }
            if (toupper(yes_no) != "YES" &
                toupper(yes_no) != "NO") {
              yes_no <- "NO"
            }
          } else {
            yes_no <- "YES"
          }
          if (toupper(yes_no) == "YES") {
            snps_to_remove <- vector("list", rep)
            if (add)
              sel_a <-
                split(selected_add_QTN$snp, selected_add_QTN$rep)
            if (dom & !same_add_dom_QTN)
              sel_d <-
                split(selected_dom_QTN$snp, selected_dom_QTN$rep)
            if (epi)
              sel_e <-
                split(selected_epi_QTN$snp, selected_epi_QTN$rep)
            if (out_geno == "BED" |
                out_geno == "gds" | output_format == "gemma") {
              genofile <- SNPRelate::snpgdsOpen(gdsfile)
              snpset <-
                gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))
              for (i in 1:rep) {
                snps_to_remove[[i]] <- unlist(c(sel_a[i], sel_d[i], sel_e[i]))
                snpset_no_QTN <-
                  setdiff(snpset, snps_to_remove[[i]])
                SNPRelate::snpgdsGDS2BED(
                  genofile,
                  bed.fn = paste0(out_name, "_noQTN_rep_", i),
                  snp.id = snpset_no_QTN,
                  verbose = F,
                  snpfirstdim = F
                )
                if (verbose)
                  cat("\nSaving genotype file ", i, "without QTNs")
              }
              gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
            } else if (out_geno == "numeric" | out_geno == "none") {
              for (i in 1:rep) {
                snps_to_remove[[i]] <- unlist(c(sel_a[i], sel_d[i], sel_e[i]))
                data.table::fwrite(
                  geno_obj[!geno_obj$snp %in% snps_to_remove[[i]], ],
                  paste0(out_name, "_noQTN_rep", i, ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA,
                  showProgress = FALSE
                )
                if (verbose)
                  cat("\nSaving numeric genotype file ",
                      i,
                      "without QTNs")
              }
            }
          }
        } else {
          snps_to_remove <- list()
          if (add)
            sel_a <-
              split(selected_add_QTN$snp, selected_add_QTN$rep)
          if (dom & !same_add_dom_QTN)
            sel_d <-
              split(selected_dom_QTN$snp, selected_dom_QTN$rep)
          if (epi)
            sel_e <-
              split(selected_epi_QTN$snp, selected_epi_QTN$rep)
          if (out_geno == "BED" |
              out_geno == "gds" | output_format == "gemma") {
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            snpset <-
              gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "snp.id"))
            snps_to_remove[[1]] <-
              unlist(c(sel_a[1], sel_d[1], sel_e[1]))
            snpset_no_QTN <-
              setdiff(snpset, snps_to_remove[[1]])
            SNPRelate::snpgdsGDS2BED(
              genofile,
              bed.fn = paste0(out_name, "_noQTN"),
              snp.id = snpset_no_QTN,
              verbose = F,
              snpfirstdim = F
            )
            if (verbose)
              cat("\nSaving genotype file without QTNs!")
            gdsfmt::showfile.gds(closeall = TRUE, verbose = F)
          } else if (out_geno == "numeric" | out_geno == "none") {
            snps_to_remove[[1]] <- unlist(c(sel_a[1], sel_d[1], sel_e[1]))
            data.table::fwrite(
              geno_obj[!geno_obj$snp %in% snps_to_remove[[1]], ],
              paste0(out_name, "_noQTN.txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA,
              showProgress = FALSE
            )
            if (verbose)
              cat("\nSaving numeric genotype file")
          }
        }
      }
      hets <- NULL
      if (dom & !is.null(dom_effect)) {
        h_num <- if (!is.null(dom_QTN_num)) {
          dom_QTN_num > 0 & any(unlist(dom_effect) > 0)
        } else if (!is.null(pleio_d) &
                   !is.null(trait_spec_d_QTN_num)) {
          all(pleio_d > 0 & trait_spec_d_QTN_num > 0)
        }
        if (h_num & any(unlist(dom_effect) > 0)) {
          if (same_add_dom_QTN & add) {
            if (rep_by == "QTN" |
                architecture == "partially" |
                architecture == "LD") {
              if (class(QTN$add_ef_trait_obj[[1]]) == "matrix") {
                hets <- lapply(QTN$add_ef_trait_obj,
                               function(x) {
                                 f <- apply(x, 2, function(b) {
                                   b == 1
                                 })
                                 hrow <- sum(apply(f, 1, sum)) > 0
                                 return(hrow)
                               })
              } else {
                hets <- lapply(QTN$add_ef_trait_obj,
                               function(x) {
                                 lapply(x, function(x2) {
                                   f <- apply(x2, 2, function(b) {
                                     b == 1
                                   })
                                   hrow <- sum(apply(f, 1, sum)) > 0
                                   return(hrow)
                                 })
                               })
              }
            } else {
              hets <- lapply(QTN$add_ef_trait_obj,
                             function(x) {
                               f <- apply(x, 2, function(b) {
                                 b == 1
                               })
                               hrow <- sum(apply(f, 1, sum)) > 0
                               return(hrow)
                             })
            }
          } else {
            if (rep_by == "QTN" |
                architecture == "partially" |
                architecture == "LD") {
              if (class(QTN$dom_ef_trait_obj[[1]]) == "matrix") {
                hets <- lapply(QTN$dom_ef_trait_obj,
                               function(x) {
                                 f <- apply(x, 2, function(b) {
                                   b == 1
                                 })
                                 hrow <- sum(apply(f, 1, sum)) > 0
                                 return(hrow)
                               })
              } else {
                hets <- lapply(QTN$dom_ef_trait_obj,
                               function(x) {
                                 lapply(x, function(x2) {
                                   f <- apply(x2, 2, function(b) {
                                     b == 1
                                   })
                                   hrow <- sum(apply(f, 1, sum)) > 0
                                   return(hrow)
                                 })
                               })
              }
            } else {
              hets <- lapply(QTN$dom_ef_trait_obj,
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
              stop(
                "All individuals are homozygote for the selected QTNs. Dominance effect will be zero! Consider using a different seed number to select new QTNs.",
                call. = F
              )
            } else {
              warning(
                "Most individuals are homozygote for the selected QTNs. Dominance effect will be zero! Consider using a different seed number to select new QTNs.",
                call. = F,
                immediate. = T
              )
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
          } else{
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
              epi_interaction = epi_interaction,
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
        if (dom) {
          if (add & same_add_dom_QTN) {
            genetic_value <-
              base_line_multi_traits(
                ntraits = ntraits,
                cor = cor,
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$add_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                epi_interaction = epi_interaction,
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
          } else{
            genetic_value <-
              base_line_multi_traits(
                ntraits = ntraits,
                cor = cor,
                add_obj = QTN$add_ef_trait_obj,
                dom_obj = QTN$dom_ef_trait_obj,
                epi_obj = QTN$epi_ef_trait_obj,
                epi_interaction = epi_interaction,
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
        } else{
          genetic_value <-
            base_line_multi_traits(
              ntraits = ntraits,
              cor = cor,
              add_obj = QTN$add_ef_trait_obj,
              epi_obj = QTN$epi_ef_trait_obj,
              epi_interaction = epi_interaction,
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
      }
      if (verbose)
        message("* Creating phenotypes")
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
        verbose = verbose,
        QTN_variance = QTN_variance,
        add  = add,
        dom = dom,
        epi = epi,
        cor_res = cor_res,
        mean = mean,
        cor = cor
      )
      if (ntraits > 1) {
        if (!is.null(cor)) {
          cat("Population Genetic Correlation \n")
          colnames(cor) <- paste0("Trait_", 1:ntraits)
          rownames(cor) <- paste0("Trait_", 1:ntraits)
          print(cor)
        }
        if (all(h2 != 1)) {
          if (rep_by == "QTN") {
            sample_cor <- matrix(0, ntraits, ntraits)
            for (v in 1:rep) {
              sample_cor <- (sample_cor + genetic_value[[v]]$sample_cor)
            }
            sample_cor <- sample_cor / rep
          } else {
            sample_cor <- genetic_value[[1]]$sample_cor
          }
          if (!is.null(sample_cor) & length(sample_cor) > 0) {
            cat("\nSample Genetic Correlation \n")
            colnames(sample_cor) <- paste0("Trait_", 1:ntraits)
            rownames(sample_cor) <- paste0("Trait_", 1:ntraits)
            print(sample_cor)
          }
          if (!is.null(cor_res)) {
            cat("\nPopulation Residual Correlation \n")
            colnames(cor_res) <- paste0("Trait_", 1:ntraits)
            rownames(cor_res) <- paste0("Trait_", 1:ntraits)
            print(cor_res)
          }
          cat("\nSample Residual Correlation \n")
          colnames(results$sample_cor) <-
            paste0("Trait_", 1:ntraits)
          rownames(results$sample_cor) <-
            paste0("Trait_", 1:ntraits)
          print(results$sample_cor)
        }
      }
      cat("\n\nResults are saved at:", path_out)
      sink()
      close(zz)
      sunk <- FALSE
      if (!quiet) {
        file.show(paste0(path_out, "/Log_Sim.txt"))
      }
      if (!is.null(gdsfile)) {
        if (out_geno != "gds" & file.exists(gdsfile)) {
          unlink(gdsfile, force = TRUE)
        } else if (file.exists(gdsfile)) {
          tempfile <- paste0(path_out,
                             "/",
                             out_name,
                             ".gds")
          if (file.exists(tempfile)) {
            while (file.exists(tempfile)) {
              tempfile <- paste0(path_out,
                                 "/",
                                 out_name, "(", j,
                                 ").gds")
              j <- j + 1
            }
            message(
              "A file named ",
              paste0(home_dir, "/", out_name, ".gds"),
              " is already present in this folder, creating ",
              paste0(out_name, "(", j, ").gds")
            )
          }
          j <- 1
          invisible(file.rename(gdsfile, tempfile))
        }
      }
      message("Simulation completed!")
      if (to_r) {
        return(results$simulated_data)
      }
    },
    error = function(cnd) {
      message(cnd)
      dir <- dir(home_dir, full.names = T)
      unlink(dir[!dir %in% files_in_dir], force = TRUE, recursive = TRUE) 
      if (!is.null(gdsfile)) {
        if (out_geno != "gds" & file.exists(gdsfile)) {
          unlink(gdsfile, force = TRUE)
        }
      }
      if (!is.null(gdsfile)) {
        if (out_geno != "gds" & file.exists(gdsfile)) {
          unlink(gdsfile, force = TRUE)
        }
      }
    },
    interrupt = function(int) {
      dir <- dir(home_dir, full.names = T)
      unlink(dir[!dir %in% files_in_dir], force = TRUE, recursive = TRUE) 
      if (!is.null(gdsfile)) {
        if (out_geno != "gds" & file.exists(gdsfile)) {
          unlink(gdsfile, force = TRUE)
        }
      }
      if (!is.null(gdsfile)) {
        if (out_geno != "gds" & file.exists(gdsfile)) {
          unlink(gdsfile, force = TRUE)
        }
      }
    })
  }
