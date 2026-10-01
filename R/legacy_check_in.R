#' Generate environmental effects based on a given heritability
#' @keywords internal
#' @param geno_obj = NULL,
#' @param geno_file = NULL,
#' @param geno_path = NULL,
#' @param QTN_list = list(add = list(NULL), dom = list(NULL), epi = list(NULL), var = list(NULL)),
#' @param prefix = NULL,
#' @param rep = NULL,
#' @param ntraits = 1,
#' @param h2 = NULL,
#' @param mean = NULL,
#' @param model = NULL,
#' @param architecture = "pleiotropic",
#' @param add_QTN_num = NULL,
#' @param dom_QTN_num = NULL,
#' @param epi_QTN_num = NULL,
#' @param var_QTN_num = NULL,
#' @param epi_type = NULL,
#' @param epi_interaction = 2,
#' @param pleio_a = NULL,
#' @param pleio_d = NULL,
#' @param pleio_e = NULL,
#' @param trait_spec_a_QTN_num = NULL,
#' @param trait_spec_d_QTN_num = NULL,
#' @param trait_spec_e_QTN_num = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param var_effect = NULL,
#' @param remove_add_effect = FALSE,
#' @param same_add_dom_QTN = FALSE,
#' @param same_mv_QTN = FALSE,
#' @param big_add_QTN_effect = NULL,
#' @param degree_of_dom = 1,
#' @param type_of_ld = "indirect",
#' @param ld_min = 0.2,
#' @param ld_max = 0.8,
#' @param ld_method = "composite",
#' @param sim_method = "geometric",
#' @param vary_QTN = FALSE,
#' @param cor = NULL,
#' @param cor_res = NULL,
#' @param QTN_variance = FALSE,
#' @param seed = NULL,
#' @param home_dir = NULL,
#' @param output_dir = NULL,
#' @param export_gt = FALSE,
#' @param output_format = "long",
#' @param to_r = FALSE,
#' @param out_geno = NULL,
#' @param chr_prefix = "chr",
#' @param remove_QTN = FALSE,
#' @param warning_file_saver = TRUE,
#' @param constraints = list(maf_above = NULL,  maf_below = NULL, hets = NULL),
#' @param maf_cutoff = NULL,
#' @param nrows = Inf,
#' @param na_string = "NA",
#' @param SNP_effect = "Add",
#' @param SNP_impute = "Middle",
#' @param quiet = FALSE,
#' @param verbose = TRUE,
#' @param RNGversion = '3.5.1'
#' @param .rng internal environment created by `create_phenotypes()` that
#'   records the caller's RNG kind and state (see `.v1_rng_capture()`).
#' @return Phenotypes for ntraits traits
#' @author Samuel Fernandes and Alexander Lipka. Last update: Apr 20, 2020
#'
check_in <-
  function(geno_obj = NULL,
           geno_file = NULL,
           geno_path = NULL,
           QTN_list = list(
             add = list(NULL),
             dom = list(NULL),
             epi = list(NULL),
             var = list(NULL)
           ),
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
           var_QTN_num = NULL,
           epi_type = NULL,
           epi_interaction = 2,
           pleio_a = NULL,
           pleio_d = NULL,
           pleio_e = NULL,
           trait_spec_a_QTN_num = NULL,
           trait_spec_d_QTN_num = NULL,
           trait_spec_e_QTN_num = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           var_effect = NULL,
           remove_add_effect = FALSE,
           same_add_dom_QTN = FALSE,
           same_mv_QTN = FALSE,
           big_add_QTN_effect = NULL,
           degree_of_dom = 1,
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
           RNGversion = '3.5.1',
           .rng = NULL) {
    #--- basic argument validation (audit v1-core, D1: reject bad inputs) ----
    .v1_validate_basic(
      rep = rep, ntraits = ntraits, model = model, architecture = architecture,
      seed = seed, RNGversion = RNGversion, output_format = output_format,
      to_r = to_r, SNP_effect = SNP_effect, SNP_impute = SNP_impute
    )
    if (is.null(seed)) {
      # The RNG version must be in force before anything is drawn, and the
      # draw is the only use of the caller's RNG stream: create_phenotypes()
      # restores the caller's RNG kind on exit and leaves the caller's stream
      # advanced by exactly this one draw (see .v1_rng_capture()).
      suppressWarnings(RNGversion(RNGversion))
      seed <- as.integer(runif(1, 0, 1000000))
      if (!is.null(.rng)) .rng$advanced <- .Random.seed_safe()
    }
    #--- home_dir ----
    if (warning_file_saver & remove_QTN & vary_QTN) {
      yes_no <- "NO"
      yes_no <-
        toupper(readline(prompt = "Are you sure that you want to save one genotypic file per replicate (remove_QTN = TRUE and vary_QTN = TRUE) [type yes or no] ?\n"))
      if (yes_no == "Y") yes_no <- "YES"
      if (yes_no == "N") yes_no <- "NO"
      if (yes_no != "YES" &
          yes_no != "NO") {
        yes_no <- toupper(readline(prompt = "Please answer yes or no: \n"))
        if (yes_no == "Y") yes_no <- "YES"
        if (yes_no == "N") yes_no <- "NO"
      }
      if (yes_no != "YES" &
          yes_no != "NO") {
        yes_no <- "NO"
        warning(
          "Setting remove_QTN = FALSE!",
          call. = F,
          immediate. = T
        )
      }
    } else {
      yes_no <- "YES"
    }
    if (is.null(home_dir)) {
      stop("Please provide a path to output results (It may be getwd())!.",
           call. = F)
    } else if (!dir.exists(home_dir)) {
      stop("Directory provided in \'home_dir\' does not exist!.",
           call. = F)
    }
    # A run writes several files (phenotypes, QTNs, genetic values, log). They
    # always go in a folder of their own so a directory is never littered with
    # loose output; `output_dir = ""` restores the pre-2.0 behavior of writing
    # straight into `home_dir`.
    if (is.null(output_dir)) {
      output_dir <- "simplePHENOTYPES_output"
    }
    if (nzchar(output_dir)) {
      tempdir <- paste0(home_dir, "/", output_dir)
      if (dir.exists(tempdir)) {
        j <- 1
        while (dir.exists(tempdir)) {
          tempdir <- paste0(home_dir, "/", output_dir, "(", j, ")")
          j <- j + 1
        }
        message("Directory name provided by \'output_dir\' alredy exists! \nCreating: ",
                tempdir)
      }
      }

    #---- model -----
    if (!is.null(model)) {
      if (any(!toupper(unlist(strsplit(model, ""))) %in% c("A", "D", "E", "V")) |
          nchar(model) > 4) {
        stop(
          "Please assign a \'model\'. Options:\'A\', \'D\', \'E\', \'V\'  or combinations such as \'ADE\'.",
          call. = F
        )
      }
    } else{
      stop(
        "Please assign a \'model\'. Options:\'A\', \'D\', \'E\', \'V\' or combinations such as \'ADE\'.",
        call. = F
      )
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
    if (grepl("V", model)) {
      var <- TRUE
      if (!add) {
        add <- TRUE
        warning(
          "Simulation of traits controlled by variance QTNs only has not being implemented. Including additive QTNs!",
          call. = F,
          immediate. = T
        )
      }
    } else {
      var <- FALSE
    }
    #----- QTN num ------
    if (!add) {
      add_QTN_num <- NULL
      add_effect <- NULL
      pleio_a <- NULL
      trait_spec_a_QTN_num <- NULL
      big_add_QTN_effect <- NULL
      same_add_dom_QTN <- FALSE
      same_mv_QTN <- FALSE
      QTN_list$add <- NULL
    }
    if (!dom) {
      dom_QTN_num <- NULL
      dom_effect <- NULL
      pleio_d <- NULL
      trait_spec_d_QTN_num <- NULL
      same_add_dom_QTN <- FALSE
      QTN_list$dom <- NULL
    } else if (same_add_dom_QTN) {
      if (!is.null(unlist(QTN_list))) {
        stop(
          "`same_add_dom_QTN = TRUE` cannot be combined with `QTN_list` (the dominance markers would silently replace the ones you listed). ",
          "To use the same markers for additive and dominance effects, set `same_add_dom_QTN = FALSE` and list the same markers in both `QTN_list$add` and `QTN_list$dom`.",
          call. = F
        )
      }
      dom_QTN_num <- add_QTN_num
      dom_effect <- lapply(add_effect, function(x) x * degree_of_dom)
      pleio_d <- pleio_a
      trait_spec_d_QTN_num <- trait_spec_a_QTN_num
      QTN_list$dom <- QTN_list$add
    }
    if (!epi) {
      epi_QTN_num <- NULL
      epi_effect <- NULL
      pleio_e <- NULL
      trait_spec_e_QTN_num <- NULL
      QTN_list$epi <- NULL
    }
    if (!var) {
      var_QTN_num <- NULL
      var_effect <- NULL
      same_mv_QTN <- FALSE
      QTN_list$var <- NULL
    } else if (same_mv_QTN) {
      if (!is.null(unlist(QTN_list))) {
        stop(
          "`same_mv_QTN = TRUE` cannot be combined with `QTN_list`. List the variance markers in `QTN_list$var` instead.",
          call. = F
        )
      }
      var_QTN_num <- add_QTN_num
      var_effect <- add_effect
      QTN_list$var <- QTN_list$add
    }
    #---- architecture -----
    mm <- ifelse(
      !is.null(unlist(QTN_list)),
      "User-Defined",
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
    )
    if (architecture == "LD") {
      if (!is.null(unlist(QTN_list))) {
        stop(
          "`QTN_list` cannot be combined with `architecture = \"LD\"`: user-specified markers are simulated as a plain user-defined architecture and `ld_min`, `ld_max` and `type_of_ld` would be silently ignored. ",
          "Remove `QTN_list` to let create_phenotypes() select markers in LD, or use `architecture = \"pleiotropic\"` with `QTN_list`.",
          call. = F
        )
      }
      if (ntraits > 2) {
        stop(
          "`architecture = \"LD\"` simulates exactly two traits (ntraits = 2); `ntraits = ", ntraits, "` is not supported.",
          call. = F
        )
      }
      ntraits <- 2
      if (type_of_ld != "indirect" & type_of_ld != "direct") {
        stop("Parameter \'type_of_ld\' should be either \'direct\' or \'indirect\'.",
             call. = F)
      }
      .v1_validate_ld(ld_min, ld_max, ld_method)
      if (type_of_ld == "indirect" && dom && !add) {
        stop(
          "`model = \"D\"` (dominance without additive effects) is not supported with `architecture = \"LD\"` and `type_of_ld = \"indirect\"`. ",
          "Use `type_of_ld = \"direct\"`, or include additive effects (e.g. `model = \"AD\"`).",
          call. = F
        )
      }
    }
    if (ntraits == 1 && dom && epi) {
      stop(
        "A single trait (`ntraits = 1`) cannot be simulated with both dominance and epistatic effects (models \"DE\" / \"ADE\") by create_phenotypes(): this combination is not supported. ",
        "Use `model = \"AE\"` or `model = \"AD\"` for one trait, simulate two or more traits (`ntraits >= 2`), or use simulate_phenotype() with additive(), dominance() and epistasis().",
        call. = F
      )
    }
    if (ntraits == 1 && dom && same_add_dom_QTN) {
      stop(
        "`same_add_dom_QTN = TRUE` is not supported for a single trait (`ntraits = 1`) in create_phenotypes(). ",
        "Simulate two or more traits (`ntraits >= 2`), select additive and dominance QTNs separately (`same_add_dom_QTN = FALSE`), or use simulate_phenotype() with additive() and dominance() on the same loci.",
        call. = F
      )
    }
    if (var && ntraits > 1) {
      stop(
        "Variance QTL (a model containing \"V\") is only implemented for a single trait (`ntraits = 1`).",
        call. = F
      )
    }
    
    #---- genotype ----
    if (is.null(out_geno)) {
      out_geno <- "none"
    }
    if (out_geno != "none" & out_geno != "numeric" & out_geno != "BED" & out_geno !=  "gds") {
      stop("Parameter \'out_geno\' should be either \'numeric\', \'BED\' or \'gds\'.",
           call. = F)
    }
    if (sum(c(!is.null(geno_obj),!is.null(geno_file),!is.null(geno_path))) != 1) {
      stop("Please provide (only) one of `geno_obj`, `geno_file` or `geno_path`.",
           call. = F)
    }
    if (!is.null(geno_obj)) {
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
    if (!nonnumeric && SNP_effect != "Add") {
      stop(
        "`SNP_effect = \"", SNP_effect, "\"` has no effect on a numeric `geno_obj` (already coded aa = -1, Aa = 0, AA = 1); ",
        "it is only used to numericalize HapMap/VCF/PLINK/GDS input. Recode the numeric genotypes yourself or use the default \"Add\".",
        call. = F
      )
    }
    path_out <- NULL
    
    #---- vary QTN ----
    if (vary_QTN) {
      rep_by <-  "QTN"
    } else {
      rep_by <- "experiment"
    }
    #---- mean ----
    if (is.null(mean)) {
      mean <- rep(0, ntraits)
    } else if (!is.numeric(mean) || anyNA(mean) || any(!is.finite(mean))) {
      stop("Parameter \'mean\' should be a finite numeric vector with one value per trait.",
           call. = F)
    } else if (length(mean) != ntraits) {
      stop("Parameter \'mean\' should have length = \'ntraits\'.",
           call. = F)
    }
    #---- h2 ----
    if (is.null(h2)) {
      stop("Please provide the heritability \'h2\' (a number between 0 and 1 for each trait).",
           call. = F)
    }
    h2 <- as.matrix(h2)
    if (!is.numeric(h2) || anyNA(h2) || any(!is.finite(h2)) ||
        any(h2 < 0) || any(h2 > 1)) {
      stop("Parameter \'h2\' should contain finite heritabilities between 0 and 1 (h2 = 0 simulates traits without genetic effects).",
           call. = F)
    }
    if (ntraits > 1) {
      if (sum(dim(h2)) == 2) {
        h2 <- rep(h2, ntraits)
        h2 <- matrix(h2, nrow = 1)
      } else {
        if (any(dim(h2) == 1)) {
          h2 <- matrix(h2, nrow = 1)
        }
        if (ntraits != ncol(h2)) {
          stop(
            "Parameter \'h2\' should either be a vector of length 1 or a matrix of ncol = \'ntraits\'.",
            call. = F
          )
        }
      }
    } else {
      if (ncol(h2) != 1) {
        stop(
          "When ntraits = 1, the parameter \'h2\' should either be a vector of length 1 or a matrix of ncol = 1.",
          call. = F
        )
      }
    }
    colnames(h2) <- paste0("Trait_", 1:ntraits)
    h2_0 <- apply(h2, 1, function(x) {
      tx <- table(x)
      (length(tx) > 1 & any(names(tx) == "0"))
    })
    if (any(h2_0)){
      warning(
        "Setting all h2 to zero for traits. Either one of none of the traits should have h2=0.",
        call. = F,
        immediate. = T
      )
      h2[h2_0,] <- 0
    }
    null_setting <- FALSE
    if (length(unique(c(h2))) == 1) {
      if (unique(c(h2)) == 0) {
        null_setting <- TRUE
      }
    }
    if (var && any(h2 == 0)) {
      stop("Variance QTL (a model containing \"V\") requires h2 > 0 for the simulated trait.",
           call. = F)
    }
    .v1_validate_seed_arith(seed = seed, rep = rep, h2 = h2,
                            null_setting = null_setting,
                            wide = (dom || epi || var),
                            ld = identical(architecture, "LD"),
                            n_qtn = c(add_QTN_num, dom_QTN_num, epi_QTN_num,
                                      var_QTN_num))
    if (to_r && nrow(h2) > 1 && (ntraits > 1 || vary_QTN)) {
      stop(
        "`to_r = TRUE` returns the simulated data of every row of `h2` only when a single trait is simulated with `vary_QTN = FALSE`; with several rows in `h2`, ",
        if (ntraits > 1) "`ntraits > 1`" else "`vary_QTN = TRUE`",
        " would return only the last row. Use `to_r = FALSE` and read the output files, or call create_phenotypes() once per row of `h2`.",
        call. = F
      )
    }
    #----- QTN_list and QTN number -----
    if (!is.null(unlist(QTN_list))) {
      if (ntraits == 1) {
        stop(
          "`QTN_list` cannot be used with `ntraits = 1`: a single trait cannot be simulated from a user-specified marker list by create_phenotypes() (this combination is not supported). ",
          "Either select the QTNs at random (`add_QTN_num`, `dom_QTN_num`, ...), or simulate two or more traits and set `ntraits` to the number of marker vectors in each element of `QTN_list`, ",
          "or use simulate_phenotype() with `qtn =` in additive()/dominance()/epistasis().",
          call. = F
        )
      }
      if(is.null(names(QTN_list))){
        if (length(QTN_list) == 4){
          names(QTN_list) <- c("add", "dom", "epi", "var")
          warning(
            "Each list inside QTN_list should be named (one of: \"add\", \"dom\", \"epi\", \"var\"). The following order will be assumed: add = QTN_list[[1]], dom = QTN_list[[2]], epi = QTN_list[[3]], var = QTN_list[[4]]",
            call. = F,
            immediate. = T
          )
        } else if (length(QTN_list) == 3) {
          names(QTN_list) <- c("add", "dom", "epi")
          warning(
            "Each list inside QTN_list should be named (one of: \"add\", \"dom\", \"epi\", \"var\"). The following order will be assumed: add = QTN_list[[1]], dom = QTN_list[[2]], epi = QTN_list[[3]]",
            call. = F,
            immediate. = T
          )
        } else if (length(QTN_list) == 2) {
          names(QTN_list) <- c("add", "dom")
          warning(
            "Each list inside QTN_list should be named (one of: \"add\", \"dom\", \"epi\", \"var\"). The following order will be assumed: add = QTN_list[[1]], dom = QTN_list[[2]]",
            call. = F,
            immediate. = T
          )
        } else if (length(QTN_list) == 1) {
          names(QTN_list) <- c("add")
          warning(
            "Each list inside QTN_list should be named (one of: \"add\", \"dom\", \"epi\", \"var\"). The following order will be assumed: add = QTN_list[[1]]",
            call. = F,
            immediate. = T
          )
        } else {
          stop(
            "QTN_list should have maximum length = 4. E.g., QTN_list = list(add = list(marker_name), dom = list(marker_name), epi = list(marker_name), var = list(marker_name))!",
            call. = F
          )
        }
      }
      if ((
        length(QTN_list$add) +
        length(QTN_list$dom) +
        length(QTN_list$epi) +
        length(QTN_list$var)
      ) / sum(add + dom + epi + var) != ntraits &  architecture != "pleiotropic") {
        stop(
          "QTN_list should contain one list of markers for each trait (i.e., if ntraits = 2, QTN_list$add should be composed of 2 lists)",
          call. = F
        )
      } else if ((
        length(QTN_list$add) +
        length(QTN_list$dom) +
        length(QTN_list$epi) +
        length(QTN_list$var)
      ) / sum(add + dom + epi + var) == 1 & ntraits != 1) {
        if (add) {
          QTN_list$add[1:ntraits] <- QTN_list$add
          }
        if (dom) {
          QTN_list$dom[1:ntraits] <- QTN_list$dom
          }
        if (epi) {
          QTN_list$epi[1:ntraits] <- QTN_list$epi
          }
        if (var) {
          QTN_list$var[1:ntraits] <- QTN_list$var
          }
      }
     architecture <- "User-Defined"
      if (vary_QTN) {
        stop(
          "The option for using user inputted QTNs is only valid if \'vary_QTN = FALSE\'.",
          call. = F
        )
      }
      # `same_add_dom_QTN` / `same_mv_QTN` together with `QTN_list` are rejected
      # above, so no QTN_list$dom / QTN_list$var back-fill is needed here.
     if (is.null(QTN_list$add)) {
       add <- FALSE
     }
     if (is.null(QTN_list$dom) & !same_add_dom_QTN) {
       dom <- FALSE
     }
     if (is.null(QTN_list$epi)) {
       epi <- FALSE
     }
     if (is.null(QTN_list$var) & !same_mv_QTN) {
       var <- FALSE
     }
      if (add) {
        #TODO
        # order based on trait name
        # if (!is.null(names(QTN_list$add))) {
        #   QTN_list$add <- 
        #     QTN_list$add[order(as.numeric(gsub("[[:alpha:]]","",names(QTN_list$add))))]
        # }

        
        if (length(QTN_list$add) != ntraits) {
          stop(
            "`ntraits` (", ntraits, ") does not match the number of trait-specific marker vectors in `QTN_list$add` (", length(QTN_list$add),
            "). Set `ntraits = ", length(QTN_list$add), "` (and provide `h2`/`mean` for that many traits) or provide one marker vector per trait in `QTN_list$add`.",
            call. = F
          )
        }
        dupa <- unlist(lapply(lapply(QTN_list$add, duplicated), any))
        if (any(dupa)) {
          stop(paste0("QTN_list$add contain duplicated Markers for trait", which(dupa), ". Please remove it."),
               call. = F)
        }
        add_QTN_num <- NULL
        pleio_a <- NULL
        trait_spec_a_QTN_num <- NULL
      }
      if (dom) {
        if (length(QTN_list$dom) != ntraits) {
          stop(
            "`ntraits` (", ntraits, ") does not match the number of trait-specific marker vectors in `QTN_list$dom` (", length(QTN_list$dom),
            "). Set `ntraits = ", length(QTN_list$dom), "` (and provide `h2`/`mean` for that many traits) or provide one marker vector per trait in `QTN_list$dom`.",
            call. = F
          )
        }
        dupd <- unlist(lapply(lapply(QTN_list$dom, duplicated), any))
        if (any(dupd)) {
          stop(paste0("QTN_list$dom contain duplicated Markers for trait", which(dupd), ". Please remove it."),
               call. = F)
        }
        dom_QTN_num <- NULL
        pleio_d <- NULL
        trait_spec_d_QTN_num <- NULL
      }
      if (epi) {
        if (any(lengths(QTN_list$epi) %% epi_interaction != 0)) {
          stop(paste("epi_interaction =", epi_interaction, "Please provide",epi_interaction, "Markers should be provided for each epistatic QTN."),
               call. = F)
        }
        
        if (length(QTN_list$epi) != ntraits) {
          stop(
            "`ntraits` (", ntraits, ") does not match the number of trait-specific marker vectors in `QTN_list$epi` (", length(QTN_list$epi),
            "). Set `ntraits = ", length(QTN_list$epi), "` (and provide `h2`/`mean` for that many traits) or provide one marker vector per trait in `QTN_list$epi`.",
            call. = F
          )
        }
        dupe <- unlist(lapply(lapply(QTN_list$epi, duplicated), any))
        if (any(dupe)) {
          stop(paste0("QTN_list$epi contain duplicated Markers for trait", which(dupe), ". Please remove it."),
               call. = F)
        }
        epi_QTN_num <- NULL
        pleio_e <- NULL
        trait_spec_e_QTN_num <- NULL
      }
      if (var) {
        if (length(QTN_list$var) != 1) {
          stop(
            "Currently, variance QTL is only implemented for single trait simulations!",
            call. = F
          )
        }
        dupv <- unlist(lapply(lapply(QTN_list$var, duplicated), any))
        if (any(dupv)) {
          stop(paste0("QTN_list$var contain duplicated Markers for trait", which(dupv), ". Please remove it."),
               call. = F)
        }
        var_QTN_num <- NULL
      }
      len_a <- lengths(QTN_list$add)
      len_d <- lengths(QTN_list$dom)
      len_e <- lengths(QTN_list$epi)/epi_interaction
      len_v <- lengths(QTN_list$var)
    } else {
      if (architecture == "partially") {
        len_a <- (trait_spec_a_QTN_num + pleio_a)
        len_d <- (trait_spec_d_QTN_num + pleio_d)
        len_e <- (trait_spec_e_QTN_num + pleio_e)
        if (add & length(len_a) != ntraits) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set values for additive QTN number (\'trait_spec_a_QTN_num\' and \'pleio_a\')",
            call. = F
          )
        }
        if (dom & length(len_d) != ntraits) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set values for dominance QTN number (\'trait_spec_d_QTN_num\' and \'pleio_d\')",
            call. = F
          )
        }
        if (epi & length(len_e) != ntraits) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set values for epistatic QTN number (\'trait_spec_e_QTN_num\' and \'pleio_e\')",
            call. = F
          )
        }
      } else {
        len_a <- add_QTN_num
        len_d <- dom_QTN_num
        len_e <- epi_QTN_num
        if (add & length(len_a) != 1) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set one value for the number of additive QTNs (\'add_QTN_num\')",
            call. = F
          )
        }
        if (dom & length(len_d) != 1) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set one value for the number of dominance QTNs (\'dom_QTN_num\')",
            call. = F
          ) 
        }
        if (epi & length(len_e) != 1) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set one value for the number of epistatic QTNs (\'epi_QTN_num\')",
            call. = F
          )
        }
        len_a <- rep(len_a, ntraits)
        len_d <- rep(len_d, ntraits)
        len_e <- rep(len_e, ntraits)
      }
      if (var) {
        len_v <- var_QTN_num
        if (length(var_QTN_num) == 0) {
          stop(
            "Please provide a list of SNPs to be used as QTNs (\'QTN_list\') or set the number of variance QTNs (\'var_QTN_num\')",
            call. = F
          )
        }
      }
    }
    #---- correlation matrices and output format ----
    if (!is.null(cor)) {
      if (ntraits == 1) {
        warning("`cor` is ignored when a single trait is simulated (ntraits = 1).",
                call. = F, immediate. = T)
      } else {
        .v1_validate_cor(cor, ntraits, "cor", positive_definite = TRUE)
      }
    }
    if (!is.null(cor_res)) {
      if (ntraits == 1) {
        warning("`cor_res` is ignored when a single trait is simulated (ntraits = 1).",
                call. = F, immediate. = T)
      } else {
        .v1_validate_cor(cor_res, ntraits, "cor_res", positive_definite = FALSE)
      }
    }
    if (output_format == "wide" && ntraits > 1 && rep < 2) {
      stop(
        "`output_format = \"wide\"` needs at least two replicates (`rep >= 2`) when ntraits > 1; use `output_format = \"long\"` (or \"multi-file\") for a single replicate.",
        call. = F
      )
    }
    #---- allelic effects ----
    if (!is.null(big_add_QTN_effect)) {
      if (length(big_add_QTN_effect) != ntraits) {
        stop("Parameter \'big_add_QTN_effect\' should be a vector of length ntraits",
             call. = F)
      }
    }
    if (add) {
      if (is.null(add_effect)) {
        stop(
          "Please provide either a vector or a list of additive allelic effects \'add_effect\'.",
          call. = F
        )
      } else if (is.vector(add_effect)) {
        if (ntraits > 1) {
          add_effect <- as.list(add_effect)
        } else {
          add_effect <- list(add_effect)
        }
      } else if (!is.list(add_effect)) {
        stop("\'add_effect\' should be either a vector or a list of length = ntraits.",
             call. = F)
      }
    }
    if (dom) {
      if (is.null(dom_effect)) {
        stop(
          "Please set \'same_add_dom_QTN\'=TRUE and \'degree_of_dom\' to a value between -2 and 2, or provide either a vector or a list of dominance allelic effects \'dom_effect\'. ",
          call. = F
        )
      } else if (is.vector(dom_effect)) {
        if (ntraits > 1) {
          dom_effect <- as.list(dom_effect)
        } else {
          dom_effect <- list(dom_effect)
        }
      } else if (!is.list(dom_effect)) {
        stop("\'dom_effect\' should be either a vector or a list of length = ntraits.",
             call. = F)
      }
    }
    if (epi) {
      if (is.null(epi_effect)) {
        stop(
          "Please provide either a vector or a list of epistatic allelic effects \'epi_effect\'.",
          call. = F
        )
      } else if (is.vector(epi_effect)) {
        if (ntraits > 1) {
          epi_effect <- as.list(epi_effect)
        } else {
          epi_effect <- list(epi_effect)
        }
      } else if (!is.list(epi_effect)) {
        stop("\'epi_effect\' should be either a vector or a list of length = ntraits.",
             call. = F)
      }
    }
    if (var) {
      if (is.null(var_effect)) {
        stop(
          "Please set \'same_mv_QTN\'=TRUE, or provide either a vector or a list of variance allelic effects \'var_effect\'. ",
          call. = F
        )
      } else if (is.vector(var_effect)) {
        var_effect <- list(var_effect)
      } else if (!is.list(var_effect)) {
        stop("\'var_effect\' should be either a vector or a list of length == 1 or length == var_QTN_num.",
             call. = F)
      }
    }
    if (ntraits > 1) {
      if (add & length(add_effect) != ntraits) {
        stop("Parameter \'add_effect\' should be of length ntraits",
             call. = F)
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
    
    #----- geometric method -----    
    if (!is.character(sim_method) || length(sim_method) != 1L ||
        (sim_method != "geometric" & sim_method != "custom")) {
      stop("Parameter \'sim_method\' should be either \'geometric\' or \'custom\'!",
           call. = F)
    }
    if (sim_method == "geometric") {
      .v1_validate_effects(
        add = add, dom = dom, epi = epi, var = var,
        add_effect = add_effect, dom_effect = dom_effect,
        epi_effect = epi_effect, var_effect = var_effect,
        len_a = if (add) len_a else NULL,
        len_d = if (dom) len_d else NULL,
        len_e = if (epi) len_e else NULL,
        len_v = if (var) len_v else NULL,
        big = !is.null(big_add_QTN_effect)
      )
    }
    sm <- sim_method
      s1 = s2 = s3 = s4 <- NULL
      if (add) {
        if (!is.null(big_add_QTN_effect)) {
            if (all(lengths(add_effect) == (len_a -1))) {
              s1 <- "custom"
            } else {
              s1 <- "geometric"
            }
        } else {
          if (all(lengths(add_effect) == len_a)) {
            s1 <- "custom"
          } else {
            s1 <- "geometric"
          }
        }
      }
      if (dom) {
        if (all(lengths(dom_effect) == len_d)) {
          s2 <- "custom"
        } else {
          s2 <- "geometric"
        }
      }
      if (epi) {
        if (all(lengths(epi_effect) == len_e)) {
          s3 <- "custom"
        } else {
          s3 <- "geometric"
        }
      }
      if (var) {
        if (all(lengths(var_effect) == var_QTN_num)) {
          s4 <- "custom"
        } else {
          s4 <- "geometric"
        }
      }
      
      if (all(unique(c(s1, s2, s3, s4)) == "custom") & sm != "custom") {
        sim_method <- "custom"
        if (verbose)
          message("One effect size has been provided for each QTN. Setting sim_method = \"custom\"! ")
      }
      if (sim_method == "geometric") {
        if (add) {
          temp_add <- add_effect
          add_effect <- vector("list", ntraits)
          if (!is.null(big_add_QTN_effect)) {
            len_a <- ifelse(len_a == 0, 0,  len_a - 1)
            for (i in 1:ntraits) {
              if (len_a[i] == 0) {
                add_effect[[i]] <- 0
              } else {
                add_effect[[i]] <- c(big_add_QTN_effect[i],
                                     rep(temp_add[[i]], len_a[i]) ^
                                       (1:len_a[i]))
              }
            }
          } else {
            for (i in 1:ntraits) {
              if (len_a[i] == 0) {
                add_effect[[i]] <- 0
              } else {
                add_effect[[i]] <-
                  rep(temp_add[[i]], len_a[i]) ^
                  (1:len_a[i])
              }
            }
          }
        }
        if (dom) {
          temp_dom <- dom_effect
          dom_effect <- vector("list", ntraits)
          for (i in 1:ntraits) {
            if (len_d[i] == 0) {
              dom_effect[[i]] <- 0
            } else {
              dom_effect[[i]] <-
                rep(temp_dom[[i]], len_d[i]) ^
                (1:len_d[i])
            }
          }
        }
        if (epi) {
          temp_epi <- epi_effect
          epi_effect <- vector("list", ntraits)
          for (i in 1:ntraits) {
            if (len_e[i] == 0) {
              epi_effect[[i]] <- 0
            } else {
              epi_effect[[i]] <-
                rep(temp_epi[[i]], len_e[i]) ^
                (1:len_e[i])
            }
          }
        }
        if (var) {
          if (var_QTN_num == 0) {
            var_effect[[1]] <- 0
          } else {
            var_effect[[1]] <-
              rep(var_effect[[1]], var_QTN_num[1]) ^
              (1:var_QTN_num[1])
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
            if (any(lengths(add_effect) != len_a)) {
              stop(
                "When simulating big effect QTNs, \'add_effect\' must be of length = \'add_QTN_num\'-1 (or (\'trait_spec_a_QTN_num\' + \'pleio_a\') -1).",
                call. = F
              )
            }
            } else if (any(lengths(add_effect) != len_a)) {
            stop(
              "Please provide an \'add_effect\' object of length = \'add_QTN_num\' (or \'trait_spec_a_QTN_num\' + \'pleio_a\' if architecture = \'partially\').",
              call. = F
            )
          }
        }
        if (dom) {
          if (any(lengths(dom_effect) !=  len_d))
            stop("Please provide a \'dom_effect\' object of length  = \'dom_QTN_num\' (or \'trait_spec_d_QTN_num\' + \'pleio_d\' if architecture = \'partially\').",
                 call. = F)
        }
        if (epi) {
          if (any(lengths(epi_effect) !=  len_e))
            stop("Please provide an \'epi_effect\' object of length = \'epi_QTN_num\' (or \'trait_spec_e_QTN_num\' + \'pleio_e\' if architecture = \'partially\').",
                 call. = F)
        }
        if (var) {
          if (any(lengths(var_effect) != len_v)) {
            stop(
              "Please provide a \'var_effect\' object of length = \'var_QTN_num\'.",
              call. = F
            )
          }
        }
      }
    #----- print ------
      a1 <- NULL
      d1 <- NULL
      e1 <- NULL
      v1 <- NULL
      if (add) a1 <- "Additive"
      if (dom) d1 <- "Dominance"
      if (epi) e1 <- "Epistatic"
      if (var) v1 <- "Variance"
      adev <- c(a1, d1, e1, v1)
      if (length(adev) > 2) {
        adev <- c(paste0(adev[-length(adev)], ", ", collapse = " "), adev[length(adev)])
      }
      title <- ifelse(ntraits > 1,
                      paste("\nSimulation of a", mm,"Genetic Architecture with", 
                            paste(adev, collapse = " and "), "Effects"),
                      paste("\nSimulation of a Single Trait Genetic Architecture with", 
                            paste(adev, collapse = " and "), "Effects"))
      p1 = p2 = p3 = p4 = p5 = p6 = p7 = p8 = p9 = p10 = p11 = p12 = p13 = p14 <- NULL
      p1 <- paste(paste0(rep("*", nchar(title) - 1), collapse = ""), title,
                  paste0("\n", paste0(rep("*", nchar(title) - 1), collapse = "")))
      p2 <- paste("\n\n",paste0(rep(" ", (nchar(title)-21)/2), collapse = ""),
        "SIMULATION PARAMETERS",
        paste0(rep(" ", (nchar(title)-21)/2), collapse = ""),
        paste0("\n",paste0(rep("_", nchar(title)-1), collapse = ""),"\n"),
        "\nDate/Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
        "\nMaster Seed:", seed, "\nNumber of traits:", ntraits)
      if (!is.null(unlist(QTN_list))) {
        if (add) p3 <- paste("\nNumber of additive QTNs:", paste0(len_a, collapse = ", "))
        if (dom) p4 <- paste("\nNumber of dominance QTNs:", paste0(len_d, collapse = ", "))
        if (epi) p5 <- paste("\nNumber of epistatic QTNs:", paste0(len_e, collapse = ", "))
        if (var) p6 <- paste("\nNumber of variance QTNs:", paste0(len_v, collapse = ", "))
      } else if (architecture == "partially") {
        if (add) {
          p3 <- paste("\nNumber of pleiotropic additive QTNs:", pleio_a,
            "\nNumber of trait specific additive QTNs:",
            paste(trait_spec_a_QTN_num, collapse = ", "))
        }
        if (dom) {
          if (same_add_dom_QTN) {
            p4 <- "\nNumber of pleiotropic and trait specific dominance QTNs: Same as for the additive model (same_add_dom_QTN = TRUE)!"
          } else {
            p5 <- paste("\nNumber of pleiotropic dominance QTNs:", pleio_d,
              "\nNumber of trait specific dominance QTNs:",
              paste(trait_spec_d_QTN_num, collapse = ", ")
            )
          }
        }
        if (epi) {
          p7 <- paste("\nNumber of pleiotropic epistatic QTNs:", pleio_e,
            "\nNumber of trait specific epistatic QTNs:",
            paste(trait_spec_e_QTN_num, collapse = ", ")
          )
        }
      } else if (architecture == "pleiotropic" | architecture == "LD" |  ntraits == 1) {
        if (add)
          p3 <- paste("\nNumber of additive QTNs:", add_QTN_num)
        if (dom) {
          if (add & same_add_dom_QTN) {
            p4 <- paste("\nNumber of dominance QTNs: Same QTNs used for the additive model!")
            if (!is.null(degree_of_dom)) {
              p5 <- paste("\nDegree of dominance:", degree_of_dom)
              if (degree_of_dom < -2 | degree_of_dom > +2) {
                p6 <- "Note: suggested values should range between -2 and 2."
              }
            }
          } else {
            p4 <- paste("\nNumber of dominance QTNs:", dom_QTN_num)
          }
        }
        if (epi)
          p7 <- paste("\nNumber of epistatic QTNs:", epi_QTN_num)
        if (var) {
          if (add & same_mv_QTN) {
            p8 <- paste("\nNumber of variance QTNs: Same QTNs used for the additive model!")
          } else {
            p9 <- paste("\nNumber of variance QTNs:", var_QTN_num)
          }
        }
      }
      if (vary_QTN) p10 <- "\nReplicating set of QTNs at each simulation (vary_QTN = T)!"
      print1 <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
      
      if (same_mv_QTN) p11 <- paste("\nAdditive and Variance QTNs are the same (same_mv_QTN = T) and are saved as \'Additive_QTNs.txt\'\n")
      if (dom & same_add_dom_QTN) p12 <- paste("\nAdditive and Dominance QTNs are the same (same_add_dom_QTN = T) and are saved as \'Additive_QTNs.txt\' \n")
      p13 <- paste0("\nOutput file format: \'", output_format, "\'\n")
      
      p14 <- paste("\n",paste0(rep(" ", (nchar(title)-11)/2), collapse = ""),
                      "DIAGNOSTICS", 
                      paste0(rep(" ", (nchar(title)-11)/2), collapse = ""),
                      paste0("\n",paste0(rep("_", nchar(title)-1), collapse = ""), "\n"))
      print2 <- c(p11, p12, p13, p14)
    #---- clean parent environment ----
      args <- c("QTN_list", "ntraits", "h2", "mean", "model", "architecture", "add_QTN_num", "dom_QTN_num", "epi_QTN_num", "var_QTN_num", "pleio_a", "pleio_d", "pleio_e", "trait_spec_a_QTN_num", "trait_spec_d_QTN_num", "trait_spec_e_QTN_num", "add_effect", "dom_effect", "epi_effect", "var_effect", "same_add_dom_QTN", "same_mv_QTN", "big_add_QTN_effect", "degree_of_dom", "sim_method", "vary_QTN", "seed", "output_dir",  "out_geno")
      rm(list = args[args %in% ls(envir = parent.frame())], envir = parent.frame())
      if (!interactive()){
        quiet <- TRUE
        assign("quiet", quiet, envir = parent.frame())
      }
   #---- output variables -----
    assign("seed", seed, envir = parent.frame())
    assign("add", add, envir = parent.frame())
    assign("dom", dom, envir = parent.frame())
    assign("epi", epi, envir = parent.frame())
    assign("var", var, envir = parent.frame())
    if (dom) {
      assign("len_d", len_d, envir = parent.frame())
    }
    if (architecture == "pleiotropic" | architecture == "LD"){
      assign("add_QTN_num", add_QTN_num, envir = parent.frame())
      assign("dom_QTN_num", dom_QTN_num, envir = parent.frame())
      assign("epi_QTN_num", epi_QTN_num, envir = parent.frame())
      assign("var_QTN_num", var_QTN_num, envir = parent.frame())
    }
    if (architecture == "partially"){
      assign("pleio_a", pleio_a, envir = parent.frame())
      assign("pleio_e", pleio_e, envir = parent.frame())
      assign("pleio_d", pleio_d, envir = parent.frame())
      assign("trait_spec_d_QTN_num", trait_spec_d_QTN_num, envir = parent.frame())
      assign("trait_spec_a_QTN_num", trait_spec_a_QTN_num, envir = parent.frame())
      assign("trait_spec_e_QTN_num", trait_spec_e_QTN_num, envir = parent.frame())
    }
    assign("add_effect", add_effect, envir = parent.frame())
    assign("dom_effect", dom_effect, envir = parent.frame())
    assign("epi_effect", epi_effect, envir = parent.frame())
    assign("var_effect", var_effect, envir = parent.frame())
    assign("QTN_list", QTN_list, envir = parent.frame())
    assign("same_add_dom_QTN", same_add_dom_QTN, envir = parent.frame())
    assign("same_mv_QTN", same_mv_QTN, envir = parent.frame())
    assign("rep_by", rep_by, envir = parent.frame())
    assign("yes_no", yes_no, envir = parent.frame())
    assign("mm", mm, envir = parent.frame())
    assign("architecture", architecture, envir = parent.frame())
    assign("ntraits", ntraits, envir = parent.frame())
    assign("out_geno", out_geno, envir = parent.frame())
    assign("output_dir", output_dir, envir = parent.frame())
    assign("nonnumeric", nonnumeric, envir = parent.frame())
    assign("mean", mean, envir = parent.frame())
    assign("h2", h2, envir = parent.frame())
    assign("null_setting", null_setting, envir = parent.frame())
    assign("tempdir", tempdir, envir = parent.frame())
    assign("path_out", path_out, envir = parent.frame())
    assign("print1", print1, envir = parent.frame())
    assign("print2", print2, envir = parent.frame())
  
  }
# ---------------------------------------------------------------------------
# Argument validation helpers for the frozen v1 engine (audit v1-core, D1:
# bad inputs are rejected up front with a message that names the argument and
# the remedy; valid inputs are never altered).
# ---------------------------------------------------------------------------

#' Scalar / choice validation shared by check_in()
#' @keywords internal
#' @noRd
.v1_validate_basic <- function(rep, ntraits, model, architecture, seed,
                               RNGversion, output_format, to_r,
                               SNP_effect, SNP_impute) {
  whole1 <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == round(x)
  }
  if (is.null(rep)) {
    stop("Please provide the number of replicates `rep` (a whole number >= 1).",
         call. = FALSE)
  }
  if (!whole1(rep)) {
    stop("`rep` must be a single whole number >= 1.", call. = FALSE)
  }
  if (!whole1(ntraits)) {
    stop("`ntraits` must be a single whole number >= 1.", call. = FALSE)
  }
  if (!is.null(model) &&
      (!is.character(model) || length(model) != 1L || is.na(model) ||
       !nzchar(model) || nchar(model) > 4L ||
       any(!strsplit(model, "")[[1]] %in% c("A", "D", "E", "V")))) {
    stop("Please assign a \'model\'. Options:\'A\', \'D\', \'E\', \'V\' or combinations such as \'ADE\' (upper case).",
         call. = FALSE)
  }
  if (!is.character(architecture) || length(architecture) != 1L ||
      !architecture %in% c("pleiotropic", "partially", "LD")) {
    stop("The genetic architecture used is not valid! Please choose one of: \'pleiotropic\', \'partially\' or \'LD\' ",
         call. = FALSE)
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)
  }
  if (!is.character(RNGversion) || length(RNGversion) != 1L ||
      is.na(RNGversion)) {
    stop("`RNGversion` must be a single character string such as \'3.5.1\'.",
         call. = FALSE)
  }
  if (!is.character(output_format) || length(output_format) != 1L ||
      !output_format %in% c("multi-file", "long", "wide", "gemma")) {
    stop("`output_format` must be one of \'multi-file\', \'long\', \'wide\' or \'gemma\'.",
         call. = FALSE)
  }
  if (!(isTRUE(to_r) || isFALSE(to_r))) {
    stop("`to_r` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(SNP_effect) || length(SNP_effect) != 1L ||
      !SNP_effect %in% c("Add", "Dom", "Left", "Right")) {
    stop("`SNP_effect` must be one of \'Add\', \'Dom\', \'Left\' or \'Right\'.",
         call. = FALSE)
  }
  if (!is.character(SNP_impute) || length(SNP_impute) != 1L ||
      !SNP_impute %in% c("Major", "Middle", "Minor")) {
    stop("`SNP_impute` must be one of \'Major\', \'Middle\' or \'Minor\'.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' LD window validation
#' @keywords internal
#' @noRd
.v1_validate_ld <- function(ld_min, ld_max, ld_method) {
  ok <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!ok(ld_min) || !ok(ld_max)) {
    stop("`ld_min` and `ld_max` must be single finite numbers.", call. = FALSE)
  }
  if (ld_min > ld_max) {
    stop("`ld_min` (", ld_min, ") must not be larger than `ld_max` (", ld_max,
         ").", call. = FALSE)
  }
  if (ld_max >= 1) {
    stop("`ld_max` must be smaller than 1 (got ", ld_max, "): an absolute LD ",
         "of 1 is met by a marker paired with itself, so the search would ",
         "return the same marker for both traits. Use e.g. `ld_max = 0.99`.",
         call. = FALSE)
  }
  if (!is.character(ld_method) || length(ld_method) != 1L ||
      !ld_method %in% c("composite", "r", "dprime", "corr")) {
    stop("`ld_method` must be one of \'composite\', \'r\', \'dprime\' or \'corr\'.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Seed arithmetic validation (residual seed is (seed + rep) * round(10 * h2))
#'
#' Besides the residual seed, the LD architecture derives retry seeds
#' `seed * s + z (+ rep) (+ x)` with the retry counter `s` running up to 10
#' (`z` <= `rep`, `x` <= the number of QTNs of one class), so its bound is
#' about `10 * seed`. `ld = TRUE` adds that bound; `n_qtn` is the total number
#' of QTNs requested (an upper bound for `x`).
#' @keywords internal
#' @noRd
.v1_validate_seed_arith <- function(seed, rep, h2, null_setting, wide,
                                    ld = FALSE, n_qtn = 0) {
  mult <- 1
  if (!null_setting) {
    h1 <- h2[, 1]
    bad <- h1 > 0 & round(h1 * 10) == 0
    if (rep > 1 && any(bad)) {
      stop(
        "With `rep > 1`, every replicate would use the same residual seed (the ",
        "residual seed is `(seed + replicate) * round(10 * h2[, 1])`, which is 0 ",
        "for h2 <= 0.05, because round(10 * 0.05) is 0 under R\'s round-half-to-even ",
        "rule): identical replicates would be returned. Use h2 > 0.05 ",
        "(or `rep = 1`) for this simulation.",
        call. = FALSE
      )
    }
    mult <- max(1, round(h1 * 10))
  }
  if (!is.null(seed)) {
    ld_mult <- if (ld) 10 else 1   # linkage retry seeds: seed * s, s <= 10
    ld_extra <- if (ld) 2 * rep + max(0, sum(n_qtn, na.rm = TRUE)) else 0
    worst_of <- function(sd) {
      max(mult * abs(sd + rep),
          (if (wide) 2 else 1) * abs(sd) + rep + 100,
          ld_mult * abs(sd) + ld_extra)
    }
    if (worst_of(seed) > .Machine$integer.max) {
      # exact inclusive magnitude bound: the largest a >= 0 such that both
      # seed = a and seed = -a stay within R's integer range (the derived
      # seeds grow with a, so a binary search on the predicate is exact)
      ok <- function(a) {
        worst_of(a) <= .Machine$integer.max &&
          worst_of(-a) <= .Machine$integer.max
      }
      lo <- 0
      hi <- .Machine$integer.max
      centred <- ok(0)
      if (centred) {
        while (lo < hi) {
          mid <- floor((lo + hi + 1) / 2)
          if (ok(mid)) lo <- mid else hi <- mid - 1
        }
      }
      head_msg <- paste0(
        "`seed` (", format(seed, scientific = FALSE), ") is too large: the seeds derived from it ",
        "(e.g. `(seed + rep) * round(10 * h2)`",
        if (ld) paste0(" and, for `architecture = \"LD\"`, the marker-search ",
                       "retry seeds `seed * s + ...` with `s` up to 10") else "",
        ") would overflow R\'s integer range. ")
      if (centred) {
        stop(head_msg, "Use a seed with abs(seed) <= ",
             format(lo, scientific = FALSE), " for this call.",
             call. = FALSE)
      }
      # No acceptance interval centred at zero (a large `rep` shifts it, since
      # the residual seed is (seed + rep) * round(10 * h2)): the accepted
      # seeds are an interval [a, b] (each derived-seed bound is convex in the
      # seed, so the feasible set is an interval). Solve the three bounds in
      # closed form, then nudge to the exact integer boundaries of the real
      # acceptance rule.
      M <- .Machine$integer.max
      a <- max(-rep - floor(M / mult),
               -floor((M - rep - 100) / (if (wide) 2 else 1)),
               -floor((M - ld_extra) / ld_mult))
      b <- min(-rep + floor(M / mult),
               floor((M - rep - 100) / (if (wide) 2 else 1)),
               floor((M - ld_extra) / ld_mult))
      accepted <- function(sd) worst_of(sd) <= M
      if (a <= b) {
        while (a <= b && !accepted(a)) a <- a + 1
        while (a <= b && !accepted(b)) b <- b - 1
        while (a > -M && accepted(a - 1)) a <- a - 1
        while (b < M && accepted(b + 1)) b <- b + 1
      }
      if (a > b) {
        stop(head_msg, "No seed is accepted for this call: the derived seeds ",
             "overflow for every integer seed. Reduce `rep` (or `n_qtn`).",
             call. = FALSE)
      }
      stop(head_msg, "For this call the accepted seeds are the integers in [",
           format(a, scientific = FALSE), ", ", format(b, scientific = FALSE),
           "] (this interval is not centred at 0, so a symmetric magnitude ",
           "bound does not apply; `rep` is large). Use a seed in that ",
           "range, or reduce `rep` / `n_qtn`.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Correlation-matrix validation for `cor` / `cor_res`
#' @keywords internal
#' @noRd
.v1_validate_cor <- function(m, ntraits, name, positive_definite) {
  if (is.data.frame(m)) m <- as.matrix(m)
  if (!is.matrix(m) || !is.numeric(m) || nrow(m) != ntraits ||
      ncol(m) != ntraits) {
    stop("`", name, "` must be a numeric ", ntraits, " x ", ntraits,
         " correlation matrix (one row and column per trait).", call. = FALSE)
  }
  if (anyNA(m) || any(!is.finite(m))) {
    stop("`", name, "` contains missing or non-finite values.", call. = FALSE)
  }
  if (!isSymmetric(unname(m), tol = 1e-8)) {
    stop("`", name, "` is not symmetric: the upper and lower triangles must ",
         "agree.", call. = FALSE)
  }
  if (name == "cor_res" && any(abs(diag(m) - 1) > 1e-8)) {
    stop("`cor_res` must be a correlation matrix (unit diagonal): a diagonal ",
         "different from 1 would change the residual variances and hence the ",
         "heritability.", call. = FALSE)
  }
  if (name == "cor" && any(abs(diag(m) - 1) > 1e-8)) {
    warning("`cor` has a diagonal different from 1; it is treated as a ",
            "covariance-like matrix and the realized correlation is ",
            "cov2cor(cor).", call. = FALSE, immediate. = TRUE)
  }
  ev <- eigen(unname(m), symmetric = TRUE, only.values = TRUE)$values
  if (positive_definite) {
    ok <- all(ev > 0) && !inherits(try(chol(unname(m)), silent = TRUE),
                                   "try-error")
    if (!ok) {
      stop("`", name, "` is not positive definite (smallest eigenvalue ",
           format(min(ev), digits = 3), "), so the requested correlation ",
           "cannot be realized. Supply a positive-definite correlation matrix ",
           "(the v1 engine no longer repairs it silently).", call. = FALSE)
    }
  } else if (min(ev) < -1e-8) {
    stop("`", name, "` is not positive semi-definite (smallest eigenvalue ",
         format(min(ev), digits = 3), "). Supply a valid correlation matrix.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Effect-size specification validation (geometric method)
#'
#' Each effect vector must be either one value (base of a geometric series) or
#' exactly one value per QTN (custom); the two styles cannot be mixed because
#' one `sim_method` decision is taken for all effect classes.
#' @keywords internal
#' @noRd
.v1_validate_effects <- function(add, dom, epi, var, add_effect, dom_effect,
                                 epi_effect, var_effect, len_a, len_d, len_e,
                                 len_v, big) {
  kinds <- character(0)
  chk <- function(eff, len, arg, count_arg, offset = 0) {
    len <- pmax(len - offset, 0)
    for (i in seq_along(eff)) {
      n_exp <- len[min(i, length(len))]
      n_eff <- length(eff[[i]])
      custom_ok <- n_eff == n_exp
      geom_ok <- n_eff == 1L
      if (!custom_ok && !geom_ok) {
        stop("`", arg, "` has ", n_eff, " value(s)",
             if (length(eff) > 1L) paste0(" for trait ", i) else "",
             " but ", n_exp, " QTN(s) need an effect. Provide either a single ",
             "value (the base of a geometric series) or exactly one effect per ",
             "QTN (", count_arg, if (offset > 0) " - 1, because `big_add_QTN_effect` fills the first QTN", ").",
             call. = FALSE)
      }
      kinds <<- c(kinds, arg = if (n_exp == 0) "either" else
        if (custom_ok && geom_ok) "either" else
          if (custom_ok) "custom" else "geometric")
      names(kinds)[length(kinds)] <<- arg
    }
  }
  if (add) chk(add_effect, len_a, "add_effect", "`add_QTN_num`", if (big) 1 else 0)
  if (dom) chk(dom_effect, len_d, "dom_effect", "`dom_QTN_num`")
  if (epi) chk(epi_effect, len_e, "epi_effect", "`epi_QTN_num`")
  if (var) chk(var_effect, len_v, "var_effect", "`var_QTN_num`")
  if (any(kinds == "custom") && any(kinds == "geometric")) {
    stop("Effect sizes are specified inconsistently: ",
         paste(unique(names(kinds)[kinds == "custom"]), collapse = ", "),
         " give one effect per QTN (custom) while ",
         paste(unique(names(kinds)[kinds == "geometric"]), collapse = ", "),
         " give a single geometric-series base. create_phenotypes() applies ",
         "one `sim_method` to all effect classes, so a custom vector would be ",
         "silently exponentiated. Give every effect class one value per QTN ",
         "(with `sim_method = \"custom\"`), or a single base value each.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Record the caller's RNG kind and state (restored by .v1_rng_restore())
#' @keywords internal
#' @noRd
.v1_rng_capture <- function() {
  e <- new.env(parent = emptyenv())
  e$kind <- RNGkind()
  e$seed <- .Random.seed_safe()
  e$advanced <- NULL
  e
}

#' Restore the caller's RNG kind and state on exit from create_phenotypes()
#'
#' The legacy engine switches to `RNGversion("3.5.1")` (sample.kind
#' "Rounding") and calls `set.seed()`. On exit the caller's kind (including
#' `sample.kind`) and `.Random.seed` are put back. When the master seed was
#' drawn (seed = NULL) the caller's stream is left advanced by that one draw,
#' so successive calls with `seed = NULL` still differ. The engine's draw is
#' made under `RNGversion()`, i.e. with the Mersenne-Twister generator; when the
#' caller's generator is the same kind, its state after the draw is spliced back
#' in. When the caller uses another generator (for example L'Ecuyer-CMRG, whose
#' state has a different length) the engine's state cannot be reused, so the
#' caller's own state is restored and then advanced by one `runif()` draw of the
#' caller's own generator.
#' @keywords internal
#' @noRd
.v1_rng_restore <- function(e) {
  suppressWarnings(RNGkind(e$kind[1L], e$kind[2L], e$kind[3L]))
  st <- e$seed
  adv <- e$advanced
  spliced <- FALSE
  if (!is.null(st) && !is.null(adv) && length(adv) == length(st) &&
      (adv[1L] %% 100L) == (st[1L] %% 100L)) {
    st <- c(st[1L], adv[-1L])
    spliced <- TRUE
  }
  .restore_seed(st)
  if (!is.null(st) && !is.null(adv) && !spliced) {
    # foreign generator: advance the caller's own stream by one draw
    invisible(stats::runif(1L))
  }
  invisible(NULL)
}

#' Fail early when more QTNs are requested than the data set has markers
#'
#' Reads the (already validated) QTN counts that `check_in()` injected into the
#' `create_phenotypes()` frame. Constraint-based filtering is checked later by
#' the QTN-selection functions.
#' @keywords internal
#' @noRd
.v1_check_qtn_capacity <- function(env, n_markers) {
  g <- function(x) get0(x, envir = env, inherits = FALSE)
  if (isTRUE(g("null_setting")) || !is.null(unlist(g("QTN_list")))) {
    return(invisible(TRUE))
  }
  epi_w <- g("epi_interaction")
  if (is.null(epi_w)) epi_w <- 2
  need <- c(add = 0, dom = 0, epi = 0, var = 0)
  if (identical(g("architecture"), "partially")) {
    need["add"] <- sum(g("pleio_a"), g("trait_spec_a_QTN_num"))
    need["dom"] <- sum(g("pleio_d"), g("trait_spec_d_QTN_num"))
    need["epi"] <- epi_w * sum(g("pleio_e"), g("trait_spec_e_QTN_num"))
  } else {
    need["add"] <- sum(g("add_QTN_num"))
    need["dom"] <- sum(g("dom_QTN_num"))
    need["epi"] <- epi_w * sum(g("epi_QTN_num"))
    need["var"] <- sum(g("var_QTN_num"))
  }
  worst <- which.max(need)
  if (need[worst] > n_markers) {
    stop("Not enough markers: ", need[worst], " ",
         switch(names(worst), add = "additive", dom = "dominance",
                epi = "epistatic (number of QTNs x `epi_interaction`)",
                var = "variance"),
         " QTN marker(s) were requested but the marker data has only ",
         n_markers, ". Lower the number of QTNs or use a larger marker set.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Variance-QTL standard-deviation multiplier must not be negative
#'
#' The vQTL residual standard deviation is proportional to
#' `1 + sum_i var_effect[i] * (dosage_i + 1)`; a negative value makes
#' `rnorm(sd < 0)` return NaN for those individuals.
#' @keywords internal
#' @noRd
.v1_check_vqtl_sd <- function(QTN, var_effect, var_QTN_num) {
  if (is.null(QTN) || is.null(var_QTN_num) || length(var_QTN_num) != 1L) {
    return(invisible(TRUE))
  }
  QTN <- as.matrix(QTN) + 1
  sigma <- matrix(1, nrow(QTN), 1)
  for (i in seq_len(min(var_QTN_num, ncol(QTN)))) {
    sigma <- sigma + var_effect[i] * QTN[, i]
  }
  if (anyNA(sigma) || any(sigma < 0)) {
    stop("Variance QTL: the standard-deviation multiplier 1 + sum(var_effect * (dosage + 1)) is negative (or missing) for some individuals, ",
         "which would produce missing phenotypes. Use smaller or positive `var_effect` values.",
         call. = FALSE)
  }
  invisible(TRUE)
}
