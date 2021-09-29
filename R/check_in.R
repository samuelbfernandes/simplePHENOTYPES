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
#' @return Phenotypes for ntraits traits
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Apr 20, 2020
#'
#'----------------------------phenotypes---------------------------------------
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
           RNGversion = '3.5.1') {
    if (is.null(seed)) {
      seed <- as.integer(runif(1, 0, 1000000))
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
      ntraits <- 2
      if (type_of_ld != "indirect" & type_of_ld != "direct") {
        stop("Parameter \'type_of_ld\' should be either \'direct\' or \'indirect\'.",
             call. = F)
      }
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
    } else if (length(mean) != ntraits) {
      stop("Parameter \'mean\' should have length = \'ntraits\'.",
           call. = F)
    }
    #---- h2 ----
    h2 <- as.matrix(h2)
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
    #----- QTN_list and QTN number -----
    if (!is.null(unlist(QTN_list))) {
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
      if (same_add_dom_QTN) {
        if (is.null(unlist(QTN_list$dom))) {
          QTN_list$dom <- QTN_list$add
        } else if (!all.equal(QTN_list$dom, QTN_list$add)) {
          stop(
            "If \'same_add_dom_QTN = TRUE\', \'QTN_list$dom\'  should be NULL or identical to \'QTN_list$add\'. Instead, the QTNs in \'QTN_list$add\' will be used.",
            call. = F
          )
        }
      }
      if (same_mv_QTN) {
        if (is.null(unlist(QTN_list$var))) {
          QTN_list$var <- QTN_list$add
        } else if (!all.equal(QTN_list$var, QTN_list$add)) {
          stop(
            "If \'same_mv_QTN = TRUE\', \'QTN_list$var\' should be NULL or identical to \'QTN_list$add\'. Instead, the QTNs in \'QTN_list$add\' will be used.",
            call. = F
          )
        }
      }
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
          ntraits <- length(QTN_list$add)
          warning(
            "Setting ntraits = length(QTN_list$add)!",
            call. = F,
            immediate. = T
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
          ntraits <- length(QTN_list$dom)
          warning(
            "Setting ntraits = length(QTN_list$dom)!",
            call. = F,
            immediate. = T
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
        if (unique(lengths(QTN_list$epi) %% epi_interaction) != 0) {
          stop(paste("epi_interaction =", epi_interaction, "Please provide",epi_interaction, "Markers should be provided for each epistatic QTN."),
               call. = F)
        }
        
        if (length(QTN_list$epi) != ntraits) {
          ntraits <- length(QTN_list$epi)
          warning(
            "Setting ntraits = length(QTN_list$epi)!",
            call. = F,
            immediate. = T
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
    if (sim_method != "geometric" & sim_method != "custom") {
      stop("Parameter \'sim_method\' should be either \'geometric\' or \'custom\'!",
           call. = F)
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
        paste0("\n",paste0(rep("_", nchar(title)-1), collapse = ""),"\n"), "\nMaster Seed:", seed, "\nNumber of traits:", ntraits)
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