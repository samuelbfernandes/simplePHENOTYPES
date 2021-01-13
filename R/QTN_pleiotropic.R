#' Select SNPs to be assigned as QTNs.
#' @keywords internal
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add_QTN_num = NULL,
#' @param dom_QTN_num = NULL,
#' @param epi_QTN_num = NULL
#' @param epi_type = NULL,
#' @param epi_interaction = 2,
#' @param same_add_dom_QTN = NULL,
#' @param constraints = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @param verbose = verbose
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Apr 20, 2020
#'
#------------------------------  QTN_pleiotropic -------------------------------
qtn_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           same_add_dom_QTN = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           epi_QTN_num = NULL,
           epi_type = NULL,
           epi_interaction = 2,
           constraints = list(maf_above = NULL, maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL,
           verbose = verbose) {
    #---------------------------------------------------------------------------
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    epi_ef_trait_obj <-  NULL
    add_QTN <- TRUE
    dom_QTN <- TRUE
    epi_QTN <- TRUE
    if (!is.null(add_QTN_num)) {
      if (add_QTN_num == 0) {
        add_QTN <- FALSE
        add_QTN_num <- 1
      }
    }
    if (!is.null(dom_QTN_num)) {
      if (dom_QTN_num == 0) {
        dom_QTN <- FALSE
        dom_QTN_num <- 1
      }
    }
    if (!is.null(epi_QTN_num)) {
      if (epi_QTN_num == 0) {
        epi_QTN <- FALSE
        epi_QTN_num <- 1
      }
    }
    if (any(lengths(constraints) > 0)) {
      index <- constraint(
        genotypes = genotypes,
        maf_above = constraints$maf_above,
        maf_below = constraints$maf_below,
        hets = constraints$hets,
        verbose = verbose
      )
    } else {
      index <- seq_len(nrow(genotypes))
    }
    if (verbose)
      message("* Selecting QTNs")
    if (rep_by != "QTN") {
      rep <- 1
    }
    if (same_add_dom_QTN & add) {
      add_QTN_geno_info <- vector("list", rep)
      add_ef_trait_obj <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + i)
        }
        vector_of_add_QTN <-
          sample(index, add_QTN_num, replace = FALSE)
        times <- 1
        dif <- c()
        while (!any(genotypes[vector_of_add_QTN, - (1:5)] == 0) &
               times <= 10 & dom) {
          if (!is.null(seed)) {
            set.seed(seed + i)
          }
          dif <- c(dif, vector_of_add_QTN)
          vector_of_add_QTN <-
            sample(setdiff(index, dif), add_QTN_num, replace = FALSE)
          times <- times + 1
        }
        add_QTN_geno_info[[i]] <-
          as.data.frame(genotypes[vector_of_add_QTN, ],
                        check.names = FALSE,
                        fix.empty.names = FALSE)
        add_ef_trait_obj[[i]] <-
          t(add_QTN_geno_info[[i]][, - c(1:5)])
        colnames(add_ef_trait_obj[[i]]) <-
          paste0("Chr_",
                 unlist(add_QTN_geno_info[[i]][, 3]),
                 "_",
                 unlist(add_QTN_geno_info[[i]][, 4]))
      }
      add_QTN_geno_info <-
        do.call(rbind, add_QTN_geno_info)
      ns <- ncol(genotypes) - 5
      maf <-
        round(apply(add_QTN_geno_info[, -c(1:5)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
      names(maf) <- add_QTN_geno_info[, 1]
      add_QTN_geno_info <- data.frame(
        add_QTN_geno_info[, 1:5],
        maf = maf,
        add_QTN_geno_info[, - c(1:5)],
        check.names = FALSE,
        fix.empty.names = FALSE
      )
      add_QTN_geno_info <-
        data.frame(
          rep = rep(1:rep, each = add_QTN_num),
          add_QTN_geno_info,
          check.names = FALSE,
          fix.empty.names = FALSE
        )
      if (!export_gt) {
        add_QTN_geno_info <- add_QTN_geno_info[, 1:7]
      }
      if (!is.null(seed)) {
        s <- as.matrix(seed + 1:rep)
      } else {
        s <- "set.seed not assigned"
      }
      if (add_QTN) {
        if (verbose){
          write.table(
          s,
          paste0("Seed_num_for_", add_QTN_num,
                 "_Add_and_Dom_QTN.txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
      }
        data.table::fwrite(
          add_QTN_geno_info,
          "Additive_and_Dominance_Selected_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
    } else {
      if (add) {
        add_QTN_geno_info <- vector("list", rep)
        add_ef_trait_obj <- vector("list", rep)
        for (i in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + i)
          }
          vector_of_add_QTN <-
            sample(index, add_QTN_num, replace = FALSE)
          add_QTN_geno_info[[i]] <-
            as.data.frame(genotypes[vector_of_add_QTN, ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
          add_ef_trait_obj[[i]] <-
            t(add_QTN_geno_info[[i]][, - c(1:5)])
          colnames(add_ef_trait_obj[[i]]) <-
            paste0("Chr_",
                   unlist(add_QTN_geno_info[[i]][, 3]),
                   "_",
                   unlist(add_QTN_geno_info[[i]][, 4]))
        }
        add_QTN_geno_info <-
          do.call(rbind, add_QTN_geno_info)
        ns <- ncol(genotypes) - 5
        maf <-
          round(apply(add_QTN_geno_info[, -c(1:5)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
        names(maf) <- add_QTN_geno_info[, 1]
        add_QTN_geno_info <- data.frame(
          add_QTN_geno_info[, 1:5],
          maf = maf,
          add_QTN_geno_info[, - c(1:5)],
          check.names = FALSE,
          fix.empty.names = FALSE
        )
        add_QTN_geno_info <-
          data.frame(
            rep = rep(1:rep, each = add_QTN_num),
            add_QTN_geno_info,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        if (!export_gt) {
          add_QTN_geno_info <- add_QTN_geno_info[, 1:7]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed + 1:rep)
        } else {
          s <- "set.seed not assigned"
        }
        if (add_QTN) {
          if (verbose){
          write.table(
            s,
            paste0("Seed_num_for_", add_QTN_num,
                   "_Add_QTN.txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          }
          data.table::fwrite(
            add_QTN_geno_info,
            "Additive_Selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
      if (dom) {
        dom_QTN_geno_info <- vector("list", rep)
        dom_ef_trait_obj <- vector("list", rep)
        for (i in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + i + rep)
          }
          vector_of_dom_QTN <-
            sample(index, dom_QTN_num, replace = FALSE)
          times <- 1
          dif <- c()
          while (!any(genotypes[vector_of_dom_QTN, - (1:5)] == 0) &
                 times <= 10) {
            if (!is.null(seed)) {
              set.seed(seed + i + rep)
            }
            dif <- c(dif, vector_of_dom_QTN)
            vector_of_dom_QTN <-
              sample(setdiff(index, dif), dom_QTN_num, replace = FALSE)
            times <- times + 1
          }
          dom_QTN_geno_info[[i]] <-
            as.data.frame(genotypes[vector_of_dom_QTN, ],
                          check.names = FALSE,
                          fix.empty.names = FALSE)
          dom_ef_trait_obj[[i]] <-
            t(dom_QTN_geno_info[[i]][, - c(1:5)])
          colnames(dom_ef_trait_obj[[i]]) <-
            paste0("Chr_",
                   unlist(dom_QTN_geno_info[[i]][, 3]),
                   "_",
                   unlist(dom_QTN_geno_info[[i]][, 4]))
        }
        dom_QTN_geno_info <-
          do.call(rbind, dom_QTN_geno_info)
        ns <- ncol(genotypes) - 5
        maf <-
          round(apply(dom_QTN_geno_info[, -c(1:5)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
        names(maf) <- dom_QTN_geno_info[, 1]
        dom_QTN_geno_info <- data.frame(
          dom_QTN_geno_info[, 1:5],
          maf = maf,
          dom_QTN_geno_info[, - c(1:5)],
          check.names = FALSE,
          fix.empty.names = FALSE
        )
        dom_QTN_geno_info <-
          data.frame(
            rep = rep(1:rep, each = dom_QTN_num),
            dom_QTN_geno_info,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        if (!export_gt) {
          dom_QTN_geno_info <- dom_QTN_geno_info[, 1:7]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed + 1:rep) + rep
        } else {
          s <- "set.seed not assigned"
        }
        if (dom_QTN) {
          if (verbose){
          write.table(
            s,
            paste0("Seed_num_for_", dom_QTN_num,
                   "_Dom_QTN.txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          }
          data.table::fwrite(
            dom_QTN_geno_info,
            "Dominance_Selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
    }
    if (epi) {
      epi_QTN_gen_infor <- vector("list", rep)
      epi_ef_trait_obj <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + seed + i + rep)
        }
        vector_of_epi_QTN <-
          sample(index, (epi_interaction * epi_QTN_num), replace = FALSE)
        epi_QTN_gen_infor[[i]] <-
          as.data.frame(genotypes[vector_of_epi_QTN, ],
                        check.names = FALSE,
                        fix.empty.names = FALSE)
        epi_ef_trait_obj[[i]] <-
          t(epi_QTN_gen_infor[[i]][, - c(1:5)])
        colnames(epi_ef_trait_obj[[i]]) <-
          paste0("Chr_",
                 unlist(epi_QTN_gen_infor[[i]][, 3]),
                 "_",
                 unlist(epi_QTN_gen_infor[[i]][, 4]))
      }
      epi_QTN_gen_infor <-
        do.call(rbind, epi_QTN_gen_infor)
      ns <- ncol(genotypes) - 5
      maf <-
        round(apply(epi_QTN_gen_infor[, -c(1:5)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
      names(maf) <- epi_QTN_gen_infor[, 1]
      epi_QTN_gen_infor <- data.frame(
        epi_QTN_gen_infor[, 1:5],
        maf = maf,
        epi_QTN_gen_infor[, - c(1:5)],
        check.names = FALSE,
        fix.empty.names = FALSE
      )
      epi_QTN_gen_infor <-
        data.frame(
          rep = rep(rep(1:rep, each = epi_QTN_num), each = epi_interaction),
          QTN = rep(1:epi_QTN_num, each = epi_interaction),
          epi_QTN_gen_infor,
          check.names = FALSE,
          fix.empty.names = FALSE
        )
      if (!export_gt) {
        epi_QTN_gen_infor <- epi_QTN_gen_infor[, 1:8]
      }
      if (!is.null(seed)) {
        ss <- as.matrix((seed + seed) + 1:rep + rep)
      } else {
        ss <- "set.seed not assigned"
      }
      if (epi_QTN) {
        if (verbose){
        write.table(
          ss,
          paste0("Seed_num_for_", epi_QTN_num,
                 "_Epi_QTN.txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        }
        data.table::fwrite(
          epi_QTN_gen_infor,
          "Epistatic_Selected_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
    }
    if (!is.null(add_ef_trait_obj) & !add_QTN) {
      add_ef_trait_obj <- lapply(add_ef_trait_obj, function(x) {
        rnames <- rownames(x)
        x <- as.matrix(rep(0, nrow(x)))
        rownames(x)  <- rnames
        return(x)
      })
    }
    if (!is.null(dom_ef_trait_obj) & !dom_QTN) {
      dom_ef_trait_obj <- lapply(dom_ef_trait_obj, function(x) {
        rnames <- rownames(x)
        x <- as.matrix(rep(0, nrow(x)))
        rownames(x)  <- rnames
        return(x)
      })
    }
    if (!is.null(epi_ef_trait_obj) & !epi_QTN) {
      epi_ef_trait_obj <- lapply(epi_ef_trait_obj, function(x) {
        rnames <- rownames(x)
        x <- as.matrix(rep(0, nrow(x)))
        rownames(x)  <- rnames
        return(x)
      })
    }
    if (!is.null(add_ef_trait_obj)) {
      biallelic <- any(unlist(lapply(add_ef_trait_obj, function(x) {
        apply(x, 2, function(y) {
          length(unique(y)) > 3
        })
      }),
      recursive = T))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    if (!is.null(dom_ef_trait_obj)) {
      biallelic <- any(unlist(lapply(dom_ef_trait_obj, function(x) {
        apply(x, 2, function(y) {
          length(unique(y)) > 3
        })
      }),
      recursive = T))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    if (!is.null(epi_ef_trait_obj)) {
      biallelic <- any(sapply(epi_ef_trait_obj, function(x) {
        apply(x, 2, function(y) {
          length(unique(y)) > 3
        })
      }))
      if (biallelic) {
        stop("Please use only biallelic markers.",
             call. = F)
      }
    }
    return(
      list(
        add_ef_trait_obj = add_ef_trait_obj,
        dom_ef_trait_obj = dom_ef_trait_obj,
        epi_ef_trait_obj = epi_ef_trait_obj
      )
    )
  }
