#' Select SNPs to be assigned as QTNs
#' @keywords internal
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add_QTN_num = NULL,
#' @param dom_QTN_num = NULL,
#' @param same_add_dom_QTN = NULL,
#' @param dom = NULL,
#' @param add = NULL,
#' @param ld = NULL,
#' @param gdsfile NULL
#' @param constraints = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @param type_of_ld = NULL
#' @param verbose = verbose
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes
#' Last update: Apr 20, 2020
#'
#'----------------------------- QTN_linkage ------------------------------------
qtn_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           ld = NULL,
           gdsfile = NULL,
           constraints = list(maf_above = NULL, maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           same_add_dom_QTN = NULL,
           add = NULL,
           dom = NULL,
           type_of_ld = NULL,
           verbose = verbose) {
    #---------------------------------------------------------------------------
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    add_QTN <- TRUE
    dom_QTN <- TRUE
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
    if (rep_by != "QTN") {
      rep <- 1
    }
    if (any(lengths(constraints) > 0)) {
      index <- constraint(
        genotypes = genotypes,
        maf_above = constraints$maf_above,
        maf_below = constraints$maf_below,
        hets = constraints$hets,
        verbose = verbose
      )
      if (add) {
        if (length(index) < add_QTN_num) {
          stop("Not enough SNP left after applying the selected constrain!",
               call. = F)
        }
      }
      if (dom) {
        if (length(index) < dom_QTN_num) {
          stop("Not enough SNP left after applying the selected constrain!",
               call. = F)
        }
      }
    } else {
      index <- seq_len(nrow(genotypes))
    }
    n <- max(index)
    if (verbose)
      message("* Selecting QTNs")
    if (type_of_ld == "indirect") {
      if (same_add_dom_QTN & add) {
        sup <- vector("list", rep)
        inf <- vector("list", rep)
        add_gen_info_sup <- vector("list", rep)
        add_gen_info_inf <- vector("list", rep)
        QTN_causing_ld <- vector("list", rep)
        results <- vector("list", rep)
        LD_summary <- vector("list", rep)
        seed_num <- c()
        for (z in 1:rep) {
          s <- 1
          border <- TRUE
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          while (s <= 10 & border) {
            if (!is.null(seed)) {
              seed_num[z] <- (seed * s) + z
              set.seed(seed_num[z])
            }
            vector_of_add_QTN <-
              sample(index, add_QTN_num, replace = FALSE)
            sup_temp <- c()
            inf_temp <- c()
            ld_between_QTNs_temp <- c()
            actual_ld_sup <- c()
            actual_ld_inf <- c()
            again <- TRUE
            times <- 1
            dif <- c()
            while (again) {
              x <- 1
              for (j in vector_of_add_QTN) {
                ldsup <- 1
                i <- j + 1
                while (ldsup > ld) {
                  if (i > n) {
                    if (verbose)
                      warning(
                        "There are no SNPs downstream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp1 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, j),
                      count = c(-1, 1)
                    )
                  snp2 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i),
                      count = c(-1, 1)
                    )
                  ldsup <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                  if (is.nan(ldsup)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i <- i + 1
                }
                actual_ld_sup[x] <- ldsup
                sup_temp[x] <- i - 1
                ldinf <- 1
                i2 <- j - 1
                while (ldinf > ld) {
                  if (i2 < 1) {
                    if (verbose)
                      warning(
                        "There are no SNPs uptream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp3 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i2),
                      count = c(-1, 1)
                    )
                  ldinf <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = "composite"))
                  if (is.nan(ldinf)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i2 <- i2 - 1
                }
                actual_ld_inf[x]  <- ldinf
                inf_temp[x] <- i2 + 1
                snp_sup <-
                  gdsfmt::read.gdsn(
                    gdsfmt::index.gdsn(genofile, "genotype"),
                    start = c(1, sup_temp[x]),
                    count = c(-1, 1)
                  )
                snp_inf <-
                  gdsfmt::read.gdsn(
                    gdsfmt::index.gdsn(genofile, "genotype"),
                    start = c(1, inf_temp[x]),
                    count = c(-1, 1)
                  )
                ld_between_QTNs_temp[x] <-
                  SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = "composite")
                x <- x + 1
              }
              if ((!any(genotypes[sup_temp, - (1:5)] == 0) |
                   !any(genotypes[inf_temp, - (1:5)] == 0)) &
                  times <= 10 & dom) {
                if (!is.null(seed)) {
                  seed_num[z] <- (seed * s) + z
                  set.seed(seed_num[z])
                }
                dif <- c(dif, vector_of_add_QTN)
                vector_of_add_QTN <-
                  sample(setdiff(index, dif), add_QTN_num, replace = FALSE)
                again <- TRUE
              } else {
                again <- FALSE
              }
              times <- times + 1
            }
            if (i > n | i2 < 1) {
              border <- TRUE
            } else {
              border <- FALSE
            }
            s <- s + 1
          }
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- inf_temp
          QTN_causing_ld[[z]] <-
            data.frame(
              snp_type = "cause_of_LD",
              genotypes[vector_of_add_QTN, ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup[[z]] <-
            data.frame(
              snp_type = "QTN_upstream",
              genotypes[sup[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_inf[[z]] <-
            data.frame(
              snp_type = "QTN_downstream",
              genotypes[inf[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                add_gen_info_sup[[z]],
                                add_gen_info_inf[[z]])
          LD_summary[[z]] <- data.frame(
            z,
            QTN_causing_ld[[z]][, 2],
            ld,
            actual_ld_inf,
            actual_ld_sup,
            add_gen_info_inf[[z]][, 2],
            add_gen_info_sup[[z]][, 2],
            ld_between_QTNs_temp,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary[[z]]) <-
            c(
              "rep",
              "SNP_causing_LD",
              "input_LD (absolute value)",
              "Actual_LD_with_QTN_of_Trait_1",
              "Actual_LD_with_QTN_of_Trait_2",
              "QTN_for_trait_1",
              "QTN_for_trait_2",
              "LD_between_QTNs"
            )
        }
        LD_summary <- do.call(rbind, LD_summary)
        data.table::fwrite(
          LD_summary,
          "LD_summary.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        results <- do.call(rbind, results)
        ns <- nrow(genotypes) - 5
        maf <- round(apply(results[, - c(1:6)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx, (1 - sumx))
        }), 4)
        names(maf) <- results[, 2]
        results <- data.frame(
          results[, 1:6],
          maf = maf,
          results[, - c(1:6)],
          check.names = FALSE,
          fix.empty.names = FALSE
        )
        results <-
          data.frame(
            rep = rep(1:rep, each = add_QTN_num * 3),
            results,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        if (!export_gt) {
          results <- results[, 1:8]
        }
        if (add_QTN) {
          write.table(
            seed_num,
            paste0("seed_num_for_", add_QTN_num,
                   "_Add_and_Dom_QTN.txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          data.table::fwrite(
            results,
            "Additive_and_Dominance_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
        add_ef_trait_obj <- mapply(function(x, y) {
          rownames(x) <-
            paste0("Chr_", x$chr, "_", x$pos)
          rownames(y) <-
            paste0("Chr_", y$chr, "_", y$pos)
          b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
          return(b)
        },
        x = add_gen_info_sup,
        y = add_gen_info_inf,
        SIMPLIFY = F)
      } else {
        if (add) {
          sup <- vector("list", rep)
          inf <- vector("list", rep)
          add_gen_info_sup <- vector("list", rep)
          add_gen_info_inf <- vector("list", rep)
          QTN_causing_ld <- vector("list", rep)
          results_add <- vector("list", rep)
          LD_summary_add <- vector("list", rep)
          seed_num <- c()
          for (z in 1:rep) {
            s <- 1
            border <- TRUE
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            while (s <= 10 & border) {
              if (!is.null(seed)) {
                seed_num[z] <-  (seed * s) + z
                set.seed(seed_num[z])
              }
              vector_of_add_QTN <-
                sample(index, add_QTN_num, replace = FALSE)
              x <- 1
              sup_temp <- c()
              inf_temp <- c()
              ld_between_QTNs_temp <- c()
              actual_ld_sup <- c()
              actual_ld_inf <- c()
              for (j in vector_of_add_QTN) {
                ldsup <- 1
                i <- j + 1
                while (ldsup > ld) {
                  if (i > n) {
                    if (verbose)
                      warning(
                        "There are no SNPs downstream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp1 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, j),
                      count = c(-1, 1)
                    )
                  snp2 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i),
                      count = c(-1, 1)
                    )
                  ldsup <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                  if (is.nan(ldsup)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i <- i + 1
                }
                actual_ld_sup[x] <- ldsup
                sup_temp[x] <- i - 1
                ldinf <- 1
                i2 <- j - 1
                while (ldinf > ld) {
                  if (i2 < 1) {
                    if (verbose)
                      warning(
                        "There are no SNPs upstream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp3 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i2),
                      count = c(-1, 1)
                    )
                  ldinf <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = "composite"))
                  if (is.nan(ldinf)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i2 <- i2 - 1
                }
                actual_ld_inf[x] <- ldinf
                inf_temp[x] <- i2 + 1
                snp_sup <-
                  gdsfmt::read.gdsn(
                    gdsfmt::index.gdsn(genofile, "genotype"),
                    start = c(1, sup_temp[x]),
                    count = c(-1, 1)
                  )
                snp_inf <-
                  gdsfmt::read.gdsn(
                    gdsfmt::index.gdsn(genofile, "genotype"),
                    start = c(1, inf_temp[x]),
                    count = c(-1, 1)
                  )
                ld_between_QTNs_temp[x] <-
                  SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = "composite")
                x <- x + 1
              }
              if (i > n | i2 < 1) {
                border <- TRUE
              } else {
                border <- FALSE
              }
              s <- s + 1
            }
            SNPRelate::snpgdsClose(genofile)
            sup[[z]] <- sup_temp
            inf[[z]] <- inf_temp
            QTN_causing_ld[[z]] <-
              data.frame(
                snp_type = "cause_of_LD",
                genotypes[vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup[[z]] <-
              data.frame(
                snp_type = "QTN_upstream",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_inf[[z]] <-
              data.frame(
                snp_type = "QTN_downstream",
                genotypes[inf[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add[[z]] <- rbind(QTN_causing_ld[[z]],
                                      add_gen_info_sup[[z]],
                                      add_gen_info_inf[[z]])
            LD_summary_add[[z]] <- data.frame(
              z,
              QTN_causing_ld[[z]][, 2],
              ld,
              actual_ld_inf,
              actual_ld_sup,
              add_gen_info_inf[[z]][, 2],
              add_gen_info_sup[[z]][, 2],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add[[z]]) <-
              c(
                "rep",
                "SNP_causing_LD",
                "input_LD (absolute value)",
                "Actual_LD_with_QTN_of_Trait_1",
                "Actual_LD_with_QTN_of_Trait_2",
                "QTN_for_trait_1",
                "QTN_for_trait_2",
                "LD_between_QTNs"
              )
          }
          LD_summary_add <- do.call(rbind, LD_summary_add)
          data.table::fwrite(
            LD_summary_add,
            "LD_summary_additive.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_add <- do.call(rbind, results_add)
          ns <- nrow(genotypes) - 5
          maf <- round(apply(results_add[, -c(1:6)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_add[, 2]
          results_add <-
            data.frame(
              results_add[, 1:6],
              maf = maf,
              results_add[, -c(1:6)],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results_add <-
            data.frame(
              rep = rep(1:rep, each = add_QTN_num * 3),
              results_add,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results_add <- results_add[, 1:8]
          }
          if (add_QTN) {
            write.table(
              seed_num,
              paste0("seed_num_for_", add_QTN_num,
                     "_Add_QTN.txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            data.table::fwrite(
              results_add,
              "Additive_selected_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          add_ef_trait_obj <- mapply(function(x, y) {
            rownames(x) <-
              paste0("Chr_", x$chr, "_", x$pos)
            rownames(y) <-
              paste0("Chr_",  y$chr, "_", y$pos)
            b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
            return(b)
          },
          x = add_gen_info_sup,
          y = add_gen_info_inf,
          SIMPLIFY = F)
        }
        if (dom) {
          sup <- vector("list", rep)
          inf <- vector("list", rep)
          dom_gen_info_sup <- vector("list", rep)
          dom_gen_info_inf <- vector("list", rep)
          QTN_causing_ld <- vector("list", rep)
          results_dom <- vector("list", rep)
          LD_summary_dom <- vector("list", rep)
          seed_num <- c()
          for (z in 1:rep) {
            s <- 1
            border <- TRUE
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            while (s <= 10 & border) {
              if (!is.null(seed)) {
                seed_num[z] <- (seed * s) + z + rep
                set.seed(seed_num[z])
              }
              vector_of_dom_QTN <-
                sample(index, dom_QTN_num, replace = FALSE)
              sup_temp <- c()
              inf_temp <- c()
              ld_between_QTNs_temp <- c()
              actual_ld_sup <- c()
              actual_ld_inf <- c()
              again <- TRUE
              times <- 1
              dif <- c()
              while (again) {
                x <- 1
                for (j in vector_of_dom_QTN) {
                  ldsup <- 1
                  i <- j + 1
                  while (ldsup > ld) {
                    if (i > n) {
                      if (verbose)
                        warning(
                          "There are no SNPs downstream. Selecting a different seed number",
                          call. = F,
                          immediate. = T
                        )
                      break
                    }
                    snp1 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, j),
                        count = c(-1, 1)
                      )
                    snp2 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, i),
                        count = c(-1, 1)
                      )
                    ldsup <-
                      abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                    if (is.nan(ldsup)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    i <- i + 1
                  }
                  actual_ld_sup[x] <- ldsup
                  sup_temp[x] <- i - 1
                  ldinf <- 1
                  i2 <- j - 1
                  while (ldinf > ld) {
                    if (i2 < 1) {
                      if (verbose)
                        warning(
                          "There are no SNPs upstream. Selecting a different seed number",
                          call. = F,
                          immediate. = T
                        )
                      break
                    }
                    snp3 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, i2),
                        count = c(-1, 1)
                      )
                    ldinf <-
                      abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = "composite"))
                    if (is.nan(ldinf)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    i2 <- i2 - 1
                  }
                  actual_ld_inf[x] <- ldinf
                  inf_temp[x] <- i2 + 1
                  snp_sup <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, sup_temp[x]),
                      count = c(-1, 1)
                    )
                  snp_inf <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, inf_temp[x]),
                      count = c(-1, 1)
                    )
                  ld_between_QTNs_temp[x] <-
                    SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = "composite")
                  x <- x + 1
                }
                if ((!any(genotypes[sup_temp, - (1:5)] == 0) |
                     !any(genotypes[inf_temp, - (1:5)] == 0)) &
                    times <= 10) {
                  if (!is.null(seed)) {
                    seed_num[z] <- (seed * s) + z + rep
                    set.seed(seed_num[z])
                  }
                  dif <- c(dif, vector_of_dom_QTN)
                  vector_of_dom_QTN <-
                    sample(setdiff(index, dif), dom_QTN_num, replace = FALSE)
                  again <- TRUE
                } else {
                  again <- FALSE
                }
                times <- times + 1
              }
              if (i > n | i2 < 1) {
                border <- TRUE
              } else {
                border <- FALSE
              }
              s <- s + 1
            }
            SNPRelate::snpgdsClose(genofile)
            sup[[z]] <- sup_temp
            inf[[z]] <- inf_temp
            QTN_causing_ld[[z]] <-
              data.frame(
                snp_type = "cause_of_LD",
                genotypes[vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup[[z]] <-
              data.frame(
                snp_type = "QTN_upstream",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_inf[[z]] <-
              data.frame(
                snp_type = "QTN_downstream",
                genotypes[inf[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom[[z]] <- rbind(QTN_causing_ld[[z]],
                                      dom_gen_info_sup[[z]],
                                      dom_gen_info_inf[[z]])
            LD_summary_dom[[z]] <- data.frame(
              z,
              QTN_causing_ld[[z]][, 2],
              ld,
              actual_ld_inf,
              actual_ld_sup,
              dom_gen_info_inf[[z]][, 2],
              dom_gen_info_sup[[z]][, 2],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom[[z]]) <-
              c(
                "rep",
                "SNP_causing_LD",
                "input_LD (absolute value)",
                "Actual_LD_with_QTN_of_Trait_1",
                "Actual_LD_with_QTN_of_Trait_2",
                "QTN_for_trait_1",
                "QTN_for_trait_2",
                "LD_between_QTNs"
              )
          }
          LD_summary_dom <- do.call(rbind, LD_summary_dom)
          data.table::fwrite(
            LD_summary_dom,
            "LD_summary_dominance.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_dom <- do.call(rbind, results_dom)
          ns <- nrow(genotypes) - 5
          maf <- round(apply(results_dom[, -c(1:6)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_dom[, 2]
          results_dom <-
            data.frame(
              results_dom[, 1:6],
              maf = maf,
              results_dom[, - c(1:6)],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results_dom <-
            data.frame(
              rep = rep(1:rep, each = dom_QTN_num * 3),
              results_dom,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results_dom <- results_dom[, 1:8]
          }
          if (add_QTN) {
            write.table(
              seed_num,
              paste0("seed_num_for_", dom_QTN_num,
                     "_Dom_QTN.txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            data.table::fwrite(
              results_dom,
              "Dominance_selected_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          dom_ef_trait_obj <- mapply(function(x, y) {
            rownames(x) <-
              paste0("Chr_", x$chr, "_", x$pos)
            rownames(y) <-
              paste0("Chr_", y$chr, "_", y$pos)
            b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
            return(b)
          },
          x = dom_gen_info_sup,
          y = dom_gen_info_inf,
          SIMPLIFY = F)
        }
      }
      if (!is.null(add_ef_trait_obj) & !add_QTN) {
        add_ef_trait_obj <- lapply(add_ef_trait_obj, function(x) {
          lapply(x, function(y) {
            rnames <- rownames(y)
            y <- matrix(0, nrow = nrow(y), ncol =  1)
            rownames(y)  <- rnames
            return(y)
          })
        })
      }
      if (!is.null(dom_ef_trait_obj) & !dom_QTN) {
        dom_ef_trait_obj <- lapply(dom_ef_trait_obj, function(x) {
          lapply(x, function(y) {
            rnames <- rownames(y)
            y <- matrix(0, nrow = nrow(y), ncol =  1)
            rownames(y)  <- rnames
            return(y)
          })
        })
      }
      if (!is.null(add_ef_trait_obj)) {
        biallelic <- any(unlist(lapply(add_ef_trait_obj, function(x) {
          sapply(x, function(x2) {
            apply(x2, 2, function(y) {
              length(unique(y)) > 3
            })
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
          sapply(x, function(x2) {
            apply(x2, 2, function(y) {
              length(unique(y)) > 3
            })
          })
        }),
        recursive = T))
        if (biallelic) {
          stop("Please use only biallelic markers.",
               call. = F)
        }
      }
      return(list(
        add_ef_trait_obj = add_ef_trait_obj,
        dom_ef_trait_obj = dom_ef_trait_obj
      ))
    } else {
      if (same_add_dom_QTN & add) {
        sup <- vector("list", rep)
        inf <- vector("list", rep)
        add_gen_info_sup <- vector("list", rep)
        add_gen_info_inf <- vector("list", rep)
        results <- vector("list", rep)
        LD_summary <- vector("list", rep)
        seed_num <- c()
        for (z in 1:rep) {
          s <- 1
          border <- TRUE
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          while (s <= 10 & border) {
            if (!is.null(seed)) {
              seed_num[z] <- (seed * s) + z
              set.seed(seed_num[z])
            }
            vector_of_add_QTN <-
              sample(index, add_QTN_num, replace = FALSE)
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            x <- 1
            sup_temp <- c()
            ld_between_QTNs_temp <- c()
            again <- TRUE
            times <- 1
            dif <- c()
            while (again) {
              for (j in vector_of_add_QTN) {
                ldsup <- 1
                i <- j + 1
                while (ldsup > ld) {
                  if (i > n) {
                    if (verbose)
                      warning(
                        "There are no SNPs downstream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp1 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, j),
                      count = c(-1, 1)
                    )
                  snp2 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i),
                      count = c(-1, 1)
                    )
                  ldsup <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                  if (is.nan(ldsup)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i <- i + 1
                }
                ld_between_QTNs_temp[x] <- ldsup
                sup_temp[x] <- i - 1
                x <- x + 1
              }
              if ((!any(genotypes[sup_temp, - (1:5)] == 0) |
                   !any(genotypes[vector_of_dom_QTN, - (1:5)] == 0)) &
                  times <= 10 & dom) {
                if (!is.null(seed)) {
                  seed_num[z] <- (seed * s) + z
                  set.seed(seed_num[z])
                }
                dif <- c(dif, vector_of_add_QTN)
                vector_of_add_QTN <-
                  sample(setdiff(index, dif), add_QTN_num, replace = FALSE)
                again <- TRUE
              } else {
                again <- FALSE
              }
              times <- times + 1
            }
            if (i > n | i2 < 1) {
              border <- TRUE
            } else {
              border <- FALSE
            }
            s <- s + 1
          }
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- vector_of_add_QTN
          add_gen_info_inf[[z]] <-
            data.frame(
              snp_type = "QTN_for_trait_1",
              genotypes[vector_of_add_QTN, ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup[[z]] <-
            data.frame(
              snp_type = "QTN_for_trait_2",
              genotypes[sup[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results[[z]] <- rbind(add_gen_info_inf[[z]],
                                add_gen_info_sup[[z]])
          LD_summary[[z]] <- data.frame(
            z,
            ld,
            ld_between_QTNs_temp,
            add_gen_info_inf[[z]][, 2],
            add_gen_info_sup[[z]][, 2],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary[[z]]) <-
            c(
              "rep",
              "Aimed_LD (absolute value)",
              "Actual_LD ",
              "QTN_for_trait_1",
              "QTN_for_trait_2"
            )
        }
        LD_summary <- do.call(rbind, LD_summary)
        data.table::fwrite(
          LD_summary,
          "LD_summary.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        results <- do.call(rbind, results)
        ns <- nrow(genotypes) - 5
        maf <- round(apply(results[, -c(1:6)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        names(maf) <- results[, 2]
        results <-
          data.frame(
            results[, 1:6],
            maf = maf,
            results[, - c(1:6)],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        results <-
          data.frame(
            rep = rep(1:rep, each = add_QTN_num * 2),
            results,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        if (!export_gt) {
          results <- results[, 1:8]
        }
        if (add_QTN) {
          write.table(
            seed_num,
            paste0(
              "seed_num_for_",
              add_QTN_num,
              "_Add_and_Dom_QTN",
              ".txt"
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          data.table::fwrite(
            results,
            "Additive_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
        add_ef_trait_obj <- mapply(function(x, y) {
          rownames(x) <-
            paste0("Chr_", x$chr, "_", x$pos)
          rownames(y) <-
            paste0("Chr_", y$chr, "_", y$pos)
          b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
          return(b)
        },
        x = add_gen_info_sup,
        y = add_gen_info_inf,
        SIMPLIFY = F)
      } else {
        if (add) {
          sup <- vector("list", rep)
          inf <- vector("list", rep)
          add_gen_info_sup <- vector("list", rep)
          add_gen_info_inf <- vector("list", rep)
          results_add <- vector("list", rep)
          LD_summary_add <- vector("list", rep)
          seed_num <- c()
          for (z in 1:rep) {
            s <- 1
            border <- TRUE
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            while (s <= 10 & border) {
              if (!is.null(seed)) {
                seed_num[z] <- (seed * s) + z
                set.seed(seed_num[z])
              }
              vector_of_add_QTN <-
                sample(index, add_QTN_num, replace = FALSE)
              genofile <- SNPRelate::snpgdsOpen(gdsfile)
              x <- 1
              sup_temp <- c()
              ld_between_QTNs_temp <- c()
              for (j in vector_of_add_QTN) {
                ldsup <- 1
                i <- j + 1
                while (ldsup > ld) {
                  if (i > n) {
                    if (verbose)
                      warning(
                        "There are no SNPs downstream. Selecting a different seed number",
                        call. = F,
                        immediate. = T
                      )
                    break
                  }
                  snp1 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, j),
                      count = c(-1, 1)
                    )
                  snp2 <-
                    gdsfmt::read.gdsn(
                      gdsfmt::index.gdsn(genofile, "genotype"),
                      start = c(1, i),
                      count = c(-1, 1)
                    )
                  ldsup <-
                    abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                  if (is.nan(ldsup)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i <- i + 1
                }
                sup_temp[x] <- i - 1
                ld_between_QTNs_temp[x] <- ldsup
                x <- x + 1
              }
              if (i > n | i2 < 1) {
                border <- TRUE
              } else {
                border <- FALSE
              }
              s <- s + 1
            }
            SNPRelate::snpgdsClose(genofile)
            sup[[z]] <- sup_temp
            inf[[z]] <- vector_of_add_QTN
            add_gen_info_inf[[z]] <-
              data.frame(
                snp_type = "QTN_for_trait_1",
                genotypes[vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup[[z]] <-
              data.frame(
                snp_type = "QTN_for_trait_2",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add[[z]] <- rbind(add_gen_info_inf[[z]],
                                      add_gen_info_sup[[z]])
            LD_summary_add[[z]] <- data.frame(
              z,
              ld,
              ld_between_QTNs_temp,
              add_gen_info_inf[[z]][, 2],
              add_gen_info_sup[[z]][, 2],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add[[z]]) <-
              c(
                "rep",
                "Aimed_LD (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
          }
          LD_summary_add <- do.call(rbind, LD_summary_add)
          data.table::fwrite(
            LD_summary_add,
            "LD_summary_Additive.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_add <- do.call(rbind, results_add)
          ns <- nrow(genotypes) - 5
          maf <- round(apply(results_add[, -c(1:6)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_add[, 2]
          results_add <-
            data.frame(
              results_add[, 1:6],
              maf = maf,
              results_add[, - c(1:6)],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results_add <-
            data.frame(
              rep = rep(1:rep, each = add_QTN_num * 2),
              results_add,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results_add <- results_add[, 1:8]
          }
          if (add_QTN) {
            write.table(
              seed_num,
              paste0("seed_num_for_", add_QTN_num,
                     "_Add_QTN",
                     ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            data.table::fwrite(
              results_add,
              "Additive_selected_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          add_ef_trait_obj <- mapply(function(x, y) {
            rownames(x) <-
              paste0("Chr_", x$chr, "_", x$pos)
            rownames(y) <-
              paste0("Chr_",  y$chr, "_", y$pos)
            b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
            return(b)
          },
          x = add_gen_info_sup,
          y = add_gen_info_inf,
          SIMPLIFY = F)
        }
        if (dom) {
          sup <- vector("list", rep)
          inf <- vector("list", rep)
          dom_gen_info_sup <- vector("list", rep)
          dom_gen_info_inf <- vector("list", rep)
          results_dom <- vector("list", rep)
          LD_summary_dom <- vector("list", rep)
          seed_num <- c()
          for (z in 1:rep) {
            s <- 1
            border <- TRUE
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            while (s <= 10 & border) {
              if (!is.null(seed)) {
                seed_num[z] <- (seed * s) + z + rep
                set.seed(seed_num[z])
              }
              vector_of_dom_QTN <-
                sample(index, dom_QTN_num, replace = FALSE)
              genofile <- SNPRelate::snpgdsOpen(gdsfile)
              x <- 1
              sup_temp <- c()
              ld_between_QTNs_temp <- c()
              again <- TRUE
              times <- 1
              dif <- c()
              while (again) {
                for (j in vector_of_dom_QTN) {
                  ldsup <- 1
                  i <- j + 1
                  while (ldsup > ld) {
                    if (i > n) {
                      if (verbose)
                        warning(
                          "There are no SNPs downstream. Selecting a different seed number",
                          call. = F,
                          immediate. = T
                        )
                      break
                    }
                    snp1 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, j),
                        count = c(-1, 1)
                      )
                    snp2 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, i),
                        count = c(-1, 1)
                      )
                    ldsup <-
                      abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = "composite"))
                    if (is.nan(ldsup)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    i <- i + 1
                  }
                  sup_temp[x] <- i - 1
                  ld_between_QTNs_temp[x] <- ldsup
                  x <- x + 1
                }
                if ((!any(genotypes[sup_temp, - (1:5)] == 0) |
                     !any(genotypes[vector_of_dom_QTN, - (1:5)] == 0)) &
                    times <= 10) {
                  if (!is.null(seed)) {
                    seed_num[z] <- (seed * s) + z + rep
                    set.seed(seed_num[z])
                  }
                  dif <- c(dif, vector_of_dom_QTN)
                  vector_of_dom_QTN <-
                    sample(setdiff(index, dif), dom_QTN_num, replace = FALSE)
                  again <- TRUE
                } else {
                  again <- FALSE
                }
                times <- times + 1
              }
              if (i > n | i2 < 1) {
                border <- TRUE
              } else {
                border <- FALSE
              }
              s <- s + 1
            }
            SNPRelate::snpgdsClose(genofile)
            sup[[z]] <- sup_temp
            inf[[z]] <- vector_of_dom_QTN
            dom_gen_info_inf[[z]] <-
              data.frame(
                snp_type = "QTN_for_trait_1",
                genotypes[vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup[[z]] <-
              data.frame(
                snp_type = "QTN_for_trait_2",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom[[z]] <- rbind(dom_gen_info_inf[[z]],
                                      dom_gen_info_sup[[z]])
            LD_summary_dom[[z]] <- data.frame(
              z,
              ld,
              ld_between_QTNs_temp,
              dom_gen_info_inf[[z]][, 2],
              dom_gen_info_sup[[z]][, 2],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom[[z]]) <-
              c(
                "rep",
                "Aimed_LD (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
          }
          LD_summary_dom <- do.call(rbind, LD_summary_dom)
          data.table::fwrite(
            LD_summary_dom,
            "LD_summary_Dominance.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_dom <- do.call(rbind, results_dom)
          ns <- nrow(genotypes) - 5
          maf <- round(apply(results_dom[, -c(1:6)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_dom[, 2]
          results_dom <-
            data.frame(
              results_dom[, 1:6],
              maf,
              results_dom[, - c(1:6)],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results_dom <-
            data.frame(
              rep = rep(1:rep, each = dom_QTN_num * 2),
              results_dom,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results_dom <- results_dom[, 1:8]
          }
          if (add_QTN) {
            write.table(
              seed_num,
              paste0("seed_num_for_", dom_QTN_num,
                     "_Dom_QTN",
                     ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            data.table::fwrite(
              results_dom,
              "Dominance_selected_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          dom_ef_trait_obj <- mapply(function(x, y) {
            rownames(x) <-
              paste0("Chr_", x$chr, "_", x$pos)
            rownames(y) <-
              paste0("Chr_", y$chr, "_", y$pos)
            b <- list(t(x[, - (1:6)]), t(y[, - (1:6)]))
            return(b)
          },
          x = dom_gen_info_sup,
          y = dom_gen_info_inf,
          SIMPLIFY = F)
        }
      }
      if (!is.null(add_ef_trait_obj) & !add_QTN) {
        add_ef_trait_obj <- lapply(add_ef_trait_obj, function(x) {
          lapply(x, function(y) {
            rnames <- rownames(y)
            y <- matrix(0, nrow = nrow(y), ncol =  1)
            rownames(y)  <- rnames
            return(y)
          })
        })
      }
      if (!is.null(dom_ef_trait_obj) & !dom_QTN) {
        dom_ef_trait_obj <- lapply(dom_ef_trait_obj, function(x) {
          lapply(x, function(y) {
            rnames <- rownames(y)
            y <- matrix(0, nrow = nrow(y), ncol =  1)
            rownames(y)  <- rnames
            return(y)
          })
        })
      }
      if (!is.null(add_ef_trait_obj)) {
        biallelic <- any(unlist(lapply(add_ef_trait_obj, function(x) {
          sapply(x, function(x2) {
            apply(x2, 2, function(y) {
              length(unique(y)) > 3
            })
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
          sapply(x, function(x2) {
            apply(x2, 2, function(y) {
              length(unique(y)) > 3
            })
          })
        }),
        recursive = T))
        if (biallelic) {
          stop("Please use only biallelic markers.",
               call. = F)
        }
      }
      return(list(
        add_ef_trait_obj = add_ef_trait_obj,
        dom_ef_trait_obj = dom_ef_trait_obj
      ))
    }
  }
