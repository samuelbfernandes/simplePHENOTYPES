#' Select SNPs to be assigned as QTNs
#' @keywords internal
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add_QTN_num = NULL,
#' @param dom_QTN_num = NULL,
#' @param same_add_dom_QTN = NULL,
#' @param dom = NULL,
#' @param add = NULL,
#' @param ld_min = NULL,
#' @param ld_max = NULL,
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param ld_method Four methods can be used to calculate linkage disequilibrium values: "composite" for LD composite measure, "r" for R coefficient (by EM algorithm assuming HWE, it could be negative), "dprime" for D', and "corr" for correlation coefficient.
#' @param gdsfile NULL
#' @param constraints = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @param type_of_ld = NULL
#' @param verbose = verbose
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes. Last update: Apr 20, 2020
#'
qtn_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           ld_max = NULL,
           ld_min = NULL,
           ld_method = "composite",
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
    if (!requireNamespace("SNPRelate", quietly = TRUE) ||
        !requireNamespace("gdsfmt", quietly = TRUE)) {
      stop(.gds_needed("LD-architecture"), call. = FALSE)
    }
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
              seed_num[z] <- (seed * s) + z
              set.seed(seed_num[z])
            vector_of_add_QTN <-
              sample(index, add_QTN_num, replace = FALSE)
            sup_temp <- c()
            inf_temp <- c()
            ld_between_QTNs_temp <- c()
            actual_ld_sup <- c()
            actual_ld_inf <- c()
            x <- 1
            for (j in vector_of_add_QTN) {
              times <- 1
              dif <- c()
              ldsup <- 1
              ldinf <- 1
              while (times <= 100 & (ldsup < ld_min | ldinf < ld_min | ldsup > ld_max | ldinf > ld_max)) {
                i <- j + 1
                i2 <- j - 1
                while (ldsup > ld_max) {
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
                    abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                  if (is.nan(ldsup)) {
                    SNPRelate::snpgdsClose(genofile)
                    stop("Monomorphic SNPs are not accepted", call. = F)
                  }
                  i <- i + 1
                }
                actual_ld_sup[x] <- ldsup
                sup_temp[x] <- i - 1
                while (ldinf > ld_max) {
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
                    abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
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
                  SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = ld_method)[1]
                if (((!any(genotypes[sup_temp, - (1:5)] == 0) |
                      !any(genotypes[inf_temp, - (1:5)] == 0)) & dom) |
                    ldsup > ld_max | ldinf > ld_max | ldsup < ld_min | ldinf < ld_min) {
                    seed_num[z] <- (seed * s) + z
                    set.seed(seed_num[z])
                  dif <- c(dif, vector_of_add_QTN, j)
                  j <-
                    sample(setdiff(index, dif), 1, replace = FALSE)
                  ldsup <- 1
                  ldinf <- 1
                }
                times <- times + 1
              }
              vector_of_add_QTN[x] <- j
              x <- x + 1
              if (((!any(genotypes[sup_temp, - (1:5)] == 0) |
                    !any(genotypes[inf_temp, - (1:5)] == 0)) &
                   dom)) {
                warning(
                  "No heterozygote found to simulate dominance.",
                  call. = F,
                  immediate. = T
                )
              }
              if (s == 10 & (ldsup < ld_min | ldinf < ld_min)) {
                stop(
                  "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
                  call. = F
                )
              }
              if (s == 10 & (ldsup > ld_max | ldinf > ld_max)) {
                stop(
                  "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
                  call. = F
                )
              }
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
          .ld_check_indirect(genotypes, vector_of_add_QTN, sup_temp, inf_temp, z,
                            ld_method = ld_method, ld_min = ld_min, ld_max = ld_max,
                            reported_sup = actual_ld_sup, reported_inf = actual_ld_inf)
          QTN_causing_ld[[z]] <-
            data.frame(
              type = "cause_of_LD",
              trait = "none",
              genotypes[vector_of_add_QTN, ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup[[z]] <-
            data.frame(
              type = "QTN_downstream",
              trait = "trait_1",
              genotypes[sup[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_inf[[z]] <-
            data.frame(
              type = "QTN_upstream",
              trait = "trait_2",
              genotypes[inf[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                add_gen_info_sup[[z]],
                                add_gen_info_inf[[z]])
          LD_summary[[z]] <- data.frame(
            z,
            QTN_causing_ld[[z]][, "snp"],
            ld_min,
            ld_max,
            actual_ld_sup,
            actual_ld_inf,
            add_gen_info_sup[[z]][, "snp"],
            add_gen_info_inf[[z]][, "snp"],
            ld_between_QTNs_temp,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary[[z]]) <-
            c(
              "rep",
              "SNP_causing_LD",
              "ld_min (absolute value)",
              "ld_max (absolute value)",
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
          "LD_Summary.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        results <- do.call(rbind, results)
        ns <- ncol(genotypes) - 5
        maf <- round(apply(results[, - c(1:7)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx, (1 - sumx))
        }), 4)
        names(maf) <- results[, "snp"]
        results <- data.frame(
          results[, 1:2],
          additive_effect = c(rep("-", add_QTN_num), unlist(add_effect)),
          dominance_effect = c(rep("-", add_QTN_num), unlist(dom_effect)),
          results[, 3:7],
          maf = maf,
          results[, - c(1:7)],
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
          results <- results[, 1:11]
        }
        if (add_QTN) {
          if (verbose){
            write.table(
              seed_num,
              paste0("Seed_num_for_", add_QTN_num,
                     "_Add_QTN.txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          }
          data.table::fwrite(
            results,
            "Additive_QTNs.txt",
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
          b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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
                seed_num[z] <-  (seed * s) + z
                set.seed(seed_num[z])
              vector_of_add_QTN <-
                sample(index, add_QTN_num, replace = FALSE)
              sup_temp <- c()
              inf_temp <- c()
              ld_between_QTNs_temp <- c()
              actual_ld_sup <- c()
              actual_ld_inf <- c()
              x <- 1
              for (j in vector_of_add_QTN) {
                times <- 1
                dif <- c()
                ldsup <- 1
                ldinf <- 1
                while (times <= 100 & (ldsup < ld_min | ldinf < ld_min | ldsup > ld_max | ldinf > ld_max)) {
                  i <- j + 1
                  i2 <- j - 1
                  while (ldsup > ld_max) {
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
                      abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                    if (is.nan(ldsup)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    i <- i + 1
                  }
                  actual_ld_sup[x] <- ldsup
                  sup_temp[x] <- i - 1
                  while (ldinf > ld_max) {
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
                      abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
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
                    SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = ld_method)[1]
                  if (ldsup > ld_max | ldinf > ld_max | ldsup < ld_min | ldinf < ld_min) {
                      seed_num[z] <- (seed * s) + z
                      set.seed(seed_num[z])
                    dif <- c(dif, vector_of_add_QTN, j)
                    j <-
                      sample(setdiff(index, dif), 1, replace = FALSE)
                    ldsup <- 1
                    ldinf <- 1
                  }
                  times <- times + 1
                }
                vector_of_add_QTN[x] <- j
                x <- x + 1
                if (s == 10 & (ldsup < ld_min | ldinf < ld_min)) {
                  stop(
                    "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
                    call. = F
                  )
                }
                if (s == 10 & (ldsup > ld_max | ldinf > ld_max)) {
                  stop(
                    "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
                    call. = F
                  )
                }
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
            .ld_check_indirect(genotypes, vector_of_add_QTN, sup_temp, inf_temp, z,
                            ld_method = ld_method, ld_min = ld_min, ld_max = ld_max,
                            reported_sup = actual_ld_sup, reported_inf = actual_ld_inf)
            QTN_causing_ld[[z]] <-
              data.frame(
                type = "cause_of_LD",
                trait = "none",
                genotypes[vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup[[z]] <-
              data.frame(
                type = "QTN_downstream",
                trait = "trait_1",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_inf[[z]] <-
              data.frame(
                type = "QTN_upstream",
                trait = "trait_2",
                genotypes[inf[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add[[z]] <- rbind(QTN_causing_ld[[z]],
                                      add_gen_info_sup[[z]],
                                      add_gen_info_inf[[z]])
            LD_summary_add[[z]] <- data.frame(
              z,
              QTN_causing_ld[[z]][, "snp"],
              ld_min,
              ld_max,
              actual_ld_sup,
              actual_ld_inf,
              add_gen_info_sup[[z]][, "snp"],
              add_gen_info_inf[[z]][, "snp"],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add[[z]]) <-
              c(
                "rep",
                "SNP_causing_LD",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
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
            "LD_Summary_Additive.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_add <- do.call(rbind, results_add)
          ns <- ncol(genotypes) - 5
          maf <- round(apply(results_add[, -c(1:7)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_add[, "snp"]
          results_add <-
            data.frame(
              results_add[, 1:2],
              additive_effect = c(rep("-", add_QTN_num), unlist(add_effect)),
              results_add[, 3:7],
              maf = maf,
              results_add[, -c(1:7)],
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
            results_add <- results_add[, 1:10]
          }
          if (add_QTN) {
            if (verbose){
              write.table(
                seed_num,
                paste0("Seed_num_for_", add_QTN_num,
                       "_Add_QTN.txt"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE
              )
            }
            data.table::fwrite(
              results_add,
              "Additive_QTNs.txt",
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
            b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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
                seed_num[z] <- (seed * s) + z + rep
                set.seed(seed_num[z])
              vector_of_dom_QTN <-
                sample(index, dom_QTN_num, replace = FALSE)
              sup_temp <- c()
              inf_temp <- c()
              ld_between_QTNs_temp <- c()
              actual_ld_sup <- c()
              actual_ld_inf <- c()
              x <- 1
              for (j in vector_of_dom_QTN) {
                times <- 1
                dif <- c()
                ldsup <- 1
                ldinf <- 1
                while (times <= 100 & (ldsup < ld_min | ldinf < ld_min | ldsup > ld_max | ldinf > ld_max)) {
                  i <- j + 1
                  i2 <- j - 1
                  while (ldsup > ld_max) {
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
                      abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                    if (is.nan(ldsup)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    i <- i + 1
                  }
                  actual_ld_sup[x] <- ldsup
                  sup_temp[x] <- i - 1
                  while (ldinf > ld_max) {
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
                      abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
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
                    SNPRelate::snpgdsLDpair(snp_sup, snp_inf, method = ld_method)[1]
                  if (((!any(genotypes[sup_temp, - (1:5)] == 0) |
                        !any(genotypes[inf_temp, - (1:5)] == 0)) & dom) |
                      ldsup > ld_max | ldinf > ld_max | ldsup < ld_min | ldinf < ld_min) {
                    SNPRelate::snpgdsClose(genofile)
                    stop(
                      "Indirect LD with dominance QTNs had to re-sample an intermediate marker (no heterozygote among the selected QTNs, or LD outside [ld_min, ld_max]); ",
                      "this path is not supported (it would overwrite the additive marker vector and report wrong intermediate markers). ",
                      "Try a different `seed`, a different LD window, `constraints = list(hets = \"include\")`, or `type_of_ld = \"direct\"`.",
                      call. = F
                    )
                    dif <- c(dif, vector_of_add_QTN, j)
                    j <-
                      sample(setdiff(index, dif), 1, replace = FALSE)
                    ldsup <- 1
                    ldinf <- 1
                  }
                  times <- times + 1
                }
                vector_of_add_QTN[x] <- j
                x <- x + 1
                if (((!any(genotypes[sup_temp, - (1:5)] == 0) |
                      !any(genotypes[inf_temp, - (1:5)] == 0)) &
                     dom)) {
                  warning(
                    "No heterozygote found to simulate dominance.",
                    call. = F,
                    immediate. = T
                  )
                }
                if (s == 10 & (ldsup < ld_min | ldinf < ld_min)) {
                  stop(
                    "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
                    call. = F
                  )
                }
                if (s == 10 & (ldsup > ld_max | ldinf > ld_max)) {
                  stop(
                    "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
                    call. = F
                  )
                }
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
            .ld_check_indirect(genotypes, vector_of_dom_QTN, sup_temp, inf_temp, z,
                            ld_method = ld_method, ld_min = ld_min, ld_max = ld_max,
                            reported_sup = actual_ld_sup, reported_inf = actual_ld_inf)
            QTN_causing_ld[[z]] <-
              data.frame(
                type = "cause_of_LD",
                trait = "none",
                genotypes[vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup[[z]] <-
              data.frame(
                type = "QTN_downstream",
                trait = "trait_1",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_inf[[z]] <-
              data.frame(
                type = "QTN_upstream",
                trait = "trait_2",
                genotypes[inf[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom[[z]] <- rbind(QTN_causing_ld[[z]],
                                      dom_gen_info_sup[[z]],
                                      dom_gen_info_inf[[z]])
            LD_summary_dom[[z]] <- data.frame(
              z,
              QTN_causing_ld[[z]][, "snp"],
              ld_min,
              ld_max,
              actual_ld_sup,
              actual_ld_inf,
              dom_gen_info_sup[[z]][, "snp"],
              dom_gen_info_inf[[z]][, "snp"],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom[[z]]) <-
              c(
                "rep",
                "SNP_causing_LD",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
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
            "LD_Summary_Dominance.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_dom <- do.call(rbind, results_dom)
          ns <- ncol(genotypes) - 5
          maf <- round(apply(results_dom[, -c(1:7)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_dom[, "snp"]
          results_dom <-
            data.frame(
              results_dom[, 1:2],
              dominance_effect = c(rep("-", dom_QTN_num), unlist(dom_effect)),
              results_dom[, 3:7],
              maf = maf,
              results_dom[, - c(1:7)],
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
            results_dom <- results_dom[, 1:10]
          }
          if (dom_QTN) {
            if (verbose){
              write.table(
                seed_num,
                paste0("Seed_num_for_", dom_QTN_num,
                       "_Dom_QTN.txt"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE
              )
            }
            data.table::fwrite(
              results_dom,
              "Dominance_QTNs.txt",
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
            b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          attempt <- 1L
          repeat {
            seed_a <- .ld_attempt_seed(seed, attempt)
            repaired <- attempt > 1L
            s <- 1
            border <- TRUE
            run <- .ld_run_attempt({
              while (s <= 10 & border) {
                  seed_num[z] <- (seed_a * s) + z
                  set.seed(seed_num[z])
                vector_of_add_QTN <-
                  sample(index, add_QTN_num, replace = FALSE)
                sup_temp <- c()
                anchor_temp <- c()
                ld_between_QTNs_temp <- c()
                x <- 1
                for (j in vector_of_add_QTN) {
                  times <- 1
                  dif <- c()
                  ldsup <- 1
                  ldinf <- 1
                  i <- j + 1
                  i2 <- j - 1
                  while (times <= 100 & (ldsup > ld_max | ldsup < ld_min) &
                         (ldinf > ld_max | ldinf < ld_min)) {
                    if (i > n & i2 < 1) {
                      if (i > n) {
                        if (verbose)
                          warning(
                            "There are no SNPs downstream. Selecting a different seed number",
                            call. = F,
                            immediate. = T
                          )
                        break
                      }
                      if (i2 < 1) {
                        if (verbose)
                          warning(
                            "There are no SNPs upstream. Selecting a different seed number",
                            call. = F,
                            immediate. = T
                          )
                        break
                      }
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
                    snp3 <-
                      gdsfmt::read.gdsn(
                        gdsfmt::index.gdsn(genofile, "genotype"),
                        start = c(1, i2),
                        count = c(-1, 1)
                      )
                    ldsup <-
                      abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                    ldinf <-
                      abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
                    if (is.nan(ldinf) | is.nan(ldsup)) {
                      SNPRelate::snpgdsClose(genofile)
                      stop("Monomorphic SNPs are not accepted", call. = F)
                    }
                    if ((!any(genotypes[i, - (1:5)] == 0) &
                         !any(genotypes[i2, - (1:5)] == 0) & dom) |
                        (ldsup > ld_max | ldsup < ld_min) &
                        (ldinf > ld_max | ldinf < ld_min)) {
                        seed_num[z] <- (seed_a * s) + z + (if (repaired) x else 0)
                        set.seed(seed_num[z])
                      dif <- c(dif, vector_of_add_QTN, j)
                      j <-
                        sample(setdiff(index, dif), 1, replace = FALSE)
                      ldsup <- 1
                      ldinf <- 1
                      if (repaired) {
                        i <- j
                        i2 <- j
                      }
                    }
                    i <- i + 1
                    i2 <- i2 - 1
                    times <- times + 1
                  }
                  if (ldsup > ld_max) {
                    ldsup <- 100
                  } else if (ldinf > ld_max) {
                    ldinf <- 100
                  } else if (ldsup < ld_min) {
                    ldsup <- 100
                  } else if (ldinf < ld_min) {
                    ldinf <- 100
                  }
                  closest <- which.min(abs(c(ld_max - ldinf,  ld_max -ldsup)))
                  ld_between_QTNs_temp[x] <- 
                    ifelse(closest == 1 , ldinf, ldsup)
                  sup_temp[x] <-
                    ifelse(closest == 1 , i2 + 1, i - 1)
                  anchor_temp[x] <- j
                  x <- x + 1
                  if ((!any(genotypes[i - 1, - (1:5)] == 0) &
                       !any(genotypes[i2 + 1, - (1:5)] == 0) & dom)) {
                    warning(
                      "No heterozygote found to simulate dominance.",
                      call. = F,
                      immediate. = T
                    )
                  }
                  if (s == 10 & (ldsup < ld_min & ldinf < ld_min)) {
                    .ld_search_stop(
                      "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!"
                    )
                  }
                  if (s == 10 & (ldsup > ld_max & ldinf > ld_max)) {
                    .ld_search_stop(
                      "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!"
                    )
                  }
                }
                if (i > n | i < 1) {
                  border <- TRUE
                } else {
                  border <- FALSE
                }
                s <- s + 1
              }
            }, quiet = attempt > 1L)
            failure <- run$failure
            if (is.null(failure)) {
              sup[[z]] <- sup_temp
              inf[[z]] <- if (attempt == 1L) vector_of_add_QTN else anchor_temp
              failure <- .ld_direct_violation(genotypes, inf[[z]], sup_temp,
                                              ld_between_QTNs_temp, ld_method,
                                              ld_min, ld_max)
            }
            if (is.null(failure)) {
              .ld_emit_warnings(run$warnings)
              if (verbose && attempt > 1L) {
                message("* LD search, replicate ", z, ": contract met on attempt ", attempt,
                        " (derived seed ", format(seed_a, scientific = FALSE), ")")
              }
              break
            }
            if (.ld_give_up_now(failure, attempt)) {
              SNPRelate::snpgdsClose(genofile)
              .ld_search_giveup(failure, z, attempt)
            }
            attempt <- attempt + 1L
          }
          SNPRelate::snpgdsClose(genofile)
          add_gen_info_inf[[z]] <-
            data.frame(
              type = "QTN_selected",
              trait = "trait_2",
              genotypes[inf[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup[[z]] <-
            data.frame(
              type = "QTN_in_LD",
              trait = "trait_1",
              genotypes[sup[[z]], ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results[[z]] <- rbind(add_gen_info_inf[[z]],
                                add_gen_info_sup[[z]])
          LD_summary[[z]] <- data.frame(
            z,
            ld_min,
            ld_max,
            ld_between_QTNs_temp,
            add_gen_info_sup[[z]][, "snp"],
            add_gen_info_inf[[z]][, "snp"],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary[[z]]) <-
            c(
              "rep",
              "ld_min (absolute value)",
              "ld_max (absolute value)",
              "Actual_LD ",
              "QTN_for_trait_1",
              "QTN_for_trait_2"
            )
        }
        LD_summary <- do.call(rbind, LD_summary)
        data.table::fwrite(
          LD_summary,
          "LD_Summary.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        results <- do.call(rbind, results)
        ns <- ncol(genotypes) - 5
        maf <- round(apply(results[, -c(1:7)], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        names(maf) <- results[, "snp"]
        results <-
          data.frame(
            results[, 1:2],
            additive_effect = unlist(add_effect[2:1]),
            dominance_effect = unlist(dom_effect[2:1]),
            results[, 3:7],
            maf = maf,
            results[, - c(1:7)],
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
          results <- results[, 1:11]
        }
        if (add_QTN) {
          if (verbose){
            write.table(
              seed_num,
              paste0(
                "Seed_num_for_",
                add_QTN_num,
                "_Add_and_Dom_QTN",
                ".txt"
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          }
          data.table::fwrite(
            results,
            "Additive_QTNs.txt",
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
          b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            attempt <- 1L
            repeat {
              seed_a <- .ld_attempt_seed(seed, attempt)
              repaired <- attempt > 1L
              s <- 1
              border <- TRUE
              run <- .ld_run_attempt({
                while (s <= 10 & border) {
                    seed_num[z] <- (seed_a * s) + z
                    set.seed(seed_num[z])
                  vector_of_add_QTN <-
                    sample(index, add_QTN_num, replace = FALSE)
                  x <- 1
                  sup_temp <- c()
                  ld_between_QTNs_temp <- c()
                  dif <- c()
                  new_vector_of_add_QTN <- c()
                  for (j in vector_of_add_QTN) {
                    times <- 1
                    ldsup <- 1
                    ldinf <- 1
                    i <- j + 1
                    i2 <- j - 1
                    while (times <= length(index)/2 & (ldsup > ld_max | ldsup < ld_min) &
                           (ldinf > ld_max | ldinf < ld_min)) {
                      if (i > n | i2 < 1) {
                        warning(
                          "Trying to find SNPs that match the \'ld_max\' and \'ld_min\' criteria.",
                          call. = F,
                          immediate. = T
                        )
                        seed_num[z] <- (seed_a * s) + z + x
                        set.seed(seed_num[z])
                        j <-
                          sample(setdiff(index, dif), 1, replace = FALSE)
                        ldsup <- 1
                        ldinf <- 1
                        i <- j + 1
                        i2 <- j - 1
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
                      snp3 <-
                        gdsfmt::read.gdsn(
                          gdsfmt::index.gdsn(genofile, "genotype"),
                          start = c(1, i2),
                          count = c(-1, 1)
                        )
                      ldsup <-
                        abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                      ldinf <-
                        abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
                      if (is.nan(ldinf) | is.nan(ldsup)) {
                        SNPRelate::snpgdsClose(genofile)
                        stop("Monomorphic SNPs are not accepted", call. = F)
                      }
                      if ((ldsup > ld_max | ldsup < ld_min) &
                          (ldinf > ld_max | ldinf < ld_min)) {
                        seed_num[z] <- (seed_a * s) + z + x
                        set.seed(seed_num[z])
                        j <-
                          sample(setdiff(index, dif), 1, replace = FALSE)
                        ldsup <- 1
                        ldinf <- 1
                        i <- j + 1
                        i2 <- j - 1
                      } else {
                        i <- i + 1
                        i2 <- i2 - 1 
                      }
                      dif <- c(dif, j)                  
                      times <- times + 1
                    }
                    if (ldsup > ld_max) {
                      ldsup <- 100
                    } else if (ldinf > ld_max) {
                      ldinf <- 100
                    } else if (ldsup < ld_min) {
                      ldsup <- 100
                    } else if (ldinf < ld_min) {
                      ldinf <- 100
                    }
                    closest <- which.min(abs(c(ld_max - ldinf,  ld_max -ldsup)))
                    ld_between_QTNs_temp[x] <- 
                      ifelse(closest == 1 , ldinf, ldsup)
                    sup_temp[x] <-
                      ifelse(closest == 1 , i2 + 1, i - 1)
                    new_vector_of_add_QTN[x] <- j
                    x <- x + 1
                    if (s == 10 & (ldsup < ld_min & ldinf < ld_min)) {
                      .ld_search_stop(
                        "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!"
                      )
                    }
                    if (s == 10 & (ldsup > ld_max & ldinf > ld_max)) {
                      .ld_search_stop(
                        "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!"
                      )
                    }
                  }
                  if (i > n | i < 1) {
                    border <- TRUE
                  } else {
                    border <- FALSE
                  }
                  s <- s + 1
                }
              }, quiet = attempt > 1L)
              failure <- run$failure
              if (is.null(failure)) {
                sup[[z]] <- sup_temp
                inf[[z]] <- new_vector_of_add_QTN
                failure <- .ld_direct_violation(genotypes, inf[[z]], sup_temp,
                                                ld_between_QTNs_temp, ld_method,
                                                ld_min, ld_max)
              }
              if (is.null(failure)) {
                .ld_emit_warnings(run$warnings)
                if (verbose && attempt > 1L) {
                  message("* LD search, replicate ", z, ": contract met on attempt ", attempt,
                          " (derived seed ", format(seed_a, scientific = FALSE), ")")
                }
                break
              }
              if (.ld_give_up_now(failure, attempt)) {
                SNPRelate::snpgdsClose(genofile)
                .ld_search_giveup(failure, z, attempt)
              }
              attempt <- attempt + 1L
            }
            SNPRelate::snpgdsClose(genofile)
            add_gen_info_inf[[z]] <-
              data.frame(
                type = "QTN_selected",
                trait = "trait_2",
                genotypes[new_vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup[[z]] <-
              data.frame(
                type = "QTN_in_LD",
                trait = "trait_1",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add[[z]] <- rbind(add_gen_info_inf[[z]],
                                      add_gen_info_sup[[z]])
            LD_summary_add[[z]] <- data.frame(
              z,
              ld_min,
              ld_max,
              ld_between_QTNs_temp,
              add_gen_info_sup[[z]][, "snp"],
              add_gen_info_inf[[z]][, "snp"],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add[[z]]) <-
              c(
                "rep",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
          }
          LD_summary_add <- do.call(rbind, LD_summary_add)
          data.table::fwrite(
            LD_summary_add,
            "LD_Summary_Additive.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_add <- do.call(rbind, results_add)
          ns <- ncol(genotypes) - 5
          maf <- round(apply(results_add[, -c(1:7)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_add[, "snp"]
          results_add <-
            data.frame(
              results_add[, 1:2],
              additive_effect = unlist(add_effect[2:1]),
              results_add[, 3:7],
              maf = maf,
              results_add[, - c(1:7)],
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
            results_add <- results_add[, 1:10]
          }
          if (add_QTN) {
            if (verbose){
              write.table(
                seed_num,
                paste0("Seed_num_for_", add_QTN_num,
                       "_Add_QTN",
                       ".txt"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE
              )
            }
            data.table::fwrite(
              results_add,
              "Additive_QTNs.txt",
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
            b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            attempt <- 1L
            repeat {
              seed_a <- .ld_attempt_seed(seed, attempt)
              repaired <- attempt > 1L
              s <- 1
              border <- TRUE
              run <- .ld_run_attempt({
                while (s <= 10 & border) {
                    seed_num[z] <- (seed_a * s) + z + rep
                    set.seed(seed_num[z])
                  vector_of_dom_QTN <-
                    sample(index, dom_QTN_num, replace = FALSE)
                  sup_temp <- c()
                  ld_between_QTNs_temp <- c()
                  x <- 1
                  ld_between_QTNs_temp <- c()
                  dif <- c()
                  new_vector_of_dom_QTN <- c()
                  for (j in vector_of_dom_QTN) {
                    times <- 1
                    ldsup <- 1
                    ldinf <- 1
                    i <- j + 1
                    i2 <- j - 1
                    while (times <= length(index)/2 & (ldsup > ld_max | ldsup < ld_min) &
                           (ldinf > ld_max | ldinf < ld_min)) {
                      if (i > n | i2 < 1) {
                        warning(
                          "Trying to find SNPs that match the \'ld_max\' and \'ld_min\' criteria.",
                          call. = F,
                          immediate. = T
                        )
                        seed_num[z] <- (seed_a * s) + z + rep + x
                        set.seed(seed_num[z])
                        j <-
                          sample(setdiff(index, dif), 1, replace = FALSE)
                        ldsup <- 1
                        ldinf <- 1
                        i <- j + 1
                        i2 <- j - 1
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
                      snp3 <-
                        gdsfmt::read.gdsn(
                          gdsfmt::index.gdsn(genofile, "genotype"),
                          start = c(1, i2),
                          count = c(-1, 1)
                        )
                      ldsup <-
                        abs(SNPRelate::snpgdsLDpair(snp1, snp2, method = ld_method))[1]
                      ldinf <-
                        abs(SNPRelate::snpgdsLDpair(snp1, snp3, method = ld_method))[1]
                      if (is.nan(ldinf) | is.nan(ldsup)) {
                        SNPRelate::snpgdsClose(genofile)
                        stop("Monomorphic SNPs are not accepted", call. = F)
                      }
                      if ((!any(genotypes[i, - (1:5)] == 0) &
                           !any(genotypes[i2, - (1:5)] == 0) & dom) |
                          (ldsup > ld_max | ldsup < ld_min) &
                          (ldinf > ld_max | ldinf < ld_min)) {
                          seed_num[z] <- (seed_a * s) + z + rep + x
                          set.seed(seed_num[z])
                        j <-
                          sample(setdiff(index, dif), 1, replace = FALSE)
                        ldsup <- 1
                        ldinf <- 1
                        if (repaired) {
                          i <- j + 1
                          i2 <- j - 1
                        } else {
                          i <- i + 1
                          i2 <- i2 - 1
                        }
                      } else {
                        i <- i + 1
                        i2 <- i2 - 1
                      }
                      dif <- c(dif, j)      
                      times <- times + 1
                    }
                    if (ldsup > ld_max) {
                      ldsup <- 100
                    } else if (ldinf > ld_max) {
                      ldinf <- 100
                    } else if (ldsup < ld_min) {
                      ldsup <- 100
                    } else if (ldinf < ld_min) {
                      ldinf <- 100
                    }
                    closest <- which.min(abs(c(ld_max - ldinf,  ld_max -ldsup)))
                    ld_between_QTNs_temp[x] <- 
                      ifelse(closest == 1 , ldinf, ldsup)
                    sup_temp[x] <-
                      ifelse(closest == 1 , i2 + 1, i - 1)
                    new_vector_of_dom_QTN[x] <- j
                    x <- x + 1
                    if ((!any(genotypes[i - 1, - (1:5)] == 0) &
                         !any(genotypes[i2 + 1, - (1:5)] == 0) & dom)) {
                      warning(
                        "No heterozygote found to simulate dominance.",
                        call. = F,
                        immediate. = T
                      )
                    }
                    if (s == 10 & (ldsup < ld_min & ldinf < ld_min)) {
                      .ld_search_stop(
                        "None of the selected SNPs met the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!"
                      )
                    }
                    if (s == 10 & (ldsup > ld_max & ldinf > ld_max)) {
                      .ld_search_stop(
                        "None of the selected SNPs met the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!"
                      )
                    }
                  }
                  if (i > n | i < 1) {
                    border <- TRUE
                  } else {
                    border <- FALSE
                  }
                  s <- s + 1
                }
              }, quiet = attempt > 1L)
              failure <- run$failure
              if (is.null(failure)) {
                sup[[z]] <- sup_temp
                inf[[z]] <- new_vector_of_dom_QTN
                failure <- .ld_direct_violation(genotypes, inf[[z]], sup_temp,
                                                ld_between_QTNs_temp, ld_method,
                                                ld_min, ld_max)
              }
              if (is.null(failure)) {
                .ld_emit_warnings(run$warnings)
                if (verbose && attempt > 1L) {
                  message("* LD search, replicate ", z, ": contract met on attempt ", attempt,
                          " (derived seed ", format(seed_a, scientific = FALSE), ")")
                }
                break
              }
              if (.ld_give_up_now(failure, attempt)) {
                SNPRelate::snpgdsClose(genofile)
                .ld_search_giveup(failure, z, attempt)
              }
              attempt <- attempt + 1L
            }
            SNPRelate::snpgdsClose(genofile)
            dom_gen_info_inf[[z]] <-
              data.frame(
                type = "QTN_selected",
                trait = "trait_2",
                genotypes[new_vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup[[z]] <-
              data.frame(
                type = "QTN_in_LD",
                trait = "trait_1",
                genotypes[sup[[z]], ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom[[z]] <- rbind(dom_gen_info_inf[[z]],
                                      dom_gen_info_sup[[z]])
            LD_summary_dom[[z]] <- data.frame(
              z,
              ld_min,
              ld_max,
              ld_between_QTNs_temp,
              dom_gen_info_sup[[z]][, "snp"],
              dom_gen_info_inf[[z]][, "snp"],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom[[z]]) <-
              c(
                "rep",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
          }
          LD_summary_dom <- do.call(rbind, LD_summary_dom)
          data.table::fwrite(
            LD_summary_dom,
            "LD_Summary_Dominance.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          results_dom <- do.call(rbind, results_dom)
          ns <- ncol(genotypes) - 5
          maf <- round(apply(results_dom[, -c(1:7)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results_dom[, "snp"]
          results_dom <-
            data.frame(
              results_dom[, 1:2],
              dominance_effect = unlist(dom_effect[2:1]),
              results_dom[, 3:7],
              maf,
              results_dom[, - c(1:7)],
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
            results_dom <- results_dom[, 1:10]
          }
          if (dom_QTN) {
            if (verbose){
              write.table(
                seed_num,
                paste0("Seed_num_for_", dom_QTN_num,
                       "_Dom_QTN",
                       ".txt"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE
              )
            }
            data.table::fwrite(
              results_dom,
              "Dominance_QTNs.txt",
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
            b <- list(t(x[, - (1:7)]), t(y[, - (1:7)]))
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

# ---------------------------------------------------------------------------
# LD-contract checks (audit v1-core, D1). The frozen selection walks are not
# altered; the pairs they return are verified and rejected with an informative
# message when the "distinct causal markers in LD" contract is not met.
# ---------------------------------------------------------------------------

#' Absolute LD between two marker rows of a numeric genotype frame
#' @keywords internal
#' @noRd
.ld_pair <- function(genotypes, a, b, ld_method) {
  g1 <- as.numeric(unlist(genotypes[a, -(1:5)])) + 1
  g2 <- as.numeric(unlist(genotypes[b, -(1:5)])) + 1
  abs(SNPRelate::snpgdsLDpair(g1, g2, method = ld_method))[1]
}

#' Marker key used to name genotype rows (Chr_<chr>_<pos>)
#' @keywords internal
#' @noRd
.ld_key <- function(genotypes, idx) {
  paste0("Chr_", genotypes$chr[idx], "_", genotypes$pos[idx])
}

#' Stop with the LD-contract message
#'
#' `attempts` > 1 states that the direct-LD search was repeated with derived
#' seeds (see `.ld_attempt_seed()`) before giving up.
#' @keywords internal
#' @noRd
.ld_contract_stop <- function(reason, z, type_of_ld, attempts = 1L) {
  stop(
    "The LD contract could not be met for this seed and LD window (",
    type_of_ld, " LD, replicate ", z, "): ", reason, ". ",
    if (attempts > 1L) {
      paste0("The marker search was repeated ", attempts - 1L,
             " more time(s) with derived seeds and still failed. ")
    } else "",
    "Try a different `seed`, a different LD window (`ld_min`/`ld_max`), ",
    if (type_of_ld == "indirect") "`type_of_ld = \"direct\"`, " else "",
    "or a denser marker set.",
    call. = FALSE
  )
}

#' Verify indirect-LD selections: one distinct marker per trait and per pair
#'
#' Structural checks (distinct, non-shared QTNs on the intermediate marker's
#' chromosome) always run. When `ld_method`, `ld_min` and `ld_max` are given
#' the absolute LD of every (cause, downstream) and (cause, upstream) pair is
#' recomputed and must lie in the inclusive window `[ld_min, ld_max]`, exactly
#' as `.ld_check_direct()` does for direct LD; when `reported_sup` /
#' `reported_inf` are given they must equal the recomputed LD. A cause marker
#' that coincides with another triple's QTN is deliberately not rejected.
#' @keywords internal
#' @noRd
.ld_check_indirect <- function(genotypes, cause, sup, inf, z,
                               ld_method = NULL, ld_min = NULL, ld_max = NULL,
                               reported_sup = NULL, reported_inf = NULL) {
  if (anyDuplicated(sup) || anyDuplicated(inf) ||
      anyDuplicated(.ld_key(genotypes, sup)) ||
      anyDuplicated(.ld_key(genotypes, inf))) {
    .ld_contract_stop("two intermediate markers resolved to the same neighbouring marker (or to markers sharing a chromosome position), so a trait would have a duplicated QTN",
                      z, "indirect")
  }
  if (length(intersect(sup, inf))) {
    .ld_contract_stop("the same marker was selected as a QTN for both traits",
                      z, "indirect")
  }
  if (any(cause == sup | cause == inf)) {
    .ld_contract_stop("an intermediate marker coincides with one of its own QTNs",
                      z, "indirect")
  }
  chr_cause <- genotypes$chr[cause]
  if (any(genotypes$chr[sup] != chr_cause | genotypes$chr[inf] != chr_cause)) {
    .ld_contract_stop("a QTN lies on a different chromosome than its intermediate marker",
                      z, "indirect")
  }
  if (!is.null(ld_method) && !is.null(ld_min) && !is.null(ld_max)) {
    if (length(sup) != length(cause) || length(inf) != length(cause)) {
      .ld_contract_stop("the numbers of intermediate and linked markers differ",
                        z, "indirect")
    }
    pair_ld <- function(other) {
      vapply(seq_along(cause), function(k) {
        .ld_pair(genotypes, cause[k], other[k], ld_method)
      }, numeric(1))
    }
    ld_sup <- pair_ld(sup)
    ld_inf <- pair_ld(inf)
    tol <- 1e-9
    out <- function(l) anyNA(l) || any(l < ld_min - tol | l > ld_max + tol)
    if (out(ld_sup) || out(ld_inf)) {
      .ld_contract_stop("a selected marker has an absolute LD with its intermediate marker outside [ld_min, ld_max]",
                        z, "indirect")
    }
    differs <- function(l, r) {
      !is.null(r) && (length(r) != length(l) || any(abs(l - r) > 1e-6))
    }
    if (differs(ld_sup, reported_sup) || differs(ld_inf, reported_inf)) {
      .ld_contract_stop("the LD reported for a marker pair differs from its actual LD",
                        z, "indirect")
    }
  }
  invisible(TRUE)
}

#' Direct-LD contract violation (reason string) or NULL when the pairs are valid
#'
#' Distinct markers, no marker shared between the two traits, both members of
#' a pair on the same chromosome, absolute LD of every pair inside the
#' inclusive window `[ld_min, ld_max]`, and the reported LD equal to the
#' recomputed one.
#' @keywords internal
#' @noRd
.ld_direct_violation <- function(genotypes, anchors, partners, reported,
                                 ld_method, ld_min, ld_max) {
  if (length(anchors) != length(partners)) {
    return("the numbers of selected and linked markers differ")
  }
  if (any(anchors == partners)) {
    return("a marker was paired with itself")
  }
  if (anyDuplicated(anchors) || anyDuplicated(partners) ||
      anyDuplicated(.ld_key(genotypes, anchors)) ||
      anyDuplicated(.ld_key(genotypes, partners))) {
    return("a trait would have a duplicated QTN")
  }
  if (length(intersect(anchors, partners))) {
    return("the same marker was selected as a QTN for both traits")
  }
  if (any(genotypes$chr[anchors] != genotypes$chr[partners])) {
    return("a pair spans two chromosomes")
  }
  ld <- vapply(seq_along(anchors), function(k) {
    .ld_pair(genotypes, anchors[k], partners[k], ld_method)
  }, numeric(1))
  tol <- 1e-9
  if (anyNA(ld) || any(ld < ld_min - tol | ld > ld_max + tol)) {
    return("a selected pair has an absolute LD outside [ld_min, ld_max]")
  }
  if (length(reported) != length(ld) || any(abs(ld - reported) > 1e-6)) {
    return("the LD reported for a pair differs from its actual LD")
  }
  NULL
}

#' Verify direct-LD selections: distinct markers, same chromosome, LD in window
#' @keywords internal
#' @noRd
.ld_check_direct <- function(genotypes, anchors, partners, reported,
                             ld_method, ld_min, ld_max, z) {
  reason <- .ld_direct_violation(genotypes, anchors, partners, reported,
                                 ld_method, ld_min, ld_max)
  if (!is.null(reason)) .ld_contract_stop(reason, z, "direct")
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Bounded retry of the direct-LD marker search (frozen v1 engine).
#
# The first attempt for every replicate is the original search with the
# original seeds (bit-identical output whenever it meets the LD contract).
# Only when it does not -- the contract check fails, or the search stops with
# "None of the selected SNPs met the minimum/maximum LD threshold" or runs off
# the marker set -- the replicate is searched again, up to
# `.ld_max_attempts()` attempts in total (fewer if a retry attempt uses up every
# candidate marker, see `.ld_give_up_now()`), from the seed
# `.ld_attempt_seed(seed, attempt)`; attempts >= 2 also reset the neighbour
# pointers after a re-sample (`repaired`), which the frozen dominance walks do
# not (see the "LD architecture" section of ?create_phenotypes).
# ---------------------------------------------------------------------------

#' Maximum number of search attempts per replicate (the first is the frozen one)
#' @keywords internal
#' @noRd
.ld_max_attempts <- function() 50L

#' Distance between the seeds of consecutive attempts (a prime far larger than
#' any replicate/QTN offset added to a seed)
#' @keywords internal
#' @noRd
.ld_retry_stride <- 1000003

#' Seed used by search attempt `attempt` (1 = the original seed, unchanged)
#'
#' Attempt `a` >= 2 moves the seed by `(a - 1) * .ld_retry_stride` towards
#' zero (away from zero for `seed <= 0`), so `abs()` of the result never
#' exceeds `max(abs(seed), (.ld_max_attempts() - 1) * .ld_retry_stride)` and the
#' integer-range bound of `.v1_validate_seed_arith()` still holds.
#' @keywords internal
#' @noRd
.ld_attempt_seed <- function(seed, attempt) {
  if (attempt <= 1L) return(seed)
  step <- (attempt - 1L) * .ld_retry_stride
  if (seed > 0) seed - step else seed + step
}

#' Signal a (retriable) marker-search failure with the frozen message
#' @keywords internal
#' @noRd
.ld_search_stop <- function(msg) {
  stop(structure(class = c("ld_search_failed", "error", "condition"),
                 list(message = msg, call = NULL)))
}

#' Evaluate one search attempt; return its retriable failure (or NULL)
#'
#' `expr` is evaluated in the caller's frame (lazy argument), so its
#' assignments are visible to the caller. A search failure -- the frozen
#' "None of the selected SNPs ..." stops, the out-of-range GDS read, and the
#' exhaustion of every candidate marker (the frozen walk then ends in
#' `sample.int()`'s "invalid first argument") -- is returned as a condition;
#' every other error propagates unchanged. With `quiet`, warnings raised during
#' the attempt are held back and returned.
#' @keywords internal
#' @noRd
.ld_run_attempt <- function(expr, quiet = FALSE) {
  held <- list()
  failure <- withCallingHandlers(
    tryCatch({
      expr
      NULL
    },
    ld_search_failed = function(e) e,
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("'start' is invalid", msg, fixed = TRUE)) {
        e
      } else if (identical(msg, "invalid first argument") &&
                 grepl("sample", paste(deparse(conditionCall(e)), collapse = ""),
                       fixed = TRUE)) {
        structure(class = c("ld_search_exhausted", "ld_search_failed", "error",
                            "condition"),
                  list(message = paste0(
                    "LD architecture: the marker search used up every candidate marker without ",
                    "finding a pair inside the LD window [ld_min, ld_max]. Widen the window, use a ",
                    "different `ld_method`, or provide a denser marker set."),
                    call = NULL))
      } else {
        stop(e)
      }
    }),
    warning = function(w) {
      if (quiet) {
        held[[length(held) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    }
  )
  list(failure = failure, warnings = held)
}

#' Stop retrying? After `.ld_max_attempts()` attempts, or when a retry attempt
#' (>= 2) has used up every candidate marker (a full scan again would repeat)
#' @keywords internal
#' @noRd
.ld_give_up_now <- function(failure, attempt) {
  attempt >= .ld_max_attempts() ||
    (attempt > 1L && inherits(failure, "ld_search_exhausted"))
}

#' Re-signal the warnings held back during a successful retry attempt
#' @keywords internal
#' @noRd
.ld_emit_warnings <- function(held) {
  for (w in held) {
    warning(conditionMessage(w), call. = FALSE, immediate. = TRUE)
  }
  invisible(NULL)
}

#' Give up after `attempts` failed searches (condition: re-raise; else contract)
#' @keywords internal
#' @noRd
.ld_search_giveup <- function(failure, z, attempts) {
  if (inherits(failure, "condition")) stop(failure)
  .ld_contract_stop(failure, z, "direct", attempts = attempts)
}

#' qtn_linkage() with an informative message when the search runs off the
#' marker set
#'
#' The neighbour search is row-index based: when an intermediate/selected
#' marker sits at the first or last row (or the data set is tiny), the frozen
#' code reads outside the GDS file and fails with `'start' is invalid`.
#' @keywords internal
#' @noRd
.qtn_linkage_checked <- function(...) {
  tryCatch(
    qtn_linkage(...),
    error = function(e) {
      if (grepl("'start' is invalid", conditionMessage(e), fixed = TRUE)) {
        stop(
          "LD architecture: the search for a marker in LD ran past the first or last marker of the data set, ",
          "so no partner within [ld_min, ld_max] was found. Use a larger and denser marker set, a wider LD window, ",
          "or a different `seed`.",
          call. = FALSE
        )
      }
      stop(e)
    }
  )
}
