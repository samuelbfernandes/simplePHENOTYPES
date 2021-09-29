#' Select SNPs to be assigned as QTNs
#' @keywords internal
#' @param genotypes = NULL, 
#' @param QTN_list ...
#' @param export_gt = NULL,
#' @param architecture = NULL,
#' @param same_add_dom_QTN = NULL,
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL, 
#' @param add_effect = NULL,
#' @param dom_effect = NULL,
#' @param epi_effect = NULL,
#' @param var_effect = NULL,
#' @param epi_type = NULL,
#' @param epi_interaction = 2,
#' @param ntraits = NULL
#' @param  type_of_ld = NULL,
#' @param ld_method = "composite",
#' @param gdsfile = NULL,
#' @param ld_min = NULL,
#' @param ld_max = NULL,
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Jan 7, 2021
#'
#'----------------------------- qtn_from_user ----------------------
qtn_from_user <-
  function(genotypes = NULL, 
           QTN_list,
           export_gt = NULL,
           architecture = NULL,
           same_add_dom_QTN = NULL,
           same_mv_QTN = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL, 
           var = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           var_effect = NULL,
           epi_type = NULL,
           epi_interaction = 2,
           ntraits = NULL,
           type_of_ld = NULL,
           ld_method = "composite",
           gdsfile = NULL,
           ld_min = NULL,
           ld_max = NULL,
           verbose = NULL){
    ns <- ncol(genotypes) - 5
    n <- nrow(genotypes)
    if (architecture == "LD") {
      add_ef_trait_obj <- NULL
      dom_ef_trait_obj <- NULL
      if (type_of_ld == "indirect") {
        if (same_add_dom_QTN & add) {
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          selected_snps <- genotypes$snp %in% unlist(QTN_list$add)
          if (sum(selected_snps) == 0) {
            stop(
              "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$add contain markers that are present in the marker dataset.",
              call. = F
            )
          }
          vector_of_add_QTN <-
            which(selected_snps)
          names(vector_of_add_QTN) <- genotypes$snp[vector_of_add_QTN]
          vector_of_add_QTN <- vector_of_add_QTN[unlist(QTN_list$add)]
          sup_temp <- c()
          inf_temp <- c()
          ld_between_QTNs_temp <- c()
          actual_ld_sup <- c()
          actual_ld_inf <- c()
          x <- 1
          for (j in vector_of_add_QTN) {
            ldsup <- 1
            ldinf <- 1
            i <- j + 1
            i2 <- j - 1
            while (ldsup > ld_max) {
              if (i > n) {
                if (verbose)
                  warning(
                    "There are no SNPs downstream. Selecting a different SNP.",
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
                    "There are no SNPs uptream. Selecting a different SNP.",
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
            x <- x + 1
          }
          if ((!any(genotypes[sup_temp,-(1:5)] == 0) |
               !any(genotypes[inf_temp,-(1:5)] == 0))) {
            warning("There are no heterozygote markers.",
                    call. = F,
                    immediate. = T)
          }
          if (ldsup < ld_min | ldinf < ld_min) {
            stop(
              "The selected SNP didn't meet the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
              call. = F
            )
          }
          if (ldsup > ld_max | ldinf > ld_max) {
            stop(
              "The selected SNP didn't meet the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
              call. = F
            )
          }
          SNPRelate::snpgdsClose(genofile)
          QTN_causing_ld <-
            data.frame(
              type = "cause_of_LD",
              trait = "none",
              genotypes[vector_of_add_QTN,],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup <-
            data.frame(
              type = "QTN_upstream",
              trait = "trait_1",
              genotypes[sup_temp,],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_inf <-
            data.frame(
              type = "QTN_downstream",
              trait = "trait_2",
              genotypes[inf_temp,],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results <- rbind(QTN_causing_ld,
                           add_gen_info_sup,
                           add_gen_info_inf)
          LD_summary <- data.frame(
            "ALL",
            QTN_causing_ld[, "snp"],
            ld_min,
            ld_max,
            actual_ld_inf,
            actual_ld_sup,
            add_gen_info_inf[, "snp"],
            add_gen_info_sup[, "snp"],
            ld_between_QTNs_temp,
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary) <-
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
          data.table::fwrite(
            LD_summary,
            "LD_Summary.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          maf <- round(apply(results[,-c(1:7)], 1, function(x) {
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
            results[,-c(1:7)],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          results <-
            data.frame(
              rep = "ALL",
              results,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results <- results[, 1:11]
          }
          data.table::fwrite(
            results,
            "Additive_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
          rownames(add_gen_info_sup) <-
            paste0("Chr_", add_gen_info_sup$chr, "_", add_gen_info_sup$pos)
          rownames(add_gen_info_inf) <-
            paste0("Chr_",  add_gen_info_inf$chr, "_", add_gen_info_inf$pos)
          add_ef_trait_obj <- list(list(t(add_gen_info_sup[,-(1:7)]),
                                        t(add_gen_info_inf[,-(1:7)])))
        } else {
          if (add) {
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            selected_snps <- genotypes$snp %in% unlist(QTN_list$add)
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$add contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            vector_of_add_QTN <-
              which(selected_snps)
            names(vector_of_add_QTN) <- genotypes$snp[vector_of_add_QTN]
            vector_of_add_QTN <- vector_of_add_QTN[unlist(QTN_list$add)]
            x <- 1
            sup_temp <- c()
            inf_temp <- c()
            ld_between_QTNs_temp <- c()
            actual_ld_sup <- c()
            actual_ld_inf <- c()
            for (j in vector_of_add_QTN) {
              ldsup <- 1
              ldinf <- 1
              i <- j + 1
              i2 <- j - 1
              while (ldsup > ld_max) {
                if (i > n) {
                  if (verbose)
                    warning(
                      "There are no SNPs downstream. Selecting a different SNP.",
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
                      "There are no SNPs upstream. Selecting a different SNP.",
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
              x <- x + 1
            }
            SNPRelate::snpgdsClose(genofile)
            if (ldsup < ld_min | ldinf < ld_min) {
              stop(
                "The selected SNP didn't meet the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
                call. = F
              )
            }
            if (ldsup > ld_max | ldinf > ld_max) {
              stop(
                "The selected SNP didn't meet the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
                call. = F
              )
            }
            QTN_causing_ld <-
              data.frame(
                type = "cause_of_LD",
                trait = "none",
                genotypes[vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup <-
              data.frame(
                type = "QTN_upstream",
                trait = "trait_1",
                genotypes[sup_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_inf <-
              data.frame(
                type = "QTN_downstream",
                trait = "trait_2",
                genotypes[inf_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add <- rbind(QTN_causing_ld,
                                 add_gen_info_sup,
                                 add_gen_info_inf)
            LD_summary_add <- data.frame(
              "ALL",
              QTN_causing_ld[, "snp"],
              ld_min,
              ld_max,
              actual_ld_inf,
              actual_ld_sup,
              add_gen_info_inf[, "snp"],
              add_gen_info_sup[, "snp"],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add) <-
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
            data.table::fwrite(
              LD_summary_add,
              "LD_Summary_Additive.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            maf <-
              round(apply(results_add[,-c(1:7)], 1, function(x) {
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
                results_add[,-c(1:7)],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add <-
              data.frame(
                rep = "ALL",
                results_add,
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            if (!export_gt) {
              results_add <- results_add[, 1:10]
            }
            data.table::fwrite(
              results_add,
              "Additive_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            rownames(add_gen_info_sup) <-
              paste0("Chr_", add_gen_info_sup$chr, "_", add_gen_info_sup$pos)
            rownames(add_gen_info_inf) <-
              paste0("Chr_",  add_gen_info_inf$chr, "_", add_gen_info_inf$pos)
            add_ef_trait_obj <- list(list(t(add_gen_info_sup[, -(1:7)]),
                                          t(add_gen_info_inf[, -(1:7)])))
          }
          if (dom) {
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            selected_snps <- genotypes$snp %in% unlist(QTN_list$dom)
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$dom contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            vector_of_dom_QTN <-
              which(selected_snps)
            names(vector_of_dom_QTN) <- genotypes$snp[vector_of_dom_QTN]
            vector_of_dom_QTN <- vector_of_dom_QTN[unlist(QTN_list$dom)]
            sup_temp <- c()
            inf_temp <- c()
            ld_between_QTNs_temp <- c()
            actual_ld_sup <- c()
            actual_ld_inf <- c()
            x <- 1
            for (j in vector_of_dom_QTN) {
              ldsup <- 1
              i <- j + 1
              while (ldsup > ld_max) {
                if (i > n) {
                  if (verbose)
                    warning(
                      "There are no SNPs downstream. Please select a different SNP.",
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
              ldinf <- 1
              i2 <- j - 1
              while (ldinf > ld_max) {
                if (i2 < 1) {
                  if (verbose)
                    warning(
                      "There are no SNPs upstream. Please select a different SNP.",
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
              x <- x + 1
            }
            if ((!any(genotypes[sup_temp, - (1:5)] == 0) | !any(genotypes[inf_temp, - (1:5)] == 0))) {
              warning(
                "There are no heterozygote markers.",
                call. = F,
                immediate. = T
              )
            }
            SNPRelate::snpgdsClose(genofile)
            if (ldsup < ld_min | ldinf < ld_min) {
              stop(
                "The selected SNP didn't meet the minimum LD threshold. Try another seed number or provide a genotypic file with enough LD!",
                call. = F
              )
            }
            if (ldsup > ld_max | ldinf > ld_max) {
              stop(
                "The selected SNP didn't meet the maximum LD threshold. Try another seed number or conduct an LD pruning in your genotypic file!",
                call. = F
              )
            }
            QTN_causing_ld <-
              data.frame(
                type = "cause_of_LD",
                trait = "none",
                genotypes[vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup <-
              data.frame(
                type = "QTN_upstream",
                trait = "trait_1",
                genotypes[sup_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_inf <-
              data.frame(
                type = "QTN_downstream",
                trait = "trait_2",
                genotypes[inf_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom <- rbind(QTN_causing_ld,
                                 dom_gen_info_sup,
                                 dom_gen_info_inf)
            LD_summary_dom <- data.frame(
              "ALL",
              QTN_causing_ld[, "snp"],
              ld_min,
              ld_max,
              actual_ld_inf,
              actual_ld_sup,
              dom_gen_info_inf[, "snp"],
              dom_gen_info_sup[, "snp"],
              ld_between_QTNs_temp,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom) <-
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
            data.table::fwrite(
              LD_summary_dom,
              "LD_Summary_Dominance.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
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
                rep = "ALL",
                results_dom,
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            if (!export_gt) {
              results_dom <- results_dom[, 1:10]
            }
            data.table::fwrite(
              results_dom,
              "Dominance_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            rownames(dom_gen_info_sup) <-
              paste0("Chr_", dom_gen_info_sup$chr, "_", dom_gen_info_sup$pos)
            rownames(dom_gen_info_inf) <-
              paste0("Chr_",  dom_gen_info_inf$chr, "_", dom_gen_info_inf$pos)
            dom_ef_trait_obj <- list(list(t(dom_gen_info_sup[, -(1:7)]),
                                          t(dom_gen_info_inf[, -(1:7)])))
          }
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
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          selected_snps <- genotypes$snp %in% unlist(QTN_list$add)
          if (sum(selected_snps) == 0) {
            stop(
              "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$add contain markers that are present in the marker dataset.",
              call. = F
            )
          }
          vector_of_add_QTN <-
            which(selected_snps)
          names(vector_of_add_QTN) <- genotypes$snp[vector_of_add_QTN]
          vector_of_add_QTN <- vector_of_add_QTN[unlist(QTN_list$add)]
          x <- 1
          sup_temp <- c()
          ld_between_QTNs_temp <- c()
          for (j in vector_of_add_QTN) {
            ldsup <- 1
            ldinf <- 1
            i <- j + 1
            i2 <- j - 1
            while ((ldsup > ld_max | ldsup < ld_min) &
                   (ldinf > ld_max | ldinf < ld_min))  {
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
              i <- i + 1
              i2 <- i2 - 1
            }
            if (((!any(genotypes[i - 1, - (1:5)] == 0) &
                  !any(genotypes[i2 + 1, - (1:5)] == 0) & dom))) {
              if (ldsup > ld_max) {
                ldsup <- 100
              } else if (ldinf > ld_max) {
                ldinf <- 100
              } else if (ldsup < ld_min) {
                ldsup <- 100
              } else if (ldinf < ld_min) {
                ldinf <- 100
              }
              warning(
                "There are no heterozygote markers.",
                call. = F,
                immediate. = T
              )
            } else if ((!any(genotypes[i - 1, - (1:5)] == 0) & dom) &
                       ldinf < ld_max & ldinf > ld_min) {
              ldsup <- 100
            } else if ((!any(genotypes[i2 + 1, - (1:5)] == 0) & dom) &
                       ldsup < ld_max & ldsup > ld_min) {
              ldinf <- 100
            } else {
              if (ldsup > ld_max) {
                ldsup <- 100
              } else if (ldinf > ld_max) {
                ldinf <- 100
              } else if (ldsup < ld_min) {
                ldsup <- 100
              } else if (ldinf < ld_min) {
                ldinf <- 100
              }
            }
            closest <- which.min(abs(c(ld_max - ldinf,  ld_max -ldsup)))
            ld_between_QTNs_temp[x] <- 
              ifelse(closest == 1 , ldinf, ldsup)
            sup_temp[x] <-
              ifelse(closest == 1 , i2 + 1, i - 1)
            x <- x + 1
          }
          add_gen_info_inf <-
            data.frame(
              type = "QTN_selected",
              trait = "trait_1",
              genotypes[vector_of_add_QTN, ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_gen_info_sup <-
            data.frame(
              type = "QTN_in_LD",
              trait = "trait_2",
              genotypes[sup_temp, ],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results <- rbind(add_gen_info_inf,
                           add_gen_info_sup)
          LD_summary <- data.frame(
            "ALL",
            ld_min,
            ld_max,
            ld_between_QTNs_temp,
            add_gen_info_inf[, "snp"],
            add_gen_info_sup[, "snp"],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
          colnames(LD_summary) <-
            c(
              "rep",
              "ld_min (absolute value)",
              "ld_max (absolute value)",
              "Actual_LD ",
              "QTN_for_trait_1",
              "QTN_for_trait_2"
            )
          SNPRelate::snpgdsClose(genofile)
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
          maf <- round(apply(results[, -c(1:7)], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          names(maf) <- results[, "snp"]
          results <-
            data.frame(
              results[, 1:2],
              additive_effect = unlist(add_effect),
              dominance_effect = unlist(dom_effect),
              results[, 3:7],
              maf = maf,
              results[, - c(1:7)],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          results <-
            data.frame(
              rep = "ALL",
              results,
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          if (!export_gt) {
            results <- results[, 1:11]
          }
          data.table::fwrite(
            results,
            "Additive_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
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
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            selected_snps <- genotypes$snp %in% unlist(QTN_list$add)
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$add contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            vector_of_add_QTN <-
              which(selected_snps)
            names(vector_of_add_QTN) <- genotypes$snp[vector_of_add_QTN]
            vector_of_add_QTN <- vector_of_add_QTN[unlist(QTN_list$add)]
            x <- 1
            sup_temp <- c()
            ld_between_QTNs_temp <- c()
            for (j in vector_of_add_QTN) {
              ldsup <- 1
              ldinf <- 1
              i <- j + 1
              i2 <- j - 1
              while ((ldsup > ld_max | ldsup < ld_min) &
                     (ldinf > ld_max | ldinf < ld_min))  {
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
                i <- i + 1
                i2 <- i2 - 1
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
              x <- x + 1
            }
            add_gen_info_inf <-
              data.frame(
                type = "QTN_selected",
                trait = "trait_1",
                genotypes[vector_of_add_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            add_gen_info_sup <-
              data.frame(
                type = "QTN_in_LD",
                trait = "trait_2",
                genotypes[sup_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add <- rbind(add_gen_info_inf,
                                 add_gen_info_sup)
            LD_summary_add <- data.frame(
              "ALL",
              ld_min,
              ld_max,
              ld_between_QTNs_temp,
              add_gen_info_inf[, "snp"],
              add_gen_info_sup[, "snp"],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_add) <-
              c(
                "rep",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
            SNPRelate::snpgdsClose(genofile)
            data.table::fwrite(
              LD_summary_add,
              "LD_Summary_Additive.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            maf <- round(apply(results_add[, -c(1:7)], 1, function(x) {
              sumx <- ((sum(x) + ns) / ns * 0.5)
              min(sumx,  (1 - sumx))
            }), 4)
            names(maf) <- results_add[, "snp"]
            results_add <-
              data.frame(
                results_add[, 1:2],
                additive_effect = unlist(add_effect),
                results_add[, 3:7],
                maf = maf,
                results_add[, - c(1:7)],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_add <-
              data.frame(
                rep = "ALL",
                results_add,
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            if (!export_gt) {
              results_add <- results_add[, 1:10]
            }
            data.table::fwrite(
              results_add,
              "Additive_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            rownames(add_gen_info_sup) <-
              paste0("Chr_", add_gen_info_sup$chr, "_", add_gen_info_sup$pos)
            rownames(add_gen_info_inf) <-
              paste0("Chr_",  add_gen_info_inf$chr, "_", add_gen_info_inf$pos)
            add_ef_trait_obj <- list(list(t(add_gen_info_sup[,-(1:7)]),
                                          t(add_gen_info_inf[,-(1:7)])))
          }
          if (dom) {
            genofile <- SNPRelate::snpgdsOpen(gdsfile)
            selected_snps <- genotypes$snp %in% unlist(QTN_list$dom)
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$dom contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            vector_of_dom_QTN <-
              which(selected_snps)
            names(vector_of_dom_QTN) <- genotypes$snp[vector_of_dom_QTN]
            vector_of_dom_QTN <- vector_of_dom_QTN[unlist(QTN_list$dom)]
            x <- 1
            sup_temp <- c()
            ld_between_QTNs_temp <- c()
            for (j in vector_of_dom_QTN) {
              ldsup <- 1
              ldinf <- 1
              i <- j + 1
              i2 <- j - 1
              while ((ldsup > ld_max | ldsup < ld_min) &
                     (ldinf > ld_max | ldinf < ld_min))  {
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
                i <- i + 1
                i2 <- i2 - 1
              }
              if (!any(genotypes[i - 1, - (1:5)] == 0) &
                  !any(genotypes[i2 + 1, - (1:5)] == 0) ) {
                if (ldsup > ld_max) {
                  ldsup <- 100
                } else if (ldinf > ld_max) {
                  ldinf <- 100
                } else if (ldsup < ld_min) {
                  ldsup <- 100
                } else if (ldinf < ld_min) {
                  ldinf <- 100
                }
                warning(
                  "There are no heterozygote markers.",
                  call. = F,
                  immediate. = T
                )
              } else if (!any(genotypes[i - 1, - (1:5)] == 0) &
                         ldinf < ld_max & ldinf > ld_min) {
                ldsup <- 100
              } else if (!any(genotypes[i2 + 1, - (1:5)] == 0) &
                         ldsup < ld_max & ldsup > ld_min) {
                ldinf <- 100
              } else {
                if (ldsup > ld_max) {
                  ldsup <- 100
                } else if (ldinf > ld_max) {
                  ldinf <- 100
                } else if (ldsup < ld_min) {
                  ldsup <- 100
                } else if (ldinf < ld_min) {
                  ldinf <- 100
                }
              }
              closest <- which.min(abs(c(ld_max - ldinf,  ld_max -ldsup)))
              ld_between_QTNs_temp[x] <- 
                ifelse(closest == 1 , ldinf, ldsup)
              sup_temp[x] <-
                ifelse(closest == 1 , i2 + 1, i - 1)
              x <- x + 1
            }
            dom_gen_info_inf <-
              data.frame(
                type = "QTN_selected",
                trait = "trait_1",
                genotypes[vector_of_dom_QTN, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            dom_gen_info_sup <-
              data.frame(
                type = "QTN_in_LD",
                trait = "trait_2",
                genotypes[sup_temp, ],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom <- rbind(dom_gen_info_inf,
                                 dom_gen_info_sup)
            LD_summary_dom <- data.frame(
              "ALL",
              ld_min,
              ld_max,
              ld_between_QTNs_temp,
              dom_gen_info_inf[, "snp"],
              dom_gen_info_sup[, "snp"],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
            colnames(LD_summary_dom) <-
              c(
                "rep",
                "ld_min (absolute value)",
                "ld_max (absolute value)",
                "Actual_LD ",
                "QTN_for_trait_1",
                "QTN_for_trait_2"
              )
            SNPRelate::snpgdsClose(genofile)
            data.table::fwrite(
              LD_summary_dom,
              "LD_Summary_Dominance.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            maf <- round(apply(results_dom[, -c(1:7)], 1, function(x) {
              sumx <- ((sum(x) + ns) / ns * 0.5)
              min(sumx,  (1 - sumx))
            }), 4)
            names(maf) <- results_dom[, "snp"]
            results_dom <-
              data.frame(
                results_dom[, 1:2],
                dominance_effect = unlist(dom_effect),
                results_dom[, 3:7],
                maf,
                results_dom[, - c(1:7)],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            results_dom <-
              data.frame(
                rep = "ALL",
                results_dom,
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            if (!export_gt) {
              results_dom <- results_dom[, 1:10]
            }
            data.table::fwrite(
              results_dom,
              "Dominance_QTNs.txt",
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
            rownames(dom_gen_info_sup) <-
              paste0("Chr_", dom_gen_info_sup$chr, "_", dom_gen_info_sup$pos)
            rownames(dom_gen_info_inf) <-
              paste0("Chr_",  dom_gen_info_inf$chr, "_", dom_gen_info_inf$pos)
            dom_ef_trait_obj <- list(list(t(dom_gen_info_sup[,-(1:7)]),
                                          t(dom_gen_info_inf[,-(1:7)])))
          }
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
    } else {
      add_ef_trait_obj <- NULL
      dom_ef_trait_obj <- NULL
      epi_ef_trait_obj <- NULL
      var_ef_trait_obj <- NULL
      if (add) {
        add_ef_trait_obj <-
          list(lapply(QTN_list$add, function(i) {
            selected_snps <- genotypes$snp %in% i
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$add contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            a <- genotypes[selected_snps,]
            rownames(a) <- a$snp
            a <- a[i,]
          }))
        add_object <- list()
        for(i in 1:lengths(add_ef_trait_obj)){
          maf <- round(apply(add_ef_trait_obj[[1]][[i]][, -1:-5], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          if (same_add_dom_QTN) {
            add_object[[i]] <-
              data.frame(
                rep = "ALL",
                type = "user_specified",
                trait = paste0("trait_", i),
                additive_effect = add_effect[[i]],
                dominance_effect = dom_effect[[i]],
                add_ef_trait_obj[[1]][[i]][, 1:5],
                maf = maf,
                add_ef_trait_obj[[1]][[i]][, -1:-5],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
          } else {
            add_object[[i]] <-
              data.frame(
                rep = "ALL",
                type = "user_specified",
                trait = paste0("trait_", i),
                additive_effect = add_effect[[i]],
                add_ef_trait_obj[[1]][[i]][, 1:5],
                maf = maf,
                add_ef_trait_obj[[1]][[i]][, -1:-5],
                check.names = FALSE,
                fix.empty.names = FALSE
              )
          }
          add_ef_trait_obj[[1]][[i]] <- t(add_ef_trait_obj[[1]][[i]][, -1:-5])
          colnames(add_ef_trait_obj[[1]][[i]]) <-
            paste0("Chr_",  add_object[[i]]$chr, "_", add_object[[i]]$pos)
        }
        add_object <- do.call("rbind", add_object)
        if (!export_gt) {
          if (same_add_dom_QTN) {
            add_object <- add_object[, 1:11]
          } else {
            add_object <- add_object[, 1:10]
          }
        }
        pleio <-
          add_object$snp[duplicated(add_object$snp) |
                           duplicated(add_object$snp, fromLast = TRUE)]
        tab_p <- table(pleio) 
        if (any(tab_p < length(QTN_list$add))) {
          for (i in 1:length(tab_p)) {
            add_object$type[add_object$snp %in% names(tab_p[i])] <-
              ifelse(tab_p[i] < length(QTN_list$add),
                     paste0("Pleio_traits_", paste(gsub("trait_", "", add_object$trait[add_object$snp %in% names(tab_p[i])]), collapse = "_")),
                     "Pleiotropic")
          }
          add_object$type[!add_object$snp %in% names(tab_p)] <- "trait_specific"
        } else {
          add_object$type[duplicated(add_object$snp) |
                            duplicated(add_object$snp, fromLast = TRUE)] <-
            "Pleiotropic"
          add_object$type[!(duplicated(add_object$snp) |
                              duplicated(add_object$snp, fromLast = TRUE))] <-
            "trait_specific"
        }
          data.table::fwrite(
            add_object,
            "Additive_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
      }
      if (dom & !same_add_dom_QTN) {
        dom_ef_trait_obj <-
          list(lapply(QTN_list$dom, function(i) {
            selected_snps <- genotypes$snp %in% i
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$dom contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            a <- genotypes[selected_snps,]
            rownames(a) <- a$snp
            a <- a[i,]
          }))
        dom_object <- list()
        for(i in 1:lengths(dom_ef_trait_obj)){
          maf <- round(apply(dom_ef_trait_obj[[1]][[i]][, -1:-5], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          dom_object[[i]] <-
            data.frame(
              rep = "ALL",
              type = "user_specified",
              trait = paste0("trait_", i),
              dominance_effect = dom_effect[[i]],
              dom_ef_trait_obj[[1]][[i]][, 1:5],
              maf = maf,
              dom_ef_trait_obj[[1]][[i]][, -1:-5],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          dom_ef_trait_obj[[1]][[i]] <- t(dom_ef_trait_obj[[1]][[i]][, -1:-5])
          colnames(dom_ef_trait_obj[[1]][[i]]) <-
            paste0("Chr_",  dom_object[[i]]$chr, "_", dom_object[[i]]$pos)
        }
        dom_object <- do.call("rbind", dom_object)
        if (!export_gt) {
          dom_object <- dom_object[, 1:10]
        }
        pleio <-
          dom_object$snp[duplicated(dom_object$snp) |
                           duplicated(dom_object$snp, fromLast = TRUE)]
        tab_p <- table(pleio) 
        if (any(tab_p < length(QTN_list$dom))) {
          for (i in 1:length(tab_p)) {
            dom_object$type[dom_object$snp %in% names(tab_p[i])] <-
              ifelse(tab_p[i] < length(QTN_list$dom),
                     paste0("Pleio_traits_", paste(
                       gsub("trait_", "", dom_object$trait[dom_object$snp %in% names(tab_p[i])]), collapse = "_")),
                     "Pleiotropic")
          }
          dom_object$type[!dom_object$snp %in% names(tab_p)] <- "trait_specific"
        } else {
          dom_object$type[duplicated(dom_object$snp) |
                            duplicated(dom_object$snp, fromLast = TRUE)] <-
            "Pleiotropic"
          dom_object$type[!(duplicated(dom_object$snp) |
                              duplicated(dom_object$snp, fromLast = TRUE))] <-
            "trait_specific"
        }
        data.table::fwrite(
          dom_object,
          "Dominance_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      if (epi)  {
        epi_ef_trait_obj <-
          list(lapply(QTN_list$epi, function(i) {
            selected_snps <- genotypes$snp %in% i
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$epi contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            a <- genotypes[selected_snps,]
            rownames(a) <- a$snp
            a <- a[i,]
          }))
        epi_object <- list()
        e_len <- lengths(QTN_list$epi)/epi_interaction
        for(i in 1:lengths(epi_ef_trait_obj)){
          maf <- round(apply(epi_ef_trait_obj[[1]][[i]][, -1:-5], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          epi_object[[i]] <-
            data.frame(
              rep = "ALL",
              QTN = rep(1:e_len[i], each = epi_interaction),
              type = "user_specified",
              trait = paste0("trait_", i),
              epistatic_effect = rep(epi_effect[[i]], each = e_len[i]),
              epi_ef_trait_obj[[1]][[i]][, 1:5],
              maf = maf,
              epi_ef_trait_obj[[1]][[i]][, -1:-5],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          epi_ef_trait_obj[[1]][[i]] <- t(epi_ef_trait_obj[[1]][[i]][, -1:-5])
          colnames(epi_ef_trait_obj[[1]][[i]]) <-
            paste0("Chr_",  epi_object[[i]]$chr, "_", epi_object[[i]]$pos)
        }
        epi_object <- do.call("rbind", epi_object)
        if (!export_gt) {
          epi_object <- epi_object[, 1:11]
        }
        pleio <-
          epi_object$snp[duplicated(epi_object$snp) |
                           duplicated(epi_object$snp, fromLast = TRUE)]
        tab_p <- table(pleio) 
        if (any(tab_p < length(QTN_list$epi))) {
          for (i in 1:length(tab_p)) {
            epi_object$type[epi_object$snp %in% names(tab_p[i])] <-
              ifelse(tab_p[i] < length(QTN_list$epi),
                     paste0("Pleio_traits_", paste(
                       gsub("trait_", "", epi_object$trait[epi_object$snp %in% names(tab_p[i])]), collapse = "_")),
                     "Pleiotropic")
          }
          epi_object$type[!epi_object$snp %in% names(tab_p)] <- "trait_specific"
        } else {
          epi_object$type[duplicated(epi_object$snp) |
                            duplicated(epi_object$snp, fromLast = TRUE)] <-
            "Pleiotropic"
          epi_object$type[!(duplicated(epi_object$snp) |
                              duplicated(epi_object$snp, fromLast = TRUE))] <-
            "trait_specific"
        }
        data.table::fwrite(
          epi_object,
          "Epistatic_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      if (var & !same_mv_QTN) {
        var_ef_trait_obj <-
          lapply(QTN_list$var[[1]], function(i) {
            selected_snps <- genotypes$snp %in% i
            if (sum(selected_snps) == 0) {
              stop(
                "None of the markers provided where found in the genotypic dataset (for at least one of the traits). \nPlease verify that QTN_list$var contain markers that are present in the marker dataset.",
                call. = F
              )
            }
            a <- genotypes[selected_snps,]
            rownames(a) <- a$snp
            a <- a[i,]
          })
        maf <- round(apply(var_ef_trait_obj[[1]][, -1:-5], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        var_object <-
          data.frame(
            rep = "ALL",
            variance_effect = var_effect[[1]],
            var_ef_trait_obj[[1]][, 1:5],
            maf = maf,
            var_ef_trait_obj[[1]][, -1:-5],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        var_ef_trait_obj[[1]] <- t(var_ef_trait_obj[[1]][, -1:-5])
        colnames(var_ef_trait_obj[[1]]) <-
          paste0("Chr_",  var_object$chr, "_", var_object$pos)
        if (!export_gt) {
          var_object <- var_object[, 1:10]
        }
        data.table::fwrite(
          var_object,
          "Variance_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      return(
        list(
          add_ef_trait_obj = add_ef_trait_obj,
          dom_ef_trait_obj = dom_ef_trait_obj,
          epi_ef_trait_obj = epi_ef_trait_obj,
          var_ef_trait_obj = var_ef_trait_obj
        )
      )
    }
  }
