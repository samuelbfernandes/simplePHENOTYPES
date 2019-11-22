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
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes
#' Last update: Nov 05, 2019
#'
#'----------------------------- QTN_linkage ------------------------------------
QTN_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           ld = NULL,
           gdsfile = NULL,
           constrains = list(maf_above = NULL, maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           same_add_dom_QTN = NULL,
           add = NULL,
           dom = NULL) {
    #---------------------------------------------------------------------------
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    add_QTN <- TRUE 
    dom_QTN <- TRUE 
    if (!is.null(add_QTN_num)) {
      if (add_QTN_num == 0 ) {
        add_QTN <- FALSE
        add_QTN_num <- 1
      }
    }
    if (!is.null(dom_QTN_num)) {
      if (dom_QTN_num == 0 ) {
        dom_QTN <- FALSE
        dom_QTN_num <- 1
      }
    }
    if (rep_by != "QTN") {
      rep <- 1
      }
    if (any(lengths(constrains) > 0)) {
      index <- constrain(genotypes = genotypes,
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      index <- 6:nrow(genotypes)
    }
    if (same_add_dom_QTN & add) {
      sup <- vector("list", rep)
      inf <- vector("list", rep)
      add_gen_info_sup <- vector("list", rep)
      add_gen_info_inf <- vector("list", rep)
      QTN_causing_ld <- vector("list", rep)
      results <- vector("list", rep)
      for (z in 1:rep){
        if (!is.null(seed)) {
          set.seed(seed + z)
        }
        vector_of_add_QTN <-
          sample(index, add_QTN_num, replace = FALSE)
        genofile <- SNPRelate::snpgdsOpen(gdsfile)
        x <- 1
        sup_temp <- c()
        inf_temp <- c()
        for (j in vector_of_add_QTN) {
          ldsup <- 1
          i <- j + 1
          while (ldsup >= ld) {
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
          sup_temp[x] <- i
          ldinf <- 1
          i2 <- j - 1
          while (ldinf >= ld) {
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
          inf_temp[x] <- i2
          x <- x + 1
        }
        SNPRelate::snpgdsClose(genofile)
        sup[[z]] <- sup_temp
        inf[[z]] <- inf_temp
        QTN_causing_ld[[z]] <-
          data.frame(SNP = "cause_of_LD", genotypes[vector_of_add_QTN, ])
        add_gen_info_sup[[z]] <-
          data.frame(SNP = "QTN_upstream", genotypes[sup[[z]], ])
        add_gen_info_inf[[z]] <-
          data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
        results[[z]] <- rbind(QTN_causing_ld[[z]],
                              add_gen_info_sup[[z]],
                              add_gen_info_inf[[z]])
      }
      results <- do.call(rbind, results)
      results <-
        data.frame(rep = rep(1:rep, each = add_QTN_num * 3), results)
      if (!export_gt){
        results <- results[, 1:6]
      }
      if (add_QTN) {
      write.table(
        c(seed + 1:rep),
        paste0(
          "seed_num_for_", add_QTN_num,
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
        paste0(
          "Genotypic_info_for_",
          add_QTN_num,
          "Add_and_Dom_QTN_with_LD_of_",
          paste(ld, collapse = "_"),
          ".txt"
        ),
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
        results <- vector("list", rep)
        for (z in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + z)
          }
          vector_of_add_QTN <-
            sample(index, add_QTN_num, replace = FALSE)
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          x <- 1
          sup_temp <- c()
          inf_temp <- c()
          for (j in vector_of_add_QTN) {
            ldsup <- 1
            i <- j + 1
            while (ldsup >= ld) {
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
            sup_temp[x] <- i
            ldinf <- 1
            i2 <- j - 1
            while (ldinf >= ld) {
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
            inf_temp[x] <- i2
            x <- x + 1
          }
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- inf_temp
          QTN_causing_ld[[z]] <-
            data.frame(SNP = "cause_of_LD", genotypes[vector_of_add_QTN, ])
          add_gen_info_sup[[z]] <-
            data.frame(SNP = "QTN_upstream", genotypes[sup[[z]], ])
          add_gen_info_inf[[z]] <-
            data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                add_gen_info_sup[[z]],
                                add_gen_info_inf[[z]])
        }
        results <- do.call(rbind, results)
        results <-
          data.frame(rep = rep(1:rep, each = add_QTN_num * 3), results)
        if (!export_gt) {
          results <- results[, 1:6]
        }
        if (add_QTN) {
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_num_for_", add_QTN_num,
            "_Add_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          results,
          paste0(
            "Genotypic_info_for_",
            add_QTN_num,
            "Add_QTN_with_LD_of_",
            paste(ld, collapse = "_"),
            ".txt"
          ),
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
        results <- vector("list", rep)
        for (z in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + z + rep)
          }
          vector_of_dom_QTN <-
            sample(index, dom_QTN_num, replace = FALSE)
          genofile <- SNPRelate::snpgdsOpen(gdsfile)
          x <- 1
          sup_temp <- c()
          inf_temp <- c()
          for (j in vector_of_dom_QTN) {
            ldsup <- 1
            i <- j + 1
            while (ldsup >= ld) {
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
            sup_temp[x] <- i
            ldinf <- 1
            i2 <- j - 1
            while (ldinf >= ld) {
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
            inf_temp[x] <- i2
            x <- x + 1
          }
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- inf_temp
          QTN_causing_ld[[z]] <-
            data.frame(SNP = "cause_of_LD", genotypes[vector_of_dom_QTN, ])
          dom_gen_info_sup[[z]] <-
            data.frame(SNP = "QTN_upstream", genotypes[sup[[z]], ])
          dom_gen_info_inf[[z]] <-
            data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                dom_gen_info_sup[[z]],
                                dom_gen_info_inf[[z]])
        }
        results <- do.call(rbind, results)
        results <-
          data.frame(rep = rep(1:rep, each = dom_QTN_num * 3), results)
        if (!export_gt) {
          results <- results[, 1:6]
        }
        if (add_QTN) {
        write.table(
          c(seed + 1:rep + rep),
          paste0(
            "seed_num_for_", dom_QTN_num,
            "_Dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          results,
          paste0(
            "Genotypic_info_for_",
            dom_QTN_num,
            "Dom_QTN_with_LD_of_",
            paste(ld, collapse = "_"),
            ".txt"
          ),
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
        lapply(x, function(y){
          rnames <- rownames(y) 
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    if (!is.null(dom_ef_trait_obj) & !dom_QTN) {
      dom_ef_trait_obj <- lapply(dom_ef_trait_obj, function(x) {
        lapply(x, function(y){
          rnames <- rownames(y) 
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    return(list(add_ef_trait_obj = add_ef_trait_obj,
                dom_ef_trait_obj = dom_ef_trait_obj))
  }
