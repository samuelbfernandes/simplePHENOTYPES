#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param additive_QTN_number = NULL,
#' @param dominance_QTN_number = NULL,
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
#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_linkage ---------------------------------
QTN_linkage <-
  function(genotypes = NULL,
           seed = NULL,
           additive_QTN_number = NULL,
           dominance_QTN_number = NULL,
           ld = NULL,
           gdsfile = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           same_add_dom_QTN = NULL,
           add = NULL,
           dom = NULL) {
    #---------------------------------------------------------------------------
    additive_effect_trait_object <- NULL
    dominance_effect_trait_object <- NULL
    if (rep_by != 'QTN'){rep <- 1}
    if (any(lengths(constrains)>0)) { 
      index <- constrain(genotypes = genotypes, 
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      # First SNP at column 6
      index <- 6:nrow(genotypes)
    }
    if (same_add_dom_QTN) {
      sup <- vector("list", rep)
      inf <- vector("list", rep)
      add_QTN_genotypic_information_sup <- vector("list", rep)
      add_QTN_genotypic_information_inf <- vector("list", rep)
      QTN_causing_ld <- vector("list", rep)
      results <- vector("list", rep)
      for(z in 1:rep){
        if (!is.null(seed)) {
          set.seed(seed+z)
        }
        vector_of_add_QTN <-
          sample(index, additive_QTN_number, replace = FALSE)
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
                start = c(j, 1),
                count = c(1, -1)
              )
            snp2 <-
              gdsfmt::read.gdsn(
                gdsfmt::index.gdsn(genofile, "genotype"),
                start = c(i, 1),
                count = c(1, -1)
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
                start = c(i2, 1),
                count = c(1, -1)
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
        # close the genotype file
        SNPRelate::snpgdsClose(genofile)
        sup[[z]] <- sup_temp
        inf[[z]] <- inf_temp
        QTN_causing_ld[[z]] <- 
          data.frame(SNP = "cause_of_LD", genotypes[vector_of_add_QTN, ])
        add_QTN_genotypic_information_sup[[z]] <- 
          data.frame(SNP = "QTN_upstream",genotypes[sup[[z]], ])
        add_QTN_genotypic_information_inf[[z]] <- 
          data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
        results[[z]] <- rbind(QTN_causing_ld[[z]],
                              add_QTN_genotypic_information_sup[[z]],
                              add_QTN_genotypic_information_inf[[z]])
      }
      results <- do.call(rbind, results)
      results <- 
        data.frame(rep = rep(1:rep, each=additive_QTN_number*3),
                   results)
      if(!export_gt){
        results <- results[,1:6]
      }
      write.table(
        c(seed + 1:rep),
        paste0(
          "seed_number_for_",additive_QTN_number,
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
          "Genotypic_information_for_",
          additive_QTN_number,
          "Add_and_Dom_QTN_with_LD_of_",
          paste(ld, collapse = "_"),
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      additive_effect_trait_object <- mapply(function(x,y) {
        rownames(x) <-
          paste0("Chr_",  x$chr, "_", x$pos)
        rownames(y) <-
          paste0("Chr_",  y$chr, "_", y$pos)
        b<- list(t(x[,-(1:6)]),t(y[,-(1:6)]))
      },
      x=add_QTN_genotypic_information_sup,
      y=add_QTN_genotypic_information_inf,
      SIMPLIFY = F)
      
    } else {
      if (add) {
        sup <- vector("list", rep)
        inf <- vector("list", rep)
        add_QTN_genotypic_information_sup <- vector("list", rep)
        add_QTN_genotypic_information_inf <- vector("list", rep)
        QTN_causing_ld <- vector("list", rep)
        results <- vector("list", rep)
        for(z in 1:rep){
          if (!is.null(seed)) {
            set.seed(seed+z)
          }
          vector_of_add_QTN <-
            sample(index, additive_QTN_number, replace = FALSE)
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
                  start = c(j, 1),
                  count = c(1, -1)
                )
              snp2 <-
                gdsfmt::read.gdsn(
                  gdsfmt::index.gdsn(genofile, "genotype"),
                  start = c(i, 1),
                  count = c(1, -1)
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
                  start = c(i2, 1),
                  count = c(1, -1)
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
          # close the genotype file
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- inf_temp
          QTN_causing_ld[[z]] <- 
            data.frame(SNP = "cause_of_LD", genotypes[vector_of_add_QTN, ])
          add_QTN_genotypic_information_sup[[z]] <- 
            data.frame(SNP = "QTN_upstream",genotypes[sup[[z]], ])
          add_QTN_genotypic_information_inf[[z]] <- 
            data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                add_QTN_genotypic_information_sup[[z]],
                                add_QTN_genotypic_information_inf[[z]])
        }
        results <- do.call(rbind, results)
        results <- 
          data.frame(rep = rep(1:rep, each=additive_QTN_number*3),
                     results)
        if(!export_gt){
          results <- results[,1:6]
        }
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",additive_QTN_number,
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
            "Genotypic_information_for_",
            additive_QTN_number,
            "Add_QTN_with_LD_of_",
            paste(ld, collapse = "_"),
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        additive_effect_trait_object <- mapply(function(x,y) {
          rownames(x) <-
            paste0("Chr_",  x$chr, "_", x$pos)
          rownames(y) <-
            paste0("Chr_",  y$chr, "_", y$pos)
          b<- list(t(x[,-(1:6)]),t(y[,-(1:6)]))
        },
        x=add_QTN_genotypic_information_sup,
        y=add_QTN_genotypic_information_inf,
        SIMPLIFY = F)
        
      }
      if (dom) {
        sup <- vector("list", rep)
        inf <- vector("list", rep)
        dom_QTN_genotypic_information_sup <- vector("list", rep)
        dom_QTN_genotypic_information_inf <- vector("list", rep)
        QTN_causing_ld <- vector("list", rep)
        results <- vector("list", rep)
        for(z in 1:rep){
          if (!is.null(seed)) {
            set.seed(seed + z + 10)
          }
          vector_of_dom_QTN <-
            sample(index, dominance_QTN_number, replace = FALSE)
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
                  start = c(j, 1),
                  count = c(1, -1)
                )
              snp2 <-
                gdsfmt::read.gdsn(
                  gdsfmt::index.gdsn(genofile, "genotype"),
                  start = c(i, 1),
                  count = c(1, -1)
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
                  start = c(i2, 1),
                  count = c(1, -1)
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
          # close the genotype file
          SNPRelate::snpgdsClose(genofile)
          sup[[z]] <- sup_temp
          inf[[z]] <- inf_temp
          QTN_causing_ld[[z]] <- 
            data.frame(SNP = "cause_of_LD", genotypes[vector_of_dom_QTN, ])
          dom_QTN_genotypic_information_sup[[z]] <- 
            data.frame(SNP = "QTN_upstream",genotypes[sup[[z]], ])
          dom_QTN_genotypic_information_inf[[z]] <- 
            data.frame(SNP = "QTN_downstream", genotypes[inf[[z]], ])
          results[[z]] <- rbind(QTN_causing_ld[[z]],
                                dom_QTN_genotypic_information_sup[[z]],
                                dom_QTN_genotypic_information_inf[[z]])
        }
        results <- do.call(rbind, results)
        results <- 
          data.frame(rep = rep(1:rep, each=dominance_QTN_number*3),
                     results)
        if(!export_gt){
          results <- results[,1:6]
        }
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",dominance_QTN_number,
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
            "Genotypic_information_for_",
            dominance_QTN_number,
            "Dom_QTN_with_LD_of_",
            paste(ld, collapse = "_"),
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        dominance_effect_trait_object <- mapply(function(x,y) {
          rownames(x) <-
            paste0("Chr_",  x$chr, "_", x$pos)
          rownames(y) <-
            paste0("Chr_",  y$chr, "_", y$pos)
          b<- list(t(x[,-(1:6)]),t(y[,-(1:6)]))
        },
        x=dom_QTN_genotypic_information_sup,
        y=dom_QTN_genotypic_information_inf,
        SIMPLIFY = F)
        
      }
    }
    return(list(additive_effect_trait_object = additive_effect_trait_object,
                dominance_effect_trait_object = dominance_effect_trait_object))
  }
