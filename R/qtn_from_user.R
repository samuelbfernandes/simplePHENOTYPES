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
#' @param ntraits = NULL
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
           add = NULL,
           dom = NULL,
           epi = NULL, 
           ntraits = NULL){
    ns <- ncol(genotypes) - 5
    if (ntraits == 1 | architecture == "pleiotropic") {
      add_ef_trait_obj <- NULL
      dom_ef_trait_obj <- NULL
      epi_ef_trait_obj <- NULL
      if (!is.null(unlist(QTN_list$add)) & add) {
        add_ef_trait_obj <-
          lapply(QTN_list$add[[1]], function(i) {
            a <- genotypes[genotypes$snp %in% i,]
            rownames(a) <- a$snp
            a <- a[i,]
          })
          maf <- round(apply(add_ef_trait_obj[[1]][, -1:-5], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          add_object <-
            data.frame(
              rep = 1,
              add_ef_trait_obj[[1]][, 1:5],
              maf = maf,
              add_ef_trait_obj[[1]][, -1:-5],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_ef_trait_obj[[1]] <- t(add_ef_trait_obj[[1]][, -1:-5])
          colnames(add_ef_trait_obj[[1]]) <-
            paste0("Chr_",  add_object$chr, "_", add_object$pos)
        if (!export_gt) {
          add_object <- add_object[, 1:7]
        }
        if (same_add_dom_QTN) {
          data.table::fwrite(
            add_object,
            "Additive_and_Dominance_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        } else {
          data.table::fwrite(
            add_object,
            "Additive_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
      if (!is.null(unlist(QTN_list$dom)) & dom) {
        dom_ef_trait_obj <-
          lapply(QTN_list$dom, function(i) {
            a <- genotypes[genotypes$snp %in% i,]
            rownames(a) <- a$snp
            a <- a[i,]
          })
        maf <- round(apply(dom_ef_trait_obj[[1]][, -1:-5], 1, function(x) {
          sumx <- ((sum(x) + ns) / ns * 0.5)
          min(sumx,  (1 - sumx))
        }), 4)
        dom_object <-
          data.frame(
            rep = 1,
            dom_ef_trait_obj[[1]][, 1:5],
            maf = maf,
            dom_ef_trait_obj[[1]][, -1:-5],
            check.names = FALSE,
            fix.empty.names = FALSE
          )
        dom_ef_trait_obj[[1]] <- t(dom_ef_trait_obj[[1]][, -1:-5])
        colnames(dom_ef_trait_obj[[1]]) <-
          paste0("Chr_",  dom_object$chr, "_", dom_object$pos)
        if (!export_gt) {
          dom_object <- dom_object[, 1:7]
        }
          data.table::fwrite(
            dom_object,
            "Dominance_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
      }
      if (!is.null(unlist(QTN_list$epi)) & epi) {
        epi_ef_trait_obj <- 
          list(lapply(QTN_list$epi, function(i) {
            e <- genotypes[genotypes$snp %in% i,]
            rownames(e) <- e$snp
            e <- t(e[i,-c(1:5)])
          }))
      }
      return(
        list(
          add_ef_trait_obj = add_ef_trait_obj,
          dom_ef_trait_obj = dom_ef_trait_obj,
          epi_ef_trait_obj = epi_ef_trait_obj
        )
      )
    } else {
      add_ef_trait_obj <- NULL
      dom_ef_trait_obj <- NULL
      epi_ef_trait_obj <- NULL
      if (!is.null(unlist(QTN_list$add)) & add) {
        add_ef_trait_obj <-
          list(lapply(QTN_list$add, function(i) {
            a <- genotypes[genotypes$snp %in% i,]
            rownames(a) <- a$snp
            a <- a[i,]
          }))
        add_object <- list()
        for(i in 1:lengths(add_ef_trait_obj)){
          maf <- round(apply(add_ef_trait_obj[[1]][[i]][, -1:-5], 1, function(x) {
            sumx <- ((sum(x) + ns) / ns * 0.5)
            min(sumx,  (1 - sumx))
          }), 4)
          add_object[[i]] <-
            data.frame(
              rep = 1,
              type = "user_specified",
              trait = paste0("trait_", i),
              add_ef_trait_obj[[1]][[i]][, 1:5],
              maf = maf,
              add_ef_trait_obj[[1]][[i]][, -1:-5],
              check.names = FALSE,
              fix.empty.names = FALSE
            )
          add_ef_trait_obj[[1]][[i]] <- t(add_ef_trait_obj[[1]][[i]][, -1:-5])
          colnames(add_ef_trait_obj[[1]][[i]]) <-
            paste0("Chr_",  add_object[[i]]$chr, "_", add_object[[i]]$pos)
        }
        add_object <- do.call("rbind", add_object)
        if (!export_gt) {
          add_object <- add_object[, 1:9]
        }
        pleio <-
          add_object$snp[duplicated(add_object$snp) |
                           duplicated(add_object$snp, fromLast = TRUE)]
        tab_p <- table(pleio) 
        if (any(tab_p < length(QTN_list$add))) {
          for (i in 1:length(tab_p)) {
            add_object$type[add_object$snp %in% names(tab_p[i])] <-
              ifelse(tab_p[i] < length(QTN_list$add),
                     paste0("Pleio_traits_", paste(
                       gsub("trait_", "", add_object$trait[add_object$snp %in% names(tab_p[i])]), collapse = "_")),
                     "Pleiotropic")
          }
        } else {
          add_object$type[duplicated(add_object$snp) |
                            duplicated(add_object$snp, fromLast = TRUE)] <-
            "Pleiotropic"
          add_object$type[!(duplicated(add_object$snp) |
                              duplicated(add_object$snp, fromLast = TRUE))] <-
            "trait_specific"
        }
        if (same_add_dom_QTN) {
          data.table::fwrite(
            add_object,
            "Additive_and_Dominance_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        } else {
          data.table::fwrite(
            add_object,
            "Additive_selected_QTNs.txt",
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }
      if (!is.null(unlist(QTN_list$dom)) & dom) {
        dom_ef_trait_obj <-
          list(lapply(QTN_list$dom, function(i) {
            a <- genotypes[genotypes$snp %in% i,]
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
              rep = 1,
              type = "user_specified",
              trait = paste0("trait_", i),
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
          dom_object <- dom_object[, 1:9]
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
          "Dominance_selected_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      if (!is.null(unlist(QTN_list$epi)) & epi) {
        epi_ef_trait_obj <- 
          list(lapply(QTN_list$epi, function(i) {
            e <- genotypes[genotypes$snp %in% i,]
            rownames(e) <- e$snp
            e <- t(e[i,-c(1:5)])
          }))
      }
      return(
        list(
          add_ef_trait_obj = add_ef_trait_obj,
          dom_ef_trait_obj = dom_ef_trait_obj,
          epi_ef_trait_obj = epi_ef_trait_obj
        )
      )
    }
   }
