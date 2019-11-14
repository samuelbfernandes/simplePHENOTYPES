#' Select SNPs to be assigned as QTNs
#' @keywords internal
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add  null
#' @param dom = NULL,
#' @param epi = NULL
#' @param same_add_dom_QTN = NULL,
#' @param pleio_a = NULL,
#' @param pleio_d = NULL,
#' @param pleio_e = NULL,
#' @param trait_spec_a_QTN_num = NULL,
#' @param trait_spec_d_QTN_num = NULL,
#' @param trait_spec_e_QTN_num = NULL,
#' @param ntraits = NULL
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @return Genotype of selected SNPs
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Nov 05, 2019
#'
#'----------------------------- QTN_partially_pleiotropic ----------------------
QTN_partially_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           pleio_a = NULL,
           pleio_d = NULL,
           pleio_e = NULL,
           trait_spec_a_QTN_num = NULL,
           trait_spec_d_QTN_num = NULL,
           trait_spec_e_QTN_num = NULL,
           ntraits = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           same_add_dom_QTN = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL) {
    #---------------------------------------------------------------------------
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    epi_ef_trait_obj <-  NULL
    add_QTN <- TRUE 
    dom_QTN <- TRUE 
    epi_QTN <- TRUE 
    if (!is.null(pleio_a) &
        !is.null(trait_spec_a_QTN_num)) {
      if (all((pleio_a + trait_spec_a_QTN_num) == 0 )) {
        add_QTN <- FALSE
        pleio_a <- 1
        trait_spec_a_QTN_num <- rep(1, ntraits)
      }
    }
    if (!is.null(pleio_d) &
        !is.null(trait_spec_d_QTN_num)) {
      if (all((pleio_d + trait_spec_d_QTN_num) == 0 )) {
        dom_QTN <- FALSE
        pleio_d <- 1
        trait_spec_d_QTN_num <- rep(1, ntraits)
      }
    }
    if (!is.null(pleio_e) &
        !is.null(trait_spec_e_QTN_num)) {
      if (all((pleio_e + trait_spec_e_QTN_num) == 0 )) {
        epi_QTN <- FALSE
        pleio_e <- 1
        trait_spec_e_QTN_num <- rep(1, ntraits)
      }
    }
    if (any(lengths(constrains) > 0)) {
      index <- constrain(genotypes = genotypes,
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      index <- 6:nrow(genotypes)
    }
    if (rep_by != "QTN"){
      rep <- 1
      }
    if (same_add_dom_QTN & add) {
      add_pleio_gen_info <- vector("list", rep)
      add_specific_gen_info <- vector("list", rep)
      for (j in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + j)
        }
        vec_of_pleio_add_QTN <-
          sample(index, pleio_a, replace = FALSE)
        add_pleio_gen_info[[j]] <-
          as.data.frame(genotypes[vec_of_pleio_add_QTN, ])
        snps <-
          setdiff(index, vec_of_pleio_add_QTN)
        vec_spec_add_QTN_temp <- vector("list", ntraits)
        add_specific_gen_info_temp <- vector("list", ntraits)
        ss <- c()
        for (i in 1:ntraits) {
          if (!is.null(seed)) {
            ss[i] <- seed + i + j
            set.seed(seed + i + j)
          }
          vec_spec_add_QTN_temp[[i]] <-
            sample(snps, trait_spec_a_QTN_num[i], replace = FALSE)
          snps <- setdiff(snps, vec_spec_add_QTN_temp[[i]])
          add_specific_gen_info_temp[[i]] <-
            as.data.frame(genotypes[vec_spec_add_QTN_temp[[i]], ])
        }
        add_specific_gen_info_temp <-
          do.call(rbind, add_specific_gen_info_temp)
        add_specific_gen_info[[j]] <-
          data.frame(trait = paste0("trait_", rep(1:ntraits,
                                                 trait_spec_a_QTN_num)),
                     add_specific_gen_info_temp)
      }
      add_object <- mapply(function(x, y) {
        p <- split(y, y[, 1])
        names(p) <- NULL
        lapply(p, function(z) {
          x <- data.frame(type = "Pleiotropic", trait = unique(z[, 1]), x)
          z <- data.frame(type = "trait_specific", z)
          rbind(x, z)
        }
        )
      },
      x = add_pleio_gen_info,
      y = add_specific_gen_info,
      SIMPLIFY = F)
      add_ef_trait_obj <- add_object
      add_object <- unlist(add_object, recursive = FALSE)
      add_object <- do.call(rbind, add_object)
      add_object <-
        data.frame(rep = sort(c(rep(1:rep,
                                    each = pleio_a * ntraits),
                                rep(1:rep,
                                    each = sum(trait_spec_a_QTN_num)))),
                   add_object)
      add_ef_trait_obj <-
        lapply(add_ef_trait_obj, function(x) {
          lapply(x, function(b){
            rownames(b) <-
              paste0("Chr_", b$chr, "_", b$pos)
            b <- b[, - (1:7)]
            return(t(b))
          })
        })
      if (!export_gt) {
        add_object <- add_object[, 1:7]
      }
      if (add_QTN) {
      write.table(
        c(seed + 1:rep),
        paste0(
          "seed_number_for_",
          paste0(pleio_a, collapse = "_"),
          "pleiotropic_Add_and_Dom_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      write.table(
        ss,
        paste0(
          "seed_number_for_",
          paste0(trait_spec_a_QTN_num, collapse = "_"),
          "trait_specific_Add_and_Dom_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        add_object,
        "Additive_and_Dominance_selected_QTNs.txt",
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      }
    } else {
      if (add) {
        add_pleio_gen_info <- vector("list", rep)
        add_specific_gen_info <- vector("list", rep)
        for (j in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + j)
          }
          vec_of_pleio_add_QTN <-
            sample(index, pleio_a, replace = FALSE)
          add_pleio_gen_info[[j]] <-
            as.data.frame(genotypes[vec_of_pleio_add_QTN, ])
          snps <-
            setdiff(index, vec_of_pleio_add_QTN)
          vec_spec_add_QTN_temp <- vector("list", ntraits)
          add_specific_gen_info_temp <- vector("list", ntraits)
          ss <- c()
          for (i in 1:ntraits) {
            if (!is.null(seed)) {
              ss[i] <- seed + i + j
              set.seed(seed + i + j)
            }
            vec_spec_add_QTN_temp[[i]] <-
              sample(snps, trait_spec_a_QTN_num[i], replace = FALSE)
            snps <- setdiff(snps, vec_spec_add_QTN_temp[[i]])
            add_specific_gen_info_temp[[i]] <-
              as.data.frame(genotypes[vec_spec_add_QTN_temp[[i]], ])
          }
          add_specific_gen_info_temp <-
            do.call(rbind, add_specific_gen_info_temp)
          add_specific_gen_info[[j]] <-
            data.frame(trait = paste0("trait_",
                                      rep(1:ntraits,
                                          trait_spec_a_QTN_num)),
                       add_specific_gen_info_temp)
        }
        add_object <- mapply(function(x, y) {
          p <- split(y, y[, 1])
          names(p) <- NULL
          lapply(p, function(z) {
            x <- data.frame(type = "Pleiotropic", trait = unique(z[, 1]), x)
            z <- data.frame(type = "trait_specific", z)
            rbind(x, z)
          }
          )
        },
        x = add_pleio_gen_info,
        y = add_specific_gen_info,
        SIMPLIFY = F)
        add_ef_trait_obj <- add_object
        add_object <- unlist(add_object, recursive = FALSE)
        add_object <- do.call(rbind, add_object)
        add_object <-
          data.frame(rep = sort(c(rep(1:rep,
                                      each = pleio_a * ntraits),
                                  rep(1:rep,
                                      each = sum(trait_spec_a_QTN_num)))),
                     add_object)
        add_ef_trait_obj <-
          lapply(add_ef_trait_obj, function(x) {
            lapply(x, function(b) {
              rownames(b) <-
                paste0("Chr_",  b$chr, "_", b$pos)
              b <- b[, - (1:7)]
              return(t(b))
            })
          })
        if (!export_gt) {
          add_object <- add_object[, 1:7]
        }
        if (add_QTN) {
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",
            paste0(trait_spec_a_QTN_num + pleio_a, collapse = "_"),
            "pleiotropic_Add_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        write.table(
          ss,
          paste0(
            "seed_number_for_",
            paste0(trait_spec_a_QTN_num, collapse = "_"),
            "trait_specific_Add_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
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
      if (dom) {
        dom_pleio_gen_info <- vector("list", rep)
        dom_spec_gen_info <- vector("list", rep)
        for (j in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + j + rep)
          }
          vec_pleio_dom_QTN <-
            sample(index, pleio_d, replace = FALSE)
          dom_pleio_gen_info[[j]] <-
            as.data.frame(genotypes[vec_pleio_dom_QTN, ])
          snpsd <-
            setdiff(index, vec_pleio_dom_QTN)
          vec_spec_dom_QTN_temp <- vector("list", ntraits)
          dom_spec_gen_info_temp <- vector("list", ntraits)
          ssd <- c()
          for (i in 1:ntraits) {
            if (!is.null(seed)) {
              ssd[i] <- seed + i + j + rep
              set.seed(seed + i + j + rep)
            }
            vec_spec_dom_QTN_temp[[i]] <-
              sample(snpsd, trait_spec_d_QTN_num[i], replace = FALSE)
            snpsd <- setdiff(snpsd, vec_spec_dom_QTN_temp[[i]])
            dom_spec_gen_info_temp[[i]] <-
              as.data.frame(genotypes[vec_spec_dom_QTN_temp[[i]], ])
          }
          dom_spec_gen_info_temp <-
            do.call(rbind, dom_spec_gen_info_temp)
          dom_spec_gen_info[[j]] <-
            data.frame(trait = paste0("trait_",
                                      rep(1:ntraits,
                                          trait_spec_d_QTN_num)),
                       dom_spec_gen_info_temp)
        }
        dom_object <- mapply(function(x, y) {
          p <- split(y, y[, 1])
          names(p) <- NULL
          lapply(p, function(z) {
            x <- data.frame(type = "Pleiotropic", trait = unique(z[, 1]), x)
            z <- data.frame(type = "trait_specific", z)
            rbind(x, z)
          }
          )
        },
        x = dom_pleio_gen_info,
        y = dom_spec_gen_info,
        SIMPLIFY = F)
        dom_ef_trait_obj <- dom_object
        dom_object <- unlist(dom_object, recursive = FALSE)
        dom_object <- do.call(rbind, dom_object)
        dom_object <-
          data.frame(rep = sort(c(rep(1:rep,
                                      each = pleio_d * ntraits),
                                  rep(1:rep,
                                      each = sum(trait_spec_d_QTN_num)))),
                     dom_object)
        dom_ef_trait_obj <-
          lapply(dom_ef_trait_obj, function(x) {
            lapply(x, function(b){
              rownames(b) <-
                paste0("Chr_",  b$chr, "_", b$pos)
              b <- b[, - (1:7)]
              return(t(b))
            })
          })
        if (!export_gt){
          dom_object <- dom_object[, 1:7]
        }
        if (dom_QTN) {
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",
            paste0(pleio_a, collapse = "_"),
            "pleiotropic_Dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        write.table(
          ssd,
          paste0(
            "seed_number_for_",
            paste0(trait_spec_a_QTN_num, collapse = "_"),
            "trait_specific_Dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          dom_object,
          "Dominance_selected_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        }
      }
    }
    if (epi) {
      epi_pleio_QTN_gen_info <- vector("list", rep)
      epi_spec_QTN_gen_info <- vector("list", rep)
      for (j in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + seed + j)
        }
        vec_pleio_epi_QTN <-
          sample(index, (2 * pleio_e), replace = FALSE)
        epi_pleio_QTN_gen_info[[j]] <-
          as.data.frame(genotypes[vec_pleio_epi_QTN, ])
        snps_e <-
          setdiff(index, vec_pleio_epi_QTN)
        vec_spec_epi_QTN_temp <- vector("list", ntraits)
        epi_spec_QTN_gen_info_temp <- vector("list", ntraits)
        sse <- c()
        for (i in 1:ntraits) {
          if (!is.null(seed)) {
            sse[i] <- seed + i + seed + j
            set.seed(seed + i + seed + j)
          }
          vec_spec_epi_QTN_temp[[i]] <-
            sample(snps_e, (2 * trait_spec_e_QTN_num[i]), replace = FALSE)
          snps_e <- setdiff(snps_e, vec_spec_epi_QTN_temp[[i]])
          epi_spec_QTN_gen_info_temp[[i]] <-
            as.data.frame(genotypes[vec_spec_epi_QTN_temp[[i]], ])
        }
        epi_spec_QTN_gen_info_temp <-
          do.call(rbind, epi_spec_QTN_gen_info_temp)
        epi_spec_QTN_gen_info[[j]] <-
          data.frame(trait = paste0("trait_",
                                    rep(1:ntraits,
                                        (2 * trait_spec_e_QTN_num))),
                     epi_spec_QTN_gen_info_temp)
      }
      epi_object <- mapply(function(x, y) {
        p <- split(y, y[, 1])
        names(p) <- NULL
        lapply(p, function(z) {
          x <- data.frame(type = "Pleiotropic", trait = unique(z[, 1]), x)
          z <- data.frame(type = "trait_specific", z)
          rbind(x, z)
        }
        )
      },
      x = epi_pleio_QTN_gen_info,
      y = epi_spec_QTN_gen_info,
      SIMPLIFY = F)
      epi_ef_trait_obj <- epi_object
      epi_object <- unlist(epi_object, recursive = FALSE)
      epi_object <- do.call(rbind, epi_object)
      epi_object <-
        data.frame(rep = sort(c(rep(1:rep,
                                    each = pleio_e * ntraits * 2),
                                rep(1:rep,
                                    each = sum(trait_spec_e_QTN_num) * 2))),
                   epi_object)
      epi_ef_trait_obj <-
        lapply(epi_ef_trait_obj, function(x) {
          lapply(x, function(b) {
            rownames(b) <-
              paste0("Chr_", b$chr, "_", b$pos)
            b <- b[, - (1:7)]
            return(t(b))
          })
        })
      if (!export_gt) {
        epi_object <- epi_object[, 1:7]
      }
      if (epi_QTN) {
      write.table(
        c(seed + seed + 1:rep),
        paste0(
          "seed_number_for_",
          paste0(pleio_e, collapse = "_"),
          "pleiotropic_Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      write.table(
        sse,
        paste0(
          "seed_number_for_",
          paste0(trait_spec_e_QTN_num, collapse = "_"),
          "trait_specific_Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        epi_object,
        "Epistatic_selected_QTNs.txt",
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
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
    if (!is.null(epi_ef_trait_obj) & !epi_QTN) {
      epi_ef_trait_obj <- lapply(epi_ef_trait_obj, function(x) {
        lapply(x, function(y){
          rnames <- rownames(y) 
          y <- matrix(0, nrow = nrow(y), ncol =  1)
          rownames(y)  <- rnames
          return(y)
        })
      })
    }
    return(list(
      add_ef_trait_obj = add_ef_trait_obj,
      dom_ef_trait_obj = dom_ef_trait_obj,
      epi_ef_trait_obj = epi_ef_trait_obj
    ))
  }
