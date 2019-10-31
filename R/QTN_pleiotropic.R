#' Select SNPs to be assigned as QTNs.
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param add_QTN_num = NULL,
#' @param dom_QTN_num = NULL,
#' @param epi_QTN_num = NULL
#' @param same_add_dom_QTN = NULL,
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#------------------------------  QTN_pleiotropic -------------------------------
QTN_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           same_add_dom_QTN = NULL,
           add_QTN_num = NULL,
           dom_QTN_num = NULL,
           epi_QTN_num = NULL,
           constrains = list(maf_above = NULL, maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL
  ) {
    #---------------------------------------------------------------------------
    add_ef_trait_obj <- NULL
    dom_ef_trait_obj <- NULL
    epi_ef_trait_obj <-  NULL
    if (any(lengths(constrains) > 0)) {
      index <- constrain(genotypes = genotypes,
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      index <- 6:nrow(genotypes)
    }
    if (rep_by != "QTN") {
      rep <- 1
      }
    if (same_add_dom_QTN) {
      add_QTN_geno_info <- vector("list", rep)
      add_ef_trait_obj <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + i)
        }
        vector_of_add_QTN <-
          sample(index, add_QTN_num, replace = FALSE)
        add_QTN_geno_info[[i]] <-
          as.data.frame(genotypes[vector_of_add_QTN, ])
        add_ef_trait_obj[[i]] <-
          t(add_QTN_geno_info[[i]][, -c(1:5)])
        colnames(add_ef_trait_obj[[i]]) <-
          paste0(
            "Chr_",
            unlist(add_QTN_geno_info[[i]][, 3]),
            "_",
            unlist(add_QTN_geno_info[[i]][, 4])
          )
      }
      add_QTN_geno_info <-
        do.call(rbind, add_QTN_geno_info)
      add_QTN_geno_info <-
        data.frame(rep = rep(1:rep, each = add_QTN_num),
                   add_QTN_geno_info)
      if (!export_gt){
        add_QTN_geno_info <- add_QTN_geno_info[, 1:5]
      }
      if (!is.null(seed)) {
        s <- as.matrix(seed + 1:rep)
      } else {
        s <- "set.seed not assigned"
      }
      write.table(
        s,
        paste0("seed_num_for_", rep, "_reps_and_", add_QTN_num,
               "_Add_and_Dom_QTN", ".txt"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        add_QTN_geno_info,
        paste0(
          "geno_info_for_", rep, "_reps_and_",
          add_QTN_num,
          "_Add_and_Dom_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
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
            as.data.frame(genotypes[vector_of_add_QTN, ])
          add_ef_trait_obj[[i]] <-
            t(add_QTN_geno_info[[i]][, -c(1:5)])
          colnames(add_ef_trait_obj[[i]]) <-
            paste0(
              "Chr_",
              unlist(add_QTN_geno_info[[i]][, 3]),
              "_",
              unlist(add_QTN_geno_info[[i]][, 4])
            )
        }
        add_QTN_geno_info <-
          do.call(rbind, add_QTN_geno_info)
        add_QTN_geno_info <-
          data.frame(rep = rep(1:rep, each = add_QTN_num),
                     add_QTN_geno_info)
        if (!export_gt) {
          add_QTN_geno_info <- add_QTN_geno_info[, 1:5]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed + 1:rep)
        } else {
          s <- "set.seed not assigned"
        }
        write.table(
          s,
          paste0("seed_num_for_", rep, "_reps_and_", add_QTN_num,
                 "_Add_QTN", ".txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          add_QTN_geno_info,
          paste0(
            "geno_info_for_", rep, "_reps_and_",
            add_QTN_num,
            "_Add_QTN",
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
      if (dom) {
        dom_QTN_geno_info <- vector("list", rep)
        dom_ef_trait_obj <- vector("list", rep)
        for (i in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + i + 10)
          }
          vector_of_dom_QTN <-
            sample(index, dom_QTN_num, replace = FALSE)
          dom_QTN_geno_info[[i]] <-
            as.data.frame(genotypes[vector_of_dom_QTN, ])
          dom_ef_trait_obj[[i]] <-
            t(dom_QTN_geno_info[[i]][, -c(1:5)])
          colnames(dom_ef_trait_obj[[i]]) <-
            paste0(
              "Chr_",
              unlist(dom_QTN_geno_info[[i]][, 3]),
              "_",
              unlist(dom_QTN_geno_info[[i]][, 4])
            )
        }
        dom_QTN_geno_info <-
          do.call(rbind, dom_QTN_geno_info)
        dom_QTN_geno_info <-
          data.frame(rep = rep(1:rep, each = dom_QTN_num),
                     dom_QTN_geno_info)
        if (!export_gt) {
          dom_QTN_geno_info <- dom_QTN_geno_info[, 1:5]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed + 1:rep) + 10
        } else {
          s <- "set.seed not assigned"
        }
        write.table(
          s,
          paste0("seed_num_for_", rep, "_reps_and_", dom_QTN_num,
                 "Dom_QTN", ".txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          dom_QTN_geno_info,
          paste0(
            "geno_info_for_", rep, "_reps_and_",
            dom_QTN_num,
            "_dom_QTN",
            ".txt"
          ),
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
      }
    }
    if (epi) {
      epi_QTN_gen_infor <- vector("list", rep)
      epi_ef_trait_obj <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed * seed + i + 20)
        }
        vector_of_epi_QTN <-
          sample(index, (2 * epi_QTN_num), replace = FALSE)
        epi_QTN_gen_infor[[i]] <-
          as.data.frame(genotypes[vector_of_epi_QTN, ])
        epi_ef_trait_obj[[i]] <-
          t(epi_QTN_gen_infor[[i]][, -c(1:5)])
        colnames(epi_ef_trait_obj[[i]]) <-
          paste0(
            "Chr_",
            unlist(epi_QTN_gen_infor[[i]][, 3]),
            "_",
            unlist(epi_QTN_gen_infor[[i]][, 4])
          )
      }
      epi_QTN_gen_infor <-
        do.call(rbind, epi_QTN_gen_infor)
      epi_QTN_gen_infor <-
        data.frame(rep = rep(rep(1:rep, each = epi_QTN_num), each = 2),
                   epi_QTN_gen_infor)
      if (!export_gt) {
        epi_QTN_gen_infor <- epi_QTN_gen_infor[, 1:5]
      }
      if (!is.null(seed)) {
        ss <- as.matrix( (seed * seed) + 1:rep + 20)
      } else {
        ss <- "set.seed not assigned"
      }
      write.table(
        ss,
        paste0(
          "seed_num_for_", rep, "_reps_and_",
          epi_QTN_num,
          "Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        epi_QTN_gen_infor,
        paste0(
          "geno_info_for_",
          epi_QTN_num,
          "_epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    }
    return(list(
      add_ef_trait_obj = add_ef_trait_obj,
      dom_ef_trait_obj = dom_ef_trait_obj,
      epi_ef_trait_obj = epi_ef_trait_obj
    ))
  }
