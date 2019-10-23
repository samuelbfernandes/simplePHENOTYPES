#' Select SNPs to be assigned as QTNs.
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param additive_QTN_number = NULL,
#' @param dominance_QTN_number = NULL,
#' @param epistatic_QTN_number = NULL
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
           additive_QTN_number = NULL,
           dominance_QTN_number = NULL,
           epistatic_QTN_number = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL
  ) {
    #---------------------------------------------------------------------------
    additive_effect_trait_object <- NULL
    dominance_effect_trait_object <- NULL
    epistatic_effect_trait_object <-  NULL
    
    if (any(lengths(constrains) > 0)) {
      index <- constrain(genotypes = genotypes, 
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      # First SNP at column 6
      index <- 6:nrow(genotypes)
    }
    if (rep_by != 'QTN'){rep <- 1}
    if (same_add_dom_QTN) {
      add_QTN_genotypic_information <- vector("list", rep)
      additive_effect_trait_object <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + i)
        }
        vector_of_add_QTN <-
          sample(index, additive_QTN_number, replace = FALSE)
        add_QTN_genotypic_information[[i]] <- 
          as.data.frame(genotypes[vector_of_add_QTN, ])
        additive_effect_trait_object[[i]] <-
          t(add_QTN_genotypic_information[[i]][, -c(1:5)])
        colnames(additive_effect_trait_object[[i]]) <-
          paste0(
            "Chr_",
            unlist(add_QTN_genotypic_information[[i]][, 3]),
            "_",
            unlist(add_QTN_genotypic_information[[i]][, 4])
          )
      }
      add_QTN_genotypic_information <- 
        do.call(rbind, add_QTN_genotypic_information)
      add_QTN_genotypic_information <- 
        data.frame(rep = rep(1:rep, each=additive_QTN_number),
                   add_QTN_genotypic_information)
      if(!export_gt){
        add_QTN_genotypic_information <- add_QTN_genotypic_information[,1:5]
      }
      if (!is.null(seed)) {
        s <- as.matrix(seed+1:rep)
      } else {
        s <- "set.seed not assigned"
      }
      write.table(
        s,
        paste0("seed_number_for_",rep,"_reps_and_", additive_QTN_number,
               "_Add_and_Dom_QTN", ".txt"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        add_QTN_genotypic_information,
        paste0(
          "Genotypic_information_for_",rep, "_reps_and_",
          additive_QTN_number,
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
        add_QTN_genotypic_information <- vector("list", rep)
        additive_effect_trait_object <- vector("list", rep)
        for (i in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + i)
          }
          vector_of_add_QTN <-
            sample(index, additive_QTN_number, replace = FALSE)
          add_QTN_genotypic_information[[i]] <- 
            as.data.frame(genotypes[vector_of_add_QTN, ])
          additive_effect_trait_object[[i]] <-
            t(add_QTN_genotypic_information[[i]][, -c(1:5)])
          colnames(additive_effect_trait_object[[i]]) <-
            paste0(
              "Chr_",
              unlist(add_QTN_genotypic_information[[i]][, 3]),
              "_",
              unlist(add_QTN_genotypic_information[[i]][, 4])
            )
        }
        add_QTN_genotypic_information <- 
          do.call(rbind, add_QTN_genotypic_information)
        add_QTN_genotypic_information <- 
          data.frame(rep = rep(1:rep, each=additive_QTN_number),
                     add_QTN_genotypic_information)
        if(!export_gt){
          add_QTN_genotypic_information <- add_QTN_genotypic_information[,1:5]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed+1:rep)
        } else {
          s <- "set.seed not assigned"
        }
        write.table(
          s,
          paste0("seed_number_for_",rep,"_reps_and_", additive_QTN_number,
                 "_Add_QTN", ".txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          add_QTN_genotypic_information,
          paste0(
            "Genotypic_information_for_",rep, "_reps_and_",
            additive_QTN_number,
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
        dom_QTN_genotypic_information <- vector("list", rep)
        dominance_effect_trait_object <- vector("list", rep)
        for (i in 1:rep) {
          if (!is.null(seed)) {
            set.seed(seed + i + 10)
          }
          vector_of_dom_QTN <-
            sample(index, dominance_QTN_number, replace = FALSE)
          dom_QTN_genotypic_information[[i]] <- 
            as.data.frame(genotypes[vector_of_dom_QTN, ])
          dominance_effect_trait_object[[i]] <-
            t(dom_QTN_genotypic_information[[i]][, -c(1:5)])
          colnames(dominance_effect_trait_object[[i]]) <-
            paste0(
              "Chr_",
              unlist(dom_QTN_genotypic_information[[i]][, 3]),
              "_",
              unlist(dom_QTN_genotypic_information[[i]][, 4])
            )
        }
        dom_QTN_genotypic_information <- 
          do.call(rbind, dom_QTN_genotypic_information)
        dom_QTN_genotypic_information <- 
          data.frame(rep = rep(1:rep, each=dominance_QTN_number),
                     dom_QTN_genotypic_information)     
        if(!export_gt){
          dom_QTN_genotypic_information <- dom_QTN_genotypic_information[, 1:5]
        }
        if (!is.null(seed)) {
          s <- as.matrix(seed + 1:rep) + 10
        } else {
          s <- "set.seed not assigned"
        }
        write.table(
          s,
          paste0("seed_number_for_",rep,"_reps_and_", dominance_QTN_number,
                 "Dom_QTN", ".txt"),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        data.table::fwrite(
          dom_QTN_genotypic_information,
          paste0(
            "Genotypic_information_for_",rep, "_reps_and_",
            dominance_QTN_number,
            "_Dominance_QTN",
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
      epi_QTN_genotypic_information <- vector("list", rep)
      epistatic_effect_trait_object <- vector("list", rep)
      for (i in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed * seed + i + 20)
        }
        vector_of_epi_QTN <-
          sample(index, (2 * epistatic_QTN_number), replace = FALSE)
        epi_QTN_genotypic_information[[i]] <-
          as.data.frame(genotypes[vector_of_epi_QTN, ])
        epistatic_effect_trait_object[[i]] <-
          t(epi_QTN_genotypic_information[[i]][, -c(1:5)])
        colnames(epistatic_effect_trait_object[[i]]) <-
          paste0(
            "Chr_",
            unlist(epi_QTN_genotypic_information[[i]][, 3]),
            "_",
            unlist(epi_QTN_genotypic_information[[i]][, 4])
          )
      }
      epi_QTN_genotypic_information <- 
        do.call(rbind, epi_QTN_genotypic_information)
      epi_QTN_genotypic_information <- 
        data.frame(rep = rep(rep(1:rep, each=epistatic_QTN_number), each=2),
                   epi_QTN_genotypic_information)
      if(!export_gt){
        epi_QTN_genotypic_information <- epi_QTN_genotypic_information[,1:5]
      }
      if (!is.null(seed)) {
        ss <- as.matrix((seed * seed) + 1:rep + 20)
      } else {
        ss <- "set.seed not assigned"
      }
      write.table(
        ss,
        paste0(
          "seed_number_for_",rep, "_reps_and_",
          epistatic_QTN_number,
          "Epi_QTN",
          ".txt"
        ),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE
      )
      data.table::fwrite(
        epi_QTN_genotypic_information,
        paste0(
          "Genotypic_information_for_",
          epistatic_QTN_number,
          "_Epistatic_QTN",
          ".txt"
        ),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
      
    }
    return(list(
      additive_effect_trait_object = additive_effect_trait_object,
      dominance_effect_trait_object = dominance_effect_trait_object,
      epistatic_effect_trait_object = epistatic_effect_trait_object
    ))
  }
