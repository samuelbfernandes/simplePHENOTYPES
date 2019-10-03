#' Select SNPs to be assigned as QTNs
#' @export
#' @param genotypes = NULL,
#' @param seed = NULL,
#' @param overlap = NULL,
#' @param overlap_e = NULL,
#' @param specific_QTN_number = NULL,
#' @param specific_e_QTN_number = NULL,
#' @param ntraits = NULL
#' @param constrains = list(maf_above = NULL, maf_below = NULL)
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param export_gt = FALSE
#' @return Genotype of selected SNPs
#' @author Alex lipka and Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'----------------------------- QTN_partially_pleiotropic ----------------------
QTN_partially_pleiotropic <-
  function(genotypes = NULL,
           seed = NULL,
           overlap = NULL,
           overlap_e = NULL,
           specific_QTN_number = NULL,
           specific_e_QTN_number = NULL,
           ntraits = NULL,
           constrains = list(maf_above = NULL,
                             maf_below = NULL),
           rep = NULL,
           rep_by = NULL,
           export_gt = NULL) {
    #---------------------------------------------------------------------------
    if (any(lengths(constrains)>0)) { 
      index <- constrain(genotypes = genotypes, 
                         maf_above = constrains$maf_above,
                         maf_below = constrains$maf_below)
    } else {
      # First SNP at column 6
      index <- 6:nrow(genotypes)
    }
    if (rep_by != 'QTN'){rep <- 1}
      add_pleiotropic_QTN_genotypic_info <- vector("list", rep)
      add_specific_QTN_genotypic_info <- vector("list", rep)
      epi_pleiotropic_QTN_genotypic_info <- vector("list", rep)
      epi_specific_QTN_genotypic_info <- vector("list", rep)
      for (j in 1:rep) {
        if (!is.null(seed)) {
          set.seed(seed + j)
        }
        vector_of_pleiotropic_add_QTN <-
          sample(index, overlap, replace = FALSE)
        add_pleiotropic_QTN_genotypic_info[[j]] <-
          as.data.frame(genotypes[vector_of_pleiotropic_add_QTN, ])
        snps <-
          setdiff(index, vector_of_pleiotropic_add_QTN)
        vector_of_specific_add_QTN_temp <- vector("list", ntraits)
        add_specific_QTN_genotypic_info_temp <- vector("list", ntraits)
        ss <- c()
        for (i in 1:ntraits) {
          if (!is.null(seed)) {
            ss[i] <- seed + i + j
            set.seed(seed + i + j)
          }
          vector_of_specific_add_QTN_temp[[i]] <-
            sample(snps, specific_QTN_number[i], replace = FALSE)
          snps <- setdiff(snps, vector_of_specific_add_QTN_temp[[i]])
          add_specific_QTN_genotypic_info_temp[[i]] <-
            as.data.frame(genotypes[vector_of_specific_add_QTN_temp[[i]], ])
        }
        add_specific_QTN_genotypic_info_temp <- 
          do.call(rbind, add_specific_QTN_genotypic_info_temp)
        add_specific_QTN_genotypic_info[[j]] <- 
          data.frame(trait = paste0("trait_",rep(1:ntraits, specific_QTN_number)),
                     add_specific_QTN_genotypic_info_temp)
        if (!is.null(overlap_e) | !is.null(specific_e_QTN_number)) {
          if (!is.null(seed)) {
            set.seed(seed + seed + j)
          }
          vector_of_pleiotropic_epi_QTN <-
            sample(index, (2 * overlap_e), replace = FALSE)
          epi_pleiotropic_QTN_genotypic_info[[j]] <-
            as.data.frame(genotypes[vector_of_pleiotropic_epi_QTN, ])
          snps_e <-
            setdiff(index, vector_of_pleiotropic_epi_QTN)
          vector_of_specific_epi_QTN_temp <- vector("list", ntraits)
          epi_specific_QTN_genotypic_info_temp <- vector("list", ntraits)
          sse <- c()
          for (i in 1:ntraits) {
            if (!is.null(seed)) {
              sse[i] <- (seed + i) + seed + j
              set.seed( (seed + i) + seed + j)
            }
            vector_of_specific_epi_QTN_temp[[i]] <-
              sample(snps_e, (2 * specific_e_QTN_number[i]), replace = FALSE)
            snps_e <- setdiff(snps_e, vector_of_specific_epi_QTN_temp[[i]])
            epi_specific_QTN_genotypic_info_temp[[i]] <-
              as.data.frame(genotypes[vector_of_specific_epi_QTN_temp[[i]], ])
          }
          epi_specific_QTN_genotypic_info_temp <- 
            do.call(rbind, epi_specific_QTN_genotypic_info_temp)
          epi_specific_QTN_genotypic_info[[j]] <- 
            data.frame(trait = paste0("trait_",rep(1:ntraits, (2 * specific_e_QTN_number))),
                       epi_specific_QTN_genotypic_info_temp)
        }
      }
      add_object <- mapply(function(x,y) {
        p <- split(y, y[,1])
        names(p) <- NULL
        lapply(p, function(z) {
          x <- data.frame(type= "Pleiotropic", trait = unique(z[,1]), x)
          z <- data.frame(type= "Specific", z)
          rbind(x,z)
        }
        )
      },
      x=add_pleiotropic_QTN_genotypic_info,
      y=add_specific_QTN_genotypic_info,
      SIMPLIFY = F)
      additive_effect_trait_object = add_object
      add_object <- unlist(add_object,recursive=FALSE)
      add_object <- do.call(rbind, add_object)
      add_object <-
        data.frame(rep = sort(c(rep(1:rep, each = overlap*ntraits),
                                rep(1:rep, each = sum(specific_QTN_number)))),
                   add_object)
      additive_effect_trait_object <-  
        lapply(additive_effect_trait_object, function(x){
          lapply(x, function(b){
            rownames(b) <-
              paste0("Chr_",  b$chr, "_", b$pos)
            b <- b[,-(1:7)]
            return(t(b))
          })
        })
      if(!export_gt){
        add_object <- add_object[,1:7]
      }
      if (!is.null(overlap_e)) {
        epi_object <- mapply(function(x,y) {
          p <- split(y, y[,1])
          names(p) <- NULL
          lapply(p, function(z) {
            x <- data.frame(type= "Pleiotropic", trait = unique(z[,1]), x)
            z <- data.frame(type= "Specific", z)
            rbind(x,z)
          }
          )
        },
        x=epi_pleiotropic_QTN_genotypic_info,
        y=epi_specific_QTN_genotypic_info,
        SIMPLIFY = F)
        epistatic_effect_trait_object <- epi_object
        epi_object <- unlist(epi_object,recursive=FALSE)
        epi_object <- do.call(rbind, epi_object)
        epi_object <-
          data.frame(rep = sort(c(rep(1:rep, each = overlap_e*ntraits*2),
                                  rep(1:rep, each = sum(specific_e_QTN_number)*2))),
                     epi_object)
        epistatic_effect_trait_object <-  
          lapply(epistatic_effect_trait_object, function(x){
            lapply(x, function(b){
              rownames(b) <-
                paste0("Chr_",  b$chr, "_", b$pos)
              b <- b[,-(1:7)]
              return(t(b))
            })
          })
        if(!export_gt){
          epi_object <- epi_object[,1:7]
        }
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",
            paste0(specific_QTN_number + overlap, collapse = "_"),
            "_Add_QTN",
            ".txt"
          ),
          row.names = FALSE,
          col.names = FALSE,
          sep = "\t",
          quote = FALSE
        )
        write.table(
          c(seed + seed + 1:rep),
          paste0(
            "seed_number_for_",
            paste0(specific_e_QTN_number + overlap_e, collapse = "_"),
            "_Epi_QTN",
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
        data.table::fwrite(
          epi_object,
          "Epistatic_selected_QTNs.txt",
          row.names = FALSE,
          sep = "\t",
          quote = FALSE,
          na = NA
        )
        return(
          list(
            additive_effect_trait_object = additive_effect_trait_object,
            epistatic_effect_trait_object = epistatic_effect_trait_object
          )
        )
      } else {
        write.table(
          c(seed + 1:rep),
          paste0(
            "seed_number_for_",
            paste0(specific_QTN_number + overlap, collapse = "_"),
            "_Add_QTN",
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
        return(list(additive_effect_trait_object = additive_effect_trait_object))
      }
  }
