#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive_object hhh
#' @param dominance_object hhh
#' @param epistatic_object hhh
#' @param additive_effect hhh
#' @param dominance_effect = NULL,
#' @param epistatic_effect kkkk
#' @param big_additive_QTN_effect jjj
#' @param big_dominance_QTN_effect = null
#' @param seed hhh
#' @param ntraits hhh
#' @param correlation hhh
#' @param architecture hhh
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @param sim_method = NULL,
#' @return A matrix of Genetic values for multiple traits
#' @author Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'-------------------------------base_line_multi_traits-------------------------
base_line_multi_traits <-
  function(additive_object = NULL,
           dominance_object = NULL,
           epistatic_object = NULL,
           additive_effect = NULL,
           dominance_effect = NULL,
           epistatic_effect = NULL,
           big_additive_QTN_effect = NULL,
           big_dominance_QTN_effect = NULL,
           seed = seed,
           ntraits = NULL,
           correlation = NULL,
           architecture = NULL,
           rep = NULL,
           rep_by = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL,
           sim_method = NULL) {
    #'---------------------------------------------------------------------------
    if (rep_by != 'QTN'){ rep <- 1}
    results <- vector("list", rep)
    for(z in 1:rep){
      if (!is.null(correlation) & architecture != "LD") {
        if (architecture == "pleiotropic") {
          if (add) {
            genetic_value <-
              matrix(NA, nrow(additive_object[[z]]), ncol = ntraits)
            rownames <- rownames(additive_object[[z]])
          } else if (dom) {
            genetic_value <-
              matrix(NA, nrow(dominance_object[[z]]), ncol = ntraits)
            rownames <- rownames(dominance_object[[z]])
          } else {
            genetic_value <-
              matrix(NA, nrow(epistatic_object[[z]]), ncol = ntraits)
            rownames <- rownames(epistatic_object[[z]])
          }
          VA <- c()
          VE <- c()
          VD <- c()
          for (j in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                additive_object = additive_object[[z]],
                dominance_object = dominance_object[[z]],
                epistatic_object = epistatic_object[[z]],
                additive_effect = additive_effect[,j],
                dominance_effect = dominance_effect[,j],
                epistatic_effect = epistatic_effect[,j],
                big_additive_QTN_effect = big_additive_QTN_effect[j],
                big_dominance_QTN_effect = big_dominance_QTN_effect[j],
                seed = seed,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            genetic_value[, j] <-
              trait_temp$base_line[[1]]
            if (add) VA[j] <- trait_temp$VA
            if (dom) VD[j] <- trait_temp$VD
            if (epi) VE[j] <- trait_temp$VE
          }
          sdg <- apply(genetic_value, 2, sd)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          #' whiting transformation
          L <- t(chol(cov(genetic_s)))
          G_white <- t(solve(L) %*% t(genetic_s))
          #' coloring transformation
          #' approximation to make it positive definite
          if (!lqmm::is.positive.definite(correlation)) {
            cat("Correlation matrix not positive definite! Using make.positive.definite \n")
            correlation <-
              lqmm::make.positive.definite(correlation)
          }
          L <- t(chol(correlation))
          traits <- t(L %*% t(G_white))
          rownames(traits) <- rownames
          #' bring it back to the original scale
          cor_original_trait <- c()
          for (i in 1:ncol(traits)) {
            traits[, i] <-
              traits[, i] * sdg[i] + meang[i]
            cor_original_trait[i] <-
              cor(traits[, i], genetic_value[,i])
          }
          sample_cor <- cor(traits)
          results[[z]] <- 
            if (add == TRUE & dom == TRUE & epi == TRUE) {
              list(
                base_line = traits,
                VA = VA,
                VE = VE,
                VD = VD,
                sample_cor = sample_cor
              )
            } else if (add == TRUE & epi == TRUE){
              list(
                base_line = traits,
                VA = VA,
                VE = VE,
                sample_cor = sample_cor
              )
            } else if (add == TRUE & dom == TRUE){
              list(
                base_line = traits,
                VA = VA,
                VD = VD,
                sample_cor = sample_cor
              )
            } else if (dom == TRUE & epi == TRUE){
              list(
                base_line = traits,
                VE = VE,
                VD = VD,
                sample_cor = sample_cor
              )
            } else if (add == TRUE){
              list(
                base_line = traits,
                VA = VA,
                sample_cor = sample_cor
              )
            } else if (dom == TRUE){
              list(
                base_line = traits,
                VD = VD,
                sample_cor = sample_cor
              )
            } else if (epi == TRUE){
              list(
                base_line = traits,
                VE = VE,
                sample_cor = sample_cor
              )
            }
        } else {
          if (add) {
            genetic_value <-
              matrix(NA, nrow(additive_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(additive_object[[z]][[1]])
          } else if (dom) {
            genetic_value <-
              matrix(NA, nrow(dominance_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(dominance_object[[z]][[1]])
          } else {
            genetic_value <-
              matrix(NA, nrow(epistatic_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(epistatic_object[[z]][[1]])
          }
          VA <- c()
          VD <- c()
          VE <- c()
          for (j in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                additive_object = additive_object[[z]][[j]],
                dominance_object = dominance_object[[z]][[j]],
                epistatic_object = epistatic_object[[z]][[j]],
                additive_effect = additive_effect[,j],
                dominance_effect = dominance_effect[,j],
                epistatic_effect = epistatic_effect[,j],
                big_additive_QTN_effect = big_additive_QTN_effect[j],
                big_dominance_QTN_effect = big_dominance_QTN_effect[j],
                seed = seed,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            genetic_value[, j] <-
              trait_temp$base_line[[1]]
            if (add) VA[j] <- trait_temp$VA
            if (dom) VD[j] <- trait_temp$VD
            if (epi) VE[j] <- trait_temp$VE
          }
          sdg <- apply(genetic_value, 2, sd)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          #' whiting transformation
          L <- t(chol(cov(genetic_s)))
          G_white <- t(solve(L) %*% t(genetic_s))
          #' coloring transformation
          if (!lqmm::is.positive.definite(correlation)) {
            cat(
              "Correlation matrix not positive definite!
              \nUsing lqmm::make.positive.definite \n"
            )
            correlation <-
              lqmm::make.positive.definite(correlation)
          }
          L <- t(chol(correlation))
          traits <- t(L %*% t(G_white))
          rownames(traits) <- rownames
          #' bring it back to the original scale
          cor_original_trait <- c()
          for (i in 1:ntraits) {
            traits[, i] <- (traits[, i] * sdg[i]) + meang[i]
            cor_original_trait[i] <-
              cor(traits[, i], genetic_value[, i])
          }
          sample_cor <- cor(traits)
          results[[z]] <- 
            if (add == TRUE & dom == TRUE & epi == TRUE) {
              list(
                base_line = traits,
                VA = VA,
                VD = VD,
                VE = VE,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (add == TRUE & epi == TRUE){
              list(
                base_line = traits,
                VA = VA,
                VE = VE,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (add == TRUE & dom == TRUE){
              list(
                base_line = traits,
                VA = VA,
                VD = VD,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (dom == TRUE & epi == TRUE){
              list(
                base_line = traits,
                VD = VD,
                VE = VE,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (add == TRUE){
              list(
                base_line = traits,
                VA = VA,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (dom == TRUE){
              list(
                base_line = traits,
                VD = VD,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            } else if (epi == TRUE){
              list(
                base_line = traits,
                VD = VD,
                sample_cor = sample_cor,
                cor_original_trait = cor_original_trait
              )
            }
        }
      } else {
        VA <- c()
        VE <- c()
        VD <- c()
        if (architecture == "pleiotropic") {
          if (add) {
            traits <- matrix(NA, nrow(additive_object[[z]]), ncol = ntraits)
            rownames <- rownames(additive_object[[1]])
          } else if (dom) {
            traits <- matrix(NA, nrow(dominance_object[[z]]), ncol = ntraits)
            rownames <- rownames(additive_object[[1]])
          } else {
            traits <- matrix(NA, nrow(epistatic_object[[z]]), ncol = ntraits)
            rownames <- rownames(additive_object[[1]])
          }
          for (i in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                additive_object = additive_object[[z]],
                dominance_object = dominance_object[[z]],
                epistatic_object = epistatic_object[[z]],
                additive_effect = additive_effect[,i],
                dominance_effect = dominance_effect[,i],
                epistatic_effect = epistatic_effect[,i],
                big_additive_QTN_effect = big_additive_QTN_effect[i],
                big_dominance_QTN_effect = big_dominance_QTN_effect[i],
                seed = seed,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            traits[, i] <- trait_temp$base_line[, 1]
            if (add) VA[i] <- trait_temp$VA
            if (dom) VD[i] <- trait_temp$VD
            if (epi) VE[i] <- trait_temp$VE
          }
          rownames(traits) <- rownames
        } else {
          if (add) {
            traits <- matrix(NA, nrow(additive_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(additive_object[[z]][[1]]) 
          } else if (dom) {
            traits <- matrix(NA, nrow(dominance_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(dominance_object[[z]][[1]]) 
          } else {
            traits <- matrix(NA, nrow(epistatic_object[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(epistatic_object[[z]][[1]]) 
          }
          for (i in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                additive_object = additive_object[[z]][[i]],
                dominance_object = dominance_object[[z]][[i]],
                epistatic_object = epistatic_object[[z]][[i]],
                additive_effect = additive_effect[,i],
                dominance_effect = dominance_effect[,i],
                epistatic_effect = epistatic_effect[,i],
                big_additive_QTN_effect = big_additive_QTN_effect[i],
                big_dominance_QTN_effect = big_dominance_QTN_effect[i],
                seed = seed,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            traits[, i] <- trait_temp$base_line[, 1]
            if (add) VA[i] <- trait_temp$VA
            if (dom) VD[i] <- trait_temp$VD
            if (epi) VE[i] <- trait_temp$VE
          }
          rownames(traits) <- rownames
        }
        sample_cor <- cor(traits)
        results[[z]] <- 
          if (add == TRUE & dom == TRUE & epi == TRUE) {
            list(
              base_line = traits,
              VA = VA,
              VD = VD,
              VE = VE,
              sample_cor = sample_cor
            )
          } else if (add == TRUE & epi == TRUE){
            list(
              base_line = traits,
              VA = VA,
              VE = VE,
              sample_cor = sample_cor
            )
          } else if (add == TRUE & dom == TRUE){
            list(
              base_line = traits,
              VA = VA,
              VD = VD,
              sample_cor = sample_cor
            )
          } else if (dom == TRUE & epi == TRUE){
            list(
              base_line = traits,
              VD = VD,
              VE = VE,
              sample_cor = sample_cor
            )
          } else if (add == TRUE){
            list(base_line = traits,
                 VA = VA,
                 sample_cor = sample_cor)
          } else if (dom == TRUE){
            list(
              base_line = traits,
              VD = VD,
              sample_cor = sample_cor
            )
          } else if (epi == TRUE){
            list(
              base_line = traits,
              VE = VE,
              sample_cor = sample_cor
            )
          }
      }
    }
    return(results)
  }
