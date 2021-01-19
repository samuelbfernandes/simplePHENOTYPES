#' Calculate genetic value based on QTN objects.
#' @keywords internal
#' @param add_obj hhh
#' @param dom_obj hhh
#' @param epi_obj hhh
#' @param add_effect hhh
#' @param dom_effect = NULL,
#' @param epi_effect kkkk
#' @param epi_interaction = NULL,
#' @param ntraits hhh
#' @param cor hhh
#' @param architecture hhh
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param add = NULL,
#' @param dom = NULL,
#' @param epi = NULL
#' @param sim_method = NULL,
#' @param verbose = TRUE
#' @return A matrix of Genetic values for multiple traits
#' @author Samuel Fernandes
#' Last update: Apr 20, 2020
#'
#'-------------------------------base_line_multi_traits-------------------------
base_line_multi_traits <-
  function(add_obj = NULL,
           dom_obj = NULL,
           epi_obj = NULL,
           add_effect = NULL,
           dom_effect = NULL,
           epi_effect = NULL,
           epi_interaction = NULL,
           ntraits = NULL,
           cor = NULL,
           architecture = NULL,
           rep = NULL,
           rep_by = NULL,
           add = NULL,
           dom = NULL,
           epi = NULL,
           sim_method = NULL,
           verbose = TRUE) {
    #'--------------------------------------------------------------------------
    traits <- NULL
    VA <- NULL
    VD <- NULL
    VE <- NULL
    sample_cor <- NULL
    QTN_var <- list(var_add = list(),
                    var_dom = list(),
                    var_epi = list())
    if (rep_by != "QTN") {
      rep <- 1
    }
    results <- vector("list", rep)
    for (z in 1:rep) {
      if (!is.null(cor) & architecture != "LD") {
        if (architecture == "pleiotropic") {
          if (add) {
            genetic_value <-
              matrix(NA, nrow(add_obj[[z]]), ncol = ntraits)
            rownames <- rownames(add_obj[[z]])
          } else if (dom) {
            genetic_value <-
              matrix(NA, nrow(dom_obj[[z]]), ncol = ntraits)
            rownames <- rownames(dom_obj[[z]])
          } else {
            genetic_value <-
              matrix(NA, nrow(epi_obj[[z]]), ncol = ntraits)
            rownames <- rownames(epi_obj[[z]])
          }
          VA <- c()
          VE <- c()
          VD <- c()
          for (j in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                add_obj = add_obj[[z]],
                dom_obj = dom_obj[[z]],
                epi_obj = epi_obj[[z]],
                add_effect = add_effect[[j]],
                dom_effect = dom_effect[[j]],
                epi_effect = epi_effect[[j]],
                epi_interaction = epi_interaction,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            genetic_value[, j] <-
              trait_temp$base_line[[1]]
            if (add) {
              VA[j] <- trait_temp$VA
              QTN_var$var_add[[j]] <- trait_temp$var_add
            }
            if (dom) {
              VD[j] <- trait_temp$VD
              QTN_var$var_dom[[j]] <- trait_temp$var_dom
            }
            if (epi) {
              VE[j] <- trait_temp$VE
              QTN_var$var_epi[[j]] <- trait_temp$var_epi
            }
          }
          sdg <- apply(genetic_value, 2, sd)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          cg <- cov(genetic_s)
          cg <- make_pd(cg, verbose = verbose)
          L <- t(chol(cg))
          G_white <- t(solve(L) %*% t(genetic_s))
          cor <- make_pd(cor, verbose = verbose)
          L <- t(chol(cor))
          traits <- t(L %*% t(G_white))
          rownames(traits) <- rownames
          cor_original_trait <- c()
          for (i in seq_len(ncol(traits))) {
            traits[, i] <-
              traits[, i] * sdg[i] + meang[i]
            cor_original_trait[i] <-
              cor(traits[, i], genetic_value[, i])
          }
          sample_cor <- cor(traits)
          results[[z]] <- list(
            base_line = traits,
            VA = VA,
            VE = VE,
            VD = VD,
            sample_cor = sample_cor,
            QTN_var = QTN_var
          )
        } else {
          if (add) {
            genetic_value <-
              matrix(NA, nrow(add_obj[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(add_obj[[z]][[1]])
          } else if (dom) {
            genetic_value <-
              matrix(NA, nrow(dom_obj[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(dom_obj[[z]][[1]])
          } else {
            genetic_value <-
              matrix(NA, nrow(epi_obj[[z]][[1]]), ncol = ntraits)
            rownames <- rownames(epi_obj[[z]][[1]])
          }
          VA <- c()
          VD <- c()
          VE <- c()
          for (j in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                add_obj = add_obj[[z]][[j]],
                dom_obj = dom_obj[[z]][[j]],
                epi_obj = epi_obj[[z]][[j]],
                add_effect = add_effect[[j]],
                dom_effect = dom_effect[[j]],
                epi_effect = epi_effect[[j]],
                epi_interaction = epi_interaction,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            genetic_value[, j] <-
              trait_temp$base_line[[1]]
            if (add) {
              VA[j] <- trait_temp$VA
              QTN_var$var_add[[j]] <- trait_temp$var_add
            }
            if (dom) {
              VD[j] <- trait_temp$VD
              QTN_var$var_dom[[j]] <- trait_temp$var_dom
            }
            if (epi) {
              VE[j] <- trait_temp$VE
              QTN_var$var_epi[[j]] <- trait_temp$var_epi
            }
          }
          sdg <- apply(genetic_value, 2, sd)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          cg <- cov(genetic_s)
          cg <- make_pd(cg, verbose = verbose)
          L <- t(chol(cg))
          G_white <- t(solve(L) %*% t(genetic_s))
          cor <- make_pd(cor, verbose = verbose)
          L <- t(chol(cor))
          traits <- t(L %*% t(G_white))
          rownames(traits) <- rownames
          cor_original_trait <- c()
          for (i in 1:ntraits) {
            traits[, i] <- (traits[, i] * sdg[i]) + meang[i]
            cor_original_trait[i] <-
              cor(traits[, i], genetic_value[, i])
          }
          sample_cor <- cor(traits)
          results[[z]] <-
            list(
              base_line = traits,
              VA = VA,
              VD = VD,
              VE = VE,
              sample_cor = sample_cor,
              cor_original_trait = cor_original_trait,
              QTN_var = QTN_var
            )
        }
      } else {
        VA <- c()
        VE <- c()
        VD <- c()
        if (architecture == "pleiotropic") {
          if (add) {
            traits <- matrix(NA, nrow(add_obj[[z]]), ncol = ntraits)
            rownames <- rownames(add_obj[[1]])
          } else if (dom) {
            traits <- matrix(NA, nrow(dom_obj[[z]]), ncol = ntraits)
            rownames <- rownames(dom_obj[[1]])
          } else {
            traits <- matrix(NA, nrow(epi_obj[[z]]), ncol = ntraits)
            rownames <- rownames(epi_obj[[1]])
          }
          for (i in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                add_obj = add_obj[[z]],
                dom_obj = dom_obj[[z]],
                epi_obj = epi_obj[[z]],
                add_effect = add_effect[[i]],
                dom_effect = dom_effect[[i]],
                epi_effect = epi_effect[[i]],
                epi_interaction = epi_interaction,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            traits[, i] <- trait_temp$base_line[, 1]
            if (add) {
              VA[i] <- trait_temp$VA
              QTN_var$var_add[[i]] <- trait_temp$var_add
            }
            if (dom) {
              VD[i] <- trait_temp$VD
              QTN_var$var_dom[[i]] <- trait_temp$var_dom
            }
            if (epi) {
              VE[i] <- trait_temp$VE
              QTN_var$var_epi[[i]] <- trait_temp$var_epi
            }
          }
          rownames(traits) <- rownames
        } else {
          if (add) {
            traits <- matrix(NA,
                             nrow(add_obj[[z]][[1]]),
                             ncol = ntraits)
            rownames <- rownames(add_obj[[z]][[1]])
          } else if (dom) {
            traits <- matrix(NA,
                             nrow(dom_obj[[z]][[1]]),
                             ncol = ntraits)
            rownames <- rownames(dom_obj[[z]][[1]])
          } else {
            traits <- matrix(NA,
                             nrow(epi_obj[[z]][[1]]),
                             ncol = ntraits)
            rownames <- rownames(epi_obj[[z]][[1]])
          }
          for (i in 1:ntraits) {
            trait_temp <-
              base_line_single_trait(
                add_obj = add_obj[[z]][[i]],
                dom_obj = dom_obj[[z]][[i]],
                epi_obj = epi_obj[[z]][[i]],
                add_effect = add_effect[[i]],
                dom_effect = dom_effect[[i]],
                epi_effect = epi_effect[[i]],
                epi_interaction = epi_interaction,
                ntraits = ntraits,
                add = add,
                dom = dom,
                epi = epi,
                sim_method = sim_method
              )
            traits[, i] <- trait_temp$base_line[, 1]
            if (add) {
              VA[i] <- trait_temp$VA
              QTN_var$var_add[[i]] <- trait_temp$var_add
            }
            if (dom) {
              VD[i] <- trait_temp$VD
              QTN_var$var_dom[[i]] <- trait_temp$var_dom
            }
            if (epi) {
              VE[i] <- trait_temp$VE
              QTN_var$var_epi[[i]] <- trait_temp$var_epi
            }
          }
          rownames(traits) <- rownames
        }
        if (!all(traits == 0)) {
          sample_cor <- cor(traits)
        }
        results[[z]] <- list(
          base_line = traits,
          VA = VA,
          VD = VD,
          VE = VE,
          sample_cor = sample_cor,
          QTN_var = QTN_var
        )
      }
    }
    return(results)
  }
