#' Calculate genetic value based on QTN objects.
#'
#' Frozen v1 engine, multi-trait genetic values.
#'
#' @details
#' \strong{How `cor` is imposed (v1 semantics, reproduced exactly).} The
#' per-trait genetic values \eqn{G} (one column per trait, built by
#' [genetic_effect()] from each trait's own effects) are standardized column
#' by column, whitened with the inverse Cholesky factor of their sample
#' covariance (after [make_pd()]), coloured with the Cholesky factor of `cor`
#' (after [make_pd()]) and finally rescaled with each trait's original
#' standard deviation and mean:
#' \code{T = scale(G) L_cg^{-T} L_cor^T * sd(G) + mean(G)}. For a
#' positive-definite `cor` with a unit diagonal this realizes the requested
#' correlation exactly in the sample and preserves each trait's variance.
#' Consequences to keep in mind:
#' \itemize{
#'   \item With a unit-diagonal `cor`, trait 1 is unchanged (its genetic value
#'     is exactly the input value) and only trait \eqn{k >= 2} is a mixture:
#'     a linear combination of the input values of traits \eqn{1, \dots, k}
#'     (it is regressed onto the earlier traits). QTNs that were declared
#'     specific to an earlier trait therefore leak into the later traits'
#'     genetic values once `cor` is supplied; the reverse does not happen.
#'   \item The per-trait effects written to the `*_QTNs.txt` files
#'     (`add_eff_t2`, ...) and the `VA`/`VD`/`VE` and per-QTN PVE files are the
#'     values \emph{before} the whitening/colouring step. For traits >= 2 they
#'     are not the effective effects and not the effective PVE; the residual
#'     variance and the reported heritability are computed from the
#'     transformed genetic values and therefore are correct.
#'   \item `cor` must be a full `ntraits x ntraits` matrix (a scalar is not
#'     accepted). A covariance-valued matrix (diagonal different from 1) is
#'     accepted and treated as covariance-like: the realized correlation is
#'     `cov2cor(cor)` and the genetic variance of trait \eqn{k} is multiplied
#'     by `cor[k, k]`, so trait 1 is \emph{rescaled} as well (for
#'     `cor = [[4, 1], [1, 1]]` the standard deviation of trait 1 is doubled,
#'     that of trait 2 is unchanged, and the realized correlation is 0.5).
#'   \item A non-positive-definite `cor` is not realized as requested: see
#'     [make_pd()] (the diagonal is raised, off-diagonals shrink).
#'     `create_phenotypes()` rejects such input.
#'   \item Traits whose genetic values are perfectly collinear (identical
#'     effect series, or a single QTN shared by all traits) or have zero
#'     variance cannot be whitened; this is reported with an informative error.
#'   \item `sample_cor` is computed per replicate (it is `NULL` for a
#'     replicate whose genetic values are all zero).
#' }
#'
#' @param add_obj additive QTN object
#' @param dom_obj dominance QTN object
#' @param epi_obj epistatic QTN object
#' @param add_effect additive effect sizes
#' @param dom_effect dominance effect sizes
#' @param epi_effect epistatic effect sizes
#' @param epi_interaction markers per epistatic interaction
#' @param ntraits number of traits
#' @param cor genetic correlation
#' @param architecture genetic architecture
#' @param rep = 1,
#' @param rep_by = 'QTN',
#' @param add additive flag
#' @param dom dominance flag
#' @param epi = NULL
#' @param sim_method effect-simulation method
#' @param verbose = TRUE
#' @return A matrix of Genetic values for multiple traits
#' @author Samuel Fernandes
#' @keywords internal
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
      # sample_cor is replicate-local: an all-zero replicate must not inherit
      # the previous replicate's matrix.
      sample_cor <- NULL
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
          .check_cor_inputs(genetic_value, sdg, cor, ntraits)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          cg <- cov(genetic_s)
          cg <- make_pd(cg, verbose = verbose,
                        what = "sample covariance of the standardized genetic values")
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
          .check_cor_inputs(genetic_value, sdg, cor, ntraits)
          meang <- apply(genetic_value, 2, mean)
          genetic_s <- apply(genetic_value, 2, scale)
          cg <- cov(genetic_s)
          cg <- make_pd(cg, verbose = verbose,
                        what = "sample covariance of the standardized genetic values")
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

#' Fail loudly on inputs the whitening/colouring step cannot handle
#' @keywords internal
#' @noRd
.check_cor_inputs <- function(genetic_value, sdg, cor, ntraits) {
  cm <- tryCatch(as.matrix(cor), error = function(e) NULL)
  if (is.null(cm) || !is.numeric(cm) || nrow(cm) != ntraits ||
      ncol(cm) != ntraits) {
    stop("`cor` must be a numeric ", ntraits, " x ", ntraits,
         " genetic correlation matrix (one row and column per trait); ",
         "a scalar or a matrix of a different size is not supported.",
         call. = FALSE)
  }
  bad <- which(!is.finite(sdg) | sdg == 0)
  if (length(bad) > 0L) {
    stop("Trait(s) ", paste(bad, collapse = ", "),
         " have zero (or undefined) genetic variance, so the genetic ",
         "correlation `cor` cannot be imposed. Check that the effects and ",
         "QTN numbers give every trait a non-zero genetic value.",
         call. = FALSE)
  }
  invisible(TRUE)
}
