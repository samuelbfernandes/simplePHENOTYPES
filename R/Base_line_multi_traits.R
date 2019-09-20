#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive_object hhh
#' @param epistatic_object hhh
#' @param additive_effect hhh
#' @param epistatic_effect kkkk
#' @param big_additive_QTN_effect jjj
#' @param seed hhh
#' @param set_cor hhh
#' @param ntraits hhh
#' @param correlation hhh
#' @param model hhh
#' @return A matrix of Genetic values for multiple traits
#' @author Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'-------------------------------base_line_multi_traits-------------------------
base_line_multi_traits <-
  function(additive_object = NULL,
           epistatic_object = NULL,
           additive_effect = NULL,
           epistatic_effect = NULL,
           big_additive_QTN_effect = NULL,
           seed = seed,
           set_cor = TRUE,
           ntraits = NULL,
           correlation = NULL,
           model = "pleiotropic") {
  #'---------------------------------------------------------------------------
    if (set_cor) {
      if (model == "pleiotropic") {
        genetic_value <-
          base_line_single_trait(
            additive_object = additive_object[[1]],
            epistatic_object = epistatic_object[[1]],
            additive_effect = additive_effect[1],
            epistatic_effect = epistatic_effect[1],
            big_additive_QTN_effect = big_additive_QTN_effect[1],
            seed = seed
          )
        n <- nrow(genetic_value$base_line)
        varg <- var(genetic_value$base_line)
        traits <- matrix(NA, n, ncol = ntraits)
        traits[, 1] <- genetic_value$base_line$additive_effect
        for (i in 2:ntraits) {
          traits[, i] <-
            genetic_value$base_line$additive_effect +
            rnorm(n, 0, sqrt(varg * 0.01))
        }
        genetic_s <- apply(traits, 2, scale)
        #' whiting transformation
        L <- t(chol(cov(genetic_s)))
        G_white <- t(solve(L) %*% t(genetic_s))
        #' coloring transformation
        #' approximation to make it positive definite
        if (!lqmm::is.positive.definite(correlation)) {
          cat("Correlation matrix not positive definite!
              Using make.positive.definite \n")
          correlation <-
            lqmm::make.positive.definite(correlation)
        }
        L <- t(chol(correlation))
        traits <- t(L %*% t(G_white))
        rownames(traits) <- rownames(additive_object[[1]])
        #' bring it back to the original scale
        cor_original_trait <- c()
        for (i in 1:ncol(traits)) {
          traits[, i] <-
            (traits[, i] * sd(genetic_value$base_line$additive_effect)) +
            mean(genetic_value$base_line$additive_effect)
          cor_original_trait[i] <-
            cor(traits[, i], genetic_value$base_line$additive_effect)
        }
        sample_cor <- cor(traits)
        results <- if (!is.null(epistatic_object)) {
          list(
            base_line = traits,
            VA = c(genetic_value$VA),
            VE = c(genetic_value$VE),
            sample_cor = sample_cor
          )
        } else {
          list(
            base_line = traits,
            VA = c(genetic_value$VA),
            sample_cor = sample_cor
          )
        }
        return(results)
      } else {
        genetic_value <-
          matrix(NA, nrow(additive_object[[1]]), ncol = ntraits)
        VA <- c()
        VE <- c()
        for (j in 1:ntraits) {
          trait_temp <-
            base_line_single_trait(
              additive_object = additive_object[[j]],
              epistatic_object = epistatic_object[[j]],
              additive_effect = additive_effect[j],
              epistatic_effect = epistatic_effect[j],
              big_additive_QTN_effect = big_additive_QTN_effect[j],
              seed = seed
            )
          genetic_value[, j] <-
            trait_temp$base_line$additive_effect
          VA[j] <- trait_temp$VA
          if (!is.null(epistatic_object)) {
            VE[j] <- trait_temp$VE
          }
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
        rownames(traits) <- rownames(additive_object[[1]])
        #' bring it back to the original scale
        cor_original_trait <- c()
        for (i in 1:ncol(traits)) {
          traits[, i] <- (traits[, i] * sdg[i]) + meang[i]
          cor_original_trait[i] <-
            cor(traits[, i], genetic_value[, i])
        }
        sample_cor <- cor(traits)
        results <- if (!is.null(epistatic_object)) {
          list(
            base_line = traits,
            VA = VA,
            VE = VE,
            sample_cor = sample_cor,
            cor_original_trait = cor_original_trait
          )
        } else {
          list(
            base_line = traits,
            VA = VA,
            sample_cor = sample_cor,
            cor_original_trait = cor_original_trait
          )
        }
        return(results)
      }
    } else {
      traits <- matrix(NA, nrow(additive_object[[1]]), ncol = ntraits)
      VA <- c()
      VE <- c()
      for (i in 1:ntraits) {
        ifelse(model == "pleiotropic", j <- 1, j <- i)
        trait_temp <-
          base_line_single_trait(
            additive_object = additive_object[[j]],
            epistatic_object = epistatic_object[[j]],
            additive_effect = additive_effect[i],
            epistatic_effect = epistatic_effect[i],
            big_additive_QTN_effect = big_additive_QTN_effect[i],
            seed = seed
          )
        traits[, i] <- trait_temp$base_line[, 1]
        VA[i] <- trait_temp$VA
        if (!is.null(epistatic_object)) {
          VE[i] <- trait_temp$VE
        }
      }
      rownames(traits) <- rownames(additive_object[[1]])
      sample_cor <- cor(traits)
      results <- if (!is.null(epistatic_object)) {
        list(
          base_line = traits,
          VA = VA,
          VE = VE,
          sample_cor = sample_cor
        )
      } else {
        list(base_line = traits,
             VA = VA,
             sample_cor = sample_cor)
      }
      return(results)
    }
  }