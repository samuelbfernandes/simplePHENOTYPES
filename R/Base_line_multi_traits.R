#' Calculate genetic value based on QTN objects.
#' @export
#' @param additive.object hhh
#' @param epistatic.object hhh
#' @param additive.effect hhh
#' @param epistatic.effect kkkk
#' @param big.additive.QTN.effect jjj
#' @param seed hhh
#' @param set.cor hhh
#' @param ntraits hhh
#' @param correlation hhh
#' @param model hhh
#' @return A matrix of Genetic values for multiple traits
#' @author Samuel Fernandes


#' Last update: Jul 22, 2019
#'
#'--------------------------------Base_line_multi_traits-------------------------
Base_line_multi_traits <-
  function(additive.object = NULL,
           epistatic.object = NULL,
           additive.effect = NULL,
           epistatic.effect = NULL,
           big.additive.QTN.effect = NULL,
           seed = seed,
           set.cor = TRUE,
           ntraits = NULL,
           correlation = NULL,
           model = "pleiotropic") {
    #'---------------------------------------------------------------------------
    if (set.cor) {
      if (model == "pleiotropic") {
        Genetic_value <-
          Base_line_single_trait(
            additive.object = additive.object[[1]],
            epistatic.object = epistatic.object[[1]],
            additive.effect = additive.effect[1],
            epistatic.effect = epistatic.effect[1],
            big.additive.QTN.effect = big.additive.QTN.effect[1],
            seed = seed
          )
        
        n <- nrow(Genetic_value$base.line)
        varg <- var(Genetic_value$base.line)
        traits <- matrix(NA, n, ncol = ntraits)
        traits[, 1] <- Genetic_value$base.line$Additive.effect
        
        for (i in 2:ntraits) {
          traits[, i] <-
            Genetic_value$base.line$Additive.effect + rnorm(n, 0, sqrt(varg * 0.01))
        }
        
        Genetic_s <- apply(traits, 2, scale)
        
        #' whiting transformation
        L <- t(chol(cov(Genetic_s)))
        G_white <- t(solve(L) %*% t(Genetic_s))
        
        #' coloring transformation
        #' approximation to make it positive definite
        if (!lqmm::is.positive.definite(correlation)) {
          cat("Correlation matrix not positive definite! Using make.positive.definite \n")
          correlation <-
            lqmm::make.positive.definite(correlation)
        }
        L <- t(chol(correlation))
        traits <- t(L %*% t(G_white))
        
        rownames(traits) <- rownames(additive.object[[1]])
        
        #' bring it back to the original scale
        cor_original_trait <- c()
        for (i in 1:ncol(traits)) {
          traits[, i] <-
            (traits[, i] * sd(Genetic_value$base.line$Additive.effect)) + mean(Genetic_value$base.line$Additive.effect)
          
          cor_original_trait[i] <-
            cor(traits[, i], Genetic_value$base.line$Additive.effect)
        }
        sample.cor <- cor(traits)
        
        results <- if (!is.null(epistatic.object)) {
          list(
            base.line = traits,
            VA = c(Genetic_value$VA),
            VE = c(Genetic_value$VE),
            sample.cor = sample.cor
          )
        } else {
          list(
            base.line = traits,
            VA = c(Genetic_value$VA),
            sample.cor = sample.cor
          )
        }
        return(results)
        
      } else {
        Genetic_value <-
          matrix(NA, nrow(additive.object[[1]]), ncol = ntraits)
        VA <- c()
        VE <- c()
        for (j in 1:ntraits) {
          trait_temp <-
            Base_line_single_trait(
              additive.object = additive.object[[j]],
              epistatic.object = epistatic.object[[j]],
              additive.effect = additive.effect[j],
              epistatic.effect = epistatic.effect[j],
              big.additive.QTN.effect = big.additive.QTN.effect[j],
              seed = seed
            )
          Genetic_value[, j] <-
            trait_temp$base.line$Additive.effect
          VA[j] <- trait_temp$VA
          if (!is.null(epistatic.object)) {
            VE[j] <- trait_temp$VE
          }
        }
        sdg <- apply(Genetic_value, 2, sd)
        meang <- apply(Genetic_value, 2, mean)
        Genetic_s <- apply(Genetic_value, 2, scale)
        
        #' whiting transformation
        L <- t(chol(cov(Genetic_s)))
        G_white <- t(solve(L) %*% t(Genetic_s))
        
        #' coloring transformation
        if (!lqmm::is.positive.definite(correlation)) {
          cat(
            "Correlation matrix not positive definite! \nUsing lqmm::make.positive.definite \n"
          )
          correlation <-
            lqmm::make.positive.definite(correlation)
        }
        L <- t(chol(correlation))
        traits <- t(L %*% t(G_white))
        
        rownames(traits) <- rownames(additive.object[[1]])
        #' bring it back to the original scale
        cor_original_trait <- c()
        for (i in 1:ncol(traits)) {
          traits[, i] <- (traits[, i] * sdg[i]) + meang[i]
          cor_original_trait[i] <-
            cor(traits[, i], Genetic_value[, i])
        }
        
        sample.cor <- cor(traits)
        results <- if (!is.null(epistatic.object)) {
          list(
            base.line = traits,
            VA = VA,
            VE = VE,
            sample.cor = sample.cor,
            cor_original_trait = cor_original_trait
          )
        } else {
          list(
            base.line = traits,
            VA = VA,
            sample.cor = sample.cor,
            cor_original_trait = cor_original_trait
          )
        }
        return(results)
      }
      
    } else {
      traits <- matrix(NA, nrow(additive.object[[1]]), ncol = ntraits)
      VA <- c()
      VE <- c()
      for (i in 1:ntraits) {
        ifelse(model == "pleiotropic", j <- 1, j <- i)
        trait_temp <-
          Base_line_single_trait(
            additive.object = additive.object[[j]],
            epistatic.object = epistatic.object[[j]],
            additive.effect = additive.effect[i],
            epistatic.effect = epistatic.effect[i],
            big.additive.QTN.effect = big.additive.QTN.effect[i],
            seed = seed
          )
        traits[, i] <- trait_temp$base.line[, 1]
        VA[i] <- trait_temp$VA
        if (!is.null(epistatic.object)) {
          VE[i] <- trait_temp$VE
        }
      }
      rownames(traits) <- rownames(additive.object[[1]])
      sample.cor <- cor(traits)
      
      results <- if (!is.null(epistatic.object)) {
        list(
          base.line = traits,
          VA = VA,
          VE = VE,
          sample.cor = sample.cor
        )
      } else {
        list(base.line = traits,
             VA = VA,
             sample.cor = sample.cor)
      }
      return(results)
    }
  }
