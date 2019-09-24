#' Generate environmental effects based on a given heritability
#' @export
#' @param base_line_trait = NULL,
#' @param h2 = NULL,
#' @param rep = NULL,
#' @param seed = NULL,
#' @param ntraits = NULL,
#' @param h2_MT = NULL,
#' @param format = 'multi-file'
#' @param fam = NULL
#' @param to_r = FALSE
#' @return Phenotypes for ntraits traits
#' @author Alex lipka and Samuel Fernandes
#' Last update: Jul 22, 2019
#'
#'----------------------------phenotypes---------------------------------------
phenotypes <-
  function(base_line_trait = NULL,
           h2 = NULL,
           rep = NULL,
           seed = NULL,
           ntraits = NULL,
           h2_MT = NULL,
           format = "multi-file",
           fam = NULL,
           to_r = FALSE) {
    #---------------------------------------------------------------------------
    if (ntraits > 1) {
      if (format == "multi-file") {
        # Create a working directory for the output results:
        dir.create("Phenotypes")
        # Set the working directory
        setwd("./Phenotypes")
      }
      # For loop through the vector of heritabilities
      for (i in h2) {
        simulated_data <- vector("list", rep)
        ss <- c()
        # If heritability is zero
        if (i == 0) {
          # Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated_data[[z]] <-
              data.frame(rownames(base_line_trait),
                         matrix(NA, nrow(base_line_trait), ntraits),
                         rep = z)
            colnames(simulated_data[[z]]) <-
              c("<Taxa>",
                "Target",
                paste0("Trait_", 2:ntraits),
                "Rep")
            if (!is.null(seed)) {
              ss <- c(ss, seed + z)
              set.seed(seed + z)
            }
            # using multivariate normal with ntrait variances and
            # 0 covariances for independent residuals
            simulated_data[[z]][, 2:(ntraits + 1)] <-
              mvtnorm::rmvnorm(
                n = nrow(base_line_trait),
                mean = rep(0, ntraits),
                sigma = diag(1,
                             ntraits)
              )
          }
          if (format == "multi-file") {
            invisible(lapply(1:rep, function(x) {
              data.table::fwrite(
                simulated_data[[x]][- (ntraits + 2)],
                paste0("Simulated_Data_", "_Rep", x, "_Herit_", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            temp_simulated_data <- do.call(rbind, simulated_data)
            data.table::fwrite(
              temp_simulated_data,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else if (format == "gemma") {
            invisible(lapply(1:rep, function(x) {
              temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)], 
                                 by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", "_Rep", x, "_Herit_", i, ".fam"),
                row.names = FALSE,
                sep = "\t",
                col.names = FALSE,
                quote = FALSE,
                na = NA
              )
            }))
          } else {
            temp <- simulated_data[[1]][- (ntraits + 2)]
            for (j in 2:rep) {
              temp <- cbind(temp, simulated_data[[j]][-c(1, (ntraits + 2))])
            }
            data.table::fwrite(
              temp,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          # Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(
              format == "multi-file",
              paste0(
                "../seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                i,
                ".txt"
              ),
              paste0("seed_number_for_",
                     rep, "_Reps", "_Herit_", i, ".txt")
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
        } else {
          # Calcualte V_e, the residual variance
          residual_cov <-
            diag( (apply(base_line_trait, 2, var) / c(i, h2_MT)) -
                    apply(base_line_trait, 2, var))
          # Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated_data[[z]] <-
              data.frame(rownames(base_line_trait),
                         matrix(NA, nrow(base_line_trait), ntraits),
                         rep = z)
            colnames(simulated_data[[z]]) <-
              c("<Taxa>", "Target",
                c(paste0("Trait_", 2:ntraits, "_H2_", i), "Rep"))
            if (!is.null(seed)) {
              set.seed(round( (seed * z * z) * i))
              ss <- c(ss, round( (seed * z * z) * i))
            }
            # using multivariate normal with ntrait variances and
            # 0 covariances for independent residuals
            simulated_data[[z]][, 2:(ntraits + 1)] <-
              base_line_trait + mvtnorm::rmvnorm(
                n = nrow(base_line_trait),
                mean = rep(0, ntraits),
                sigma = residual_cov
              )
          }
          H2 <- sapply(1:rep, function(x) {
            apply(base_line_trait, 2, var) /
              apply(simulated_data[[x]][1:ntraits + 1], 2, var)
          })
          H2 <- apply(H2, 1, mean)
          names(H2) <- NULL
          #cat("Populational heritability: \n")
          #HH <- c(i, h2_MT)
          #names(HH) <- names(H2)
          #print(HH)
          cat("\nSample heritability (Average of",
              rep,
              " replications): \n")
          print(H2)
          if (format == "multi-file") {
            invisible(lapply(1:rep, function(x) {
              data.table::fwrite(
                simulated_data[[x]][- (ntraits + 2)],
                paste0("Simulated_Data_", "_Rep", x, "_Herit_", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            temp_simulated_data <- do.call(rbind, simulated_data)
            data.table::fwrite(
              temp_simulated_data,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else if (format == "gemma") {
            invisible(lapply(1:rep, function(x) {
              temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)], 
                                 by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", "_Rep", x, "_Herit_", i, ".fam"),
                row.names = FALSE,
                sep = "\t",
                col.names = FALSE,
                quote = FALSE,
                na = NA
              )
            }))
          } else {
            temp <- simulated_data[[1]][- (ntraits + 2)]
            for (j in 2:rep) {
              temp <- cbind(temp, simulated_data[[j]][-c(1, (ntraits + 2))])
            }
            data.table::fwrite(
              temp,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
          # Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(
              format == "multi-file",
              paste0(
                "../seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                i,
                ".txt"
              ),
              paste0("seed_number_for_",
                     rep, "_Reps", "_Herit_", i, ".txt")
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
        }
      }  #End for(i in h2)
      if (to_r) {
        cat(paste("Files saved in:", getwd()))
        simulated_data <- do.call(rbind, simulated_data)
        return(simulated_data)
      } else {
        return(paste("Files saved in:", getwd())) 
      }
    } else {
      # For loop through the vector of heritabilities
      for (i in h2) {
        # If heritability is zero
        ss <- c()
        if (i == 0) {
          simulated_data <-
            data.frame(rownames(base_line_trait),
                       matrix(NA, nrow(base_line_trait), rep))
          # Format the output file for the simulated phenotypes
          colnames(simulated_data) <-
            c("<Taxa>",
              paste0("normal_random_variables_", 1:rep))
          # Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              set.seed(seed + j)
              ss <- c(ss, seed + j)
            }
            simulated_data[, j + 1] <-
              rnorm(nrow(base_line_trait),
                    mean = 0,
                    sd = 1)
          }
          # Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed_number_for_", rep, "_Reps", "_Herit_", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          
          if (format == "multi-file") {
            invisible(apply(as.matrix(1:rep), 1,function(x) {
              data.table::fwrite(
                simulated_data[, c(1, x + 1)],
                paste0("Simulated_Data", "_Rep", x, "_Herit_", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            temp <- simulated_data[, 1:2]
            colnames(temp) <- c("<Taxa>","Pheno")
            for(x in 2:rep) {
              temp2 <- simulated_data[, c(1, x + 1)]
              colnames(temp2) <- c("<Taxa>","Pheno")
              temp <- rbind(temp, temp2 )
            }
            temp$reps <- rep(1:rep, each = nrow(base_line_trait))
            data.table::fwrite(
              temp,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else if (format == "gemma") {
            temp_fam <- merge(fam, simulated_data, 
                              by.x = "V1", by.y = "<Taxa>", sort = FALSE)
            data.table::fwrite(
              temp_fam,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".fam"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            data.table::fwrite(
              simulated_data,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
        } else {
          # Calcualte V_e, the residual variance
          residual.variance <-
            (var(base_line_trait) / i) - var(base_line_trait)
          simulated_data <-
            data.frame(rownames(base_line_trait),
                       matrix(NA, nrow(base_line_trait), rep))
          # Format the output file for the simulated phenotypes
          colnames(simulated_data) <-
            c("<Taxa>", c(paste0(
              "Heritability_", i, "_Rep_", 1:rep
            )))
          # Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              ss <- c(ss, round( (seed * j * j) * i))
              set.seed(round( (seed * j * j) * i))
            }
            normal_random_variables <-
              rnorm(
                nrow(base_line_trait),
                mean = 0,
                sd = sqrt(residual.variance)
              )
            simulated_data[, j + 1] <-
              base_line_trait + normal_random_variables
          }
          H2 <-
            mean(as.vector(var(base_line_trait)) /
                   apply(simulated_data[, 1:rep + 1], 2, var))
          #cat("\n\nPopulational heritability: \n")
          #print(i)
          cat("\nSample heritability (Average of",
              rep,
              " replications): \n"
          )
          print(H2)
          # Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed_number_for_", rep, "_Reps", "_Herit_", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          if (format == "multi-file") {
            invisible(apply(as.matrix(1:rep), 1,function(x) {
              data.table::fwrite(
                simulated_data[, c(1, x + 1)],
                paste0("Simulated_Data", "_Rep", x, "_Herit_", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            temp <- simulated_data[, 1:2]
            colnames(temp) <- c("<Taxa>","Pheno")
            for(x in 2:rep) {
              temp2 <- simulated_data[, c(1, x + 1)]
              colnames(temp2) <- c("<Taxa>","Pheno")
              temp <- rbind(temp, temp2 )
            }
            temp$reps <- rep(1:rep, each = nrow(base_line_trait))
            data.table::fwrite(
              temp,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else if (format == "gemma") {
            temp_fam <- merge(fam, simulated_data, 
                              by.x = "V1", by.y = "<Taxa>", sort = FALSE)
            data.table::fwrite(
              temp_fam,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".fam"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            data.table::fwrite(
              simulated_data,
              paste0("Simulated_Data_", rep, "_Reps", "_Herit_", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }
        }
      }
      if (to_r) {
        cat(paste("Files saved in:", getwd()))
        return(simulated_data)
      } else {
        return(paste("Files saved in:", getwd())) 
      }
    }
  }
