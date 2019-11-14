#' Generate environmental effects based on a given heritability
#' @keywords internal
#' @param base_line_trait = NULL,
#' @param h2 = NULL,
#' @param rep = NULL,
#' @param seed = NULL,
#' @param ntraits = NULL,
#' @param output_format = 'multi-file'
#' @param fam = NULL
#' @param to_r = FALSE
#' @param rep_by = 'QTN'
#' @param hets = 'QTN'
#' @param verbose = TRUE
#' @return Phenotypes for ntraits traits
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Nov 05, 2019
#'
#'----------------------------phenotypes---------------------------------------
phenotypes <-
  function(base_line_trait = NULL,
           h2 = NULL,
           rep = NULL,
           seed = NULL,
           ntraits = NULL,
           output_format = NULL,
           fam = NULL,
           to_r = NULL,
           rep_by = NULL,
           hets = NULL,
           verbose = TRUE) {
    #---------------------------------------------------------------------------
    h <- 0
    n <- nrow(base_line_trait[[1]]$base_line)
    names <- if (output_format == "gemma"){
      fam[,1]
    } else {
      rownames(base_line_trait[[1]]$base_line) 
    }
    if (rep_by != "QTN") {
      if (ntraits > 1) {
        if (output_format == "multi-file") {
          dir.create("Phenotypes")
          setwd("./Phenotypes")
        }
        H2 <- matrix(NA, nrow(h2), ntraits)
        va <- apply(base_line_trait[[1]]$base_line, 2, var)
        if (any(va == 0)) {
          warning("Genetic variance = 0 for at least one trait! This will result in h2 = 0!", call. = F)
        h2[, which(va == 0)] <- 0
          }
        for (i in 1:nrow(h2)) {
          simulated_data <- vector("list", rep)
          ss <- c()
          if (any(h2[i, ] == 0)) {
            colnames <- c("<Taxa>", paste0("Trait_", 1:ntraits), "Rep")
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(names,
                           matrix(NA, n, ntraits),
                           rep = z)
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer(seed + z)
                ss <- c(ss, sss)
                set.seed(sss)
              }
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                mvtnorm::rmvnorm(
                  n = n,
                  mean = rep(0, ntraits),
                  sigma = diag(1,
                               ntraits)
                )
            }
            H2 <- matrix(0, 1, ntraits)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp_simulated_data <- do.call(rbind, simulated_data)
              data.table::fwrite(
                temp_simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)],
                                   by.x = "V1", by.y = "<Taxa>", sort = FALSE)
                data.table::fwrite(
                  temp_fam,
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".fam"),
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
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  i,
                  ".txt"
                ),
                paste0("seed_number_for_",
                       rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt")
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          } else {
            residual_cov <-  diag( (va / h2[i, ]) - va)
            colnames <- c("<Taxa>",
                          c(paste0("Trait_",
                                   1:ntraits, "_H2_", h2[i, ]), "Rep"))
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(names,
                           matrix(NA, n, ntraits),
                           rep = z)
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <-  as.integer((seed + z) * round(h2[i, 1] * 10))
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                base_line_trait[[1]]$base_line + mvtnorm::rmvnorm(
                  n = n,
                  mean = rep(0, ntraits),
                  sigma = residual_cov
                )
            }
            H2_temp <- sapply(1:rep, function(x) {
              va / apply(simulated_data[[x]][1:ntraits + 1], 2, var)
            })
                H2[i, ] <- apply(H2_temp, 1, mean)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp_simulated_data <- do.call(rbind, simulated_data)
              data.table::fwrite(
                temp_simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)],
                                   by.x = "V1", by.y = "<Taxa>", sort = FALSE)
                data.table::fwrite(
                  temp_fam,
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".fam"),
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
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  paste(h2[i, ], collapse = "_"),
                  ".txt"
                ),
                paste0("seed_number_for_",
                       rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt")
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          }
        }
        colnames(H2) <- paste0("Trait_", 1:ntraits)
          if (verbose) cat("\nSample heritability (Average of",
              rep,
              " replications): \n")
        if (verbose) print(H2)
        if (to_r) {
          simulated_data <- do.call(rbind, simulated_data)
          return(simulated_data)
        }
      } else {
        H2 <- matrix(NA, nrow = rep, ncol = ncol(h2))
        va <- var(base_line_trait[[1]]$base_line)
        if (va == 0) {
          warning("Genetic variance = 0! Please select a different set of QTNs", call. = F)
          h2[, which(va == 0)] <- 0
          }
        for (i in 1:nrow(h2)) {
          ss <- c()
          if (any(h2[i, ] == 0)) {
            simulated_data <-
              data.frame(names, matrix(NA, n, rep))
            colnames(simulated_data) <-
              c("<Taxa>", paste0("normal_random_variables_", 1:rep))
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer(seed + j)
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[, j + 1] <-
                rnorm(n, mean = 0, sd = 1)
            }
            H2 <- matrix(0, nrow = 1, ncol = ncol(h2))
            write.table(
              ss,
              paste0("seed_number_for_", rep, "_Reps", "_Herit_",
                     paste(h2[i, ], collapse = "_"), ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0("Simulated_Data", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Taxa>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Taxa>", "Pheno")
                  temp <- rbind(temp, temp2 )
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(fam, simulated_data,
                                by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".fam"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          } else {
            residual.variance <- (va / h2[i, 1]) - va
            simulated_data <-
              data.frame(names, matrix(NA, n, rep))
            colnames(simulated_data) <-
              c("<Taxa>", c(paste0( "Heritability_", h2[i, ], "_Rep_", 1:rep)))
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer((seed + j) * round(h2[i, 1] * 10))
                ss <- c(ss, sss)
                set.seed(sss)
              }
              normal_random_variables <-
                rnorm(n, mean = 0, sd = sqrt(residual.variance))
              simulated_data[, j + 1] <-
                base_line_trait[[1]]$base_line + normal_random_variables
              H2[j, i] <- va / var(simulated_data[, j + 1])
            }
            write.table(
              ss,
              paste0("seed_number_for_", rep, "_Reps", "_Herit_",
                     paste(h2[i, ], collapse = "_"), ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0("Simulated_Data", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Taxa>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Taxa>", "Pheno")
                  temp <- rbind(temp, temp2)
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(fam, simulated_data,
                                by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".fam"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          }
        }
          H2 <- apply(H2, 2, mean)
          if (verbose) cat("\nSample heritability (Average of",
              rep,
              " replications): \n"
          )
          if (verbose) print(H2)
        if (to_r) {
          return(simulated_data)
        }
      }
    } else {
      if (ntraits > 1) {
        if (output_format == "multi-file") {
          dir.create("Phenotypes")
          setwd("./Phenotypes")
        }
        H2 <- matrix(NA, nrow(h2), ntraits)
        for (i in 1:nrow(h2)) {
          simulated_data <- vector("list", rep)
          H2_temp <- matrix(NA, rep, ntraits)
          ss <- c()
          if (any(h2[i, ] == 0)) {
            colnames <- c("<Taxa>", paste0("Trait_", 1:ntraits), "Rep")
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(names,
                           matrix(NA, n, ntraits),
                           rep = z)
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer(seed + z)
                ss <- c(ss, sss)
                set.seed(sss)
              }
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                mvtnorm::rmvnorm(
                  n = n,
                  mean = rep(0, ntraits),
                  sigma = diag(1, ntraits)
                )
            }
            H2 <- matrix(0, 1, ntraits)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp_simulated_data <- do.call(rbind, simulated_data)
              data.table::fwrite(
                temp_simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)],
                                   by.x = "V1", by.y = "<Taxa>", sort = FALSE)
                data.table::fwrite(
                  temp_fam,
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".fam"),
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
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  paste(h2[i, ], collapse = "_"),
                  ".txt"
                ),
                paste0("seed_number_for_",
                       rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt")
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          } else {
            colnames <- c("<Taxa>",
                          paste0("Trait_", 1:ntraits, "_H2_", h2[i, ]), "Rep")
            for (z in 1:rep) {
              va <- apply(base_line_trait[[z]]$base_line, 2, var)
              residual_cov <- diag( (va / h2[i, ]) - va)
              simulated_data[[z]] <-
                data.frame(names,
                           matrix(NA, n, ntraits),
                           rep = z)
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer((seed + z) * round(h2[i, 1] * 10))
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                base_line_trait[[z]]$base_line + mvtnorm::rmvnorm(
                  n = n,
                  mean = rep(0, ntraits),
                  sigma = residual_cov
                )
              if (any(va == 0)) {
                warning("Genetic variance = 0 for at least one trait in rep ",z,"! Please select a different set of QTNs", call. = F)
                H2_temp[z, ] <- NA
                h <- h + 1
              } else {
                H2_temp[z, ] <-
                  va / apply(simulated_data[[z]][1:ntraits + 1], 2, var)
              }
            }
            H2[i, ] <- apply(H2_temp, 2, mean, na.rm = T)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0("Simulated_Data_", "_Rep_", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp_simulated_data <- do.call(rbind, simulated_data)
              data.table::fwrite(
                temp_simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(fam, simulated_data[[x]][- (ntraits + 2)],
                                   by.x = "V1", by.y = "<Taxa>", sort = FALSE)
                data.table::fwrite(
                  temp_fam,
                  paste0("Simulated_Data_", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".fam"),
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
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  paste(h2[i, ], collapse = "_"),
                  ".txt"
                ),
                paste0("seed_number_for_",
                       rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt")
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
          }
        }
        colnames(H2) <- paste0("Trait_", 1:ntraits)
        if (verbose) cat("\nSample heritability (Average of",
            rep - h,
            " replications): \n")
        if (verbose) print(H2)
        if (to_r) {
          simulated_data <- do.call(rbind, simulated_data)
          return(simulated_data)
        }
      } else {
        H2 <- matrix(NA, nrow = rep, ncol = ncol(h2))
        for (i in 1:nrow(h2)) {
          ss <- c()
          if (any(h2[i, ] == 0)) {
            simulated_data <-
              data.frame(names, matrix(NA, n, rep))
            colnames(simulated_data) <-
              c("<Taxa>", paste0("normal_random_variables_", 1:rep))
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer(seed + j)
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[, j + 1] <-
                rnorm(n,
                      mean = 0,
                      sd = 1)
            }
            H2 <- matrix(0, nrow = 1, ncol = ncol(h2))
            write.table(
              ss,
              paste0("seed_number_for_", rep, "_Reps", "_Herit_",
                     paste(h2[i, ], collapse = "_"), ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0("Simulated_Data", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Taxa>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Taxa>", "Pheno")
                  temp <- rbind(temp, temp2 )
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(fam, simulated_data,
                                by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".fam"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          } else {
            simulated_data <-
              data.frame(names, matrix(NA, n, rep))
            colnames(simulated_data) <-
              c("<Taxa>", c(paste0(
                "Heritability_", h2[i, ], "_Rep_", 1:rep
              )))
            for (j in 1:rep) {
              va <- var(base_line_trait[[j]]$base_line)
              residual.variance <-  (va / h2[i, 1]) - va
              if (!is.null(seed)) {
                sss <-  as.integer((seed + j) * round(h2[i, 1] * 10))
                ss <- c(ss, sss)
                set.seed(sss)
              }
              normal_random_variables <-
                rnorm(
                  n,
                  mean = 0,
                  sd = sqrt(residual.variance)
                )
              simulated_data[, j + 1] <-
                base_line_trait[[j]]$base_line + normal_random_variables
              if (va == 0) {
                warning("Genetic variance = 0 in rep ",j,"! Please select a different set of QTNs", call. = F)
                h <- h + 1
                H2[j, i] <- NA
                } else {
                H2[j, i] <- va / var(simulated_data[, j + 1]) 
              }
            }
            write.table(
              ss,
              paste0("seed_number_for_", rep, "_Reps", "_Herit_",
                     paste(h2[i, ], collapse = "_"), ".txt"),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0("Simulated_Data", "_Rep", x, "_Herit_",
                         paste(h2[i, ], collapse = "_"), ".txt"),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Taxa>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Taxa>", "Pheno")
                  temp <- rbind(temp, temp2 )
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(fam, simulated_data,
                                by.x = "V1", by.y = "<Taxa>", sort = FALSE)
              data.table::fwrite(
                temp_fam,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".fam"),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0("Simulated_Data_", rep, "_Reps", "_Herit_",
                       paste(h2[i, ], collapse = "_"), ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          }
        }
        H2 <- apply(H2, 2, mean, na.rm = T)
        if (verbose) cat("\nSample heritability (Average of",
            rep - h,
            " replications): \n")
        if (verbose) print(H2)
        if (to_r) {
          return(simulated_data)
        }
      }
    }
  }
