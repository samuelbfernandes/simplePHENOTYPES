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
#' @param QTN_variance = FALSE
#' @param cor_res = NULL
#' @return Phenotypes for ntraits traits
#' @author Samuel Fernandes and Alexander Lipka
#' Last update: Apr 20, 2020
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
           verbose = TRUE,
           QTN_variance = FALSE,
           add  = NULL,
           dom = NULL,
           epi = NULL,
           cor_res = NULL,
           mean = NULL,
           cor = NULL) {
    #---------------------------------------------------------------------------
    h <- 0
    n <- nrow(base_line_trait[[1]]$base_line)
    if (is.null(cor_res)) {
      cor_res <- diag(1, ntraits)
    }
    names <- if (output_format == "gemma") {
      fam[, 1]
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
        vg <- apply(base_line_trait[[1]]$base_line, 2, var)
        if (any(vg == 0)) {
          warning("Genetic variance = 0 for at least one trait! This will result in h2 = 0!",
                  call. = F)
          h2[, which(vg == 0)] <- 0
        }
        for (i in seq_len(nrow(h2))) {
          simulated_data <- vector("list", rep)
          simulated_cor <- vector("list", rep)
          ss <- c()
          if (any(h2[i, ] == 0)) {
            colnames <- c("<Trait>", paste0("Trait_", 1:ntraits), "Rep")
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(
                  names,
                  matrix(NA, n, ntraits),
                  rep = z,
                  check.names = FALSE,
                  fix.empty.names = FALSE
                )
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer(seed + z)
                ss <- c(ss, sss)
                set.seed(sss)
              }
              sigma <-
                sqrt(diag(1, ntraits)) %*% cor_res %*% sqrt(diag(1, ntraits))
              simulated_cor[[z]] <- sigma
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                mvtnorm::rmvnorm(n = n,
                                 mean = mean,
                                 sigma = sigma)
            }
            H2 <- matrix(0, 1, ntraits)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
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
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(
                  fam,
                  simulated_data[[x]][- (ntraits + 2)],
                  by.x = "V1",
                  by.y = "<Trait>",
                  sort = FALSE
                )
                data.table::fwrite(
                  temp_fam,
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".fam"
                  ),
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
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            if (verbose){
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  i,
                  ".txt"
                ),
                paste0(
                  "Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                )
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
          } else {
            residual_cov <-  diag((vg / h2[i, ]) - vg)
            colnames <- c("<Trait>",
                          c(paste0(
                            "Trait_",
                            1:ntraits, "_H2_", h2[i, ]
                          ), "Rep"))
            vp <- matrix(NA, rep, ntraits)
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(
                  names,
                  matrix(NA, n, ntraits),
                  rep = z,
                  check.names = FALSE,
                  fix.empty.names = FALSE
                )
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <-  as.integer((seed + z) * round(h2[i, 1] * 10))
                set.seed(sss)
                ss <- c(ss, sss)
              }
              sigma <-
                sqrt(residual_cov) %*% cor_res %*% sqrt(residual_cov)
              simulated_cor[[z]] <- sigma
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                base_line_trait[[1]]$base_line +
                mvtnorm::rmvnorm(n = n,
                                 mean = rep(0, ntraits),
                                 sigma = sigma)
              simulated_data[[z]][, 2:(ntraits + 1)] <- 
                apply(t(1:ntraits), 2, function(x) {simulated_data[[z]][, 2:(ntraits + 1)][, x] + mean[x]})
              vp[z, ] <-
                diag(var(simulated_data[[z]][, 2:(ntraits + 1)]))
            }
            H2_temp <- vg / t(vp)
            if (QTN_variance) {
              if (!is.null(cor)) {warning(
                "Please notice that when \'cor\' is provided, the PVE is only accurate for trait 1 since the Cholesky transformation changes the allelic effects of all other traits to ensure the desired correlation.",
                call. = F,
                immediate. = T
              )}
              if (add) {
                for (u in 1:ntraits) {
                  lqtna <- length(base_line_trait[[1]]$QTN_var$var_add[[u]])
                  add_var_per_QTN_temp <-
                    data.frame(matrix(NA, rep, lqtna))
                  colnames(add_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtna, sep = "_")
                  for (p in 1:rep) {
                    add_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_add[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          add_var_per_QTN_temp),
                    paste0("PVE_of_ADD_QTNs_Trait_", u , ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
              if (dom) {
                for (u in 1:ntraits) {
                  lqtnd <- length(base_line_trait[[1]]$QTN_var$var_dom[[u]])
                  dom_var_per_QTN_temp <-
                    data.frame(matrix(
                      NA,
                      rep,
                      lqtnd
                    ))
                  colnames(dom_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtnd, sep = "_")
                  for (p in 1:rep) {
                    dom_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_dom[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          dom_var_per_QTN_temp),
                    paste0("PVE_of_DOM_QTNs_Trait_", u , ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
              if (epi) {
                for (u in 1:ntraits) {
                  lqtne <- length(base_line_trait[[1]]$QTN_var$var_epi[[u]])
                  epi_var_per_QTN_temp <-
                    data.frame(matrix(
                      NA,
                      rep,
                      lqtne
                    ))
                  colnames(epi_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtne, sep = "_")
                  for (p in 1:rep) {
                    epi_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_epi[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          epi_var_per_QTN_temp),
                    paste0("PVE_of_EPI_QTNs_trait_", u , ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
            }
            H2[i, ] <- apply(H2_temp, 1, mean)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
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
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(
                  fam,
                  simulated_data[[x]][- (ntraits + 2)],
                  by.x = "V1",
                  by.y = "<Trait>",
                  sort = FALSE
                )
                data.table::fwrite(
                  temp_fam,
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".fam"
                  ),
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
                temp <- cbind(temp, simulated_data[[j]][- c(1, (ntraits + 2))])
              }
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            if (verbose){
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                paste0(
                  "Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                )
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
          }
        }
        colnames(H2) <- paste0("Trait_", 1:ntraits)
        cat("\nPopulation Heritability:\n")
        print(h2)
        cat("\nSample heritability (Average of",
            rep,
            " replications): \n")
        print(H2)
        if (all(H2 != 1)) {
          sample_cor <- matrix(0, ntraits, ntraits)
          for (v in 1:rep) {
            sample_cor <- (sample_cor + simulated_cor[[v]])
          }
          sample_cor <- stats::cov2cor(sample_cor / rep)
        } else {
          sample_cor <- NULL
        }
        if (to_r) {
          simulated_data <- do.call(rbind, simulated_data)
          return(list(simulated_data = simulated_data,
                      sample_cor = sample_cor))
        } else {
          return(list(sample_cor = sample_cor))
        }
      } else {
        H2 <- matrix(NA, nrow = rep, ncol = length(h2))
        vg <- var(base_line_trait[[1]]$base_line)
        if (vg == 0) {
          warning("Genetic variance = 0! Please select a different set of QTNs",
                  call. = F)
          h2[, which(vg == 0)] <- 0
        }
        for (i in seq_len(nrow(h2))) {
          ss <- c()
          if (any(h2[i, ] == 0)) {
            simulated_data <-
              data.frame(
                names,
                matrix(NA, n, rep),
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            colnames(simulated_data) <-
              c("<Trait>",
                paste0("normal_random_variables_", 1:rep))
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer(seed + j)
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[, j + 1] <-
                rnorm(n, mean = mean, sd = 1)
              
            }
            H2 <- matrix(0, nrow = 1, ncol = length(h2))
            if (verbose){
              write.table(
              ss,
              paste0(
                "Seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                       paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                       paste(h2[i,], collapse = "_")),
                ".txt"
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0(
                    "Simulated_Data",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Trait>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Trait>", "Pheno")
                  temp <- rbind(temp, temp2)
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(
                fam,
                simulated_data,
                by.x = "V1",
                by.y = "<Trait>",
                sort = FALSE
              )
              data.table::fwrite(
                temp_fam,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".fam"
                ),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          } else {
            residual.variance <- (vg / h2[i, 1]) - vg
            simulated_data <-
              data.frame(
                names,
                matrix(NA, n, rep),
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            colnames(simulated_data) <-
              c("<Trait>", c(paste0(
                "Heritability_", h2[i, ], "_Rep_", 1:rep
              )))
            vp <- c()
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer((seed + j) * round(h2[i, 1] * 10))
                ss <- c(ss, sss)
                set.seed(sss)
              }
              normal_random_variables <-
                rnorm(n,
                      mean = 0,
                      sd = sqrt(residual.variance))
              simulated_data[, j + 1] <-
                base_line_trait[[1]]$base_line + normal_random_variables + mean
              vp[j] <- var(simulated_data[, j + 1])
              H2[j, i] <- vg / vp[j]
            }
            if (QTN_variance) {
              if (!is.null(cor)) {warning(
                "Please notice that when \'cor\' is provided, the PVE is only accurate for trait 1 since the Cholesky transformation changes the allelic effects of all other traits to ensure the desired correlation.",
                call. = F,
                immediate. = T
              )}
              if (add) {
                lqtna <- length(base_line_trait[[1]]$var_add)
                add_var_per_QTN <-
                  data.frame(matrix(NA, rep, lqtna))
                colnames(add_var_per_QTN) <-
                  paste("QTN", 1:lqtna, sep = "_")
                for (p in 1:rep) {
                  add_var_per_QTN[p, ] <-
                    base_line_trait[[1]]$var_add /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        add_var_per_QTN),
                  "PVE_of_ADD_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
              if (dom) {
                lqtnd <- length(base_line_trait[[1]]$var_dom)
                dom_var_per_QTN <-
                  data.frame(matrix(
                    NA,
                    rep,
                    lqtnd
                  ))
                colnames(dom_var_per_QTN) <-
                  paste("QTN", 1:lqtnd, sep = "_")
                for (p in 1:rep) {
                  dom_var_per_QTN[p, ] <-
                    base_line_trait[[1]]$var_dom /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        dom_var_per_QTN),
                  "PVE_of_DOM_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
              if (epi) {
                lqtne <- length(base_line_trait[[1]]$var_epi)
                epi_var_per_QTN <-
                  data.frame(matrix(
                    NA,
                    rep,
                    lqtne
                  ))
                colnames(epi_var_per_QTN) <-
                  paste("QTN", 1:lqtne, sep = "_")
                for (p in 1:rep) {
                  epi_var_per_QTN[p, ] <-
                    base_line_trait[[1]]$var_epi /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        epi_var_per_QTN),
                  "PVE_of_EPI_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
            }
            if (verbose){
            write.table(
              ss,
              paste0(
                "Seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                       paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                       paste(h2[i,], collapse = "_")),
                ".txt"
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0(
                    "Simulated_Data",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Trait>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Trait>", "Pheno")
                  temp <- rbind(temp, temp2)
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(
                fam,
                simulated_data,
                by.x = "V1",
                by.y = "<Trait>",
                sort = FALSE
              )
              data.table::fwrite(
                temp_fam,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".fam"
                ),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          }
        }
        H2 <- apply(H2, 2, mean)
        cat("\nPopulation Heritability:", h2)
        cat("\nSample heritability (Average of",
            rep,
            " replications): \n")
        print(H2)
        if (to_r) {
          return(list(simulated_data = simulated_data))
        }
      }
    } else {
      if (ntraits > 1) {
        if (output_format == "multi-file") {
          dir.create("Phenotypes")
          setwd("./Phenotypes")
        }
        H2 <- matrix(NA, nrow(h2), ntraits)
        for (i in seq_len(nrow(h2))) {
          simulated_cor <- vector("list", rep)
          simulated_data <- vector("list", rep)
          H2_temp <- matrix(NA, rep, ntraits)
          vp <- H2_temp
          ss <- c()
          if (any(h2[i, ] == 0)) {
            colnames <- c("<Trait>", paste0("Trait_", 1:ntraits), "Rep")
            for (z in 1:rep) {
              simulated_data[[z]] <-
                data.frame(
                  names,
                  matrix(NA, n, ntraits),
                  rep = z,
                  check.names = FALSE,
                  fix.empty.names = FALSE
                )
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer(seed + z)
                ss <- c(ss, sss)
                set.seed(sss)
              }
              sigma <-
                sqrt(diag(1, ntraits)) %*% cor_res %*% sqrt(diag(1, ntraits))
              simulated_cor[[z]] <- sigma
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                mvtnorm::rmvnorm(n = n,
                                 mean = mean,
                                 sigma = sigma)
            }
            H2 <- matrix(0, 1, ntraits)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
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
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(
                  fam,
                  simulated_data[[x]][- (ntraits + 2)],
                  by.x = "V1",
                  by.y = "<Trait>",
                  sort = FALSE
                )
                data.table::fwrite(
                  temp_fam,
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".fam"
                  ),
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
                temp <- cbind(temp, simulated_data[[j]][- c(1, (ntraits + 2))])
              }
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            if (verbose){
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                paste0(
                  "Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                )
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
          } else {
            colnames <- c("<Trait>",
                          paste0("Trait_", 1:ntraits, "_H2_", h2[i, ]),
                          "Rep")
            for (z in 1:rep) {
              vg <- apply(base_line_trait[[z]]$base_line, 2, var)
              residual_cov <- diag((vg / h2[i, ]) - vg)
              simulated_data[[z]] <-
                data.frame(
                  names,
                  matrix(NA, n, ntraits),
                  rep = z,
                  check.names = FALSE,
                  fix.empty.names = FALSE
                )
              colnames(simulated_data[[z]]) <- colnames
              if (!is.null(seed)) {
                sss <- as.integer((seed + z) * round(h2[i, 1] * 10))
                set.seed(sss)
                ss <- c(ss, sss)
              }
              sigma <-
                sqrt(residual_cov) %*% cor_res %*% sqrt(residual_cov)
              simulated_cor[[z]] <- sigma
              simulated_data[[z]][, 2:(ntraits + 1)] <-
                base_line_trait[[z]]$base_line + mvtnorm::rmvnorm(n = n,
                                                                  mean = rep(0, ntraits),
                                                                  sigma = sigma)
              simulated_data[[z]][, 2:(ntraits + 1)] <- 
                apply(t(1:ntraits), 2, function(x) {simulated_data[[z]][, 2:(ntraits + 1)][, x] + mean[x]})
              if (any(vg == 0)) {
                warning(
                  "Genetic variance = 0 for at least one trait in rep ",
                  z,
                  "! Please select a different set of QTNs",
                  call. = F
                )
                H2_temp[z, ] <- NA
                h <- h + 1
              } else {
                vp[z, ] <- apply(simulated_data[[z]][1:ntraits + 1], 2, var)
                H2_temp[z, ] <-
                  vg / vp[z, ]
              }
            }
            if (QTN_variance) {
              if (!is.null(cor)) {warning(
                "Please notice that when \'cor\' is provided, the PVE is only accurate for trait 1 since the Cholesky transformation changes the allelic effects of all other traits to ensure the desired correlation.",
                call. = F,
                immediate. = T
              )}
              if (add) {
                for (u in 1:ntraits) {
                  lqtna <- length(base_line_trait[[1]]$QTN_var$var_add[[u]])
                  add_var_per_QTN_temp <-
                    data.frame(matrix(NA, rep, lqtna))
                  colnames(add_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtna, sep = "_")
                  for (p in 1:rep) {
                    add_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_add[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          add_var_per_QTN_temp),
                    paste0("PVE_of_ADD_QTNs_trait_", u, ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
              if (dom) {
                for (u in 1:ntraits) {
                  lqtnd <- length(base_line_trait[[1]]$QTN_var$var_dom[[u]])
                  dom_var_per_QTN_temp <-
                    data.frame(matrix(
                      NA,
                      rep,
                      lqtnd
                    ))
                  colnames(dom_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtnd, sep = "_")
                  for (p in 1:rep) {
                    dom_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_dom[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          dom_var_per_QTN_temp),
                    paste0("PVE_of_DOM_QTNs_trait_", u, ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
              if (epi) {
                for (u in 1:ntraits) {
                  lqtne <- length(base_line_trait[[1]]$QTN_var$var_epi[[u]])
                  epi_var_per_QTN_temp <-
                    data.frame(matrix(
                      NA,
                      rep,
                      lqtne
                    ))
                  colnames(epi_var_per_QTN_temp) <-
                    paste("QTN", 1:lqtne, sep = "_")
                  for (p in 1:rep) {
                    epi_var_per_QTN_temp[p, ] <-
                      base_line_trait[[1]]$QTN_var$var_epi[[u]] /
                      vp[p, u]
                  }
                  data.table::fwrite(
                    cbind(REP = 1:rep,
                          epi_var_per_QTN_temp),
                    paste0("PVE_of_EPI_QTNs_trait_", u, ".txt"),
                    row.names = FALSE,
                    sep = "\t",
                    quote = FALSE,
                    na = NA
                  )
                }
              }
            }
            H2[i, ] <- apply(H2_temp, 2, mean, na.rm = T)
            if (output_format == "multi-file") {
              invisible(lapply(1:rep, function(x) {
                data.table::fwrite(
                  simulated_data[[x]][- (ntraits + 2)],
                  paste0(
                    "Simulated_Data_",
                    "_Rep_",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
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
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              invisible(lapply(1:rep, function(x) {
                temp_fam <-  merge(
                  fam,
                  simulated_data[[x]][- (ntraits + 2)],
                  by.x = "V1",
                  by.y = "<Trait>",
                  sort = FALSE
                )
                data.table::fwrite(
                  temp_fam,
                  paste0(
                    "Simulated_Data_",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".fam"
                  ),
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
                temp <- cbind(temp, simulated_data[[j]][- c(1, (ntraits + 2))])
              }
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
            if (verbose){
            write.table(
              ss,
              ifelse(
                output_format == "multi-file",
                paste0(
                  "../Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                paste0(
                  "Seed_number_for_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                )
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
          }
        }
        colnames(H2) <- paste0("Trait_", 1:ntraits)
        cat("\nPopulation Heritability:\n")
        print(h2)
        cat("\nSample heritability (Average of",
            rep - h,
            " replications): \n")
        print(H2)
        if (all(H2 != 1)) {
          sample_cor <- matrix(0, ntraits, ntraits)
          for (v in 1:rep) {
            sample_cor <- (sample_cor + simulated_cor[[v]])
          }
          sample_cor <- stats::cov2cor(sample_cor / rep)
        } else {
          sample_cor <- NULL
        }
        if (to_r) {
          simulated_data <- do.call(rbind, simulated_data)
          return(list(simulated_data = simulated_data,
                      sample_cor = sample_cor))
        } else {
          return(list(sample_cor = sample_cor))
        }
      } else {
        H2 <- matrix(NA, nrow = rep, ncol = length(h2))
        for (i in seq_len(nrow(h2))) {
          ss <- c()
          if (any(h2[i, ] == 0)) {
            simulated_data <-
              data.frame(
                names,
                matrix(NA, n, rep),
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            colnames(simulated_data) <-
              c("<Trait>",
                paste0("normal_random_variables_", 1:rep))
            for (j in 1:rep) {
              if (!is.null(seed)) {
                sss <- as.integer(seed + j)
                set.seed(sss)
                ss <- c(ss, sss)
              }
              simulated_data[, j + 1] <-
                rnorm(n,
                      mean = mean,
                      sd = 1)
            }
            H2 <- matrix(0, nrow = 1, ncol = length(h2))
            if (verbose){
            write.table(
              ss,
              paste0(
                "Seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                       paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                       paste(h2[i,], collapse = "_")),
                ".txt"
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0(
                    "Simulated_Data",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Trait>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Trait>", "Pheno")
                  temp <- rbind(temp, temp2)
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(
                fam,
                simulated_data,
                by.x = "V1",
                by.y = "<Trait>",
                sort = FALSE
              )
              data.table::fwrite(
                temp_fam,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".fam"
                ),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          } else {
            simulated_data <-
              data.frame(
                names,
                matrix(NA, n, rep),
                check.names = FALSE,
                fix.empty.names = FALSE
              )
            colnames(simulated_data) <-
              c("<Trait>", c(paste0(
                "Heritability_", h2[i, ], "_Rep_", 1:rep
              )))
            vp <- c()
            for (j in 1:rep) {
              vg <- var(base_line_trait[[j]]$base_line)
              residual.variance <-  (vg / h2[i, 1]) - vg
              if (!is.null(seed)) {
                sss <-  as.integer((seed + j) * round(h2[i, 1] * 10))
                ss <- c(ss, sss)
                set.seed(sss)
              }
              normal_random_variables <-
                rnorm(n,
                      mean = 0,
                      sd = sqrt(residual.variance))
              simulated_data[, j + 1] <-
                base_line_trait[[j]]$base_line + normal_random_variables + mean
              if (vg == 0) {
                warning(
                  "Genetic variance = 0 in rep ",
                  j,
                  "! Please select a different set of QTNs",
                  call. = F
                )
                h <- h + 1
                H2[j, i] <- NA
              } else {
                vp[j] <-  var(simulated_data[, j + 1])
                H2[j, i] <- vg / vp[j]
              }
            }
            if (QTN_variance) {
              if (!is.null(cor)) {warning(
                "Please notice that when \'cor\' is provided, the PVE is only accurate for trait 1 since the Cholesky transformation changes the allelic effects of all other traits to ensure the desired correlation.",
                call. = F,
                immediate. = T
              )}
              if (add) {
                lqtna <- length(base_line_trait[[1]]$var_add)
                add_var_per_QTN <-
                  data.frame(matrix(NA, rep, lqtna))
                colnames(add_var_per_QTN) <-
                  paste("QTN", 1:lqtna, sep = "_")
                for (p in 1:rep) {
                  add_var_per_QTN[p, ] <-
                    base_line_trait[[p]]$var_add /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        add_var_per_QTN),
                  "PVE_of_ADD_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
              if (dom) {
                lqtnd <- length(base_line_trait[[1]]$var_dom)
                dom_var_per_QTN <-
                  data.frame(matrix(
                    NA,
                    rep,
                    lqtnd
                  ))
                colnames(dom_var_per_QTN) <-
                  paste("QTN", 1:lqtnd, sep = "_")
                for (p in 1:rep) {
                  dom_var_per_QTN[p, ] <-
                    base_line_trait[[p]]$var_dom /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        dom_var_per_QTN),
                  "PVE_of_DOM_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
              if (epi) {
                lqtne <- length(base_line_trait[[1]]$var_epi)
                epi_var_per_QTN <-
                  data.frame(matrix(
                    NA,
                    rep,
                    lqtne
                  ))
                colnames(epi_var_per_QTN) <-
                  paste("QTN", 1:lqtne, sep = "_")
                for (p in 1:rep) {
                  epi_var_per_QTN[p, ] <-
                    base_line_trait[[p]]$var_epi /
                    vp[p]
                }
                data.table::fwrite(
                  cbind(REP = 1:rep,
                        epi_var_per_QTN),
                  "PVE_of_EPI_QTNs_trait_1.txt",
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }
            }
            if (verbose){
            write.table(
              ss,
              paste0(
                "Seed_number_for_",
                rep,
                "_Reps",
                "_Herit_",
                ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                       paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                       paste(h2[i,], collapse = "_")),
                ".txt"
              ),
              row.names = FALSE,
              col.names = FALSE,
              sep = "\t",
              quote = FALSE
            )
            }
            if (output_format == "multi-file") {
              invisible(apply(as.matrix(1:rep), 1, function(x) {
                data.table::fwrite(
                  simulated_data[, c(1, x + 1)],
                  paste0(
                    "Simulated_Data",
                    "_Rep",
                    x,
                    "_Herit_",
                    ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                           paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                           paste(h2[i,], collapse = "_")),
                    ".txt"
                  ),
                  row.names = FALSE,
                  sep = "\t",
                  quote = FALSE,
                  na = NA
                )
              }))
            } else if (output_format == "long") {
              temp <- simulated_data[, 1:2]
              colnames(temp) <- c("<Trait>", "Pheno")
              if (rep > 1) {
                for (x in 2:rep) {
                  temp2 <- simulated_data[, c(1, x + 1)]
                  colnames(temp2) <- c("<Trait>", "Pheno")
                  temp <- rbind(temp, temp2)
                }
              }
              temp$reps <- rep(1:rep, each = n)
              data.table::fwrite(
                temp,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else if (output_format == "gemma") {
              temp_fam <- merge(
                fam,
                simulated_data,
                by.x = "V1",
                by.y = "<Trait>",
                sort = FALSE
              )
              data.table::fwrite(
                temp_fam,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".fam"
                ),
                row.names = FALSE,
                col.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            } else {
              data.table::fwrite(
                simulated_data,
                paste0(
                  "Simulated_Data_",
                  rep,
                  "_Reps",
                  "_Herit_",
                  ifelse(nchar(paste(h2[i,], collapse = "_")) > 100,
                         paste0(h2[i, 1], "...", h2[i, ntraits]) ,
                         paste(h2[i,], collapse = "_")),
                  ".txt"
                ),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }
          }
        }
        H2 <- apply(H2, 2, mean, na.rm = T)
        cat("\nPopulation Heritability:", h2)
        cat("\nSample heritability (Average of",
            rep - h,
            " replications): \n")
        print(H2)
        if (to_r) {
          return(list(simulated_data = simulated_data))
        }
      }
    }
  }
