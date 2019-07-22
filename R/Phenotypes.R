#' Generate environmental effects based on a given heritability
#' @export
#' @param base.line.trait = NULL,
#' @param h2 = NULL,
#' @param rep = NULL,
#' @param seed = NULL,
#' @param ntraits = NULL,
#' @param h2_MT = NULL,
#' @param format = "multi-file"
#' @return Phenotypes for ntraits traits
#' @author Alex lipka and Samuel Fernandes

#' Last update: Jul 22, 2019
#'
#'----------------------------Phenotypes-----------------------------------------
`Phenotypes` <-
  function(base.line.trait = NULL,
           h2 = NULL,
           rep = NULL,
           seed = NULL,
           ntraits = NULL,
           h2_MT = NULL,
           format = "multi-file") {
    #---------------------------------------------------------------------------
    if (ntraits > 1) {
      if (format == "multi-file") {
        #Create a working directory for the output results:
        dir.create("Phenotypes")

        #Set the working directory
        setwd("./Phenotypes")
      }
      #For loop through the vector of heritabilities
      for (i in h2) {
        simulated.data <- vector("list", rep)
        ss <- c()

        #If heritability is zero
        if (i == 0) {
          #Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated.data[[z]] <-
              data.frame(rownames(base.line.trait),
                         matrix(NA, nrow(base.line.trait), ntraits),
                         rep = z)
            colnames(simulated.data[[z]]) <-
              c("<Taxa>",
                "Target",
                paste0("Trait_", 2:ntraits),
                "Rep")

            if (!is.null(seed)) {
              ss <- c(ss, seed + z)
              set.seed(seed + z)
            }
            #using multivariate normal with ntrait variances and 0 covariances
            #for independent residuals
            simulated.data[[z]][, 2:(ntraits + 1)] <-
              mvtnorm::rmvnorm(
                n = nrow(base.line.trait),
                mean = rep(0, ntraits),
                sigma = diag(1, ntraits)
              )
          }

          if (format == "multi-file") {
            invisible(lapply(1:rep, function (x) {
              data.table::fwrite(
                simulated.data[[x]][-(ntraits + 2)],
                paste0("Simulated.Data.", ".Rep", x, ".Herit.", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            simulated.data <-  do.call(rbind, simulated.data)
            data.table::fwrite(
              simulated.data,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            x <- simulated.data[[1]][-(ntraits + 2)]
            for (j in 2:rep) {
              x <- cbind(x, simulated.data[[j]][-c(1, (ntraits + 2))])
            }
            data.table::fwrite(
              x,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }

          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(
              format == "multi-file",
              paste0(
                "../seed.number.for.",
                rep,
                ".Reps",
                ".Herit.",
                i,
                ".txt"
              ),
              paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt")
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )

        } else {
          #Calcualte V_e, the residual variance
          residual.cov <-
            diag((apply(base.line.trait, 2, var) / c(i, h2_MT)) - apply(base.line.trait, 2, var))
          #Simulate rep experiments with a given heritability
          for (z in 1:rep) {
            simulated.data[[z]] <-
              data.frame(rownames(base.line.trait),
                         matrix(NA, nrow(base.line.trait), ntraits),
                         rep = z)

            colnames(simulated.data[[z]]) <-
              c("<Taxa>", "Target", c(paste0("Trait_", 2:ntraits, "_H2_", i), "Rep"))

            if (!is.null(seed)) {
              set.seed(round((seed * z * z) * i))
              ss <- c(ss, round((seed * z * z) * i))
            }
            #using multivariate normal with ntrait variances and 0 covariances
            #for independent residuals
            simulated.data[[z]][, 2:(ntraits + 1)] <-
              base.line.trait +
              mvtnorm::rmvnorm(
                n = nrow(base.line.trait),
                mean = rep(0, ntraits),
                sigma = residual.cov
              )
          }

          H2 <- sapply(1:rep, function (x) {
            apply(base.line.trait, 2, var) /
              apply(simulated.data[[x]][1:ntraits + 1], 2, var)
          })

          H2 <- apply(H2, 1, mean)
          cat("\nPopulational heritability: \n")
          HH <- c(i, h2_MT)
          names(HH) <- names(H2)
          print(HH)

          cat(paste(
            "Sample heritability (Average of",
            rep,
            " replications): \n"
          ))
          print(H2)


          if (format == "multi-file") {
            invisible(lapply(1:rep, function (x) {
              data.table::fwrite(
                simulated.data[[x]][-(ntraits + 2)],
                paste0("Simulated.Data.", ".Rep", x, ".Herit.", i, ".txt"),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE,
                na = NA
              )
            }))
          } else if (format == "long") {
            simulated.data <-  do.call(rbind, simulated.data)
            data.table::fwrite(
              simulated.data,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          } else {
            x <- simulated.data[[1]][-(ntraits + 2)]
            for (j in 2:rep) {
              x <- cbind(x, simulated.data[[j]][-c(1, (ntraits + 2))])
            }
            data.table::fwrite(
              x,
              paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE,
              na = NA
            )
          }

          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            ifelse(
              format == "multi-file",
              paste0(
                "../seed.number.for.",
                rep,
                ".Reps",
                ".Herit.",
                i,
                ".txt"
              ),
              paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt")
            ),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )

        }
      }#End for(i in h2)
      return(paste("Files saved in:", getwd()))

    } else{
      #For loop through the vector of heritabilities
      for (i in h2) {
        #If heritability is zero
        ss <- c()
        if (i == 0) {
          simulated.data <-
            data.frame(rownames(base.line.trait), matrix(NA, nrow(base.line.trait), rep))
          #Format the output file for the simulated phenotypes
          colnames(simulated.data) <-
            c("<Taxa>",
              paste0("the.normal.random.variables.", 1:rep))

          #Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              set.seed(seed + j)
              ss <- c(ss, seed + j)
            }
            simulated.data[, j + 1] <-
              rnorm(nrow(base.line.trait),
                    mean = 0,
                    sd = 1)
          }

          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          data.table::fwrite(
            simulated.data,
            paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )

        } else{
          #Calcualte V_e, the residual variance
          residual.variance <-
            (var(base.line.trait) / i) - var(base.line.trait)

          simulated.data <-
            data.frame(rownames(base.line.trait), matrix(NA, nrow(base.line.trait), rep))
          #Format the output file for the simulated phenotypes
          colnames(simulated.data) <-
            c("<Taxa>", c(paste0(
              "Heritability_", i, "_Rep_", 1:rep
            )))

          #Simulate rep experiments with a given heritability
          for (j in 1:rep) {
            if (!is.null(seed)) {
              ss <- c(ss, round((seed * j * j) * i))
              set.seed(round((seed * j * j) * i))
            }
            the.normal.random.variables <-
              rnorm(
                nrow(base.line.trait),
                mean = 0,
                sd = sqrt(residual.variance)
              )

            simulated.data[, j + 1] <-
              base.line.trait + the.normal.random.variables
          }

          H2 <- mean(as.vector(var(base.line.trait)) /
                       apply(simulated.data[, 1:rep + 1], 2, var))
          cat("\nPopulational heritability: \n")
          #HH <- c(i); names(HH) <- names(H2)
          print(i)

          cat(paste(
            "Sample heritability (Average of",
            rep,
            " replications): \n"
          ))
          print(H2)

          #Output the m rep and the seed numbers, formatted for TASSEL
          write.table(
            ss,
            paste0("seed.number.for.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t",
            quote = FALSE
          )
          data.table::fwrite(
            simulated.data,
            paste0("Simulated.Data.", rep, ".Reps", ".Herit.", i, ".txt"),
            row.names = FALSE,
            sep = "\t",
            quote = FALSE,
            na = NA
          )
        }
      }#End for(i in h2)

      return(paste("Files saved in:", getwd()))
    }
  }
