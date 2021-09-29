#' Simulating vQTL
#' @export
#' @import utils
#' @import stats
#' @importFrom rlang inform
#' @importFrom data.table fwrite
#' @param QTN Object with SNPs selected to be QTNS
#' @param base_line_trait Genetic values from mean QTN
#' @param var_QTN_num The number of vQTNs that you want to select or the number of vQTNs
#' @param var_effect This is for assigning the largest vQTN effect size if sim_method is set to "geometric"
#' @param seed The seed number used to generate random numbers
#' @param mean = mean,
#' @param h2 The heritability for each traits being simulated.
#' It could be either a vector with length equals to `ntraits`,
#' or a matrix with ncol equals to `ntraits`. If the later is used, the simulation
#' will loop over the number of rows and will generate a result for each row.
#' If a single trait is being simulated and h2 is a vector,
#' one simulation of each heritability value will be conducted. Either none or
#' all traits are expected to have `h2 = 0`.
#' @param rep The number of experiments (replicates of a trait with the same
#' genetic architecture) to be simulated.
#' @param output_format output format
#' @param fam = NULL,
#' @param to_r Option for outputting the simulated results as an R data.frame in
#' addition to saving it to file. If TRUE, results need to be assigned to an
#' R object (see vignette).
#' @param remove_add_effect = F
#' @return trait simulated under a vQTL model
#' @references Fernandes, S.B., and Lipka, A.E., 2020 simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic
#' SIMulation of Pleiotropic, Linked and Epistatic PHENOTYPES. BMC Bioinformatics 21(1):491,
#' \doi{https://doi.org/10.1186/s12859-020-03804-y} \cr
#' @author Matthew Murphy, Samuel B Fernandes and Alexander E Lipka
#' Last update: APR 2, 2021

vQTL <- function(QTN,
                 base_line_trait = NULL,
                 var_QTN_num = NULL,
                 var_effect = NULL,
                 h2 = NULL,
                 rep = NULL,
                 seed = NULL,
                 mean = NULL,
                 output_format = NULL,
                 fam = NULL,
                 to_r = NULL,
                 remove_add_effect = F) {
  msg <- "Please cite Murphy et al. (2021) when simulating vQTLs!"
  rlang::inform(msg, .frequency = "once", .frequency_id = msg)
  base_line_trait <- scale(base_line_trait)
  if (output_format == "multi-file") {
    dir.create("Phenotypes")
    setwd("./Phenotypes")
  }
  n <- nrow(QTN)
  sigma <-
    matrix(1,
           nrow = n,
           ncol = 1)
  QTN <- QTN + 1
  for (i in 1:var_QTN_num) {
    sigma <-
      sigma + ((var_effect[i]) * (QTN[, i]))
  }
  rownames(sigma) <-
    rownames(QTN)
  results <- vector("list", nrow(h2))
  for (j in seq_len(nrow(h2))) {
    ksq <-
      c((var(base_line_trait) / h2[j, 1] - var(base_line_trait)) / ((median(sigma[, 1]) ^ 2)))
    k <- sqrt(ksq)
    set.seed(seed + j)
    v <-  t(apply(as.matrix(1:n), 1, function(i) {
      k * rnorm(rep, mean = 0, sd = sigma[i, ])
    }))
    trait <-  apply(as.matrix(1:rep), 1, function(i) {
      if ( remove_add_effect) {
        v[, i] + mean
      } else {
        base_line_trait + v[, i] + mean
      }
    })
    simulated_data <- data.frame(
      taxa = rownames(QTN),
      trait,
      check.names = FALSE,
      fix.empty.names = FALSE
    )
    if (nrow(h2) > 1) {
      results[[j]] <- data.frame(
        simulated_data,
        h2 = unname(h2[j, 1]),
        check.names = FALSE,
        fix.empty.names = FALSE
      )
    }
    colnames(simulated_data) <-
      c("<Trait>", c(paste0("h2_", h2[j, 1], "_Rep_", 1:rep)))
    if (output_format == "multi-file") {
      invisible(apply(as.matrix(1:rep), 1, function(x) {
        data.table::fwrite(
          simulated_data[, c(1, x + 1)],
          paste0("Simulated_Data",
                 "_Rep",
                 x,
                 "_Herit_",
                 h2[j, 1],
                 ".txt"),
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
        paste0("Simulated_Data_",
               rep,
               "_Reps",
               "_Herit_",
               h2[j, 1],
               ".txt"),
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
        paste0("Simulated_Data_",
               rep,
               "_Reps",
               "_Herit_",
               h2[j, ],
               ".fam"),
        row.names = FALSE,
        col.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    } else {
      data.table::fwrite(
        simulated_data,
        paste0("Simulated_Data_",
               rep,
               "_Reps",
               "_Herit_",
               h2[j, ],
               ".txt"),
        row.names = FALSE,
        sep = "\t",
        quote = FALSE,
        na = NA
      )
    }
  }
  if (nrow(h2) > 1) {
    H2  <- round(unlist(lapply(results, function(o) {
      mean(apply(trait, 2, function(x)
        1 / var(x)))
    })), 4)
  } else {
    H2  <-
      round(mean(apply(simulated_data[,-1], 2, function(x)
        1 / var(x))), 4)
  }
  cat("\nPopulation Heritability:\n")
  print(h2)
  if(!remove_add_effect){cat("\nSample Heritability (Average of",
      rep,
      " replications): \n[This value might be biased, please read Murphy et al. (2021)] \n")
  print(H2)
  }
  if (to_r) {
    if (nrow(h2) > 1) {
      results <-
        as.data.frame(data.table::rbindlist(results, use.names = F))
      colnames(results) <-
        c("taxa", paste0("rep", 1:rep), "h2")
      return(list(simulated_data = results))
    } else {
      colnames(simulated_data) <-
        c("taxa", paste0("rep", 1:rep))
      return(list(simulated_data = simulated_data))
    }
  }
}
