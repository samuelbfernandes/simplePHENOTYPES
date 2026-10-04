#' Simulating vQTL
#'
#' Internal helper for the frozen legacy engine. It is called only from
#' `create_phenotypes()` and takes objects built inside that call (`QTN`,
#' `base_line_trait`), so it is not usable on its own. The user-facing
#' variance-QTL function is the grammar layer [vqtl()].
#'
#' Model (frozen v1 engine, after Murphy et al. 2022): the mean-QTN genetic value
#' is standardized to unit variance, then each individual \eqn{i} draws
#' \eqn{v_{ir} = k \cdot N(0, \sigma_i^2)} with
#' \eqn{\sigma_i = 1 + \sum_q \text{var\_effect}_q (x_{iq} + 1)} (the vQTL
#' standard deviation is linear in the 0/1/2-shifted dosage) and
#' \eqn{k^2 = (1/h^2 - 1) / \mathrm{median}(\sigma)^2}. The phenotype is
#' `scale(base_line_trait) + v + mean` (or `v + mean` when
#' `remove_add_effect = TRUE`). The scale of the mean-QTN effects therefore
#' does not matter: only their shape does.
#'
#' `h2` is therefore a calibration at the median vQTL standard deviation: an
#' individual with \eqn{\sigma_i = \mathrm{median}(\sigma)} has
#' \eqn{1 / (1 + k^2\sigma_i^2) = h^2}, but the population variance ratio is
#' \eqn{1 / (1 + k^2\,\mathrm{mean}(\sigma^2))}: lower than, equal to or
#' higher than `h2` as \eqn{\mathrm{mean}(\sigma^2)} is greater than, equal
#' to or less than \eqn{\mathrm{median}(\sigma)^2}. The
#' equation is kept for v1 compatibility; [vqtl()] is the grammar alternative.
#'
#' Requirements checked here: `h2` in (0, 1] (a zero heritability has no
#' finite vQTL scale), a single trait, one effect per vQTN and
#' \eqn{\sigma_i \ge 0} for every individual, with a positive median (a
#' negative `var_effect` can make \eqn{\sigma_i} negative, for which `rnorm()`
#' would return NaN). The "Sample Heritability" printed at the end is
#' \eqn{1 / \mathrm{var}(\text{phenotype})} averaged over replicates, computed
#' separately for every row of `h2`; it is biased (Murphy et al. 2022).
#' @keywords internal
#' @noRd
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
#' @author Matthew Murphy, Samuel B Fernandes and Alexander E Lipka. Last update: APR 2, 2021
#'

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
  # Leave the caller's RNG stream untouched (every draw is preceded by its own
  # set.seed(), so restoring the snapshot cannot change any result).
  if (!is.null(seed)) {
    .rng_state <- .Random.seed_safe()
    on.exit(.restore_seed(.rng_state), add = TRUE)
  }
  msg <- .cite_main("Murphy et al. (2022), Heredity 129:93-102, ",
                    "doi:10.1038/s41437-022-00541-1, when simulating variance ",
                    "QTL (vQTLs).")
  rlang::inform(msg, .frequency = "once", .frequency_id = msg)
  if (NCOL(base_line_trait) != 1L) {
    stop("vQTL: variance QTL simulation supports a single trait; ",
         "`base_line_trait` has ", NCOL(base_line_trait), " columns. ",
         "Use ntraits = 1 with a model containing 'V'.", call. = FALSE)
  }
  if (!is.matrix(h2)) h2 <- as.matrix(h2)
  if (anyNA(h2) || any(h2[, 1] <= 0) || any(h2[, 1] > 1)) {
    stop("vQTL: `h2` must be in (0, 1] for a model with variance QTNs ",
         "(h2 = 0 leaves no genetic variance to modulate and the vQTL ",
         "scale is undefined).", call. = FALSE)
  }
  if (length(var_effect) != var_QTN_num || anyNA(var_effect) ||
      !all(is.finite(var_effect))) {
    stop("vQTL: `var_effect` must contain one finite value per vQTN (",
         var_QTN_num, " expected, ", length(var_effect), " supplied).",
         call. = FALSE)
  }
  if (NCOL(QTN) != var_QTN_num) {
    stop("vQTL: the vQTN genotype matrix has ", NCOL(QTN), " column(s) but ",
         "var_QTN_num = ", var_QTN_num, ".", call. = FALSE)
  }
  .v1_check_vqtl_baseline(base_line_trait)
  base_line_trait <- scale(base_line_trait)
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
  if (any(sigma[, 1] < 0) || !(median(sigma[, 1]) > 0)) {
    stop("vQTL: the vQTL standard deviation 1 + sum(var_effect * (dosage + 1)) ",
         "is negative for ", sum(sigma[, 1] < 0), " individual(s) (minimum ",
         format(min(sigma[, 1]), digits = 3), ", median ",
         format(median(sigma[, 1]), digits = 3), "). Use non-negative ",
         "`var_effect` values so that every individual has a non-negative ",
         "standard deviation and the median is positive.", call. = FALSE)
  }
  if (output_format == "multi-file") {
    dir.create("Phenotypes")
    setwd("./Phenotypes")
  }
  results <- vector("list", nrow(h2))
  H2 <- numeric(nrow(h2))
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
    # Sample heritability of THIS row of h2 (1 / var(phenotype), the genetic
    # value has unit variance), averaged over replicates.
    H2[j] <- round(mean(apply(as.matrix(trait), 2, function(x)
      1 / var(x))), 4)
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
  cat("\nHeritability target (calibrated at the median vQTL standard deviation):\n")
  print(h2)
  if(!remove_add_effect){cat("\nSample Heritability (Average of",
      rep,
      " replications): \n[This value might be biased, please read Murphy et al. (2022)] \n")
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
