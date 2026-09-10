#' Combine single-architecture models under a common heritability
#'
#' `complex_phenotypes()` merges two or more realized `phenotype_sim` objects
#' built on the same genotypes with the same number of traits. The inputs'
#' genetic values are summed (so each contributes in proportion to its genetic
#' variance), the sum is scaled to the requested genetic variance share `h2`,
#' and a common residual is added; the inputs' individual residuals are
#' discarded. Realized H2 is reported from the resulting sample and can differ
#' slightly because the generated genetic and residual vectors need not be
#' exactly orthogonal. This recreates
#' partial pleiotropy by combining, for example, a "pleiotropy" model with an
#' "independent" model.
#'
#' @param ... two or more `phenotype_sim` objects.
#' @param h2 requested genetic variance share (scalar or length `n_traits`).
#' @return a combined `phenotype_sim` (architecture "complex").
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pleio <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
#'                             n_traits = 2, seed = 10, cor = 0.5)
#' pleio <- additive(pleio, prop = 0.4, n_qtn = 3)
#' indep <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, seed = 11)
#' indep <- additive(indep, prop = 0.3, n_qtn = 3)
#' both <- complex_phenotypes(pleio, indep, h2 = 0.5)
complex_phenotypes <- function(..., h2) {
  models <- list(...)
  if (length(models) < 2) {
    stop("complex_phenotypes() needs at least two phenotype_sim objects.",
         call. = FALSE)
  }
  for (m in models) .check_sim(m)

  nt <- models[[1]]$n_traits
  n  <- models[[1]]$n_ind
  nr <- models[[1]]$n_reps
  ids <- models[[1]]$ids
  reference <- models[[1]]
  for (m in models) {
    if (m$n_traits != nt) {
      stop("All inputs must have the same n_traits.", call. = FALSE)
    }
    if (m$n_ind != n || !identical(m$ids, ids)) {
      stop("All inputs must be built on the same genotypes.", call. = FALSE)
    }
    if (m$n_reps != nr) {
      stop("All inputs must have the same n_reps.", call. = FALSE)
    }
    if (!identical(m$kind, reference$kind) ||
        !identical(m$map, reference$map) ||
        !identical(m$ind_idx, reference$ind_idx) ||
        !identical(m$geno, reference$geno)) {
      stop("All inputs must contain the same genotype data, marker order, and ",
           "individual subset.", call. = FALSE)
    }
  }

  means <- lapply(models, function(m) {
    vapply(seq_len(nt), function(t) .trait_mean(m, t), numeric(1))
  })
  if (!all(vapply(means[-1L], identical, logical(1), means[[1L]]))) {
    stop("All inputs must use the same per-trait `mean`; it is applied once ",
         "to the combined phenotype.", call. = FALSE)
  }

  seed_labels <- vapply(models, function(m) {
    if (is.null(m$seed)) "NULL" else as.character(m$seed)
  }, character(1))
  seed <- models[[1]]$seed
  if (length(unique(seed_labels)) > 1) {
    warning("Inputs were built with different seeds (",
            paste(unique(seed_labels), collapse = ", "),
            "); using the first input's seed (", seed_labels[[1]], ").",
            call. = FALSE)
  }

  h2v <- .validate_proportion(h2, "h2", nt)
  h2v <- .expand_prop(h2v, nt)
  combined <- array(0, dim = c(n, nt, nr))
  for (r in seq_len(nr)) {
    raw <- Reduce(`+`, lapply(models, .genetic_matrix, rep = r))
    for (t in seq_len(nt)) {
      s <- stats::sd(raw[, t])
      if (is.finite(s) && s > 0 && h2v[t] > 0) {
        combined[, t, r] <- (raw[, t] - mean(raw[, t])) / s * sqrt(h2v[t])
      } else if (h2v[t] > 0) {
        stop("The input models have no usable genetic variation for trait ", t,
             " in replication ", r, "; h2 = ", h2v[t], " cannot be realized.",
             call. = FALSE)
      }
    }
  }

  long <- vector("list", nt * nr)
  k <- 0L
  for (r in seq_len(nr)) {
    for (t in seq_len(nt)) {
      seed_r <- .layer_seed(seed, paste0("complex_resid_t", t), r - 1L)
      e <- .seeded_residual(seed_r, n, max(0, 1 - h2v[t]))
      k <- k + 1L
      long[[k]] <- data.frame(
        id = ids, trait = paste0("Trait_", t), rep = r,
        value = combined[, t, r] + e + means[[1L]][t],
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }

  out <- models[[1]]
  out$architecture <- "complex"
  out$seed <- seed
  out$h2 <- h2v
  out$n_reps <- nr
  out$layers <- list()
  out$sources <- lapply(models, function(m) m$architecture)
  out$complex_genetic <- combined
  out$pheno <- do.call(rbind, long)
  out$var_budget <- do.call(rbind, lapply(seq_len(nt), function(t) {
    data.frame(trait = paste0("Trait_", t),
               component = c("genetic", "residual"),
               prop = c(h2v[t], 1 - h2v[t]), stringsAsFactors = FALSE)
  }))
  class(out) <- "phenotype_sim"
  out
}
