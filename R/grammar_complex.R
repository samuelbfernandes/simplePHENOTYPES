#' Combine single-architecture models under a common heritability
#'
#' `complex_phenotypes()` merges two or more realized `phenotype_sim` objects
#' built on the same genotypes with the same number of traits. The inputs'
#' genetic values are summed (so each contributes in proportion to its genetic
#' variance), a common residual is applied to reach the target heritability
#' `h2`, and the inputs' individual residuals are discarded. This recreates
#' partial pleiotropy by combining, for example, a "pleiotropy" model with an
#' "independent" model.
#'
#' @param ... two or more `phenotype_sim` objects.
#' @param h2 target heritability (scalar or length `n_traits`).
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
  n  <- nrow(models[[1]]$G)
  ids <- rownames(models[[1]]$G)
  for (m in models) {
    if (m$n_traits != nt) {
      stop("All inputs must have the same n_traits.", call. = FALSE)
    }
    if (nrow(m$G) != n || !identical(rownames(m$G), ids)) {
      stop("All inputs must be built on the same genotypes.", call. = FALSE)
    }
  }

  seeds <- unique(unlist(lapply(models, function(m) m$seed)))
  seed <- models[[1]]$seed
  if (length(seeds) > 1) {
    warning("Inputs were built with different seeds (",
            paste(seeds, collapse = ", "), "); using the first input's seed (",
            seed, ").", call. = FALSE)
  }

  combined <- Reduce(`+`, lapply(models, .genetic_matrix))   # n x nt, weighted
  h2v <- .expand_prop(h2, nt)

  long <- vector("list", nt)
  for (t in seq_len(nt)) {
    g <- combined[, t]
    s <- stats::sd(g)
    if (is.finite(s) && s > 0 && h2v[t] > 0) {
      g <- (g - mean(g)) / s * sqrt(h2v[t])
    } else {
      g <- rep(0, n)
    }
    seed_r <- .layer_seed(seed, paste0("complex_resid_t", t), 0L)
    e <- .seeded_residual(seed_r, n, max(0, 1 - h2v[t]))
    long[[t]] <- data.frame(
      id = ids, trait = paste0("Trait_", t), rep = 1L,
      value = g + e, stringsAsFactors = FALSE, row.names = NULL
    )
  }

  out <- models[[1]]
  out$architecture <- "complex"
  out$seed <- seed
  out$layers <- list()
  out$pleio_cor <- NULL
  out$sources <- lapply(models, function(m) m$architecture)
  out$pheno <- do.call(rbind, long)
  out$var_budget <- do.call(rbind, lapply(seq_len(nt), function(t) {
    data.frame(trait = paste0("Trait_", t),
               component = c("genetic", "residual"),
               prop = c(h2v[t], 1 - h2v[t]), stringsAsFactors = FALSE)
  }))
  class(out) <- "phenotype_sim"
  out
}
