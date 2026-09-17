# Mimic mode for simulate_transcriptome(): calibrate the generator to a user
# expression matrix (per-gene moments, co-expression rank/strength, and -- when
# genotypes are paired -- a per-gene GREML heritability), then GENERATE new
# synthetic expression with those calibrated parameters. Mimic calibrates the
# expression generator only; eQTL effects and phenotype slopes stay de novo
# (SPEC-transcriptome.md Section 5).

#' Genomic relationship matrix from reference-centered dosages
#'
#' `K = Z Z' / m` on the reference-centered dosages `Z` (individuals x markers),
#' then normalized to a mean diagonal of 1 so a variance-component ratio maps
#' directly to heritability. Pure R; used only by the GREML calibrator.
#' @keywords internal
#' @noRd
.tx_grm <- function(Z) {
  m <- ncol(Z)
  K <- tcrossprod(Z) / m                       # n x n
  d <- mean(diag(K))
  if (is.finite(d) && d > 0) K <- K / d
  K
}

#' Per-gene GREML heritability (EMMA-style single-component REML)
#'
#' For each gene `y` (a row of `Y`) fits `y = 1 mu + u + e`, `u ~ N(0, sigma_g^2 K)`,
#' `e ~ N(0, sigma_e^2 I)`, and returns `h2 = sigma_g^2 / (sigma_g^2 + sigma_e^2)`.
#' One eigendecomposition of `K` is shared across genes; each gene is a 1-D REML
#' optimization over `delta = sigma_e^2 / sigma_g^2` (Kang et al. 2008, EMMA). The
#' GRM has mean diagonal 1, so `h2 = 1 / (1 + delta_hat)`. This calibrates the
#' target `h2` distribution; it does not fit individual eQTL effects.
#' @keywords internal
#' @noRd
.greml_h2 <- function(Y, K) {
  n <- ncol(Y)
  if (nrow(K) != n || ncol(K) != n) {
    stop(".greml_h2(): K must be n_individuals square.", call. = FALSE)
  }
  eig <- eigen(K, symmetric = TRUE)
  xi <- pmax(eig$values, 0)                     # clamp tiny negative eigenvalues
  U  <- eig$vectors
  omega <- as.numeric(crossprod(U, rep(1, n)))  # U' 1 (intercept in rotated space)
  h2 <- numeric(nrow(Y))
  for (g in seq_len(nrow(Y))) {
    y <- Y[g, ]
    if (!is.finite(stats::sd(y)) || stats::sd(y) < 1e-12) { h2[g] <- 0; next }
    eta <- as.numeric(crossprod(U, y))          # U' y
    negll <- function(ld) {                     # ld = log(delta); minimize -REML
      d <- exp(ld); w <- 1 / (xi + d)
      A <- sum(w * omega^2)
      b <- sum(w * omega * eta)
      cc <- sum(w * eta^2)
      rss <- cc - b^2 / A                        # GLS residual sum of squares
      if (!is.finite(rss) || rss <= 0 || !is.finite(A) || A <= 0) return(1e10)
      0.5 * ((n - 1) * log(rss) + sum(log(xi + d)) + log(A))
    }
    # The single-variance-component REML profile is unimodal in delta, but a wide
    # bracket can still start optimize() poorly, so bracket it around a coarse-grid
    # minimum and compare against both boundaries (delta -> 0 is h2 = 1, delta ->
    # Inf is h2 = 0).
    grid <- seq(-12, 12, length.out = 25)
    gv <- vapply(grid, negll, numeric(1))
    j  <- which.min(gv)
    lo <- grid[max(1, j - 1)]; hi <- grid[min(length(grid), j + 1)]
    opt <- stats::optimize(negll, c(lo, hi))
    cand <- c(interior = opt$objective, h2_1 = negll(-20), h2_0 = negll(20))
    best <- which.min(cand)
    h2[g] <- if (best == 2L) 1 else if (best == 3L) 0 else
      min(max(1 / (1 + exp(opt$minimum)), 0), 1)
    if (!is.finite(h2[g])) h2[g] <- 0
  }
  h2
}

#' Number of co-expression factors by a Marchenko-Pastur edge
#'
#' Counts the eigenvalues of the individual covariance of standardized expression
#' `(1/T) Es' Es` that exceed the Marchenko-Pastur upper edge `(1 + sqrt(n/T))^2`
#' expected under pure noise -- a standard signal-factor count. Capped to the
#' generator's factor range.
#' @keywords internal
#' @noRd
.tx_estimate_factors <- function(Es) {
  Tg <- nrow(Es); n <- ncol(Es)
  if (Tg < 2L || n < 3L) return(1L)
  C  <- crossprod(Es) / Tg                       # n x n
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  lambda_plus <- (1 + sqrt(n / Tg))^2            # noise upper edge (unit-variance rows)
  q <- sum(ev > lambda_plus)
  as.integer(min(max(q, 1L), min(50L, n - 2L)))
}

#' Co-expression strength: share of the standardized spectrum in the top factors
#'
#' A documented proxy for the generator's `residual_module_fraction` (`kappa`):
#' the fraction of the standardized-expression covariance spectrum carried by the
#' leading `Q` components. It is a strength summary, not an identifiable estimate
#' of the residual module fraction (the genetic and non-genetic trans structure
#' are not separated on expression alone).
#' @keywords internal
#' @noRd
.tx_estimate_kappa <- function(Es, Q) {
  Tg <- nrow(Es); n <- ncol(Es)
  if (Tg < 2L || n < 2L) return(0)
  C  <- crossprod(Es) / Tg
  ev <- sort(eigen(C, symmetric = TRUE, only.values = TRUE)$values, decreasing = TRUE)
  tot <- sum(ev)
  if (!is.finite(tot) || tot <= 0) return(0)
  top <- sum(ev[seq_len(min(Q, length(ev)))])
  min(max(top / tot, 0), 1)
}

#' Calibrate the generator to a user expression matrix (mimic mode)
#'
#' Estimates, from `mimic` (genes x individuals) aligned to the reference
#' individuals: exact per-gene mean `mu` and variance `V`; a per-gene GREML
#' heritability `h2` (REML on the GRM of `Z`); the co-expression factor count `Q`
#' (Marchenko-Pastur). It also returns the standardized matrix `Es` so the caller
#' can compute the co-expression strength `kappa` against the FINAL factor count
#' (which a user `n_factors` may override). These calibrate the *generator's
#' targets*; eQTL effects, loadings, and phenotype slopes remain de novo.
#' @keywords internal
#' @noRd
.tx_mimic_calibrate <- function(mimic, sim, Z) {
  if (!is.matrix(mimic) || !is.numeric(mimic)) {
    stop("simulate_transcriptome(): `mimic` must be a numeric genes-by-",
         "individuals matrix.", call. = FALSE)
  }
  if (nrow(mimic) < 1L) {
    stop("simulate_transcriptome(): `mimic` has no genes.", call. = FALSE)
  }
  if (any(!is.finite(mimic))) {
    stop("simulate_transcriptome(): `mimic` has non-finite values; impute or ",
         "remove them first.", call. = FALSE)
  }
  cn <- colnames(mimic)
  if (is.null(cn)) {
    if (ncol(mimic) != sim$n_ind) {
      stop("simulate_transcriptome(): `mimic` has no individual (column) names ",
           "and its ", ncol(mimic), " columns do not match the ", sim$n_ind,
           " individuals; name its columns or match the order.", call. = FALSE)
    }
    E <- mimic
    colnames(E) <- sim$ids
  } else {
    if (anyDuplicated(cn)) {
      stop("simulate_transcriptome(): `mimic` has duplicate individual (column) ",
           "names.", call. = FALSE)
    }
    m <- match(sim$ids, cn)
    if (anyNA(m)) {
      stop("simulate_transcriptome(): `mimic` is missing individual(s): ",
           paste(utils::head(sim$ids[is.na(m)], 5), collapse = ", "), ".",
           call. = FALSE)
    }
    E <- mimic[, m, drop = FALSE]
  }
  # Gene ids key the returned truth tables (cis_eqtl/loadings), so present names
  # must be unique and non-empty; absent names (NULL rownames) are auto-assigned,
  # matching how simulate_phenotype(expression=) handles gene labels.
  gene_ids <- rownames(E)
  if (is.null(gene_ids)) {
    gene_ids <- paste0("gene", seq_len(nrow(E)))
  } else if (anyNA(gene_ids) || any(!nzchar(gene_ids)) || anyDuplicated(gene_ids)) {
    stop("simulate_transcriptome(): `mimic` has empty, NA, or duplicate gene ",
         "(row) names; give unique names or omit row names to auto-assign them.",
         call. = FALSE)
  }
  mu <- rowMeans(E)
  V  <- apply(E, 1L, stats::var)
  Es <- t(apply(E, 1L, function(r) {
    s <- stats::sd(r); if (is.finite(s) && s > 0) (r - mean(r)) / s else rep(0, length(r))
  }))
  K  <- .tx_grm(Z)
  h2 <- .greml_h2(Es, K)
  Q  <- .tx_estimate_factors(Es)          # Marchenko-Pastur factor count
  list(E = E, mu = mu, V = V, h2 = h2, Q = Q, Es = Es,
       T_genes = nrow(E), gene_ids = gene_ids)
}
