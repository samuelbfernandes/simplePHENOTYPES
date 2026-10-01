# Mimic mode for simulate_transcriptome(): calibrate the generator to a user
# expression matrix (per-gene moments, the co-expression factor count and
# module strength kappa, and -- when genotypes are paired -- a per-gene GREML
# heritability *distribution*), then GENERATE new synthetic expression with
# those calibrated parameters. Mimic calibrates the expression generator only:
# it does NOT retain the reference matrix's loadings/signs (every generated gene
# loads +1 on one factor), and eQTL effects and phenotype slopes stay de novo
# (SPEC-transcriptome.md Section 5).

# Documented minimum sample sizes (see simulate_transcriptome()): below
# .TX_MIN_N_STABLE individuals the finite-sample G-R covariance makes the realized
# heritability `h2_realized` = Var(G)/Var(P) unreliable (a warning is issued); below
# .TX_MIN_N_GREML individuals mimic's per-gene GREML estimates are noise.
.TX_MIN_N_STABLE <- 30L
.TX_MIN_N_GREML  <- 100L

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
#'
#' Model reference: Kang et al. (2008) Genetics 178(3):1709-1723,
#' doi:10.1534/genetics.107.080101 (EMMA; page/equation numbers not verified).
#' Identifiability is governed by the spectrum of `K` on the subspace orthogonal to
#' the intercept (the fixed effect `1 mu`), i.e. the eigenvalues of the projected
#' matrix `P K P`, `P = I - 11'/n`, not the raw spectrum of `K`. When those `n - 1`
#' eigenvalues are all equal -- e.g. `K = 0`, `K = I`, or any `K = a (I - 11'/n) +
#' b 11'` (such as a centred GRM from a balanced marker set) -- the REML profile is
#' analytically constant and `h2` is not identifiable; every gene then gets
#' `h2 = 0` and the result carries `attr(, "identifiable") = FALSE` (never
#' optimizer noise); otherwise the plain numeric vector is returned.
#' @keywords internal
#' @noRd
.greml_h2 <- function(Y, K) {
  n <- ncol(Y)
  if (nrow(K) != n || ncol(K) != n) {
    stop(".greml_h2(): K must be n_individuals square.", call. = FALSE)
  }
  eig <- eigen(K, symmetric = TRUE)
  xi <- pmax(eig$values, 0)                     # clamp tiny negative eigenvalues
  # spectrum of K on the complement of the intercept: orthonormal basis Qc of
  # 1-perp (columns 2..n of the QR of 1), then the eigenvalues of Qc' K Qc
  Qc <- qr.Q(qr(matrix(1, n, 1L)), complete = TRUE)[, -1L, drop = FALSE]
  xp <- pmax(eigen(crossprod(Qc, K %*% Qc), symmetric = TRUE,
                   only.values = TRUE)$values, 0)
  if (length(xp) < 2L || max(xp) - min(xp) <= 1e-8 * max(1, max(xp))) {
    return(structure(rep(0, nrow(Y)), identifiable = FALSE))
  }
  U  <- eig$vectors
  omega <- as.numeric(crossprod(U, rep(1, n)))  # U' 1 (intercept in rotated space)
  h2 <- numeric(nrow(Y))
  for (g in seq_len(nrow(Y))) {
    y <- Y[g, ]
    # absolute cutoff is intentional: the sole public caller standardizes every row to unit sd first, so this only catches degenerate/non-finite rows
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
    # A wide bracket can start optimize() poorly (the profile is not guaranteed
    # to be unimodal in delta), so bracket it around a coarse-grid minimum and
    # compare against both boundaries (delta -> 0 is h2 = 1, delta -> Inf is
    # h2 = 0).
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

#' Co-expression strength (residual module fraction) from a spiked-covariance fit
#'
#' Estimates the generator's residual module fraction `kappa` from standardized
#' expression `Es` (genes x individuals) WITHOUT double counting genetic trans
#' structure and without turning pure noise into structure:
#'
#' 1. the leading `Q` eigenvalues `lambda_i` of the gene-gene correlation matrix
#'    are inverted through the Baik-Ben Arous-Peche spiked-model relation
#'    `lambda = ell (1 + y / (ell - 1))`, `y = T / n`, to population spike sizes
#'    `s_i = ell_i - 1`; eigenvalues at or below the Marchenko-Pastur edge
#'    `(1 + sqrt(y))^2` are noise and contribute 0;
#' 2. with `Q` equal-loading modules the spikes sum to `r_w (T - Q)`, so
#'    `r_w = sum(s_i) / (T - Q)` is the mean within-module correlation of the
#'    OBSERVED expression;
#' 3. the generator's within-module correlation is
#'    `kappa * mean(sqrt(1 - h2))^2 + mean(sqrt(h2 (1 - omega)))^2` (non-genetic
#'    shared module factor + the shared genetic trans factor), so
#'    `kappa = (r_w - mean(sqrt(h2 (1 - omega)))^2) / mean(sqrt(1 - h2))^2`,
#'    clamped to `[0, 1]`. `h2` are the per-gene GREML estimates and `omega` the
#'    cis fraction assumed for the generated genes (their expected value).
#'
#' Assumptions: modules of about equal size, independent individuals, and that
#' the leading `Q` spikes are the module factors. The formula is this package's
#' own estimator (the spiked-model inversion is a standard random-matrix result;
#' no page-level source is claimed). With `h2 = 0` (default) no genetic trans
#' correction is applied.
#' @keywords internal
#' @noRd
.tx_estimate_kappa <- function(Es, Q, h2 = 0, omega = 0.25) {
  Tg <- nrow(Es); n <- ncol(Es)
  if (Tg < 2L || n < 3L || Tg - Q <= 0) return(0)
  C  <- crossprod(Es) / Tg                       # n x n; same nonzero spectrum
  ev <- sort(eigen(C, symmetric = TRUE, only.values = TRUE)$values,
             decreasing = TRUE)
  lam <- ev * Tg / (n - 1)                       # T x T correlation eigenvalues
  y <- Tg / n
  edge <- (1 + sqrt(y))^2
  lam <- lam[seq_len(min(Q, length(lam)))]
  s <- vapply(lam, function(l) {
    if (!is.finite(l) || l <= edge) return(0)
    b <- l + 1 - y
    ell <- (b + sqrt(max(b^2 - 4 * l, 0))) / 2
    max(ell - 1, 0)
  }, numeric(1))
  r_w <- sum(s) / (Tg - Q)
  h2 <- pmin(pmax(as.numeric(h2), 0), 1)
  om <- mean(as.numeric(omega))
  gen_trans <- mean(sqrt(h2 * (1 - om)))^2
  non_gen   <- mean(sqrt(1 - h2))^2
  if (!is.finite(non_gen) || non_gen <= 1e-8) return(0)
  min(max((r_w - gen_trans) / non_gen, 0), 1)
}

#' Calibrate the generator to a user expression matrix (mimic mode)
#'
#' Estimates, from `mimic` (genes x individuals) aligned to the reference
#' individuals: exact per-gene mean `mu` and variance `V`; a per-gene GREML
#' heritability `h2` (REML on the GRM of `Z`); the co-expression factor count `Q`
#' (Marchenko-Pastur). It also returns the standardized matrix `Es` so the caller
#' can compute the co-expression strength `kappa` against the FINAL factor count
#' (which a user `n_factors` may override). These calibrate the *generator's
#' targets*; eQTL effects, loadings (and their signs), and phenotype slopes
#' remain de novo. The per-gene GREML values are a DISTRIBUTION-level
#' calibration: for sparse eQTL architectures they do not correspond to the
#' truth gene by gene (the estimate is noisy, with sd about 0.1 at n = 280).
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
  if (sim$n_ind < .TX_MIN_N_GREML) {
    warning("simulate_transcriptome(): `mimic` calibrates per-gene heritability ",
            "by GREML on only ", sim$n_ind, " individuals (< ", .TX_MIN_N_GREML,
            "). Per-gene GREML estimates are then mostly sampling noise: treat ",
            "`calibration$h2` as a distribution-level calibration, not a ",
            "gene-by-gene estimate (with n = 280 the per-gene sd is already ",
            "about 0.1).", call. = FALSE)
  }
  K  <- .tx_grm(Z)
  h2r <- .greml_h2(Es, K)
  if (identical(attr(h2r, "identifiable"), FALSE)) {
    warning("simulate_transcriptome(): the genomic relationship matrix has no ",
            "variation across individuals once the intercept is projected out, ",
            "so per-gene GREML heritability is not identifiable; `mimic` sets ",
            "calibration$h2 to 0 for every gene.", call. = FALSE)
  }
  h2 <- as.numeric(h2r)
  Q  <- .tx_estimate_factors(Es)          # Marchenko-Pastur factor count
  list(E = E, mu = mu, V = V, h2 = h2, Q = Q, Es = Es,
       T_genes = nrow(E), gene_ids = gene_ids)
}
