#' PleioArch pleiotropic QTN draw and effects (exact genetic correlation)
#'
#' Implements the PleioArch algorithm (DECISION-007, SPEC §13) for the
#' `architecture = "pleiotropy"` additive layer. For exactly two traits it draws
#' pleiotropic effects from a bivariate normal whose covariance yields the target
#' genetic correlation `cor` in expectation; trait-specific effects are drawn
#' from univariate normals. Effects are scaled to the per-genotype scale by
#' `1 / sqrt(2 * MAF * (1 - MAF))` (the reference `scaleQTNEffects()` step).
#'
#' RNG stays in R (DECISION-006). Returns per-trait QTN indices and effects so
#' the additive layer can carry trait-specific loci alongside the shared
#' pleiotropic ones.
#'
#' @param sim a `phenotype_sim` (must be 2-trait pleiotropy).
#' @param n_qtn total QTNs per trait.
#' @param prop_vec per-trait additive variance proportions (length 2).
#' @param sub_seed deterministic sub-seed.
#' @return list(qtn = <list len 2>, effect = <list len 2>).
#' @keywords internal
#' @noRd
.pleio_draw <- function(sim, n_qtn, prop_vec, sub_seed) {
  a <- sim$arch_args
  cor_g    <- if (is.null(a$cor)) 0 else a$cor
  piT      <- if (is.null(a$pi_target)) 1 else a$pi_target
  piS      <- if (is.null(a$pi_secondary)) 1 else a$pi_secondary
  n_major  <- if (is.null(a$n_pleio_major)) 1 else a$n_pleio_major
  propMaj  <- if (is.null(a$prop_var_major)) 0.5 else a$prop_var_major

  if (cor_g^2 > piT * piS + 1e-12) {
    stop("Biological constraint violated: cor^2 (", round(cor_g^2, 3),
         ") cannot exceed pi_target * pi_secondary (", round(piT * piS, 3),
         ").", call. = FALSE)
  }

  vgT <- prop_vec[1]
  vgS <- prop_vec[2]

  # QTN-count partition: pleiotropic (shared) vs trait-specific.
  pleio_n <- max(1L, round(mean(c(piT, piS)) * n_qtn))
  pleio_n <- min(pleio_n, n_qtn)
  spec_nT <- n_qtn - pleio_n
  spec_nS <- n_qtn - pleio_n
  n_major <- min(n_major, pleio_n)
  n_minor <- pleio_n - n_major

  cand <- .candidate_markers(sim)
  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }
  drawn <- sample(cand, pleio_n + spec_nT + spec_nS, replace = FALSE)
  pleio_idx <- drawn[seq_len(pleio_n)]
  specT_idx <- if (spec_nT > 0) drawn[pleio_n + seq_len(spec_nT)] else integer(0)
  specS_idx <- if (spec_nS > 0)
    drawn[pleio_n + spec_nT + seq_len(spec_nS)] else integer(0)

  # Variance budget (phenotypic variance assumed 1; vg* are the target props).
  vpT <- piT * vgT
  vpS <- piS * vgS
  covP <- cor_g * sqrt(vgT) * sqrt(vgS)
  vsT <- (1 - piT) * vgT
  vsS <- (1 - piS) * vgS

  sigma_major <- matrix(c(vpT, covP, covP, vpS), 2) * propMaj
  sigma_minor <- matrix(c(vpT, covP, covP, vpS), 2) * (1 - propMaj)

  eff_major <- .draw_bivariate(n_major, sigma_major)
  eff_minor <- .draw_bivariate(n_minor, sigma_minor)
  eff_specT <- .draw_univariate(spec_nT, vsT)
  eff_specS <- .draw_univariate(spec_nS, vsS)

  maf <- sim$maf
  scale_idx <- function(idx) {
    m <- maf[idx]
    s <- sqrt(2 * m * (1 - m))
    s[!is.finite(s) | s == 0] <- 1
    1 / s
  }

  pleio_scale <- if (pleio_n > 0) scale_idx(pleio_idx) else numeric(0)
  effT_pleio <- c(eff_major[, 1], eff_minor[, 1]) * pleio_scale
  effS_pleio <- c(eff_major[, 2], eff_minor[, 2]) * pleio_scale
  effT_spec <- if (spec_nT > 0) eff_specT * scale_idx(specT_idx) else numeric(0)
  effS_spec <- if (spec_nS > 0) eff_specS * scale_idx(specS_idx) else numeric(0)

  list(
    qtn = list(
      c(pleio_idx, specT_idx),
      c(pleio_idx, specS_idx)
    ),
    effect = list(
      c(effT_pleio, effT_spec),
      c(effS_pleio, effS_spec)
    )
  )
}

#' Bivariate-normal effect draw with per-SNP variance scaling
#' @keywords internal
#' @noRd
.draw_bivariate <- function(n, sigma) {
  if (n <= 0 || all(diag(sigma) <= 0)) {
    return(matrix(0, nrow = max(n, 0), ncol = 2))
  }
  sigma_per <- sigma / n
  z <- matrix(stats::rnorm(n * 2), ncol = 2)
  z %*% chol(.make_pd_2x2(sigma_per))
}

#' Univariate-normal effect draw with per-SNP variance scaling
#' @keywords internal
#' @noRd
.draw_univariate <- function(n, v) {
  if (n <= 0) {
    return(numeric(0))
  }
  stats::rnorm(n, mean = 0, sd = sqrt(max(v, 0) / n))
}

#' Target genetic-correlation matrix for the >2-trait Cholesky fallback
#'
#' Uses a single scalar `cor` for every trait pair, or a supplied n_traits x
#' n_traits matrix passed as `cor`.
#' @keywords internal
#' @noRd
.pleio_cor_matrix <- function(sim) {
  nt <- sim$n_traits
  cor_g <- sim$arch_args$cor
  if (is.null(cor_g)) cor_g <- 0
  if (is.matrix(cor_g)) {
    if (!all(dim(cor_g) == c(nt, nt))) {
      stop("`cor` matrix must be ", nt, " x ", nt, ".", call. = FALSE)
    }
    return(cor_g)
  }
  R <- matrix(cor_g, nt, nt)
  diag(R) <- 1
  R
}

#' Nudge a 2x2 covariance toward positive-definiteness for chol()
#' @keywords internal
#' @noRd
.make_pd_2x2 <- function(s) {
  eps <- 1e-12
  if (s[1, 1] <= 0) s[1, 1] <- eps
  if (s[2, 2] <= 0) s[2, 2] <- eps
  d <- s[1, 1] * s[2, 2] - s[1, 2] * s[2, 1]
  if (d <= 0) {
    s[1, 2] <- s[2, 1] <- sign(s[1, 2]) *
      sqrt(s[1, 1] * s[2, 2]) * (1 - 1e-6)
  }
  s
}
