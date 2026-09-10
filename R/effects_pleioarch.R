#' PleioArch pleiotropic QTN draw and effects (exact genetic correlation)
#'
#' Implements the PleioArch algorithm for the
#' `architecture = "pleiotropy"` additive layer, generalized to any number of
#' traits. Pleiotropic effects are drawn from a multivariate normal whose
#' covariance yields the target genetic correlations `cor` in expectation;
#' trait-specific effects are drawn from univariate normals. Effects are scaled
#' to the per-genotype scale by `1 / sqrt(2 * MAF * (1 - MAF))` (the reference
#' `scaleQTNEffects()` step).
#'
#' The covariance among pleiotropic effects is
#' `Sigma[i, i] = pi_i * V_i` and `Sigma[i, j] = cor_ij * sqrt(V_i * V_j)`.
#' Since trait-specific loci contribute `(1 - pi_i) * V_i` and are independent
#' across traits, each trait's total genetic variance is `V_i` while the whole
#' covariance comes from the shared loci, so the realized correlation is
#' `cor_ij` in expectation. For two traits this reduces exactly to the
#' bivariate reference implementation.
#'
#' RNG stays in R. Returns per-trait QTN indices and effects so
#' the additive layer can carry trait-specific loci alongside the shared
#' pleiotropic ones.
#'
#' @param sim a `phenotype_sim` with `architecture = "pleiotropy"`.
#' @param n_qtn total QTNs per trait.
#' @param prop_vec per-trait additive variance proportions (length `n_traits`).
#' @param sub_seed deterministic sub-seed.
#' @return list(qtn = <list len n_traits>, effect = <list len n_traits>).
#' @keywords internal
#' @noRd
.pleio_draw <- function(sim, n_qtn, prop_vec, sub_seed) {
  a  <- sim$arch_args
  nt <- sim$n_traits
  # Only when the user actually asked for correlation control: `cor` absent
  # means the engine runs at cor = 0 and no correlation is being controlled.
  if (!is.null(a$cor)) {
    .cite_pleioarch()
  }
  R  <- .pleio_cor_matrix(sim)
  pi_vec <- .pleio_pi_vector(sim)
  # Default: no major QTN. Concentrating variance in one locus makes the
  # realized correlation hinge on that single effect draw, so it no longer
  # converges on `cor` as n_qtn grows -- which defeats the reason for using
  # this engine. A major locus is available, but it must be asked for.
  n_major  <- if (is.null(a$n_pleio_major)) 0 else a$n_pleio_major
  propMaj  <- if (is.null(a$prop_var_major)) 0 else a$prop_var_major
  n_major <- .validate_count(n_major, "n_pleio_major", minimum = 0L)
  propMaj <- .validate_proportion(propMaj, "prop_var_major", 1L)

  vg <- prop_vec

  # Variance budget (phenotypic variance assumed 1; vg are the target props).
  # Off-diagonals carry the full genetic covariance because trait-specific loci
  # are independent; diagonals carry only the pleiotropic share.
  sigma <- outer(sqrt(vg), sqrt(vg)) * R
  diag(sigma) <- pi_vec * vg

  .check_pleio_feasible(sigma, R, pi_vec)

  # QTN-count partition: pleiotropic (shared) vs trait-specific.
  pleio_n <- max(0L, round(mean(pi_vec) * n_qtn))
  pleio_n <- min(pleio_n, n_qtn)
  spec_n  <- n_qtn - pleio_n
  if (n_major > pleio_n) {
    stop("`n_pleio_major` (", n_major, ") cannot exceed the number of shared ",
         "pleiotropic QTNs implied by `pi` (", pleio_n, ").", call. = FALSE)
  }
  n_minor <- pleio_n - n_major

  # Keep the major/minor split coherent: variance allocated to a major class
  # with no loci in it would simply vanish from the budget.
  if (n_major == 0 && propMaj > 0) {
    stop("`prop_var_major` is positive but `n_pleio_major` is zero.",
         call. = FALSE)
  }
  if (n_major > 0 && propMaj <= 0) {
    stop("`n_pleio_major` is positive but `prop_var_major` is zero.",
         call. = FALSE)
  }

  cand <- .candidate_markers(sim)
  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }
  drawn <- sample(cand, pleio_n + nt * spec_n, replace = FALSE)
  pleio_idx <- drawn[seq_len(pleio_n)]
  spec_idx <- lapply(seq_len(nt), function(t) {
    if (spec_n <= 0) {
      return(integer(0))
    }
    drawn[pleio_n + (t - 1L) * spec_n + seq_len(spec_n)]
  })

  eff_major <- .draw_mvnorm(n_major, sigma * propMaj)
  eff_minor <- .draw_mvnorm(n_minor, sigma * (1 - propMaj))
  eff_spec <- lapply(seq_len(nt),
                     function(t) .draw_univariate(spec_n,
                                                  (1 - pi_vec[t]) * vg[t]))

  maf <- sim$maf
  scale_idx <- function(idx) {
    m <- maf[idx]
    s <- sqrt(2 * m * (1 - m))
    s[!is.finite(s) | s == 0] <- 1
    1 / s
  }
  pleio_scale <- if (pleio_n > 0) scale_idx(pleio_idx) else numeric(0)

  qtn <- vector("list", nt)
  effect <- vector("list", nt)
  for (t in seq_len(nt)) {
    eff_pleio <- c(eff_major[, t], eff_minor[, t]) * pleio_scale
    eff_t <- if (spec_n > 0) eff_spec[[t]] * scale_idx(spec_idx[[t]]) else
      numeric(0)
    qtn[[t]] <- c(pleio_idx, spec_idx[[t]])
    effect[[t]] <- c(eff_pleio, eff_t)
  }

  list(qtn = qtn, effect = effect)
}

#' Citation notice for the correlation-control engine, once per session
#'
#' The algorithm behind `cor` in the pleiotropy architecture is published
#' separately, so credit it the first time it is used in a session. Emitted via
#' [rlang::inform()] so it is a message, not a warning, and can be silenced with
#' `suppressMessages()`.
#' @keywords internal
#' @noRd
.cite_pleioarch <- function() {
  rlang::inform(
    paste0(.cite_main(), " please cite Prado et al. (in preparation) when ",
           "controlling the correlation in the pleiotropic architecture."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_pleioarch_citation"
  )
}

#' The primary citation, shared by every "please also cite" notice
#' @keywords internal
#' @noRd
.cite_main <- function() {
  paste0("In addition to citing \"Fernandes, S.B., Lipka, A.E. ",
         "simplePHENOTYPES: SIMulation of pleiotropic, linked and epistatic ",
         "phenotypes. BMC Bioinformatics 21, 491 (2020). ",
         "https://doi.org/10.1186/s12859-020-03804-y\"")
}

#' Per-trait pleiotropic variance share
#'
#' `pi` sets the share for every trait; `pi_target` / `pi_secondary` are the
#' two-trait spelling kept from the reference implementation.
#' @keywords internal
#' @noRd
.pleio_pi_vector <- function(sim) {
  a <- sim$arch_args
  nt <- sim$n_traits
  if (!is.null(a$pi) &&
      (!is.null(a$pi_target) || !is.null(a$pi_secondary))) {
    stop("Use either `pi` or `pi_target`/`pi_secondary`, not both.",
         call. = FALSE)
  }
  if (nt != 2L &&
      (!is.null(a$pi_target) || !is.null(a$pi_secondary))) {
    stop("`pi_target` and `pi_secondary` are the two-trait interface; use ",
         "`pi` when n_traits is not 2.", call. = FALSE)
  }
  if (!is.null(a$pi)) {
    if (!length(a$pi) %in% c(1L, nt)) {
      stop("`pi` must have length 1 or n_traits (", nt, "); got ",
           length(a$pi), ".", call. = FALSE)
    }
    p <- rep_len(a$pi, nt)
  } else {
    piT <- if (is.null(a$pi_target)) 1 else a$pi_target
    p <- rep(piT, nt)
    if (nt >= 2 && !is.null(a$pi_secondary)) {
      p[2] <- a$pi_secondary
    }
  }
  if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("Pleiotropic shares (`pi`, `pi_target`, `pi_secondary`) must be ",
         "between 0 and 1.", call. = FALSE)
  }
  p
}

#' Check that the requested correlations are attainable
#'
#' The pleiotropic covariance matrix must be positive semi-definite. For two
#' traits this is exactly the reference constraint `cor^2 <= pi_1 * pi_2`; with
#' more traits it additionally rules out mutually inconsistent correlations
#' (for example three traits that are each strongly negatively correlated).
#' @keywords internal
#' @noRd
.check_pleio_feasible <- function(sigma, R, pi_vec) {
  nt <- nrow(sigma)
  if (nt == 2) {
    lhs <- R[1, 2]^2
    rhs <- pi_vec[1] * pi_vec[2]
    if (lhs > rhs + 1e-12) {
      stop("Biological constraint violated: cor^2 (", round(lhs, 3),
           ") cannot exceed pi_target * pi_secondary (", round(rhs, 3),
           ").", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  tol <- nt * max(abs(ev)) * .Machine$double.eps
  if (min(ev) < -max(tol, 1e-10)) {
    stop("The requested `cor` is not attainable with these pleiotropic ",
         "shares: the implied genetic covariance matrix is not positive ",
         "semi-definite (smallest eigenvalue ", format(min(ev), digits = 3),
         "). Lower the correlations or raise `pi`. For two traits the ",
         "condition is cor^2 <= pi_1 * pi_2.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Multivariate-normal effect draw with per-SNP variance scaling
#' @keywords internal
#' @noRd
.draw_mvnorm <- function(n, sigma) {
  nt <- nrow(sigma)
  if (n <= 0 || all(diag(sigma) <= 0)) {
    return(matrix(0, nrow = max(n, 0), ncol = nt))
  }
  sigma_per <- sigma / n
  z <- matrix(stats::rnorm(n * nt), ncol = nt)
  z %*% chol(.nudge_pd(sigma_per))
}

#' Nudge a covariance matrix to positive definiteness for chol()
#'
#' Only bites at the feasibility boundary (for example `cor^2 == pi_1 * pi_2`,
#' where the matrix is singular and `chol()` fails outright). Deliberately not
#' `make_pd()`, which rounds to two decimals -- per-SNP covariances here are on
#' the order of 1e-3 and would be flattened to zero.
#' @keywords internal
#' @noRd
.nudge_pd <- function(s) {
  e <- eigen(s, symmetric = TRUE)
  if (min(e$values) > 0) {
    return(s)
  }
  floor_val <- max(abs(e$values)) * 1e-8
  if (!is.finite(floor_val) || floor_val <= 0) {
    floor_val <- .Machine$double.eps
  }
  vals <- pmax(e$values, floor_val)
  out <- e$vectors %*% diag(vals, nrow(s)) %*% t(e$vectors)
  (out + t(out)) / 2
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
    if (!is.numeric(cor_g) || any(!is.finite(cor_g)) ||
        any(cor_g < -1 | cor_g > 1)) {
      stop("Every entry of the `cor` matrix must be finite and between -1 ",
           "and 1.", call. = FALSE)
    }
    if (!isTRUE(all.equal(cor_g, t(cor_g), tolerance = 1e-12))) {
      stop("The `cor` matrix must be symmetric.", call. = FALSE)
    }
    if (any(abs(diag(cor_g) - 1) > 1e-12)) {
      stop("The diagonal of the `cor` matrix must equal 1.", call. = FALSE)
    }
    return(cor_g)
  }
  if (!is.numeric(cor_g) || length(cor_g) != 1L || !is.finite(cor_g) ||
      cor_g < -1 || cor_g > 1) {
    stop("`cor` must be one finite value between -1 and 1, or a valid ",
         "correlation matrix.", call. = FALSE)
  }
  R <- matrix(cor_g, nt, nt)
  diag(R) <- 1
  R
}
