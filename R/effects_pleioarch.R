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
#' Provenance: the associated manuscript (Prado et al.) is in preparation and has
#' no public reference, so the authoritative definition of this algorithm is the
#' bundled reference implementation in `context/PleioArch-main/Functions/`
#' (`simulateEffects.R`, `scaleQTNEffects.R`), against which this port is
#' checked. The covariance construction and the `1/sqrt(2*MAF*(1-MAF))` scaling
#' documented here are self-contained; do not cite them to a published paper
#' until one exists.
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

  # A correlation with a zero-variance trait is undefined: cov = cor*sqrt(Vi*Vj)
  # is 0 and the realized correlation is 0/0 = NA, so cor cannot be attained.
  # Reject a nonzero requested off-diagonal against any trait whose additive
  # proportion is zero, rather than silently returning NA correlations.
  if (any(vg <= 0)) {
    nt_ <- length(vg)
    for (i in seq_len(nt_ - 1L)) {
      for (j in seq(i + 1L, nt_)) {
        if (vg[i] <= 0 || vg[j] <= 0) {
          zero_t <- if (vg[i] <= 0) i else j
          other  <- if (zero_t == i) j else i
          if (abs(R[i, j]) > 0) {
            stop("architecture = \"pleiotropy\": a nonzero `cor` (", signif(R[i, j], 3),
                 ") was requested between trait ", i, " and trait ", j,
                 ", but trait ", zero_t, " has zero additive ",
                 "variance (prop = 0), so that correlation is undefined and cannot ",
                 "be realized. Give every correlated trait a positive additive ",
                 "`prop`.", call. = FALSE)
          } else {
            # cor = 0 is undefined too when a trait is invariant: the realized
            # correlation is 0/0 = NA, not 0. Warn rather than report a silent NA.
            warning("architecture = \"pleiotropy\": trait ", zero_t, " has zero ",
                    "additive variance (prop = 0), so its genetic correlation with ",
                    "trait ", other, " is undefined and will be reported as NA (a ",
                    "requested cor = 0 cannot be realized as 0 here). Give every ",
                    "correlated trait a positive additive `prop`.", call. = FALSE)
          }
        }
      }
    }
  }

  # Variance budget (phenotypic variance assumed 1; vg are the target props).
  # Off-diagonals carry the full genetic covariance because trait-specific loci
  # are independent; diagonals carry only the pleiotropic share.
  sigma <- outer(sqrt(vg), sqrt(vg)) * R
  diag(sigma) <- pi_vec * vg

  .check_pleio_feasible(sigma, R, pi_vec)

  # QTN-count partition: pleiotropic (shared) vs trait-specific.
  pleio_n <- max(0L, round(mean(pi_vec) * n_qtn))
  pleio_n <- min(pleio_n, n_qtn)
  # A single shared (pleiotropic) QTN makes both traits scalar multiples of the
  # same dosage vector, so their realized genetic correlation is exactly +/-1
  # regardless of the requested `cor` -- one locus cannot carry an intermediate
  # target. Warn rather than silently realize +/-1 for an intermediate `cor`.
  if (pleio_n == 1L && any(abs(R[upper.tri(R)]) > 0 & abs(R[upper.tri(R)]) < 1)) {
    warning("architecture = \"pleiotropy\": only one shared (pleiotropic) QTN ",
            "results from n_qtn = ", n_qtn, ", pi = ",
            paste(round(pi_vec, 3), collapse = ", "),
            "; a single shared locus makes the realized genetic correlation ",
            "exactly +/-1, so an intermediate target `cor` cannot be realized. ",
            "Increase n_qtn (or pi) so several shared loci carry the covariance.",
            call. = FALSE)
  }
  # `pi` defines the shared/trait-specific variance partition whether or not a
  # correlation is requested, so its feasibility is judged against `pi`, not
  # against `cor`. If pi implies a shared class (mean(pi) > 0) but the QTN-count
  # rounding leaves zero shared loci, that shared variance -- and, when cor != 0,
  # the whole requested covariance -- would be lost; n_qtn is too small to
  # represent the partition. Error rather than force one (collinear) locus or
  # silently drop the share.
  if (mean(pi_vec) > 0 && n_qtn > 0L && pleio_n < 1L) {
    stop("architecture = \"pleiotropy\": n_qtn = ", n_qtn, " and pi = ",
         paste(round(pi_vec, 3), collapse = ", "),
         " round the shared (pleiotropic) class to zero, so the requested ",
         "shared variance", if (any(abs(R[upper.tri(R)]) > 0))
           " and correlation" else "", " cannot be represented. Increase ",
         "n_qtn (or pi).", call. = FALSE)
  }
  spec_n  <- n_qtn - pleio_n
  # Trait-specific variance ((1 - pi) V) needs at least one trait-specific locus
  # to carry it; if the partition leaves none while pi < 1, that variance would
  # silently vanish. Error rather than lose it.
  if (spec_n < 1L && any(pi_vec < 1)) {
    stop("architecture = \"pleiotropy\": n_qtn = ", n_qtn, " leaves no ",
         "trait-specific QTNs, but pi = ", paste(round(pi_vec, 3), collapse = ", "),
         " (< 1) requests trait-specific variance. Increase n_qtn so both ",
         "shared and trait-specific loci fit.", call. = FALSE)
  }
  if (n_major > pleio_n) {
    stop("`n_pleio_major` (", n_major, ") cannot exceed the number of shared ",
         "pleiotropic QTNs implied by `pi` (", pleio_n, ").", call. = FALSE)
  }
  n_minor <- pleio_n - n_major
  # With a major/minor split, prop_var_major < 1 allocates (1 - prop_var_major)
  # of the pleiotropic variance to minor loci. If there are none, that share
  # would silently vanish (attenuating the realized correlation); reject rather
  # than lose it.
  if (n_major > 0L && n_minor == 0L && propMaj < 1) {
    stop("`prop_var_major` (", propMaj, ") < 1 assigns ", round(1 - propMaj, 3),
         " of the pleiotropic variance to minor loci, but none remain ",
         "(n_pleio_major equals the shared-QTN count implied by `pi`: ",
         pleio_n, "). Increase `n_qtn`, lower `n_pleio_major`, or set ",
         "`prop_var_major = 1`.", call. = FALSE)
  }

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
#' The algorithm behind `cor` in the pleiotropy architecture is described in a
#' manuscript still in preparation (Prado et al.), so credit that the first time
#' it is used in a session; the concrete definition is the bundled reference
#' implementation in `context/PleioArch-main/`. Emitted via [rlang::inform()] so
#' it is a message, not a warning, and can be silenced with `suppressMessages()`.
#' @keywords internal
#' @noRd
.cite_pleioarch <- function() {
  rlang::inform(
    .cite_main("Prado et al. (in preparation), when controlling the ",
               "correlation in the pleiotropic architecture."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_pleioarch_citation"
  )
}

#' The primary citation notice, shared by every "please also cite" message
#'
#' Builds the multi-line notice: the lead-in, the main simplePHENOTYPES
#' reference, then the feature-specific reference(s) passed in `...`.
#' @keywords internal
#' @noRd
.cite_main <- function(...) {
  paste0("In addition to citing:\n",
         "Fernandes, S.B., Lipka, A.E. simplePHENOTYPES: SIMulation of ",
         "pleiotropic, linked and epistatic phenotypes. BMC Bioinformatics 21, ",
         "491 (2020). https://doi.org/10.1186/s12859-020-03804-y\n",
         "Please also cite:\n",
         ...)
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
    # Reject at the level of floating-point roundoff only, scaled to the operands
    # actually compared (not an absolute floor like max(1, rhs), which would let
    # cor^2 > 0 slip through when pi1*pi2 = 0, e.g. pi = c(0.25, 0), where any
    # nonzero correlation is impossible). cor^2 and pi1*pi2 are each one product,
    # so a few ULPs of their magnitude is the right slack.
    tol <- 8 * .Machine$double.eps * max(lhs, rhs)
    if (lhs - rhs > tol) {
      stop("Biological constraint violated: cor^2 (", signif(lhs, 6),
           ") cannot exceed pi_target * pi_secondary (", signif(rhs, 6),
           ").", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  # Reject a negative eigenvalue unless it is at the level of eigen-decomposition
  # roundoff, which is a small multiple of (matrix scale) x (dimension) x eps --
  # NOT sqrt(eps) (~1e-8), which is astronomically larger than real roundoff and
  # would wave a genuinely indefinite matrix (e.g. min eigenvalue -1e-8) straight
  # through to the effect draw. The tolerance is scale-relative so it also works
  # for small-variance matrices; an exactly singular feasible boundary (min
  # eigenvalue ~ 0) still passes and is drawn exactly by .draw_mvnorm().
  scale <- max(abs(ev))
  tol <- scale * nt * 64 * .Machine$double.eps
  if (scale > 0 && min(ev) < -tol) {
    stop("The requested `cor` is not attainable with these pleiotropic ",
         "shares: the implied genetic covariance matrix is not positive ",
         "semi-definite (smallest eigenvalue ", format(min(ev), digits = 3),
         "). Lower the correlations or raise `pi`. For two traits the ",
         "condition is cor^2 <= pi_1 * pi_2.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Multivariate-normal effect draw with per-SNP variance scaling
#'
#' Uses the symmetric eigendecomposition square root rather than a Cholesky of a
#' nudged matrix. Feasibility is already checked upstream, so `sigma` is PSD;
#' the eigen root reproduces it exactly even when it is singular (the feasibility
#' boundary `cor^2 = pi_1 * pi_2`, where `chol()` fails and any nudge would
#' silently perturb a valid request). Only numerical-roundoff negative
#' eigenvalues are floored at zero.
#' @keywords internal
#' @noRd
.draw_mvnorm <- function(n, sigma) {
  nt <- nrow(sigma)
  if (n <= 0 || all(diag(sigma) <= 0)) {
    return(matrix(0, nrow = max(n, 0), ncol = nt))
  }
  sigma_per <- sigma / n
  z <- matrix(stats::rnorm(n * nt), ncol = nt)
  e <- eigen(sigma_per, symmetric = TRUE)
  root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  z %*% root
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
