# Optional RNA-seq count observation layer for a transcriptome_sim. The generator
# works on a normalized ~Gaussian latent scale; observe_counts() maps that latent
# expression to negative-binomial counts. Heritability on the count scale is NOT
# equal to the latent-scale heritability (a nonlinear, mean-dependent link), which
# is the point of keeping counts a separate, optional layer (SPEC-transcriptome.md
# Section 7).

#' Draw RNA-seq counts from a simulated transcriptome
#'
#' Maps the normalized latent expression of a `transcriptome_sim` to integer
#' counts through a negative-binomial observation model. For gene `g` and
#' individual `i`, with `z_{gi}` the latent expression standardized per gene,
#' \deqn{\log \mu_{gi} = \log L_i + \alpha_g + \sigma_g z_{gi}, \qquad
#'       Y_{gi} \sim \mathrm{NB}(\mu_{gi}, \phi_g),}
#' where `L_i` is a per-individual library-size (sequencing-depth) factor,
#' `alpha_g` a per-gene log baseline count, `sigma_g` the latent-to-log-mean
#' coupling, and `phi_g` the NB dispersion (`Var = mu + phi mu^2`). Standardizing
#' the latent expression makes `sigma_g` a per-standard-deviation log-fold change,
#' comparable across genes and unaffected by any `mimic` rescaling.
#'
#' **Count-scale heritability differs from the latent heritability.** The link is
#' nonlinear and mean-dependent, so the fraction of count variance explained by the
#' genome is not the latent `h2`; treat counts as an observation layer for
#' count-based methods, not as a re-parameterization of the genetic model.
#'
#' @param tx a `transcriptome_sim` from [simulate_transcriptome()].
#' @param library_size per-individual depth factor `L_i`: one positive number, or
#'   one per individual. Default `1` (no depth variation).
#' @param baseline per-gene log baseline count `alpha_g`: one number or one per
#'   gene. Default `3` (a mean of about `exp(3)` counts at `z = 0`).
#' @param coupling latent-to-log-mean slope `sigma_g`: one non-negative number or
#'   one per gene. Default `0.5`.
#' @param dispersion NB dispersion `phi_g` (`Var = mu + phi mu^2`): one
#'   non-negative number or one per gene. `0` gives a Poisson limit. Default `0.1`.
#' @param seed optional seed; the RNG state is restored afterwards.
#' @return the `transcriptome_sim` with a `counts` matrix (genes x individuals,
#'   integer) and a `count_model` list (the per-gene/per-individual parameters and
#'   the realized mean matrix `mu`).
#' @seealso [simulate_transcriptome()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' tx <- simulate_transcriptome(SNP55K_maize282_maf04, n_genes = 50, seed = 1)
#' tx <- observe_counts(tx, seed = 2)
#' tx$counts[1:5, 1:5]
observe_counts <- function(tx, library_size = 1, baseline = 3, coupling = 0.5,
                           dispersion = 0.1, seed = NULL) {
  if (!inherits(tx, "transcriptome_sim")) {
    stop("observe_counts(): `tx` must be a transcriptome_sim from ",
         "simulate_transcriptome().", call. = FALSE)
  }
  seed <- .validate_seed(seed)
  Tg <- tx$n_genes
  n  <- tx$n_ind
  E  <- tx$expression

  L     <- .tx_count_par(library_size, n,  "library_size", positive = TRUE)
  alpha <- .tx_count_par(baseline,     Tg, "baseline")
  sigma <- .tx_count_par(coupling,     Tg, "coupling", nonneg = TRUE)
  phi   <- .tx_count_par(dispersion,   Tg, "dispersion", nonneg = TRUE)

  # standardize latent expression per gene (constant gene -> flat, sigma has no
  # effect there); mu = exp(log L_i + alpha_g + sigma_g z_gi).
  Z <- t(apply(E, 1L, function(r) {
    s <- stats::sd(r); if (is.finite(s) && s > 0) (r - mean(r)) / s else rep(0, length(r))
  }))
  logmu <- matrix(alpha, Tg, n) + Z * sigma +           # per-gene terms (recycled by row)
    matrix(log(L), Tg, n, byrow = TRUE)                 # per-individual depth
  mu <- exp(logmu)
  if (any(!is.finite(mu))) {
    stop("observe_counts(): the count mean overflowed (log-mean too large). ",
         "Lower `baseline`, `coupling`, or `library_size` so ",
         "log L_i + alpha_g + sigma_g z_gi stays within exp()'s range.",
         call. = FALSE)
  }

  draw <- function() {
    Y <- matrix(0L, Tg, n, dimnames = dimnames(E))
    for (g in seq_len(Tg)) {
      Y[g, ] <- if (phi[g] <= 0) {
        stats::rpois(n, mu[g, ])                        # Poisson limit
      } else {
        stats::rnbinom(n, mu = mu[g, ], size = 1 / phi[g])
      }
    }
    Y
  }
  counts <- if (is.null(seed)) draw() else {
    old <- .Random.seed_safe(); set.seed(seed); on.exit(.restore_seed(old)); draw()
  }
  dimnames(mu) <- dimnames(E)

  tx$counts <- counts
  tx$count_model <- list(library_size = L, baseline = alpha, coupling = sigma,
                         dispersion = phi, mu = mu, seed = seed)
  tx
}

#' Validate/recycle a per-gene or per-individual count parameter
#' @keywords internal
#' @noRd
.tx_count_par <- function(x, len, name, positive = FALSE, nonneg = FALSE) {
  if (!is.numeric(x) || !length(x) %in% c(1L, len) || any(!is.finite(x))) {
    stop("observe_counts(): `", name, "` must be a finite number, or one per ",
         if (len > 0) "element" else "", " (length 1 or ", len, ").", call. = FALSE)
  }
  if (positive && any(x <= 0)) {
    stop("observe_counts(): `", name, "` must be positive.", call. = FALSE)
  }
  if (nonneg && any(x < 0)) {
    stop("observe_counts(): `", name, "` must be non-negative.", call. = FALSE)
  }
  if (length(x) == 1L) rep(x, len) else x
}
