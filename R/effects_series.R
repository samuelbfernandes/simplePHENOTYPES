#' Within-layer effect series (geometric default or custom)
#'
#' Parity-critical effect generation for the v2 grammar. Deterministic given its
#' inputs; the absolute magnitude is normalized away by the per-layer
#' proportion scaling in `.realize_phenotype()`, so only the *relative* effects
#' across QTNs matter.
#'
#' @param n number of effects to generate (number of QTNs / pairs).
#' @param dist within-layer distribution; only "geometric" is supported.
#' @param effect optional override. A scalar is treated as the geometric base;
#'   a length-`n` vector is used verbatim as a custom series (v1
#'   `sim_method = "custom"`).
#' @return numeric vector of length `n`.
#' @keywords internal
#' @noRd
.effect_series <- function(n, dist = "geometric", effect = NULL) {
  if (n <= 0) {
    return(numeric(0))
  }
  if (!is.null(effect) && length(effect) == n && n > 1) {
    return(as.numeric(effect))
  }
  if (!is.null(effect) && length(effect) == 1) {
    base <- as.numeric(effect)
  } else if (!is.null(effect) && length(effect) == n) {
    return(as.numeric(effect))
  } else if (!is.null(effect)) {
    stop(
      "`effect` must be either a single value (geometric base) or a vector of ",
      "length n_qtn (custom series); got length ", length(effect), ".",
      call. = FALSE
    )
  } else {
    base <- 0.5
  }
  if (dist != "geometric") {
    stop("Only dist = \"geometric\" (or an explicit `effect` series) is ",
         "supported.", call. = FALSE)
  }
  base ^ seq_len(n)
}

#' Draw a residual vector to hit a target residual variance proportion
#'
#' Assumes total phenotypic variance is scaled to 1, so the residual variance
#' equals `1 - sum(genetic proportions)`. RNG stays in R (DECISION-006).
#'
#' @param n number of individuals.
#' @param resid_var target residual variance (>= 0).
#' @return numeric vector of length `n`.
#' @keywords internal
#' @noRd
.draw_residual <- function(n, resid_var) {
  if (resid_var <= 0) {
    return(rep(0, n))
  }
  e <- stats::rnorm(n, mean = 0, sd = sqrt(resid_var))
  # Rescale to the exact target variance so realized h2 tracks the requested
  # proportions without finite-sample residual-variance drift.
  s <- stats::sd(e)
  if (s > 0) {
    e <- (e - mean(e)) / s * sqrt(resid_var)
  }
  e
}
