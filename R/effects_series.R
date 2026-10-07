#' Within-layer effect series (geometric default or custom)
#'
#' Parity-critical effect generation for the v2 grammar. Deterministic given its
#' inputs; the absolute magnitude is normalized away by the per-layer
#' proportion scaling in `.realize_phenotype()`, so only the *relative* effects
#' across QTNs matter.
#'
#' @param n number of effects to generate (number of QTNs / pairs).
#' @param dist within-layer distribution; only "geometric" is supported (checked
#'   even when an explicit series is given).
#' @param effect optional override. A scalar is treated as the geometric base;
#'   a length-`n` vector is used verbatim as a custom series (v1
#'   `sim_method = "custom"`). A base of 0 gives the all-zero series (used by the
#'   orthogonal model's `a = 0`, pure dominance); a base whose power overflows or
#'   underflows to exactly 0 within `n` terms is an error.
#' @param count name of the count argument, for the length error message
#'   (`"n_qtn"`, or `"n_pairs"` for epistasis).
#' @param arg name of the user's effect argument in messages (`"effect"`, or
#'   `"a"` for the orthogonal model).
#' @return numeric vector of length `n`.
#' @keywords internal
#' @noRd
.effect_series <- function(n, dist = "geometric", effect = NULL,
                           count = "n_qtn", arg = "effect") {
  if (n <= 0) {
    return(numeric(0))
  }
  if (!is.character(dist) || length(dist) != 1L || is.na(dist)) {
    stop("`dist` must be one non-missing character value.", call. = FALSE)
  }
  # Validated first, so an invalid `dist` is rejected even when an explicit
  # `effect` series is supplied (it used to be accepted silently then).
  if (dist != "geometric") {
    stop("Only dist = \"geometric\" is supported (or supply an explicit ",
         "`", arg, "` series); got dist = \"", dist, "\".", call. = FALSE)
  }
  if (!is.null(effect) &&
      (!is.numeric(effect) || any(!is.finite(effect)))) {
    stop("`", arg, "` must contain only finite numeric values.", call. = FALSE)
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
      "`", arg, "` must be either a single value (geometric base) or a vector of ",
      "length ", count, " (custom series); got length ", length(effect), ".",
      call. = FALSE
    )
  } else {
    base <- 0.5
  }
  series <- base ^ seq_len(n)
  if (any(!is.finite(series))) {
    stop("The geometric effect series overflows: base ", base, " to the power ",
         "of ", count, " = ", n, " is not finite from position ",
         which(!is.finite(series))[1L], ". Use a base closer to 1 or fewer ",
         "effects (or supply an explicit `", arg, "` series).", call. = FALSE)
  }
  if (base != 0 && any(series == 0)) {      # base 0 is the documented all-zero series
    stop("The geometric effect series underflows to exactly 0 from position ",
         which(series == 0)[1L], " (base ", base, ", ", count, " = ", n, "), so ",
         "those effects could never carry variance. Use a base closer to 1 or ",
         "fewer effects (or supply an explicit `", arg, "` series).",
         call. = FALSE)
  }
  series
}

#' Draw a residual vector to hit a target residual variance proportion
#'
#' Assumes total phenotypic variance is scaled to 1, so the residual variance
#' equals `1 - sum(genetic proportions)`. RNG stays in R.
#'
#' @param n number of individuals (>= 2).
#' @param resid_var target residual variance (>= 0).
#' @param residual_mode `fixed` standardizes the sample; `random` keeps the normal draw.
#' @return numeric vector of length `n`.
#' @keywords internal
#' @noRd
.draw_residual <- function(n, resid_var, residual_mode = "fixed") {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2) {
    stop("A residual needs n >= 2 individuals so its variance is defined; got ",
         "n = ", paste(n, collapse = ", "), ".", call. = FALSE)
  }
  if (resid_var <= 0) {
    return(rep(0, n))
  }
  e <- stats::rnorm(n, mean = 0, sd = sqrt(resid_var))
  if (identical(residual_mode, "random")) return(e)
  # Rescale to the exact target variance so realized h2 tracks the requested
  # proportions without finite-sample residual-variance drift.
  s <- stats::sd(e)
  if (s > 0) {
    e <- (e - mean(e)) / s * sqrt(resid_var)
  }
  e
}
