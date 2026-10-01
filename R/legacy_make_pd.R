#' Eigenvalue clamp used by the frozen v1 engine (NOT a nearest-correlation-matrix fix)
#'
#' Internal helper of the frozen legacy engine (`create_phenotypes()`), applied
#' to the user's genetic correlation matrix `cor` and to the sample covariance
#' of the standardized genetic values before they are Cholesky-factorized.
#'
#' The v1 rule is reproduced exactly (and must not be "improved": the RDS
#' references and every published v1 result depend on it):
#' \enumerate{
#'   \item If every eigenvalue of `m` is strictly positive, `m` is returned
#'     unchanged.
#'   \item Otherwise, with `tol = n * max(|lambda|) * .Machine$double.eps` and
#'     `tau_k = max(0, 2 * tol - lambda_k)`, the matrix
#'     `m + V diag(tau) V'` is formed and rounded to 2 decimals.
#' }
#' Consequences that callers and users must know:
#' \itemize{
#'   \item The result is a positive-definite *covariance-like* matrix, not a
#'     correlation matrix. Adding `V diag(tau) V'` raises the diagonal (for
#'     `[[1, .9, -.9], [.9, 1, .9], [-.9, .9, 1]]` the diagonal becomes 1.27
#'     and the off-diagonals \eqn{\pm}0.63), so the correlation that is
#'     realized afterwards is `m2[i, j] / sqrt(m2[i, i] * m2[j, j])`
#'     (0.496 instead of 0.9 in that example) and the genetic variance is
#'     inflated by the diagonal (27 percent there). Callers must therefore not
#'     assume that a non-positive-definite `cor` is realized as requested;
#'     `create_phenotypes()` rejects a non-positive-definite `cor` up front.
#'   \item `round(, 2)` can push a repaired matrix back to (numerically)
#'     singular or indefinite. Previously this surfaced as a cryptic
#'     `chol()` failure ("the leading minor of order k is not positive").
#'     It is now detected here and reported with an informative error.
#'   \item A perfectly collinear `m` (for example two traits with identical
#'     effect series, or a single QTN shared by all traits, so that
#'     `cov(scale(genetic values))` is `[[1, 1], [1, 1]]`) cannot be repaired
#'     by the clamp after rounding and is rejected with an informative error.
#' }
#'
#' @keywords internal
#' @param m square symmetric matrix
#' @param verbose = TRUE
#' @param what short description of the matrix, used in error messages
#' @return the input matrix when it is already positive definite, otherwise the
#'   clamped and rounded matrix described above. It is guaranteed to be
#'   Cholesky-factorizable (an informative error is raised otherwise).
#' @author Samuel Fernandes. Last update: Jan 5, 2021
#'
make_pd <- function(m, verbose = TRUE, what = "genetic correlation matrix"){
  if (!is.matrix(m)) m <- as.matrix(m)
  if (!is.numeric(m) || nrow(m) != ncol(m)) {
    stop("make_pd(): the ", what, " must be a square numeric matrix.",
         call. = FALSE)
  }
  if (anyNA(m) || any(!is.finite(m))) {
    stop("make_pd(): the ", what, " contains NA/NaN/Inf values.",
         call. = FALSE)
  }
  if (!isSymmetric(unname(m))) {
    stop("make_pd(): the ", what, " is not symmetric. Supply a symmetric ",
         "matrix (the lower and the upper triangle must agree).",
         call. = FALSE)
  }
  e <- eigen(m)
  if(any(e$values <= 0)){
    if (verbose) cat(
      "Modifing the genetic correlation matrix to make it positive definite! \n"
    )
  n <- nrow(m)
  tol <- nrow(m) * max(abs(e$values)) * .Machine$double.eps
  delta <- 2 * tol
  tau <- pmax(0, delta - e$values)
  dm <- e$vectors %*% diag(tau, n) %*% t(e$vectors)
  m2 <- m + dm
  m2 <- round(m2, 2)
  } else {
    m2 <- m
  }
  # Post-condition (frozen v1 arithmetic above is untouched): the caller takes
  # t(chol(.)); fail here, with a message, instead of inside chol().
  ok <- tryCatch({
    chol(m2)
    TRUE
  }, error = function(err) FALSE)
  if (!ok) {
    mn <- suppressWarnings(min(eigen(m2, symmetric = TRUE,
                                     only.values = TRUE)$values))
    stop(
      "The ", what, " is not positive definite (smallest eigenvalue after the ",
      "v1 eigenvalue clamp and rounding to 2 decimals: ",
      format(mn, digits = 3), "), so it cannot be factorized. ",
      "Supply a strictly positive-definite `cor` matrix",
      if (identical(what, "genetic correlation matrix")) {
        "."
      } else {
        " and make sure the traits do not have identical (perfectly collinear) genetic values, e.g. identical effect series or a single shared QTN."
      },
      call. = FALSE
    )
  }
  m2
}
