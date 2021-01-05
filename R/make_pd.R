#' obtain a positive definite matrix
#' @keywords internal
#' @param m square matrix
#' @param verbose = TRUE
#' @return positive definite matrix
#' @author Samuel Fernandes
#' Last update: Jan 5, 2021
#'
make_pd <- function(m, verbose = TRUE){
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
  return(m2)
  } else {
    return(m)
  }
}
