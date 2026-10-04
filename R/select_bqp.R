# Relatedness-penalized selection of exactly N individuals by binary quadratic
# programming (BQP), after Montesinos-Lopez et al. (2025), Plant Methods 22:7.
# The paper maximizes (Eq. 1; unweighted there, `weights` is an extension)
#     Z = sum_j sum_i w_j s_ij x_i - k sum_i sum_k G_ik x_i x_k,   x_i in {0, 1},
# subject to sum_i x_i = s (Eq. 2) and, per trait, sum_i s_ij x_i >= l_j with
# l_j = R_j * s / 100 (Eq. 3-4, R_j in percent; verified against the published
# equations, PMC12849579). The paper
# solves it with CVXR; this implementation is dependency-free and deterministic:
# exact enumeration when choose(n, N) is small, otherwise greedy construction
# followed by best-improvement 1-swap local search (a heuristic, not a proof of
# optimality). No random numbers are drawn.

#' Citation notice for the BQP selection method, once per session
#' @keywords internal
#' @noRd
.cite_bqp <- function() {
  rlang::inform(
    .cite_main("Montesinos-Lopez, O.A., Montesinos-Lopez, A., Hernandez-Suarez, ",
               "C.M. and Alemu, A. (2025). A selection index with minimal ",
               "genetic relatedness for multi-trait data via binary quadratic ",
               "programming. Plant Methods 22:7, doi:10.1186/s13007-025-01484-4, ",
               "when using select_ind(method = \"bqp\")."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_bqp_citation"
  )
}

#' Merit matrix, merit score and relationship matrix for method = "bqp"
#'
#' Trait columns are standardized (mean 0, sd 1) as in the paper's standardized
#' values s_ij; the merit is sum_j w_j s_ij. With `weights = NULL` the single
#' criterion `on`/`trait` is used (weight 1). The relationship matrix is the
#' VanRaden G of the sim's individuals in `sim$ids` order.
#' @keywords internal
#' @noRd
.bqp_inputs <- function(sim, on, trait, weights, rep) {
  if (!inherits(sim$geno, "Population")) {
    stop("method = \"bqp\" needs a Population-backed simulation (it penalizes ",
         "genomic relatedness, g_matrix()).", call. = FALSE)
  }
  nt <- sim$n_traits
  if (is.null(weights)) {
    S <- matrix(as.numeric(.criterion_values(sim, on, trait, rep)), ncol = 1L)
    w <- 1
  } else {
    if (!is.numeric(weights) || length(weights) != nt || any(!is.finite(weights))) {
      stop("method = \"bqp\": `weights` must be finite numeric, one economic ",
           "weight per trait (", nt, "); omit it to select on the single ",
           "criterion `on`/`trait`.", call. = FALSE)
    }
    if (!is.character(on)) {
      stop("method = \"bqp\" with `weights` scores every trait on the named ",
           "criterion `on` (\"pheno\", \"gv\" or \"bv\"); a numeric/function ",
           "`on` is a single score, so omit `weights`.", call. = FALSE)
    }
    S <- do.call(cbind, lapply(seq_len(nt), function(k) {
      as.numeric(.criterion_values(sim, on, k, rep))
    }))
    w <- as.numeric(weights)
  }
  if (anyNA(S) || any(!is.finite(S))) {
    stop("method = \"bqp\": the criterion has non-finite values.", call. = FALSE)
  }
  sds <- apply(S, 2L, stats::sd)
  if (any(!is.finite(sds)) || any(sds <= 0)) {
    stop("method = \"bqp\": a trait criterion has no variation, so it cannot ",
         "be standardized.", call. = FALSE)
  }
  S <- sweep(sweep(S, 2L, colMeans(S), "-"), 2L, sds, "/")
  G <- .bqp_g(sim)
  list(S = S, w = w, G = G, score = stats::setNames(as.numeric(S %*% w), sim$ids))
}

#' Genomic relationship matrix of the BQP paper: G = W W' / p
#'
#' Montesinos-Lopez et al. (2025) compute G from the scaled marker matrix `W`
#' (individuals x markers, each marker column centered and divided by its
#' standard deviation) as `W W' / p`, `p` the number of markers (VanRaden 2008,
#' standardized form). This is not [g_matrix()]'s VanRaden method 1, and the two
#' can rank sets differently (unequal allele frequencies weigh markers
#' differently), so the paper's form is used here. Monomorphic markers have no
#' scale and are dropped (they carry no relationship information).
#' @keywords internal
#' @noRd
.bqp_g <- function(sim) {
  W <- .geno_cols(sim, seq_len(sim$n_markers))     # individuals x markers
  storage.mode(W) <- "double"
  sdv <- apply(W, 2L, stats::sd)
  keep <- which(is.finite(sdv) & sdv > 0)
  if (!length(keep)) {
    stop("method = \"bqp\": every marker is monomorphic among these ",
         "individuals, so the relationship matrix is undefined.", call. = FALSE)
  }
  W <- W[, keep, drop = FALSE]
  W <- sweep(sweep(W, 2L, colMeans(W), "-"), 2L, sdv[keep], "/")
  G <- tcrossprod(W) / ncol(W)
  dimnames(G) <- list(sim$ids, sim$ids)
  G
}

#' Validate `lambda` and `min_gain` for method = "bqp"
#' @keywords internal
#' @noRd
.bqp_check_args <- function(lambda, min_gain, n_traits_used) {
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda < 0) {
    stop("`lambda` must be one finite number >= 0 (the relatedness penalty ",
         "weight; 0 = rank on merit alone).", call. = FALSE)
  }
  if (!is.null(min_gain)) {
    if (!is.numeric(min_gain) || any(!is.finite(min_gain)) ||
        !(length(min_gain) %in% c(1L, n_traits_used)) ||
        any(min_gain < 0) || any(min_gain > 100)) {
      stop("`min_gain` must be NULL, or numbers in [0, 100]: one value or one ",
           "per trait (", n_traits_used, "), the minimum desired gain R_j in ",
           "percent (Montesinos-Lopez et al. 2025, Eq. 4).", call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Solve the BQP exactly or by greedy + 1-swap local search
#'
#' maximize c'x - lambda x'Gx over x in {0,1}^n with sum(x) = N and, when `S` and
#' `rhs` are given, colSums(S[sel, ]) >= rhs. Deterministic; ties go to the lowest
#' index (exact: first combination in lexicographic order).
#' @return list(sel, objective, method, feasible)
#' @keywords internal
#' @noRd
.bqp_solve <- function(cvec, G, N, lambda, S = NULL, rhs = NULL,
                       max_enum = 2e5) {
  n <- length(cvec)
  tol <- 1e-9
  obj_of <- function(sel) {
    sum(cvec[sel]) - lambda * sum(G[sel, sel, drop = FALSE])
  }
  feas_of <- function(sel) {
    is.null(rhs) || all(colSums(S[sel, , drop = FALSE]) >= rhs - tol)
  }
  if (N == n) {
    sel <- seq_len(n)
    return(list(sel = sel, objective = obj_of(sel), method = "exact",
                feasible = feas_of(sel)))
  }
  if (choose(n, N) <= max_enum) {
    cm <- utils::combn(n, N)                       # N x K, lexicographic
    K <- ncol(cm)
    obj <- colSums(matrix(cvec[cm], nrow = N))
    quad <- numeric(K)
    for (a in seq_len(N)) for (b in seq_len(N)) {
      quad <- quad + G[cbind(cm[a, ], cm[b, ])]
    }
    obj <- obj - lambda * quad
    ok <- rep(TRUE, K)
    if (!is.null(rhs)) {
      for (j in seq_len(ncol(S))) {
        ok <- ok & colSums(matrix(S[, j][cm], nrow = N)) >= rhs[j] - tol
      }
    }
    if (!any(ok)) {
      return(list(sel = NULL, objective = NA_real_, method = "exact",
                  feasible = FALSE))
    }
    obj[!ok] <- -Inf
    best <- max(obj)
    k <- which(obj >= best - 1e-12 * (1 + abs(best)))[1L]
    return(list(sel = cm[, k], objective = obj[k], method = "exact",
                feasible = TRUE))
  }
  # ---- greedy construction (objective only), then penalized 1-swap search ----
  dg <- diag(G)
  in_set <- rep(FALSE, n)
  q <- numeric(n)                                  # q_k = sum_{s in set} G[k, s]
  for (step in seq_len(N)) {
    gain <- cvec - lambda * (2 * q + dg)
    gain[in_set] <- -Inf
    i <- which.max(gain)
    in_set[i] <- TRUE
    q <- q + G[, i]
  }
  Mpen <- 1e6 * (1 + max(abs(cvec)) + lambda * max(abs(G)) * N)
  pen_of <- function(Tj) {
    if (is.null(rhs)) 0 else Mpen * sum(pmax(0, rhs - Tj))
  }
  Tj <- if (is.null(rhs)) NULL else colSums(S[in_set, , drop = FALSE])
  iter <- 0L
  repeat {
    iter <- iter + 1L
    ins <- which(in_set); outs <- which(!in_set)
    # change in objective for swapping o (in) -> i (out), as an |ins| x |outs| matrix
    d_lin <- outer(-cvec[ins], cvec[outs], "+")
    d_quad <- outer(-2 * q[ins] + dg[ins], 2 * q[outs] + dg[outs], "+") -
      2 * G[ins, outs, drop = FALSE]
    delta <- d_lin - lambda * d_quad
    if (!is.null(rhs)) {
      new_pen <- 0
      for (j in seq_len(ncol(S))) {
        Tn <- Tj[j] + outer(-S[ins, j], S[outs, j], "+")
        new_pen <- new_pen + pmax(0, rhs[j] - Tn)
      }
      delta <- delta - (Mpen * new_pen - pen_of(Tj))
    }
    m <- which.max(delta)
    if (!is.finite(delta[m]) || delta[m] <= tol || iter > 1e5) break
    r <- (m - 1L) %% length(ins) + 1L
    cc <- (m - 1L) %/% length(ins) + 1L
    o <- ins[r]; i <- outs[cc]
    in_set[o] <- FALSE; in_set[i] <- TRUE
    q <- q - G[, o] + G[, i]
    if (!is.null(rhs)) Tj <- Tj - S[o, ] + S[i, ]
  }
  sel <- which(in_set)
  list(sel = sel, objective = obj_of(sel), method = "local_search",
       feasible = feas_of(sel))
}

#' Run the BQP for `select_ind(method = "bqp")`
#' @keywords internal
#' @noRd
.sel_bqp <- function(bq, direction, keep_n, lambda, min_gain) {
  sgn <- if (direction == "low") -1 else 1
  S <- sgn * bq$S
  cvec <- as.numeric(S %*% bq$w)
  rhs <- NULL
  if (!is.null(min_gain)) {
    rhs <- keep_n * rep_len(as.numeric(min_gain), ncol(S)) / 100
  }
  sol <- .bqp_solve(cvec, bq$G, keep_n, lambda, S = S, rhs = rhs)
  if (!isTRUE(sol$feasible)) {
    if (identical(sol$method, "exact")) {
      stop("method = \"bqp\": no set of ", keep_n, " individuals meets the ",
           "`min_gain` constraints; lower `min_gain` or raise `n`.", call. = FALSE)
    }
    # the 1-swap search can miss a feasible set that needs several swaps
    stop("method = \"bqp\": the greedy + 1-swap search found no set of ", keep_n,
         " individuals meeting the `min_gain` constraints. This is a heuristic ",
         "(the problem is too large to enumerate): a feasible set may still ",
         "exist. Lower `min_gain`, or select among fewer candidates.",
         call. = FALSE)
  }
  sol
}
