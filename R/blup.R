# Known-variance BLUP (DECISION-030): predicted breeding values from observable
# phenotypes and a relationship matrix, the pedigree relationship matrix, accuracy
# reporting and a manifest of the engine's selection operators.

#' Numerator relationship matrix from the recorded pedigree
#'
#' The additive relationship matrix \eqn{A} of a `Population`'s individuals,
#' built from its recorded pedigree (see [parentage()]) by the tabular method
#' (Emik & Terrill 1949) on coefficients of relationship (Wright 1922):
#' processing individuals parents-first, \eqn{A_{ij} = (A_{i,m(j)} + A_{i,f(j)}) / 2}
#' and \eqn{A_{jj} = 1 + F_j} with \eqn{F_j = A_{m(j) f(j)} / 2}. A self has both
#' parents equal, so \eqn{F = (1 + F_P) / 2}; a doubled haploid is fully inbred,
#' \eqn{A_{jj} = 2}, and relates to others as its parent does. Founders are the
#' base: non-inbred and unrelated. Ancestors not in `pop` still enter through the
#' pedigree.
#'
#' @param pop a `Population`.
#' @param ids optional character ids (of `pop`'s individuals) to return;
#'   default all.
#' @return A symmetric matrix with dimnames `ids`.
#' @references
#' Wright S (1922) Coefficients of inbreeding and relationship. \emph{The American
#'   Naturalist} 56:330--338. \doi{10.1086/279872}
#'
#' Emik LO, Terrill CE (1949) Systematic procedures for calculating inbreeding
#'   coefficients. \emph{Journal of Heredity} 40:51--55.
#'   \doi{10.1093/oxfordjournals.jhered.a105986}
#' @seealso [parentage()], [predict_ebv()], [g_matrix()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:2)
#' f1 <- cross(pop[1], pop[2], n = 2, seed = 1)
#' s1 <- selfcross(f1[1], n = 2, seed = 2)
#' a_matrix(c(f1, s1))
a_matrix <- function(pop, ids = NULL) {
  .check_population(pop)
  pop <- .ensure_pedigree(pop)
  ped <- pop$pedigree
  ped <- ped[order(ped$generation), , drop = FALSE]       # parents first
  n <- nrow(ped)
  A <- matrix(0, n, n)
  mi <- match(ped$mother, ped$key)
  fi <- match(ped$father, ped$key)
  for (k in seq_len(n)) {
    m <- mi[k]; f <- fi[k]
    if (k > 1L) {
      prev <- seq_len(k - 1L)
      row <- (if (is.na(m)) 0 else A[m, prev]) / 2 +
        (if (is.na(f)) 0 else A[f, prev]) / 2
      A[k, prev] <- row
      A[prev, k] <- row
    }
    A[k, k] <- if (identical(ped$design[k], "dh")) {
      2
    } else if (is.na(m) || is.na(f)) {
      1
    } else {
      1 + A[m, f] / 2
    }
  }
  at <- match(pop$keys, ped$key)
  out <- A[at, at, drop = FALSE]
  dimnames(out) <- list(pop$ids, pop$ids)
  if (!is.null(ids)) {
    # ids are identifiers, never positions
    if (!is.character(ids) || anyNA(ids)) {
      stop("a_matrix(): `ids` must be a character vector of individual ids.",
           call. = FALSE)
    }
    if (!all(ids %in% pop$ids)) {
      stop("a_matrix(): `ids` not in the population: ",
           paste(utils::head(setdiff(ids, pop$ids), 5), collapse = ", "), ".",
           call. = FALSE)
    }
    out <- out[ids, ids, drop = FALSE]
  }
  out
}

#' Predicted breeding values by BLUP with known variance components
#'
#' Best linear unbiased prediction (Henderson 1975) of every individual's
#' breeding value from the observable phenotypes of a phenotyped subset and a
#' relationship matrix: genomic (`"gblup"`, \eqn{K = G} from [g_matrix()],
#' VanRaden 2008) or pedigree (`"pedigree"`, \eqn{K = A} from [a_matrix()]).
#' Individuals without a record -- selection candidates -- are predicted through
#' their relationships to the phenotyped ones.
#'
#' The model is \eqn{y = 1\mu + Zu + e}, \eqn{u \sim N(0, K\sigma^2_A)},
#' \eqn{e \sim N(0, I\sigma^2_e)}, with the variance components **known** -- a
#' simulation knows them -- so this is a deterministic linear solve, with no
#' variance-component estimation. It is computed in the equivalent generalized
#' least-squares form \eqn{\hat u = K Z' (Z K Z' + \lambda I)^{-1}(y - 1\hat\mu)},
#' \eqn{\lambda = \sigma^2_e / \sigma^2_A}, which needs no inverse of `K` (so a
#' singular genomic matrix is fine). Reliability is
#' \eqn{1 - PEV_i / (K_{ii}\sigma^2_A)}.
#'
#' The result is an **estimate** from phenotypes: pass observable records (e.g.
#' [phenotype_value()] or a simulation's phenotypes), never true genetic or
#' breeding values, and compare it with the truth using [prediction_accuracy()].
#' Across selection cycles keep the variance components and, for `"gblup"`, the
#' base allele frequencies (`base_freq`) fixed at the base population's; per-cycle
#' values inflate the apparent accuracy. Records of every individual that
#' selection was based on should be included, or predictions are biased
#' (Henderson 1975).
#'
#' @param x a `Population` holding every individual to predict (phenotyped and
#'   not).
#' @param pheno a numeric vector of phenotypes named by id (a subset of `x`).
#' @param method `"gblup"` or `"pedigree"`.
#' @param h2,var_a,var_e the variance components: either `h2` (then
#'   \eqn{\sigma^2_A = h^2 \sigma^2_P}, \eqn{\sigma^2_P} the variance of `ref`, default
#'   `pheno`) or both `var_a` (base additive variance) and `var_e`.
#' @param ref optional numeric vector of base-population phenotypes setting
#'   \eqn{\sigma^2_P} for `h2`.
#' @param K optional precomputed relationship matrix over `x`'s individuals
#'   overriding the one `method` builds: finite, symmetric and positive
#'   semidefinite (it is a covariance matrix up to \eqn{\sigma^2_A}; only
#'   floating-point rounding below zero is tolerated -- \eqn{n \epsilon} times the
#'   largest eigenvalue on the correlation scale -- and
#'   \eqn{K + (\sigma^2_e / \sigma^2_A) I} over the phenotyped individuals must be
#'   positive definite), with row and
#'   column names both the ids of `x`, in the same order (it is reordered to
#'   `x`'s order).
#' @param base_freq,ridge passed to [g_matrix()] for `"gblup"`.
#' @return A numeric vector of predicted breeding values named by `x`'s ids, with
#'   attributes `reliability`, `mu` (the estimated mean), `lambda`, `var_a`,
#'   `var_e`, `method` and, for `"gblup"` with the genomic matrix built here (no
#'   `K` supplied) and `ridge = 0`, `marker_effects` (back-solved, on the
#'   \eqn{M - 2p} gene-content scale of [g_matrix()]; a supplied `K` carries no
#'   marker scale, so none are returned).
#' @references
#' Henderson CR (1975) Best linear unbiased estimation and prediction under a
#'   selection model. \emph{Biometrics} 31:423--447. \doi{10.2307/2529430}
#'
#' VanRaden PM (2008) Efficient methods to compute genomic predictions.
#'   \emph{Journal of Dairy Science} 91:4414--4423. \doi{10.3168/jds.2007-0980}
#' @seealso [prediction_accuracy()], [a_matrix()], [g_matrix()], [select_ind()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:60)
#' q <- c("ss196442916", "ss196439337", "ss196480535")
#' y <- phenotype_value(pop, q, c(1, 0.5, 0.25), h2 = 0.5, seed = 1)
#' ebv <- predict_ebv(pop, y[1:40], h2 = 0.5)      # 20 unphenotyped candidates
#' head(ebv)
predict_ebv <- function(x, pheno, method = c("gblup", "pedigree"), h2 = NULL,
                        var_a = NULL, var_e = NULL, ref = NULL, K = NULL,
                        base_freq = NULL, ridge = 0) {
  .check_population(x)
  method <- match.arg(method)
  ids <- x$ids
  if (!is.numeric(pheno) || is.null(names(pheno)) || !length(pheno) ||
      any(!is.finite(pheno))) {
    stop("predict_ebv(): `pheno` must be a finite numeric vector named by id.",
         call. = FALSE)
  }
  if (anyDuplicated(names(pheno)) || !all(names(pheno) %in% ids)) {
    stop("predict_ebv(): `pheno` names must be distinct ids of `x`.",
         call. = FALSE)
  }
  vc <- .blup_variances(pheno, h2, var_a, var_e, ref)
  built <- is.null(K)
  # a supplied K is validated and used in its symmetric form (rounding-level
  # asymmetry is tolerated by the check, so it must not survive into the solve)
  if (!built) K <- .check_relationship(K)
  if (built) {
    K <- if (method == "gblup") {
      g_matrix(x, ridge = ridge, base_freq = base_freq)
    } else {
      a_matrix(x)
    }
  }
  if (!is.matrix(K) || !identical(dim(K), c(length(ids), length(ids)))) {
    stop("predict_ebv(): `K` must be a ", length(ids), " x ", length(ids),
         " relationship matrix over `x`'s individuals.", call. = FALSE)
  }
  # a supplied K must say which individual each row and column is, identically
  # on both axes, or it could be silently applied to the wrong individuals
  if (!built) {
    rn <- rownames(K)
    if (is.null(rn) || !identical(rn, colnames(K)) || anyDuplicated(rn) ||
        !setequal(rn, ids)) {
      stop("predict_ebv(): a supplied `K` needs row and column names that are ",
           "both `x`'s ids, in the same order on both axes.", call. = FALSE)
    }
  }
  if (!is.null(rownames(K))) K <- K[ids, ids, drop = FALSE]
  rec <- match(names(pheno), ids)
  y <- as.numeric(pheno)
  lambda <- vc$var_e / vc$var_a
  # both components are positive, so the ratio must be finite and positive;
  # overflow (Inf) or underflow (0) means it is not representable
  if (!is.finite(lambda) || lambda <= 0) {
    stop("predict_ebv(): the variance ratio var_e / var_a is not representable ",
         "(", signif(vc$var_e, 3), " / ", signif(vc$var_a, 3), "); rescale the ",
         "variance components.", call. = FALSE)
  }
  # V = K_rr + lambda I must be positive definite (it is the phenotypic
  # covariance over sigma^2_A); its Cholesky factor both checks that -- a K
  # indefinite beyond what lambda absorbs fails here -- and inverts it
  R <- tryCatch(chol(K[rec, rec, drop = FALSE] + diag(lambda, length(rec))),
                error = function(e) NULL)
  if (is.null(R)) {
    stop("predict_ebv(): K + (var_e / var_a) I over the phenotyped individuals ",
         "is not positive definite; `K` is not a valid relationship ",
         "(covariance) matrix at this precision.", call. = FALSE)
  }
  W <- chol2inv(R)
  one <- rep(1, length(rec))
  w1 <- W %*% one
  mu <- as.numeric(crossprod(w1, y) / crossprod(one, w1))
  alpha <- W %*% (y - mu)
  ebv <- as.numeric(K[, rec, drop = FALSE] %*% alpha)
  if (any(!is.finite(ebv)) || !is.finite(mu)) {
    stop("predict_ebv(): the BLUP solution is not finite (numerical ",
         "overflow); rescale the phenotypes or variance components.",
         call. = FALSE)
  }
  # PEV / var_a = K - K_{.,r} P K_{r,.} with P = W - W 1 (1'W1)^{-1} 1'W
  Pm <- W - (w1 %*% t(w1)) / as.numeric(crossprod(one, w1))
  KP <- K[, rec, drop = FALSE] %*% Pm
  explained <- rowSums(KP * K[, rec, drop = FALSE])
  rel <- ifelse(diag(K) > 0, explained / diag(K), NA_real_)
  # a reliability is a squared correlation: outside [0, 1] beyond rounding the
  # system is numerically unstable (K not semidefinite at the precision the
  # variance ratio needs); within rounding, clamp
  if (any(rel < -1e-8 | rel > 1 + 1e-8, na.rm = TRUE)) {
    stop("predict_ebv(): reliabilities fall outside [0, 1]; the system is ",
         "numerically unstable (`K` is not positive semidefinite at the ",
         "precision var_e / var_a requires).", call. = FALSE)
  }
  rel <- pmin(pmax(rel, 0), 1)
  names(ebv) <- names(rel) <- ids
  attr(ebv, "reliability") <- rel
  attr(ebv, "mu") <- mu
  attr(ebv, "lambda") <- lambda
  attr(ebv, "var_a") <- vc$var_a
  attr(ebv, "var_e") <- vc$var_e
  attr(ebv, "method") <- method
  if (method == "gblup" && ridge == 0 && built) {
    attr(ebv, "marker_effects") <- .gblup_marker_effects(x, rec, alpha,
                                                         base_freq)
  }
  ebv
}

#' Validate a user-supplied relationship matrix: numeric, finite, symmetric and
#' positive semidefinite (up to rounding); returns its symmetric form
#' @keywords internal
#' @noRd
.check_relationship <- function(K) {
  if (!is.matrix(K) || !is.numeric(K) || any(!is.finite(K))) {
    stop("predict_ebv(): `K` must be a finite numeric matrix.", call. = FALSE)
  }
  if (nrow(K) != ncol(K)) {
    stop("predict_ebv(): `K` must be square.", call. = FALSE)
  }
  # element-wise tolerances, so a large entry elsewhere cannot hide a local
  # defect
  if (any(abs(K - t(K)) > 1e-8 * pmax(1, abs(K), abs(t(K))))) {
    stop("predict_ebv(): `K` must be symmetric.", call. = FALSE)
  }
  K <- (K + t(K)) / 2
  dg <- diag(K)
  if (any(dg < 0)) {
    stop("predict_ebv(): `K` has a negative diagonal (a negative variance); a ",
         "relationship matrix is a covariance matrix.", call. = FALSE)
  }
  # PSD on the correlation scale (unit diagonal), so the test does not depend
  # on how large the other variances are; a zero-variance row must be all zero
  z <- dg == 0
  off <- K
  diag(off) <- 0
  bad_zero <- any(z) && any(off[z, , drop = FALSE] != 0)
  # divide by the root variances one axis at a time (never form 1/sqrt(v_i v_j),
  # which overflows for tiny variances); for a PSD K every entry is then in
  # [-1, 1], so a non-finite one means K is not PSD
  rs <- sqrt(dg[!z])
  C <- K[!z, !z, drop = FALSE] / rs
  C <- t(t(C) / rs)
  ev <- if (!any(!z)) {
    0
  } else if (any(!is.finite(C))) {
    -Inf
  } else {
    eigen(C, symmetric = TRUE, only.values = TRUE)$values
  }
  # only floating-point rounding below zero (n eps max|eigenvalue|, the usual
  # numerical-rank tolerance) is accepted as semidefinite
  tol <- length(ev) * .Machine$double.eps * max(1, abs(ev))
  if (bad_zero || min(ev) < -tol) {
    stop("predict_ebv(): `K` must be positive semidefinite; a relationship ",
         "matrix is a covariance matrix.", call. = FALSE)
  }
  K
}

#' Variance components for BLUP: h2 (with the phenotypic variance of ref/pheno)
#' or var_a + var_e
#' @keywords internal
#' @noRd
.blup_variances <- function(pheno, h2, var_a, var_e, ref) {
  if (!is.null(h2)) {
    if (!is.null(var_a) || !is.null(var_e)) {
      stop("predict_ebv(): give either `h2` or both `var_a` and `var_e`, not ",
           "both.", call. = FALSE)
    }
    if (!is.numeric(h2) || length(h2) != 1L || !is.finite(h2) || h2 <= 0 ||
        h2 >= 1) {
      stop("predict_ebv(): `h2` must be one value in (0, 1).", call. = FALSE)
    }
    vp <- stats::var(if (is.null(ref)) pheno else ref)
    if (!is.finite(vp) || vp <= 0) {
      stop("predict_ebv(): the phenotypic variance is zero or undefined; give ",
           "`var_a` and `var_e` instead.", call. = FALSE)
    }
    return(list(var_a = h2 * vp, var_e = (1 - h2) * vp))
  }
  ok <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v > 0
  if (!ok(var_a) || !ok(var_e)) {
    stop("predict_ebv(): give `h2`, or positive `var_a` and `var_e`.",
         call. = FALSE)
  }
  list(var_a = var_a, var_e = var_e)
}

#' GBLUP marker effects, u_m = M_rec' alpha / (2 sum p(1 - p)), on g_matrix()'s
#' centred gene-content scale (so Z u_m reproduces the GEBVs).
#' @keywords internal
#' @noRd
.gblup_marker_effects <- function(x, rec, alpha, base_freq) {
  M <- t(dosages(x)) + 1
  p <- if (is.null(base_freq)) colMeans(M) / 2 else base_freq
  poly <- which(p > 0 & p < 1)
  Zr <- sweep(M[rec, poly, drop = FALSE], 2, 2 * p[poly], "-")
  u <- numeric(ncol(M))
  u[poly] <- as.numeric(crossprod(Zr, alpha)) / (2 * sum(p[poly] * (1 - p[poly])))
  stats::setNames(u, x$map$snp)
}

#' Accuracy of predicted breeding values
#'
#' The realized accuracy \eqn{cor(\hat u, u)} of predictions against the true
#' breeding values a simulation knows, with the regression slope of the truth on
#' the prediction (1 for an unbiased predictor, BLUP's \eqn{E[u \mid \hat u] = \hat u}
#' property; below 1 means over-dispersed predictions).
#'
#' @param ebv predicted values (named by id).
#' @param truth true breeding values (named by id; e.g. from a simulation's
#'   `on = "bv"` criterion). Matched by name when both are named.
#' @return A named numeric vector: `accuracy`, `slope`, `n`.
#' @seealso [predict_ebv()]
#' @export
#' @examples
#' prediction_accuracy(c(a = 1, b = 2, c = 3), c(a = 1.2, b = 1.9, c = 3.3))
prediction_accuracy <- function(ebv, truth) {
  if (!is.numeric(ebv) || !is.numeric(truth)) {
    stop("prediction_accuracy(): `ebv` and `truth` must be numeric.",
         call. = FALSE)
  }
  if (!is.null(names(ebv)) && !is.null(names(truth))) {
    common <- intersect(names(ebv), names(truth))
    ebv <- ebv[common]; truth <- truth[common]
  } else if (length(ebv) != length(truth)) {
    stop("prediction_accuracy(): unnamed `ebv` and `truth` must have equal ",
         "length.", call. = FALSE)
  }
  ok <- is.finite(ebv) & is.finite(truth)
  if (sum(ok) < 3L) {
    stop("prediction_accuracy(): need at least three matched, finite pairs.",
         call. = FALSE)
  }
  e <- as.numeric(ebv[ok]); u <- as.numeric(truth[ok])
  c(accuracy = stats::cor(e, u),
    slope = if (stats::var(e) > 0) stats::cov(u, e) / stats::var(e) else NA_real_,
    n = sum(ok))
}

#' The engine's selection operators, as a manifest
#'
#' A table of the selection and scoring operators the engine exports, with their
#' parameters and sources, for front ends (e.g. breedingDesigner) that build
#' their interface from it.
#'
#' @return A data frame with columns `id`, `label`, `fn` (the function and method
#'   to call), `params` (main arguments of that function, comma-separated; `|`
#'   separates alternatives, `+` arguments given together) and `source`.
#' @seealso [select_ind()], [predict_ebv()], [combining_ability()]
#' @export
#' @examples
#' selection_methods()[, c("id", "fn")]
selection_methods <- function() {
  m <- function(id, label, fn, params, source) {
    data.frame(id = id, label = label, fn = fn, params = params,
               source = source, stringsAsFactors = FALSE)
  }
  rbind(
    m("mass", "Mass (truncation) selection", 'select_ind(method = "mass")',
      "n | prop | intensity, on, trait, direction",
      "Falconer & Mackay 1996"),
    m("within_family", "Within-family selection",
      'select_ind(method = "within_family")', "family, n | prop",
      "Falconer & Mackay 1996"),
    m("among_family", "Among-family (family) selection",
      'select_ind(method = "among_family")', "family, n | prop",
      "Falconer & Mackay 1996"),
    m("combined", "Combined (own + family) selection",
      'select_ind(method = "combined")', "family, h2, family_relationship",
      "Lush 1947"),
    m("index", "Smith-Hazel selection index", 'select_ind(method = "index")',
      "weights", "Smith 1936; Hazel 1943"),
    m("quadratic_index", "Quadratic (nonlinear) selection index",
      'select_ind(method = "quadratic_index")', "weights, quad_weights",
      "Ceron-Rojas et al. 2026"),
    m("culling", "Independent culling levels", 'select_ind(method = "culling")',
      "culling, trait, sequential, direction", "Hazel & Lush 1942"),
    m("tandem", "Tandem selection (one trait per generation, recycled)",
      "pedigree(trait = c(...)); recurrent_selection(trait = c(...))",
      "trait", "Hazel & Lush 1942"),
    m("random", "Random selection (drift control)",
      'select_ind(method = "random")', "n | prop", "package"),
    m("ocs", "Optimum contribution selection", "optimum_contribution()",
      "merit, G, lambda | target_coancestry | max_coancestry", "Meuwissen 1997"),
    m("usefulness", "Cross usefulness", "cross_usefulness()",
      "pairs, scheme, n_progeny, select_top",
      "Zhong & Jannink 2007; Lehermeier et al. 2017"),
    m("mabc", "Marker-assisted backcrossing", "mabc_select()",
      "recurrent, donor, target_markers, flanking_markers",
      "Frisch & Melchinger 2001, 2005"),
    m("mas", "Marker-assisted selection / gene pyramiding", "marker_select()",
      "markers, favorable, requirement, min_markers, rank_on",
      "Lande & Thompson 1990"),
    m("combining_ability", "GCA / SCA / testcross merit", "combining_ability()",
      "testers, qtn, a, d, design, method",
      "Sprague & Tatum 1942; Griffing 1956"),
    m("progeny_test", "Progeny test", "progeny_test()",
      "mates, qtn, a, d, n_progeny, h2 | var_e", "package derivation"),
    m("blup", "BLUP breeding values (known variances)", "predict_ebv()",
      "pheno, method, h2 | var_a + var_e", "Henderson 1975; VanRaden 2008")
  )
}
