# Known-variance BLUP (DECISION-030; multi-trait, DECISION-032): predicted
# breeding values from observable phenotypes and a relationship matrix, the pedigree relationship matrix, accuracy
# reporting and a manifest of the engine's selection operators.

#' Numerator relationship matrix from the recorded pedigree
#'
#' The additive relationship matrix \eqn{A} of a `Population`'s individuals,
#' built from its recorded pedigree (see [parentage()]) by the tabular method
#' (Emik & Terrill 1949) on coefficients of relationship (Wright 1922):
#' processing individuals parents-first, \eqn{A_{ij} = (A_{i,m(j)} + A_{i,f(j)}) / 2}
#' and \eqn{A_{jj} = 1 + F_j} with \eqn{F_j = A_{m(j) f(j)} / 2}. A self has both
#' parents equal, so \eqn{F = (1 + F_P) / 2}; a doubled haploid is fully inbred,
#' \eqn{A_{jj} = 2}, and relates to others as its parent does. Ancestors not in
#' `pop` still enter through the pedigree.
#'
#' Founders are the base of the pedigree: unrelated to each other, and by
#' default **non-inbred** (`founder_f = 0`, \eqn{A_{jj} = 1}), which is right for
#' outbred (random-mating) founders. Inbred founders -- e.g. the bundled maize
#' inbred lines -- are not: their inbreeding coefficient enters the base as
#' \eqn{A_{jj} = 1 + F_j} and their relationships with descendants, so pass
#' `founder_f = 1` for fully inbred lines (a founder \eqn{P} crossed to another
#' inbred line then has \eqn{A_{P,F1} = 1} and the two founders' \eqn{F1} sibs
#' \eqn{A = 1}, against \eqn{0.5} with the default). Leaving the default on inbred
#' founders halves those relationships and disagrees systematically with the
#' genomic matrix of [g_matrix()], whose founder diagonal is close to 2 for inbred
#' lines. Downstream individuals follow from the recursion; selfing and doubled-haploid
#' rules are not altered by `founder_f`. Note that a `Population` in which the same
#' individual appears twice (e.g. `c(pop, pop[1])`) lists that individual under
#' two ids with a relationship equal to its own diagonal; `a_matrix()` does not
#' refuse it.
#'
#' @param pop a `Population`.
#' @param ids optional character ids (of `pop`'s individuals) to return;
#'   default all.
#' @param founder_f inbreeding coefficient of the founders, in \eqn{[0, 1]}: one
#'   value for every founder (default `0`, non-inbred), or a numeric vector named
#'   by founder (founders not named keep `0`). A name is a pedigree key (the
#'   `key` column of `parentage(pop, ancestors = TRUE)`, unique by construction)
#'   or a founder's display id. Display ids may repeat across founder pools, so an
#'   id shared by several founders is ambiguous and is an error naming the
#'   pedigree keys to use instead; names must be distinct and identify pedigree
#'   founders. `1` is a fully inbred line.
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
#' # the bundled lines are inbred: declare it (A between a line and its F1 is 1)
#' a_matrix(c(pop, f1), founder_f = 1)
a_matrix <- function(pop, ids = NULL, founder_f = 0) {
  .check_population(pop)
  pop <- .ensure_pedigree(pop)
  ped <- pop$pedigree
  ped <- ped[order(ped$generation), , drop = FALSE]       # parents first
  n <- nrow(ped)
  A <- matrix(0, n, n)
  mi <- match(ped$mother, ped$key)
  fi <- match(ped$father, ped$key)
  ff <- .founder_inbreeding(founder_f, ped, mi, fi)
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
    } else if (is.na(m) && is.na(f)) {
      1 + ff[k]                                   # a founder: base inbreeding
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

#' Founder inbreeding coefficients aligned with the pedigree rows (0 for
#' non-founders); validates `founder_f`
#' @keywords internal
#' @noRd
.founder_inbreeding <- function(founder_f, ped, mi, fi) {
  bad <- function(why) {
    stop("a_matrix(): `founder_f` ", why, call. = FALSE)
  }
  if (!is.numeric(founder_f) || !length(founder_f) || anyNA(founder_f) ||
      any(!is.finite(founder_f)) || any(founder_f < 0 | founder_f > 1)) {
    bad("must be numeric, finite and in [0, 1].")
  }
  is_f <- is.na(mi) & is.na(fi)
  ff <- numeric(nrow(ped))
  if (length(founder_f) == 1L && is.null(names(founder_f))) {
    ff[is_f] <- founder_f
    return(ff)
  }
  nm <- names(founder_f)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) {
    bad("must be one value or a vector named by founder (pedigree key or id).")
  }
  if (anyDuplicated(nm)) {
    bad(paste0("has duplicated names (", paste(unique(nm[duplicated(nm)]),
                                              collapse = ", "), ")."))
  }
  # a name is a pedigree KEY (unique) or, failing that, a founder display id;
  # display ids may legally repeat across pools, so an id shared by several
  # founders is ambiguous and must be given as keys
  fk <- ped$key[is_f]
  fid <- ped$id[is_f]
  target <- vapply(nm, function(z) {
    if (z %in% fk) return(z)
    hit <- fk[fid == z]
    if (length(hit) == 1L) return(hit)
    if (length(hit) > 1L) {
      bad(paste0("name \"", z, "\" is ambiguous: ", length(hit),
                 " founders share that id (pedigree keys ",
                 paste(utils::head(hit, 5), collapse = ", "),
                 "); name them by pedigree key (see parentage())."))
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
  if (anyNA(target)) {
    bad(paste0("names are not founder ids or keys of `pop`'s pedigree: ",
               paste(utils::head(nm[is.na(target)], 5), collapse = ", "), "."))
  }
  if (anyDuplicated(target)) {
    bad(paste0("names the same founder twice (",
               paste(utils::head(nm[target %in% target[duplicated(target)]], 5),
                     collapse = ", "), ")."))
  }
  ff[is_f] <- ifelse(fk %in% target, unname(founder_f)[match(fk, target)], 0)
  ff
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
#' **Scope (single trait).** One trait, one record per individual (a `pheno` name must not be
#' repeated, and two ids of `x` that are the same individual -- e.g. after
#' `c(pop, pop[1])` -- may not both carry a record), an intercept as the only
#' fixed effect, and no missing phenotypes (they are refused, not dropped); a
#' single record gives every EBV zero and every reliability zero -- `NA` for an
#' individual whose prior variance \eqn{K_{ii}\sigma^2_A} is zero, where reliability
#' is undefined -- since one record cannot separate the mean from the breeding
#' value. Repeated-record and single-step models are not covered; for several traits see
#' *Multi-trait BLUP* below.
#' With `h2`, \eqn{\sigma^2_A = h^2 \sigma^2_P} and
#' \eqn{\sigma^2_e = (1 - h^2) \sigma^2_P}, with \eqn{\sigma^2_P} the sample
#' variance of `ref` or, by default, of the **phenotyped subset** `pheno`. Because
#' both components are proportional to \eqn{\sigma^2_P}, the ratio
#' \eqn{\lambda = (1 - h^2) / h^2}, the EBVs and the reliabilities do not depend on
#' which \eqn{\sigma^2_P} is used (multiplying it by any factor leaves them
#' unchanged); only the reported absolute `var_a` and `var_e` attributes scale
#' with it. So `ref` matters when those absolute variances do -- e.g. after
#' selection the phenotyped subset has a truncated variance, and the base
#' population's phenotypes as `ref` (or `var_a` and `var_e` directly) keep the
#' reported components at the base population's scale.
#'
#' **Multi-trait BLUP.** Give `pheno` as an individuals x traits matrix (row names
#' ids, column names traits, `NA` for a trait an individual was not recorded on)
#' and `var_a`, `var_e` as the known traits x traits base additive genetic
#' covariance matrix \eqn{G_0} and residual covariance matrix \eqn{R_0}. The model
#' (Henderson & Quaas 1976) is \eqn{y_t = 1\mu_t + Z_t u_t + e_t} for every trait
#' \eqn{t}, with \eqn{Var(u) = G_0 \otimes K} and residuals correlated only within
#' an individual, \eqn{Cov(e_{it}, e_{js}) = R_{0,ts}} if \eqn{i = j} and 0
#' otherwise. It is solved in the same generalized least-squares form over the
#' stacked records, \eqn{\hat u = Cov(u, y) V^{-1} (y - X\hat\mu)} with
#' \eqn{V = Var(y)}, so `K` and \eqn{G_0} may be singular (e.g. a genetic
#' correlation of 1). A trait is predicted for every individual, including those
#' not recorded on it, through the genetic covariances with the traits that were
#' recorded; an individual's own record on a correlated trait also enters through
#' the residual covariance. The records need not be complete, but every trait needs
#' at least one. With uncorrelated traits (diagonal \eqn{G_0} and \eqn{R_0}) the
#' result is single-trait BLUP of each trait. The BLUP of an aggregate genotype
#' \eqn{H = a'u} is \eqn{a'\hat u}, so rank on `ebv %*% a` (e.g.
#' `select_ind(method = "mass", on = ...)`), or pass the matrix to
#' `select_ind(method = "culling", on = ...)`.
#'
#' @param x a `Population` holding every individual to predict (phenotyped and
#'   not).
#' @param pheno a numeric vector of phenotypes named by id (a subset of `x`), or,
#'   for multi-trait BLUP, a numeric individuals x traits matrix (or data frame)
#'   with row names ids of `x`, column names the traits and `NA` for a missing
#'   record.
#' @param method `"gblup"` or `"pedigree"`.
#' @param h2,var_a,var_e the variance components: either `h2` (then
#'   \eqn{\sigma^2_A = h^2 \sigma^2_P}, \eqn{\sigma^2_P} the variance of `ref`, default
#'   `pheno`) or both `var_a` (base additive variance) and `var_e`. For
#'   multi-trait BLUP, `var_a` and `var_e` are required: the traits x traits base
#'   additive genetic and residual covariance matrices (named by trait, or in
#'   `pheno`'s column order), each positive semidefinite with positive variances;
#'   `h2` and `ref` are not used.
#' @param ref optional numeric vector (at least two finite values) of
#'   base-population phenotypes setting \eqn{\sigma^2_P} for `h2`; it is an error
#'   without `h2` (with `var_a` and `var_e` it would be silently unused). Not used
#'   by multi-trait BLUP.
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
#'   marker scale, so none are returned). For multi-trait BLUP, an individuals x
#'   traits matrix, with `reliability` the matching matrix, `mu` a vector named by
#'   trait, `var_a` and `var_e` the covariance matrices used, `method`, and
#'   `marker_effects` (markers x traits) under the same condition; no `lambda`.
#' @references
#' Henderson CR, Quaas RL (1976) Multiple trait evaluation using relatives'
#'   records. \emph{Journal of Animal Science} 43:1188--1197.
#'   \doi{10.2527/jas1976.4361188x}
#'
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
#'
#' # multi-trait (illustrative covariances): a second, correlated trait recorded
#' # on every individual, the first unrecorded on the 20 candidates
#' y2 <- phenotype_value(pop, q, c(0.8, 0.6, 0), h2 = 0.5, seed = 2)
#' Y <- cbind(t1 = y, t2 = y2)
#' Y[41:60, "t1"] <- NA
#' tr <- c("t1", "t2")
#' G0 <- matrix(c(1, 0.6, 0.6, 1), 2, dimnames = list(tr, tr))
#' R0 <- diag(2); dimnames(R0) <- list(tr, tr)
#' ebv2 <- predict_ebv(pop, Y, var_a = G0, var_e = R0)
#' head(ebv2)
predict_ebv <- function(x, pheno, method = c("gblup", "pedigree"), h2 = NULL,
                        var_a = NULL, var_e = NULL, ref = NULL, K = NULL,
                        base_freq = NULL, ridge = 0) {
  .check_population(x)
  method <- match.arg(method)
  ids <- x$ids
  # an individuals x traits table of records is multi-trait BLUP
  multi <- is.matrix(pheno) || is.data.frame(pheno)
  if (multi) {
    pheno <- .check_pheno_matrix(pheno, ids)
    vc <- .mt_variances(colnames(pheno), h2, var_a, var_e, ref)
  } else {
    if (!is.numeric(pheno) || is.null(names(pheno)) || !length(pheno) ||
        any(!is.finite(pheno))) {
      stop("predict_ebv(): `pheno` must be a finite numeric vector named by ",
           "id, or an individuals x traits matrix.", call. = FALSE)
    }
    if (anyDuplicated(names(pheno)) || !all(names(pheno) %in% ids)) {
      stop("predict_ebv(): `pheno` names must be distinct ids of `x`.",
           call. = FALSE)
    }
    vc <- .blup_variances(pheno, h2, var_a, var_e, ref)
  }
  # one individual listed under two ids would have its record counted twice
  rec_ids <- if (multi) rownames(pheno)[rowSums(!is.na(pheno)) > 0] else names(pheno)
  kx <- .ensure_pedigree(x)$keys[match(rec_ids, ids)]
  if (anyDuplicated(kx)) {
    stop("predict_ebv(): `pheno` gives records under two ids (",
         paste0("\"", rec_ids[kx == kx[anyDuplicated(kx)]], "\"",
                collapse = ", "),
         ") of the same individual; give each individual's record once.",
         call. = FALSE)
  }
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
  if (multi) {
    return(.predict_ebv_mt(x, pheno, method, K, built, vc$var_a, vc$var_e,
                           base_freq, ridge))
  }
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

#' Validate a multi-trait record table: numeric, rows named by distinct ids of
#' `x`, columns by distinct trait names, `NA` for a missing record, and at least
#' one record per trait; returns a numeric matrix
#' @keywords internal
#' @noRd
.check_pheno_matrix <- function(pheno, ids) {
  if (is.data.frame(pheno)) {
    if (!all(vapply(pheno, is.numeric, logical(1)))) {
      stop("predict_ebv(): a `pheno` data frame must have only numeric trait ",
           "columns.", call. = FALSE)
    }
    pheno <- as.matrix(pheno)
  }
  if (!is.numeric(pheno) || !length(pheno)) {
    stop("predict_ebv(): a `pheno` matrix must be numeric, individuals x traits.",
         call. = FALSE)
  }
  rn <- rownames(pheno)
  if (is.null(rn) || anyNA(rn) || anyDuplicated(rn) || !all(rn %in% ids)) {
    stop("predict_ebv(): a `pheno` matrix needs row names that are distinct ids ",
         "of `x`.", call. = FALSE)
  }
  tn <- colnames(pheno)
  if (is.null(tn) || anyNA(tn) || any(!nzchar(tn)) || anyDuplicated(tn)) {
    stop("predict_ebv(): a `pheno` matrix needs distinct, non-empty trait ",
         "(column) names.", call. = FALSE)
  }
  # NA is a missing record; anything else must be a finite value
  if (any(is.nan(pheno) | (!is.na(pheno) & !is.finite(pheno)))) {
    stop("predict_ebv(): `pheno` records must be finite (NA marks a missing ",
         "record).", call. = FALSE)
  }
  empty <- tn[colSums(!is.na(pheno)) == 0L]
  if (length(empty)) {
    stop("predict_ebv(): every trait needs at least one record (none for ",
         paste(empty, collapse = ", "), ").", call. = FALSE)
  }
  pheno
}

#' Multi-trait variance components: the known T x T base additive genetic
#' (`var_a`) and residual (`var_e`) covariance matrices over the traits of
#' `pheno`, validated as covariance matrices with positive variances
#' @keywords internal
#' @noRd
.mt_variances <- function(traits, h2, var_a, var_e, ref) {
  nt <- length(traits)
  if (!is.null(h2) || !is.null(ref)) {
    stop("predict_ebv(): multi-trait BLUP needs the covariance matrices ",
         "`var_a` and `var_e` (traits x traits); `h2` and `ref` cannot set ",
         "the genetic and residual covariances between traits.", call. = FALSE)
  }
  one <- function(S, arg) {
    if (is.null(S)) {
      stop("predict_ebv(): give `", arg, "` as a ", nt, " x ", nt,
           " covariance matrix over the traits of `pheno`.", call. = FALSE)
    }
    if (!is.matrix(S) && is.numeric(S) && length(S) == 1L && nt == 1L) {
      S <- matrix(S, 1L, 1L)
    }
    if (!is.matrix(S) || !identical(dim(S), c(nt, nt))) {
      stop("predict_ebv(): `", arg, "` must be a ", nt, " x ", nt,
           " covariance matrix over the traits of `pheno`.", call. = FALSE)
    }
    # named axes must be the traits, identically on both (then reordered);
    # unnamed ones are taken in `pheno`'s column order
    rn <- rownames(S); cn <- colnames(S)
    if (!is.null(rn) || !is.null(cn)) {
      if (is.null(rn) || !identical(rn, cn) || !setequal(rn, traits) ||
          anyDuplicated(rn)) {
        stop("predict_ebv(): `", arg, "` row and column names must both be ",
             "the traits of `pheno`, in the same order on both axes.",
             call. = FALSE)
      }
      S <- S[traits, traits, drop = FALSE]
    }
    S <- .check_relationship(S, arg)
    if (any(diag(S) <= 0)) {
      stop("predict_ebv(): `", arg, "` needs a positive variance for every ",
           "trait.", call. = FALSE)
    }
    dimnames(S) <- list(traits, traits)
    S
  }
  list(var_a = one(var_a, "var_a"), var_e = one(var_e, "var_e"))
}

#' Multi-trait BLUP (Henderson & Quaas 1976) in generalized least-squares form:
#' records stacked over (trait, individual) pairs, Var(u) = G0 (x) K,
#' Var(e) = R0 (x) I (residuals correlated within an individual only), a mean
#' per trait
#' @keywords internal
#' @noRd
.predict_ebv_mt <- function(x, pheno, method, K, built, G0, R0, base_freq,
                            ridge) {
  ids <- x$ids
  n <- length(ids)
  traits <- colnames(pheno)
  nt <- length(traits)
  obs <- which(!is.na(pheno), arr.ind = TRUE)
  ia <- match(rownames(pheno)[obs[, 1]], ids)     # individual of each record
  ta <- as.integer(obs[, 2])                      # trait of each record
  y <- pheno[obs]
  # V = Var(y) = [G0 (x) K + R0 (x) I] over the records; its Cholesky factor
  # checks it is positive definite and inverts it
  V <- G0[ta, ta, drop = FALSE] * K[ia, ia, drop = FALSE] +
    R0[ta, ta, drop = FALSE] * outer(ia, ia, "==")
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) {
    stop("predict_ebv(): the phenotypic covariance of the records ",
         "(var_a (x) K + var_e (x) I) is not positive definite; check `K`, ",
         "`var_a` and `var_e`.", call. = FALSE)
  }
  W <- chol2inv(R)
  X <- outer(ta, seq_len(nt), "==") * 1           # a mean per trait
  WX <- W %*% X
  XWX <- crossprod(X, WX)
  mu <- as.numeric(solve(XWX, crossprod(WX, y)))
  alpha <- as.numeric(W %*% (y - X %*% mu))
  # u_hat = Cov(u, y) V^-1 (y - X mu), Cov(u_t, y_r) = G0[t, t_r] K[, i_r]
  ebv <- K[, ia, drop = FALSE] %*% (G0[ta, , drop = FALSE] * alpha)
  if (any(!is.finite(ebv)) || any(!is.finite(mu))) {
    stop("predict_ebv(): the BLUP solution is not finite (numerical ",
         "overflow); rescale the phenotypes or variance components.",
         call. = FALSE)
  }
  # reliability 1 - PEV / (G0_tt K_ii), with G0_tt K_ii - PEV = c' P c and
  # P = W - W X (X'WX)^-1 X'W (the mean is estimated)
  Pm <- W - WX %*% solve(XWX, t(WX))
  rel <- matrix(NA_real_, n, nt)
  for (t in seq_len(nt)) {
    Ct <- K[, ia, drop = FALSE] * rep(G0[t, ta], each = n)
    explained <- rowSums((Ct %*% Pm) * Ct)
    rel[, t] <- ifelse(diag(K) > 0, explained / (G0[t, t] * diag(K)), NA_real_)
  }
  if (any(rel < -1e-8 | rel > 1 + 1e-8, na.rm = TRUE)) {
    stop("predict_ebv(): reliabilities fall outside [0, 1]; the system is ",
         "numerically unstable (`K` is not positive semidefinite at the ",
         "precision the variance components require).", call. = FALSE)
  }
  rel <- pmin(pmax(rel, 0), 1)
  dimnames(ebv) <- dimnames(rel) <- list(ids, traits)
  attr(ebv, "reliability") <- rel
  attr(ebv, "mu") <- stats::setNames(mu, traits)
  attr(ebv, "var_a") <- G0
  attr(ebv, "var_e") <- R0
  attr(ebv, "method") <- method
  if (method == "gblup" && ridge == 0 && built) {
    # u_t = K[, i_r] (G0[t_r, t] alpha): the single-trait back-solve with the
    # record weights G0[t_r, t] alpha
    # cbind, not vapply: vapply drops to a vector when there is one marker
    me <- do.call(cbind, lapply(seq_len(nt), function(t) {
      .gblup_marker_effects(x, ia, G0[ta, t] * alpha, base_freq)
    }))
    dimnames(me) <- list(x$map$snp, traits)
    attr(ebv, "marker_effects") <- me
  }
  ebv
}

#' Validate a user-supplied relationship (or trait covariance) matrix: numeric,
#' finite, symmetric and positive semidefinite (up to rounding); returns its
#' symmetric form. `arg` names the argument in messages.
#' @keywords internal
#' @noRd
.check_relationship <- function(K, arg = "K") {
  what <- if (arg == "K") "; a relationship matrix is a covariance matrix" else ""
  if (!is.matrix(K) || !is.numeric(K) || any(!is.finite(K))) {
    stop("predict_ebv(): `", arg, "` must be a finite numeric matrix.",
         call. = FALSE)
  }
  if (nrow(K) != ncol(K)) {
    stop("predict_ebv(): `", arg, "` must be square.", call. = FALSE)
  }
  # element-wise tolerances, so a large entry elsewhere cannot hide a local
  # defect
  if (any(abs(K - t(K)) > 1e-8 * pmax(1, abs(K), abs(t(K))))) {
    stop("predict_ebv(): `", arg, "` must be symmetric.", call. = FALSE)
  }
  # symmetrize only where the two triangles differ (rounding-level): an exactly
  # symmetric K is returned bit-for-bit. For a differing pair (x, y) the
  # correctly rounded mean is (x + y) / 2 whenever x + y cannot overflow (both
  # magnitudes <= xmax / 2, or opposite signs); only otherwise is x/2 + y/2
  # used, because halving first would flush subnormal entries (s and 2s would
  # average to s, not 2s)
  asym <- K != t(K)
  if (any(asym)) {
    Kt <- t(K)
    x <- K[asym]
    y <- Kt[asym]
    half <- .Machine$double.xmax / 2
    safe <- (abs(x) <= half & abs(y) <= half) | (sign(x) != sign(y))
    K[asym] <- ifelse(safe, (x + y) / 2, x / 2 + y / 2)
  }
  dg <- diag(K)
  if (any(dg < 0)) {
    stop("predict_ebv(): `", arg, "` has a negative diagonal (a negative ",
         "variance)", what, ".", call. = FALSE)
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
  # a non-finite correlation-scale entry is itself proof of indefiniteness; it
  # must fail directly (as an eigenvalue of -Inf it would also make the
  # tolerance below infinite and pass)
  nonfinite <- any(!z) && any(!is.finite(C))
  ev <- if (!any(!z) || nonfinite) {
    0
  } else {
    eigen(C, symmetric = TRUE, only.values = TRUE)$values
  }
  # only floating-point rounding below zero (n eps max|eigenvalue|, the usual
  # numerical-rank tolerance) is accepted as semidefinite
  tol <- length(ev) * .Machine$double.eps * max(1, abs(ev))
  if (bad_zero || nonfinite || min(ev) < -tol) {
    stop("predict_ebv(): `", arg, "` must be positive semidefinite", what, ".",
         call. = FALSE)
  }
  K
}

#' Variance components for BLUP: h2 (with the phenotypic variance of ref/pheno)
#' or var_a + var_e
#' @keywords internal
#' @noRd
.blup_variances <- function(pheno, h2, var_a, var_e, ref) {
  if (!is.null(ref) && is.null(h2)) {
    stop("predict_ebv(): `ref` only sets the phenotypic variance for `h2`; ",
         "it is not used with `var_a` and `var_e`. Drop `ref` or give `h2`.",
         call. = FALSE)
  }
  if (!is.null(h2)) {
    if (!is.null(var_a) || !is.null(var_e)) {
      stop("predict_ebv(): give either `h2` or both `var_a` and `var_e`, not ",
           "both.", call. = FALSE)
    }
    if (!is.numeric(h2) || length(h2) != 1L || !is.finite(h2) || h2 <= 0 ||
        h2 >= 1) {
      stop("predict_ebv(): `h2` must be one value in (0, 1).", call. = FALSE)
    }
    if (!is.null(ref) && (!is.numeric(ref) || length(ref) < 2L ||
                          any(!is.finite(ref)))) {
      stop("predict_ebv(): `ref` must be a finite numeric vector of at least ",
           "two base-population phenotypes.", call. = FALSE)
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
#' @return A named numeric vector: `accuracy`, `slope`, `n`. Named vectors are
#'   matched by name (each name may occur once; a duplicated name in either vector
#'   is an error, even when the other vector is unnamed) and
#'   only the common, finite pairs are scored. If `ebv` or `truth` is constant
#'   the correlation is undefined: `accuracy` is `NA` (and `slope` is `NA` for a
#'   constant `ebv`), with a warning.
#' @seealso [predict_ebv()]
#' @export
#' @examples
#' prediction_accuracy(c(a = 1, b = 2, c = 3), c(a = 1.2, b = 1.9, c = 3.3))
prediction_accuracy <- function(ebv, truth) {
  if (!is.numeric(ebv) || !is.numeric(truth)) {
    stop("prediction_accuracy(): `ebv` and `truth` must be numeric.",
         call. = FALSE)
  }
  # a duplicated name is ambiguous whether or not the other vector is named
  dup <- c(names(ebv)[duplicated(names(ebv))],
           names(truth)[duplicated(names(truth))])
  if (length(dup)) {
    stop("prediction_accuracy(): duplicated names (",
         paste(utils::head(unique(dup), 5), collapse = ", "),
         ") in `ebv` or `truth`; each id must occur once, or the pairing ",
         "would silently keep the first match.", call. = FALSE)
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
  ve <- stats::var(e); vu <- stats::var(u)
  if (!(ve > 0) || !(vu > 0)) {
    warning("prediction_accuracy(): `",
            if (!(ve > 0)) "ebv" else "truth",
            "` is constant over the matched individuals, so the accuracy is ",
            "undefined (NA).", call. = FALSE)
  }
  c(accuracy = if (ve > 0 && vu > 0) stats::cor(e, u) else NA_real_,
    slope = if (ve > 0) stats::cov(u, e) / ve else NA_real_,
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
    m("bqp", "Relatedness-penalized selection (binary quadratic programming)",
      'select_ind(method = "bqp")', "weights, lambda, min_gain",
      "Montesinos-Lopez et al. 2025"),
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
      "pheno, method, h2 | var_a + var_e", "Henderson 1975; VanRaden 2008"),
    m("mt_blup", "Multi-trait BLUP (known covariance matrices)",
      "predict_ebv()", "pheno, method, var_a + var_e",
      "Henderson & Quaas 1976")
  )
}
