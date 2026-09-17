# Optimum contribution selection (Meuwissen 1997) on a genomic relationship
# matrix (VanRaden 2008). The optimizer is dependency-free: a Frank-Wolfe active
# set on the simplex maximizes merit penalized by group coancestry, and a
# bisection on the penalty hits a target coancestry when one is requested. All of
# this is deterministic linear algebra and stays in R.

#' Genomic relationship matrix (VanRaden 2008)
#'
#' Builds the additive genomic relationship matrix \eqn{G = ZZ'/(2\sum p_j(1-p_j))}
#' from marker dosages, where \eqn{Z = M - 2p} centres the 0/1/2 genotype matrix by
#' twice the allele frequency (VanRaden 2008, method 1).
#'
#' By default `p` is the allele frequency of the current individuals (`base_freq =
#' NULL`), and markers monomorphic *in this set* are dropped. Then `G` is a
#' moving-base matrix: relatedness and the diagonal are expressed relative to this
#' sample, so \eqn{F_i = G_{ii} - 1} is inbreeding relative to the current set, not
#' to a fixed founder base -- across cycles of selection the base shifts, and a
#' marker fixed by selection (p -> 0 or 1) drops out rather than reading as
#' fixation. For a fixed reference (e.g. recurrent selection or founder-relative
#' inbreeding) pass `base_freq`: markers polymorphic in the base are kept even when
#' currently monomorphic, and `2p` centring uses the base frequencies.
#'
#' @param x a `Population`, a Population-backed `phenotype_sim`, or a marker x
#'   individual dosage matrix coded -1/0/1 with individual column names.
#' @param ridge optional blend toward the identity, `G* = (1 - ridge) G + ridge I`,
#'   for a positive-definite matrix when one is needed (default 0).
#' @param base_freq optional fixed base allele frequencies, one per marker (in the
#'   marker order of `x`), each in \[0, 1] and matching the allele the 0/1/2 coding
#'   counts (dosage 2 = homozygote for it). When supplied, `G` uses this fixed base
#'   instead of the current-sample frequencies, and only markers with `base_freq`
#'   strictly in (0, 1) are kept. Default `NULL` (current-sample base).
#' @return an individuals x individuals matrix with id dimnames.
#' @references VanRaden PM (2008) Efficient methods to compute genomic predictions.
#'   \emph{Journal of Dairy Science} 91(11):4414--4423. \doi{10.3168/jds.2007-0980}
#' @seealso [optimum_contribution()], [dosages()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' G <- g_matrix(pop)
#' dim(G)
#' # Genomic inbreeding of the first few lines:
#' round(diag(G)[1:3] - 1, 3)
g_matrix <- function(x, ridge = 0, base_freq = NULL) {
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge < 0 || ridge > 1) {
    stop("`ridge` must be one value in [0, 1]: G is blended as ",
         "(1 - ridge) * G + ridge * I, so a value outside [0, 1] would not be a ",
         "valid (positive semi-definite) relationship matrix.", call. = FALSE)
  }
  d <- .dosage_from(x)                       # markers x individuals, -1/0/1
  if (anyNA(d)) {
    stop("g_matrix(): the genotypes contain missing values. A missing call ",
         "makes a marker's allele frequency NA, which would silently drop that ",
         "marker from the relationship matrix; impute or remove missing ",
         "genotypes first (e.g. filter_geno() / as_numeric(impute=)).",
         call. = FALSE)
  }
  if (!all(d %in% c(-1, 0, 1))) {
    # VanRaden's M = d + 1 assumes -1/0/1 dosages (0/1/2 gene content). An
    # out-of-domain value (e.g. a raw 0/1/2 matrix, or a coding error) would build
    # a plausible-looking but biologically invalid relationship matrix; reject it.
    stop("g_matrix(): genotypes must be coded -1/0/1 (found other values). ",
         "Convert with as_numeric(code_as = \"-101\") first.", call. = FALSE)
  }
  M <- t(d) + 1                              # individuals x markers, 0/1/2
  if (is.null(base_freq)) {
    p <- colMeans(M) / 2                     # current-sample base (moving base)
  } else {
    if (!is.numeric(base_freq) || length(base_freq) != ncol(M) ||
        any(!is.finite(base_freq)) || any(base_freq < 0 | base_freq > 1)) {
      stop("`base_freq` must be a numeric vector of one allele frequency in ",
           "[0, 1] per marker (", ncol(M), "), in the marker order of `x`.",
           call. = FALSE)
    }
    p <- base_freq                           # fixed reference base
  }
  poly <- which(p > 0 & p < 1)
  if (length(poly) == 0L) {
    stop("g_matrix(): no polymorphic markers to build a relationship matrix.",
         call. = FALSE)
  }
  M <- M[, poly, drop = FALSE]
  p <- p[poly]
  Z <- sweep(M, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  G <- tcrossprod(Z) / denom
  if (ridge > 0) {
    G <- (1 - ridge) * G + ridge * diag(nrow(G))
  }
  dimnames(G) <- list(colnames(d), colnames(d))
  G
}

#' Optimum contribution selection
#'
#' Chooses how much each candidate should contribute to the next generation to
#' maximize genetic merit while holding down relatedness --- the modern answer to
#' "truncation selection erodes diversity and long-term gain" (Meuwissen 1997). It
#' maximizes \eqn{c'g - \tfrac{\lambda}{2} c'Gc} over contributions `c` on the
#' simplex (`c >= 0`, `sum(c) = 1`), where `g` is merit and `G` the genomic
#' relationship matrix; the group coancestry of the chosen parents is
#' \eqn{\tfrac{1}{2}c'Gc}. Give a `lambda` to trace the gain--diversity frontier,
#' or a `target_coancestry` / `max_coancestry` and the penalty is found by
#' bisection.
#'
#' @param x a `phenotype_sim` (supplies both merit and, via its genotypes, `G`) or
#'   a `Population` (then `merit` must be a numeric vector or `G` is built from it
#'   and `merit` supplied).
#' @param merit selection criterion: `"bv"` (default), `"gv"`, `"pheno"`, a
#'   numeric vector (named by id or in order), or a function `merit(sim)` --- the
#'   same hook as [select_ind()], so genomic/phenomic EBVs plug straight in.
#'   Meuwissen's OCS maximizes the *transmissible* merit of the parents' progeny,
#'   so the default is the additive breeding value `"bv"` (the Fisher average-effect
#'   projection), not the total genotypic value `"gv"`: with `"gv"` a
#'   non-transmissible dominance/epistasis advantage (e.g. an overdominant
#'   heterozygote) would be chased even though it is not passed on. `"bv"` is
#'   unavailable for `architecture = "complex"`; pass numeric merit there.
#' @param trait trait index when `merit` is `"gv"`/`"pheno"` (default 1).
#' @param direction `"high"` (default) maximizes merit, `"low"` minimizes it.
#' @param lambda penalty on group coancestry (>= 0). Larger spreads contributions
#'   and lowers coancestry. Give exactly one of `lambda`, `target_coancestry`,
#'   `max_coancestry`.
#' @param target_coancestry desired group coancestry \eqn{\tfrac12 c'Gc}; the
#'   penalty is tuned to meet it (or the minimum attainable, with a warning).
#' @param max_coancestry as `target_coancestry`, but only enforced if the
#'   unconstrained optimum exceeds it.
#' @param G optional precomputed relationship matrix (from [g_matrix()]); built
#'   from `x` when omitted.
#' @param rep replication to read merit from (default 1).
#' @param min_contribution contributions below this are treated as zero when
#'   listing selected parents (default 1e-4).
#' @param max_iter,tol Iteration cap (default 10000) and duality-gap tolerance
#'   for the away-step Frank-Wolfe optimizer. A warning is issued if it does not
#'   reach `tol` within `max_iter`, meaning the contributions may be sub-optimal.
#' @return an `ocs` object: `contributions` (named, summing to 1), `parents` (ids
#'   with contribution above `min_contribution`), `n_parents`, `merit` (expected
#'   \eqn{c'g}), `coancestry` (\eqn{\tfrac12 c'Gc}), and `lambda`.
#' @references Meuwissen THE (1997) Maximizing the response of selection with a
#'   predefined rate of inbreeding. \emph{Journal of Animal Science}
#'   75(4):934--940. \doi{10.2527/1997.754934x}
#' @seealso [g_matrix()], [select_ind()], [recurrent_selection()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
#' f1  <- cross(pop[1], pop[2], n = 1, seed = 1)
#' f2  <- selfcross(f1, n = 60, seed = 2)
#' ph  <- simulate_phenotype(f2, h2 = 0.5, seed = 3) |> additive(n_qtn = 30)
#'
#' # Trace the frontier: more penalty -> lower coancestry, lower gain.
#' hi <- optimum_contribution(ph, merit = "bv", lambda = 1)
#' lo <- optimum_contribution(ph, merit = "bv", lambda = 50)
#' c(gain_hi = hi$merit, coan_hi = hi$coancestry,
#'   gain_lo = lo$merit, coan_lo = lo$coancestry)
optimum_contribution <- function(x, merit = "bv", trait = 1L,
                                 direction = c("high", "low"),
                                 lambda = NULL, target_coancestry = NULL,
                                 max_coancestry = NULL, G = NULL, rep = 1L,
                                 min_contribution = 1e-4, max_iter = 10000L,
                                 tol = 1e-10) {
  direction <- match.arg(direction)
  given <- !c(is.null(lambda), is.null(target_coancestry), is.null(max_coancestry))
  if (sum(given) != 1L) {
    stop("Give exactly one of `lambda`, `target_coancestry` or `max_coancestry`.",
         call. = FALSE)
  }
  # The supplied trade-off argument must be a finite numeric scalar; a non-finite
  # (Inf/NaN) or non-numeric value otherwise slips through to the optimizer or the
  # penalty bisection and returns a meaningless result or fails cryptically.
  .tradeoff_name <- c("lambda", "target_coancestry", "max_coancestry")[given]
  .tradeoff_val  <- list(lambda, target_coancestry, max_coancestry)[given][[1L]]
  if (!is.numeric(.tradeoff_val) || length(.tradeoff_val) != 1L ||
      !is.finite(.tradeoff_val)) {
    stop("`", .tradeoff_name, "` must be a single finite numeric value.",
         call. = FALSE)
  }

  if (inherits(x, "phenotype_sim")) {
    g <- .criterion_values(x, merit, trait, rep)   # named by the sim's ids
    if (is.null(G)) G <- g_matrix(x)
    # Align merit to G's row order. When G comes from g_matrix(x) the order already
    # matches; a user-supplied G may be in a different order, and pairing merit with
    # G positionally would optimize contributions for the wrong individuals.
    if (!is.null(rownames(G)) && !is.null(names(g))) {
      if (!all(rownames(G) %in% names(g))) {
        stop("optimum_contribution(): `G` has individuals absent from the merit ",
             "values; cannot align merit to G by name.", call. = FALSE)
      }
      g <- stats::setNames(as.numeric(g[rownames(G)]), rownames(G))
    }
  } else {
    if (is.null(G)) G <- g_matrix(x)
    if (!is.numeric(merit)) {
      stop("With a Population, `merit` must be a numeric vector of merit values ",
           "(one per individual).", call. = FALSE)
    }
    # A named merit vector is aligned to G by name; it must name every individual
    # in G. Falling back to positional pairing when some names are missing would
    # optimize contributions for the wrong individuals (and silently accept
    # extraneous names). Only an unnamed vector is taken positionally.
    if (!is.null(names(merit))) {
      missing <- setdiff(rownames(G), names(merit))
      if (length(missing)) {
        stop("`merit` is named but has no value for individual(s): ",
             paste(utils::head(missing, 10L), collapse = ", "),
             if (length(missing) > 10L) ", ..." else "",
             ". Name every individual in G, or pass an unnamed vector in G's ",
             "row order.", call. = FALSE)
      }
      g <- merit[rownames(G)]
    } else {
      g <- merit
    }
    g <- stats::setNames(as.numeric(g), rownames(G))
  }
  if (length(g) != nrow(G)) {
    stop("merit has length ", length(g), " but G is ", nrow(G), "x", nrow(G), ".",
         call. = FALSE)
  }
  .validate_coancestry_matrix(G)
  if (direction == "low") g <- -g

  lam <- if (!is.null(lambda)) {
    if (lambda < 0) stop("`lambda` must be >= 0.", call. = FALSE)
    lambda
  } else {
    target <- if (!is.null(target_coancestry)) target_coancestry else max_coancestry
    .tune_lambda(g, G, target, max_iter, tol,
                 enforce_below = !is.null(max_coancestry))
  }

  c_opt <- .frank_wolfe(g, G, lam, max_iter, tol)
  # The away-step optimizer converges to `tol` on well-posed problems, so a
  # failure to converge within max_iter is a genuine signal that the returned
  # contributions are sub-optimal (rather than the spurious near-optimum that a
  # too-tight tol used to produce with plain Frank-Wolfe).
  if (!isTRUE(attr(c_opt, "converged"))) {
    warning("optimum_contribution(): the optimizer did not reach tol (",
            signif(tol, 3), ") within max_iter (", max_iter, ") -- final ",
            "duality gap ", signif(attr(c_opt, "gap"), 3), "; the returned ",
            "contributions may be sub-optimal. Raise max_iter or relax tol.",
            call. = FALSE)
  }
  attributes(c_opt) <- NULL
  merit_val <- sum(c_opt * (if (direction == "low") -g else g))
  coan <- as.numeric(0.5 * crossprod(c_opt, G %*% c_opt))
  keep <- c_opt > min_contribution
  structure(
    list(contributions = stats::setNames(c_opt, rownames(G)),
         parents = rownames(G)[keep],
         n_parents = sum(keep),
         merit = merit_val,
         coancestry = coan,
         lambda = lam),
    class = "ocs"
  )
}

#' @export
print.ocs <- function(x, ...) {
  cat("<ocs>  optimum contribution selection\n")
  cat(sprintf("  Parents with contribution > 0: %d of %d\n",
              x$n_parents, length(x$contributions)))
  cat(sprintf("  Expected merit (c'g): %.4g   Group coancestry (0.5 c'Gc): %.4g\n",
              x$merit, x$coancestry))
  cat(sprintf("  Penalty lambda: %.4g\n", x$lambda))
  top <- sort(x$contributions[x$contributions > 0], decreasing = TRUE)
  show <- utils::head(top, 6)
  cat("  Top contributions:\n")
  for (i in seq_along(show)) {
    cat(sprintf("    %-16s %.3f\n", names(show)[i], show[i]))
  }
  if (length(top) > length(show)) cat(sprintf("    ... (%d more)\n",
                                              length(top) - length(show)))
  invisible(x)
}

#' Draw parents for mating in proportion to OCS contributions
#'
#' Turns continuous [optimum_contribution()] contributions into a concrete set of
#' `n` parents by sampling with probability equal to the contributions. Feed the
#' result to the crossing primitives or [recurrent_selection()].
#'
#' @param ocs an `ocs` object from [optimum_contribution()].
#' @param pop the `Population` the contributions were computed on.
#' @param n number of parent slots to draw.
#' @param seed optional RNG seed.
#' @return a `Population` of `n` sampled parents (with replacement, weighted by
#'   contribution). Sampling with replacement can place the same parent in several
#'   slots; because the crossing and phenotyping primitives require unique
#'   individual names, repeated draws of one parent are given disambiguated ids
#'   (`"A"`, `"A_1"`, `"A_2"`, ...) while carrying identical genotypes.
#' @seealso [optimum_contribution()], [recurrent_selection()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
#' f2  <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 60, seed = 2)
#' ph  <- simulate_phenotype(f2, h2 = 0.5, seed = 3) |> additive(n_qtn = 30)
#' oc  <- optimum_contribution(ph, merit = "bv", lambda = 5)
#' mates <- sample_parents(oc, f2, n = 10, seed = 4)
#' n_individuals(mates)
sample_parents <- function(ocs, pop, n, seed = NULL) {
  if (!inherits(ocs, "ocs")) stop("`ocs` must be an ocs object.", call. = FALSE)
  .check_population(pop)
  n <- .validate_count(n, "n", minimum = 1L)
  ids <- names(ocs$contributions)
  idx <- match(ids, pop$ids)
  if (anyNA(idx)) {
    stop("The OCS contributions name individuals that are not in `pop`.",
         call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  drawn <- sample(idx, size = n, replace = TRUE, prob = ocs$contributions)
  sampled <- pop[drawn]
  # Repeated draws of the same parent carry duplicate ids, which violates the
  # unique-name invariant every downstream primitive (cross/selfcross/
  # simulate_phenotype) enforces. Keep the repeated parent slots (same genotypes)
  # but give them unique names so the documented OCS -> mating workflow runs.
  if (anyDuplicated(sampled$ids)) {
    uid <- make.unique(sampled$ids, sep = "_")
    colnames(sampled$cis) <- colnames(sampled$trans) <- uid
    sampled$ids <- uid
  }
  sampled
}

# ---- internal ---------------------------------------------------------------

#' Validate a relationship/coancestry matrix for the OCS optimizer
#'
#' The away-step Frank-Wolfe optimizer maximizes `c'g - (lambda/2) c'Gc` and its
#' line search assumes non-negative curvature `c'Gc >= 0` -- i.e. that `G` is
#' positive semi-definite. `g_matrix()` builds a PSD matrix by construction, but a
#' user may pass any matrix through `G`; a non-PSD (or non-symmetric/non-finite)
#' one makes the objective unbounded along a negative-curvature direction, so the
#' optimizer can return an interior point with a worse objective than a vertex and
#' still report convergence. Reject such matrices up front.
#' @keywords internal
#' @noRd
.validate_coancestry_matrix <- function(G) {
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("`G` must be a square relationship matrix.", call. = FALSE)
  }
  if (any(!is.finite(G))) {
    stop("`G` has non-finite entries (NA/NaN/Inf); it must be a finite ",
         "relationship matrix.", call. = FALSE)
  }
  if (!isSymmetric(unname(G), tol = 1e-8 * max(1, max(abs(G))))) {
    stop("`G` must be symmetric.", call. = FALSE)
  }
  ev <- min(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  tol_psd <- -sqrt(.Machine$double.eps) * max(1, max(abs(diag(G))))
  if (ev < tol_psd) {
    stop("`G` is not positive semi-definite (smallest eigenvalue ", signif(ev, 3),
         "); the optimizer assumes non-negative coancestry curvature and would ",
         "otherwise return a spurious optimum. Build G with g_matrix(), or add a ",
         "ridge toward the identity (g_matrix(ridge = ) / optimum_contribution ",
         "on a valid relationship matrix).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Marker x individual -1/0/1 dosage matrix (with id column names) from any input
#' @keywords internal
#' @noRd
.dosage_from <- function(x) {
  if (inherits(x, "Population")) return(dosages(x))
  if (inherits(x, "phenotype_sim")) {
    if (inherits(x$geno, "Population")) return(dosages(x$geno))
    stop("g_matrix(): this phenotype_sim is not built on a Population; pass a ",
         "Population or a dosage matrix.", call. = FALSE)
  }
  if (is.matrix(x)) {
    if (is.null(colnames(x))) {
      stop("g_matrix(): a dosage matrix needs individual ids as column names.",
           call. = FALSE)
    }
    return(x)
  }
  stop("g_matrix(): expected a Population, a Population-backed phenotype_sim, ",
       "or a dosage matrix.", call. = FALSE)
}

#' Away-step Frank-Wolfe maximization of c'g - (lambda/2) c'Gc on the simplex
#'
#' Plain Frank-Wolfe converges only at an O(1/k) rate here, so a moderate penalty
#' `lambda` leaves it materially short of the optimum within any practical
#' iteration budget (e.g. merit ~7% low at lambda = 50 after 1000 steps). The
#' away-step variant (Lacoste-Julien & Jaggi 2015) adds, at each step, the option
#' to move *away* from the worst vertex in the current support; for this
#' strongly-concave quadratic that yields linear convergence, reaching a tiny
#' duality gap in a few thousand steps. The Frank-Wolfe duality gap
#' `max_k grad_k - grad'c` is the reported convergence measure.
#' @keywords internal
#' @noRd
.frank_wolfe <- function(g, G, lambda, max_iter, tol) {
  n <- length(g)
  c_vec <- rep(1 / n, n)                    # start at the barycentre
  Gc <- as.numeric(G %*% c_vec)
  converged <- FALSE
  last_gap <- Inf
  for (it in seq_len(max_iter)) {
    grad <- g - lambda * Gc
    s <- which.max(grad)                    # Frank-Wolfe vertex (best direction)
    last_gap <- grad[s] - sum(grad * c_vec) # FW duality gap (>= 0)
    if (last_gap <= tol) {
      converged <- TRUE
      break
    }
    active <- which(c_vec > 1e-12)
    a <- active[which.min(grad[active])]    # away vertex (worst in the support)
    d_fw <- -c_vec; d_fw[s] <- d_fw[s] + 1  # toward e_s
    d_aw <- c_vec;  d_aw[a] <- d_aw[a] - 1  # away from e_a
    # Take whichever ascent direction aligns better with the gradient; the away
    # step's length is capped so the moved-from vertex weight stays non-negative.
    if (sum(grad * d_fw) >= sum(grad * d_aw)) {
      d <- d_fw; gmax <- 1
    } else {
      d <- d_aw; gmax <- c_vec[a] / (1 - c_vec[a])
    }
    Gd <- as.numeric(G %*% d)
    quad <- lambda * sum(d * Gd)            # curvature along d (>= 0 for PSD G)
    step <- if (quad <= 0) gmax else min(gmax, sum(grad * d) / quad)
    step <- max(0, step)
    c_vec <- c_vec + step * d
    Gc <- Gc + step * Gd
  }
  c_vec[c_vec < 0] <- 0
  out <- c_vec / sum(c_vec)
  # Report convergence so the caller can warn once (warning here would spam
  # during the lambda-tuning bisection, which calls this repeatedly).
  attr(out, "converged") <- converged
  attr(out, "gap") <- last_gap
  out
}

#' Bisect the penalty so group coancestry meets a target
#' @keywords internal
#' @noRd
.tune_lambda <- function(g, G, target, max_iter, tol, enforce_below) {
  coan_at <- function(lam) {
    cc <- .frank_wolfe(g, G, lam, max_iter, tol)
    as.numeric(0.5 * crossprod(cc, G %*% cc))
  }
  c0 <- coan_at(0)                          # unconstrained (merit only)
  if (enforce_below && c0 <= target) return(0)   # constraint slack
  # grow an upper penalty until coancestry drops to/below target
  hi <- 1
  for (i in seq_len(60)) {
    if (coan_at(hi) <= target) break
    hi <- hi * 2
  }
  cmin <- coan_at(hi)
  if (cmin > target) {
    warning("optimum_contribution(): requested coancestry ", signif(target, 3),
            " is below the minimum attainable (", signif(cmin, 3),
            "); using the minimum-coancestry contributions.", call. = FALSE)
    return(hi)
  }
  lo <- 0
  for (i in seq_len(100)) {                 # bisection on a monotone function
    mid <- 0.5 * (lo + hi)
    if (coan_at(mid) > target) lo <- mid else hi <- mid
    if (hi - lo < 1e-8 * (1 + hi)) break
  }
  hi
}
