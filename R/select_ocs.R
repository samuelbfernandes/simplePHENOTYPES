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
  if (!is.numeric(d)) {
    stop("g_matrix(): the dosage matrix must be numeric (coded -1/0/1); got a ",
         class(d)[1L], " matrix. Convert with as_numeric(code_as = \"-101\") ",
         "first.", call. = FALSE)
  }
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
#'   so the default is the additive breeding value `"bv"`, not the total
#'   genotypic value `"gv"`: with `"gv"` a non-transmissible dominance/epistasis
#'   advantage (e.g. an overdominant heterozygote) would be chased even though it
#'   is not passed on. `"bv"` is the *transmitting* average effect, computed
#'   analytically from the simulation's own effects as
#'   \eqn{\alpha_j = a_j + d_j(q_j - p_j)} per locus (DECISION-019; the same
#'   quantity as `select_ind(on = "bv")`) -- not a least-squares (Fisher
#'   statistical) projection of the genotypic value, which differs off
#'   Hardy-Weinberg. It is refused (an error) for a model with an epistasis layer
#'   and for `architecture = "complex"`, because those have no per-locus
#'   \eqn{a}/\eqn{d}; pass a numeric `merit` there.
#' @param trait trait index when `merit` is `"gv"`/`"pheno"` (default 1).
#' @param direction `"high"` (default) maximizes merit, `"low"` minimizes it.
#' @param lambda penalty on group coancestry (>= 0). Larger spreads contributions
#'   and lowers coancestry. Give exactly one of `lambda`, `target_coancestry`,
#'   `max_coancestry`.
#' @param target_coancestry desired group coancestry \eqn{\tfrac12 c'Gc}; the
#'   penalty is tuned to meet it (or the minimum attainable, with a warning). A
#'   target above the coancestry of the unconstrained (merit-only) optimum cannot
#'   be met by a non-negative penalty: a warning is issued and `lambda = 0` is
#'   used.
#' @param max_coancestry as `target_coancestry`, but only enforced if the
#'   unconstrained optimum exceeds it.
#' @param G optional precomputed relationship matrix (from [g_matrix()]); built
#'   from `x` when omitted. A supplied `G` must be a symmetric positive
#'   semi-definite matrix with the individual ids as identical row and column
#'   names (an unnamed matrix is an error: merit is aligned to `G` by name).
#' @param rep replication to read merit from (default 1).
#' @param min_contribution contributions below this are treated as zero when
#'   listing selected parents (default 1e-4).
#' @param max_iter,tol Iteration cap (default 10000; a whole number >= 1) and
#'   duality-gap tolerance (positive) for the away-step Frank-Wolfe optimizer. A
#'   warning is issued if it does not reach `tol` within `max_iter`, meaning the
#'   contributions may be sub-optimal. `tol` is on the merit scale; when a
#'   coancestry target is tuned and the merit spread is extreme (outside
#'   1e-6..1e6) the search and the final solve use merit rescaled to unit spread,
#'   which leaves the contributions unchanged (the penalty scales linearly with
#'   merit).
#' @return an `ocs` object: `contributions` (named, summing to 1), `parents` (ids
#'   with contribution above `min_contribution`), `n_parents`, `merit` (expected
#'   \eqn{c'g}, on the original merit scale), `coancestry`
#'   (\eqn{\tfrac12 c'Gc}), and `lambda`. The contributions do not depend on a
#'   constant added to every merit (the optimum lies on the simplex
#'   \eqn{\sum c = 1}); merit is centred internally before the penalty is tuned.
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
  # Control arguments: an NA/non-numeric value would otherwise yield an `ocs`
  # object with NA parents or a cryptic "missing value where TRUE/FALSE needed".
  if (!is.numeric(min_contribution) || length(min_contribution) != 1L ||
      !is.finite(min_contribution) || min_contribution < 0 ||
      min_contribution >= 1) {
    stop("`min_contribution` must be one finite number in [0, 1).",
         call. = FALSE)
  }
  max_iter <- .validate_count(max_iter, "max_iter", minimum = 1L)
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be one finite positive number.", call. = FALSE)
  }

  if (inherits(x, "phenotype_sim")) {
    g <- .criterion_values(x, merit, trait, rep)   # named by the sim's ids
    if (is.null(G)) G <- g_matrix(x)
    .check_G_ids(G)
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
    .check_G_ids(G)
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
  if (any(!is.finite(g))) {
    stop("`merit` must be finite for every individual (found NA/NaN/Inf); ",
         "impute or drop individuals without a merit value.", call. = FALSE)
  }
  if (length(g) != nrow(G)) {
    stop("merit has length ", length(g), " but G is ", nrow(G), "x", nrow(G), ".",
         call. = FALSE)
  }
  .validate_coancestry_matrix(G)
  if (direction == "low") g <- -g
  # The simplex constraint sum(c) = 1 makes the optimum invariant to adding a
  # constant to every merit, so centre the merit before tuning/solving: at a large
  # common offset (e.g. 1e16) the spread would otherwise be lost to floating-point
  # cancellation. `shift` restores the original scale for the reported merit only.
  # When the finite merits span more than .Machine$double.xmax (e.g. -1e308 and
  # 1e308) the subtraction itself overflows, so the merit is halved first
  # (`hf = 2`): the centred merit is then (g - min(g)) / 2, the optimum is
  # unchanged (c is invariant under g -> g/2, lambda -> lambda/2) and lambda is
  # mapped back to the original scale below. Ordinary problems take `hf = 1` and
  # are untouched.
  hf <- if (is.finite(max(g) - min(g))) 1 else 2
  shift <- min(g) / hf
  g <- g / hf - shift

  sc <- 1                                       # merit rescaling (tuned path only)
  lam_n <- NULL                                 # penalty on the solved (g / sc) scale
  gsc <- 1                                      # G is solved as G / gsc
  lam <- if (!is.null(lambda)) {
    if (lambda < 0) stop("`lambda` must be >= 0.", call. = FALSE)
    lambda / hf
  } else {
    target <- if (!is.null(target_coancestry)) target_coancestry else max_coancestry
    tuned <- .tune_lambda(g, G, target, max_iter, tol,
                          enforce_below = !is.null(max_coancestry))
    sc <- tuned$scale
    gsc <- tuned$gscale
    lam_n <- tuned$lambda_n
    tuned$lambda
  }
  if (is.null(lam_n)) lam_n <- lam / sc

  # lambda scales linearly with the merit scale (c is unchanged under
  # g -> g/s, lambda -> lambda/s), so an extreme-scale tuned problem is solved on
  # the unit-spread scale it was tuned on; sc = 1 leaves the problem untouched.
  # (the tuned penalty is used on its own scale: multiplying by sc and dividing
  # back would overflow for merits near .Machine$double.xmax)
  # (likewise a tuned penalty is used on the G / gsc scale it was tuned on: mapping
  # it back to G would overflow for a subnormal-scale G)
  c_opt <- .frank_wolfe(g / sc, if (gsc == 1) G else G / gsc, lam_n, max_iter, tol)
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
  # Merit on the original scale: c sums to 1, so c'(g + shift) = c'g + shift.
  merit_val <- hf * (sum(c_opt * g) + shift)
  if (direction == "low") merit_val <- -merit_val
  coan <- as.numeric(0.5 * crossprod(c_opt, G %*% c_opt))
  keep <- c_opt > min_contribution
  structure(
    list(contributions = stats::setNames(c_opt, rownames(G)),
         parents = rownames(G)[keep],
         n_parents = sum(keep),
         merit = merit_val,
         coancestry = coan,
         lambda = lam * hf),
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
#' `n` parent slots. Feed the result to the crossing primitives or
#' [recurrent_selection()].
#'
#' Two methods are available. `method = "allocate"` (the default) **allocates**
#' the slots so that the realized contributions match the optimized ones as
#' closely as integers allow. This controls the *marginal* count error of each
#' individual; it is an approximation, not a guarantee on the group coancestry
#' \eqn{\tfrac12 c'Gc}, which is quadratic in the contributions and is not
#' preserved exactly by integer counts. The gap between the realized and the
#' optimized coancestry depends on `n` and on \eqn{G} and shrinks as `n` grows
#' (for example with \eqn{G = I} and \eqn{c = (0.5, 0.5)} the optimum is 0.25, but
#' `n = 1` must pick one parent, whose realized coancestry is 0.5). Individual
#' \eqn{i} receives
#' \eqn{\lfloor n c_i \rfloor} slots and the remaining slots go, one each, to the
#' individuals with the largest fractional parts of \eqn{n c_i} (largest-remainder
#' rule; ties are broken by the larger contribution and then at random with
#' `seed`), so every count is within one slot of \eqn{n c_i}. The counts are
#' deterministic except for such exact ties. The order of the slots is shuffled
#' (with `seed`), so repeated copies of a parent are not adjacent.
#'
#' `method = "multinomial"` is the former behaviour: `n` independent draws with
#' replacement, weighted by the contributions. Its expected number of slots is
#' also \eqn{n c_i}, but the *realized* counts carry sampling noise (with
#' \eqn{c = (0.5, 0.5)} and `n = 10` one parent's count ranges widely around 5),
#' so the realized group coancestry is on average above the optimum, more so for
#' small `n`. Use it only when that sampling variation is itself wanted.
#'
#' @param ocs an `ocs` object from [optimum_contribution()].
#' @param pop the `Population` the contributions were computed on.
#' @param n number of parent slots to fill.
#' @param seed optional RNG seed: one non-negative whole number. It orders the
#'   slots (`"allocate"`) or drives the draws (`"multinomial"`); the caller's RNG
#'   state is restored on exit.
#' @param method `"allocate"` (default, largest-remainder allocation) or
#'   `"multinomial"` (independent weighted draws with replacement).
#' @return a `Population` of `n` parents. An individual can occupy several slots;
#'   because the crossing and phenotyping primitives require unique individual
#'   names, repeated copies of one parent are given disambiguated ids (`"A"`,
#'   `"A_1"`, `"A_2"`, ...) while carrying identical genotypes. With
#'   `"allocate"`, individuals whose `n * contribution` is below one may receive no
#'   slot when `n` is small; use a larger `n` if every contributor must be kept.
#'   The returned object carries attribute `source`, a data frame with one row per
#'   slot, in slot (draw) order: `slot` (1..`n`, the position in the returned
#'   `Population`), `id` (the parent's id in `pop`, before the copies are
#'   disambiguated), `index` (its position in `pop`) and `name` (the unique id the
#'   slot has in the returned `Population`). It is the same for both methods and is
#'   dropped by operations that rebuild the `Population` (`c()`, subsetting).
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
sample_parents <- function(ocs, pop, n, seed = NULL,
                           method = c("allocate", "multinomial")) {
  if (!inherits(ocs, "ocs")) stop("`ocs` must be an ocs object.", call. = FALSE)
  .check_population(pop)
  n <- .validate_count(n, "n", minimum = 1L)
  method <- match.arg(method)
  ids <- names(ocs$contributions)
  idx <- match(ids, pop$ids)
  if (anyNA(idx)) {
    stop("The OCS contributions name individuals that are not in `pop`.",
         call. = FALSE)
  }
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    set.seed(seed)
    on.exit(.restore_seed(old))                # the caller's RNG state is restored
  }
  if (method == "multinomial") {
    # idx[sample.int(...)], not sample(idx, ...): with a single candidate `idx` is
    # a length-1 number and sample() would draw from 1:idx.
    drawn <- idx[sample.int(length(idx), size = n, replace = TRUE,
                            prob = ocs$contributions)]
  } else {
    counts <- .allocate_slots(ocs$contributions, n,
                              tiebreak = sample.int(length(idx)))
    drawn <- rep(idx, counts)
    # Only the ORDER of the slots is random (exchangeable), so copies of one
    # parent are not adjacent for downstream pairing; the counts are fixed.
    drawn <- drawn[sample.int(length(drawn))]
  }
  source_id <- pop$ids[drawn]                    # before copies are renamed
  sampled <- pop[drawn]
  # Repeated slots of the same parent carry duplicate ids, which violates the
  # unique-name invariant every downstream primitive (cross/selfcross/
  # simulate_phenotype) enforces. Keep the repeated parent slots (same genotypes)
  # but give them unique names so the documented OCS -> mating workflow runs.
  if (anyDuplicated(sampled$ids)) {
    uid <- make.unique(sampled$ids, sep = "_")
    colnames(sampled$cis) <- colnames(sampled$trans) <- uid
    sampled$ids <- uid
  }
  # Which individual of `pop` filled each slot, in draw order (a pure record of
  # the draw above; it consumes no RNG).
  attr(sampled, "source") <- data.frame(
    slot = seq_along(drawn), id = source_id, index = as.integer(drawn),
    name = sampled$ids, stringsAsFactors = FALSE)
  sampled
}

#' Largest-remainder allocation of `n` slots to contributions
#'
#' `floor(n * c_i)` slots each, then one extra slot to the individuals with the
#' largest fractional parts (ties: larger contribution, then the `tiebreak` rank),
#' so every count is within one of `n * c_i` and the counts sum to `n`.
#' @param contributions numeric vector of non-negative contributions (any scale;
#'   normalised to sum 1).
#' @param n total number of slots.
#' @param tiebreak integer ranks used to order exact ties (default: position).
#' @return an integer vector of slot counts, one per contribution.
#' @keywords internal
#' @noRd
.allocate_slots <- function(contributions, n,
                            tiebreak = seq_along(contributions)) {
  cc <- as.numeric(contributions)
  cc[!is.finite(cc) | cc < 0] <- 0
  if (sum(cc) <= 0) stop("The OCS contributions are all zero.", call. = FALSE)
  cc <- cc / sum(cc)
  q <- n * cc
  base <- floor(q + 1e-9)                        # guard 0.3 * 10 = 2.9999999999
  frac <- pmax(q - base, 0)
  r <- n - sum(base)
  if (r > 0) {
    # Exact ties are judged on rounded keys so 0.1 * 3 vs 0.3 do not split on
    # floating-point noise.
    ord <- order(-round(frac, 9), -round(cc, 12), tiebreak)
    base[ord[seq_len(r)]] <- base[ord[seq_len(r)]] + 1
  } else if (r < 0) {                            # floating overshoot (rare)
    pos <- which(base > 0)
    ord <- pos[order(round(frac[pos], 9), round(cc[pos], 12), -tiebreak[pos])]
    base[ord[seq_len(-r)]] <- base[ord[seq_len(-r)]] - 1
  }
  as.integer(base)
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

#' A supplied G must carry the individual ids the contributions are named by
#'
#' Without dimnames the contributions vector would be unnamed and `parents`
#' empty (an `ocs` object `sample_parents()` cannot use), and merit could only be
#' paired with G by position. Require identical, unique row/column names.
#' @keywords internal
#' @noRd
.check_G_ids <- function(G) {
  if (!is.matrix(G) || nrow(G) != ncol(G)) {
    stop("`G` must be a square relationship matrix.", call. = FALSE)
  }
  if (is.null(rownames(G)) || is.null(colnames(G)) ||
      !identical(rownames(G), colnames(G))) {
    stop("`G` must carry individual ids as identical row and column names ",
         "(as g_matrix() returns); set dimnames(G) <- list(ids, ids).",
         call. = FALSE)
  }
  if (anyDuplicated(rownames(G))) {
    stop("`G` row/column names must be unique individual ids.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Marker x individual -1/0/1 dosage matrix (with id column names) from any input
#' @keywords internal
#' @noRd
.dosage_from <- function(x) {
  if (inherits(x, "Population")) return(dosages(x))
  if (inherits(x, "phenotype_sim")) {
    if (inherits(x$geno, "Population")) {
      # A subset sim (simulate_phenotype(individuals = )) keeps `ind_idx` of the
      # backing population's columns; G must cover exactly the sim's individuals.
      d <- dosages(x$geno)
      if (!is.null(x$ind_idx)) d <- d[, x$ind_idx, drop = FALSE]
      return(d)
    }
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
#'
#' The penalty scales linearly with the merit scale, so when the merit spread is
#' extreme (outside 1e-6..1e6) the search runs on merit rescaled to unit spread
#' and the penalty is returned in the original units together with the scale used
#' (`sc = 1` otherwise, leaving ordinary problems bit-identical). The doubling
#' bracket is therefore always relative to the problem's own scale instead of a
#' fixed absolute cap. Likewise, when the relationship scale (mean diagonal of
#' `G`) is extreme the search runs on `G` and the target divided by that scale,
#' and the penalty is mapped back (a penalty on `G / s` is `lambda / s` on `G`).
#' @return list(lambda, scale, gscale, lambda_n): `lambda_n` is the penalty on
#'   the scale the search ran on (merit `g / scale`, relationships
#'   `G / gscale`); `lambda = lambda_n * scale / gscale` in original units (it
#'   may overflow when the two scales differ by more than double range; the
#'   solve uses `lambda_n`).
#' @keywords internal
#' @noRd
.tune_lambda <- function(g, G, target, max_iter, tol, enforce_below) {
  spread <- diff(range(g))
  sc <- if (is.finite(spread) && spread > 0 && (spread > 1e6 || spread < 1e-6)) {
    spread
  } else {
    1
  }
  gn <- g / sc
  sg <- mean(diag(G))
  sg <- if (is.finite(sg) && sg > 0 && (sg > 1e6 || sg < 1e-6)) sg else 1
  Gn <- G / sg
  target <- target / sg
  coan_at <- function(lam) {
    cc <- .frank_wolfe(gn, Gn, lam, max_iter, tol)
    as.numeric(0.5 * crossprod(cc, Gn %*% cc))
  }
  done <- function(lam_n) {
    list(lambda = lam_n * (sc / sg), scale = sc, gscale = sg, lambda_n = lam_n)
  }
  c0 <- coan_at(0)                          # unconstrained (merit only)
  if (enforce_below && c0 <= target) return(done(0))   # constraint slack
  # Floating-point-correct comparison: c0 and the target are O(1) coancestries, so
  # equality is judged to a small multiple of machine epsilon relative to their
  # magnitude, not by an absolute band (which would hide a genuinely-too-high
  # target such as 0.5000005 against 0.5, or 2e-16 against 1e-16 when G is tiny).
  # The band is purely relative (no absolute floor, which would hide e.g.
  # target = 2*double.xmin against c0 = double.xmin); it is computed by
  # multiplication only (no division) so subnormals cannot overflow, and it
  # collapses to exact equality when both quantities are zero/underflow.
  band <- 16 * .Machine$double.eps * max(abs(target), abs(c0))
  if (!enforce_below && abs(c0 - target) <= band) {
    return(done(0))                         # target equals the unconstrained optimum
  }
  if (!enforce_below && c0 < target) {
    # A non-negative penalty can only lower coancestry, so a target above the
    # unconstrained optimum's coancestry cannot be met; say so instead of quietly
    # returning the (near-)unconstrained solution.
    warning("optimum_contribution(): requested coancestry ", signif(target * sg, 8),
            " is above the coancestry of the unconstrained (merit-only) optimum (",
            signif(c0 * sg, 8), "), which a non-negative penalty cannot raise; using ",
            "lambda = 0 (the unconstrained optimum). Use `max_coancestry` for a ",
            "ceiling that is only enforced when it binds.", call. = FALSE)
    return(done(0))
  }
  # grow an upper penalty until coancestry drops to/below target
  hi <- 1
  for (i in seq_len(60)) {
    if (coan_at(hi) <= target) break
    hi <- hi * 2
  }
  cmin <- coan_at(hi)
  if (cmin > target) {
    warning("optimum_contribution(): requested coancestry ", signif(target * sg, 3),
            " is below the minimum attainable, or needs a larger penalty than the ",
            "search range (lambda up to 2^60 on a unit merit scale): coancestry ",
            signif(cmin * sg, 3), " was reached at the largest penalty tried. Using the ",
            "minimum-coancestry contributions.", call. = FALSE)
    return(done(hi))
  }
  lo <- 0
  for (i in seq_len(100)) {                 # bisection on a monotone function
    mid <- 0.5 * (lo + hi)
    if (coan_at(mid) > target) lo <- mid else hi <- mid
    if (hi - lo < 1e-8 * (1 + hi)) break
  }
  done(hi)
}
