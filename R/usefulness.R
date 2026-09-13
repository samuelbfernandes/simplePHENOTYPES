# Cross usefulness: rank candidate biparental crosses by the expected value of
# their best progeny, U = mu + i * sigma (Zhong & Jannink 2007; Lehermeier et al.
# 2017). Because simplePHENOTYPES can simulate the family, the progeny mean and
# genetic SD come from the crossing engine, which captures linkage in the variance
# rather than approximating it. Progeny are scored on the template model's FIXED
# additive effects (a plain dosage x effect dot product, no per-family variance
# rescaling), so mu and sigma are comparable across crosses -- the property a
# re-simulated, self-standardizing genetic value would destroy. This is the
# simulation ("true effects known") counterpart of a PopVar-style prediction from
# estimated marker effects (roadmap: cross selection from real data).

#' Predict the usefulness of candidate crosses
#'
#' For each candidate biparental cross, simulates a progeny family and reports its
#' additive genetic mean, genetic standard deviation, and \emph{usefulness}
#' \eqn{U = \mu + i\,\sigma} (Zhong \& Jannink 2007; Lehermeier et al. 2017). A
#' cross with a lower mean but more genetic variance can outrank a safer one,
#' which is the information parent-mean ranking misses.
#'
#' \eqn{U = \mu + i\,\sigma} is a **normal-theory criterion**: it equals the
#' expected genetic value of the best `select_top` fraction of progeny only when
#' that family's genetic values are approximately normally distributed, with `i`
#' the standardized selection intensity for the fraction. For a large
#' polygenic family this holds well; for a small family or one driven by a few
#' large-effect QTNs the distribution is discrete/skewed and `U` is an
#' approximation that can exceed the realized top-fraction mean (and even the
#' family maximum). Use it as a *ranking* criterion under that assumption; read
#' the realized family from the returned `mean`/`sd` (or simulate the family
#' directly with [cross()]/[selfcross()]) if you need the exact tail mean.
#'
#' Progeny are scored on the additive QTN effects already realized in `sim` (the
#' template), applied on a common scale, so the means and standard deviations are
#' directly comparable between crosses. Only the additive component enters the
#' score (the breeding-value basis of the usefulness criterion); dominance and
#' epistasis layers are ignored, and doubled-haploid / inbred families carry no
#' dominance anyway.
#'
#' @param sim a realized, Population-backed `phenotype_sim` whose candidate parents
#'   are its individuals and whose additive layer supplies the scoring model.
#' @param pairs a two-column matrix or data frame of parent pairs (ids or
#'   positions). `NULL` (default) uses all pairwise combinations of the parents ---
#'   guard the count yourself for large panels.
#' @param scheme how the progeny family is derived: `"dh"` (doubled haploids, fully
#'   inbred in one step), `"selfcross"` (self to near-inbred over `generations`), or
#'   `"cross"` (single-cross progeny; only informative for non-inbred parents).
#' @param n_progeny family size simulated per cross.
#' @param generations selfing generations for `scheme = "selfcross"`.
#' @param select_top fraction of progeny whose intensity `i` sets the usefulness
#'   horizon (default 0.1; smaller = a more elite target).
#' @param trait trait index to evaluate (default 1).
#' @param direction `"high"` (default) adds `i*sigma`; `"low"` subtracts it.
#' @param seed optional RNG seed for the whole evaluation.
#' @return a data frame with `parent1`, `parent2`, `mean`, `sd`, `usefulness`,
#'   sorted best first, carrying attribute `intensity` (the `i` used).
#' @references Zhong S, Jannink J-L (2007) Using quantitative trait loci results to
#'   discriminate among crosses on the basis of their progeny mean and variance.
#'   \emph{Genetics} 177(1):567--576. \doi{10.1534/genetics.107.075358};
#'   Lehermeier C, Teyssedre S, Schon C-C (2017) Genetic gain increases by applying
#'   the usefulness criterion with improved variance prediction in selection of
#'   crosses. \emph{Genetics} 207(4):1651--1661. \doi{10.1534/genetics.117.300403}
#' @seealso [cross()], [single_seed_descent()], [select_ind()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:6)
#' sim <- simulate_phenotype(pop, h2 = 0.5, seed = 1) |> additive(n_qtn = 30)
#' u <- cross_usefulness(sim, scheme = "dh", n_progeny = 30, seed = 2)
#' head(u)   # crosses ranked by expected best-progeny value
cross_usefulness <- function(sim, pairs = NULL,
                             scheme = c("dh", "selfcross", "cross"),
                             n_progeny = 100L, generations = 5L,
                             select_top = 0.1, trait = 1L,
                             direction = c("high", "low"), seed = NULL) {
  .check_sim(sim)
  pop <- .as_founder_pop(sim)
  scheme <- match.arg(scheme)
  direction <- match.arg(direction)
  n_progeny <- .validate_count(n_progeny, "n_progeny", minimum = 2L)
  if (!is.numeric(select_top) || length(select_top) != 1L ||
      select_top <= 0 || select_top >= 1) {
    stop("`select_top` must be a single value in (0, 1).", call. = FALSE)
  }
  i_val <- .intensity_from_p(select_top)
  model <- .additive_model(sim, trait)          # fixed loci + effects (by name)
  pairs <- .resolve_pairs(pairs, pop)
  if (!is.null(seed)) set.seed(seed)

  res <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    a <- pop[pairs[k, 1]]
    b <- pop[pairs[k, 2]]
    fam <- .make_family(a, b, scheme, n_progeny, generations)
    gv <- .additive_gv(fam, model)
    mu <- mean(gv)
    sdev <- stats::sd(gv)
    u <- if (direction == "low") mu - i_val * sdev else mu + i_val * sdev
    res[[k]] <- data.frame(parent1 = a$ids, parent2 = b$ids,
                           mean = mu, sd = sdev, usefulness = u,
                           stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, res)
  out <- out[order(out$usefulness, decreasing = direction == "high"), ]
  rownames(out) <- NULL
  attr(out, "intensity") <- i_val
  out
}

# ---- internal ---------------------------------------------------------------

#' Standardized selection intensity for the top proportion p under normality
#' @keywords internal
#' @noRd
.intensity_from_p <- function(p) {
  stats::dnorm(stats::qnorm(1 - p)) / p
}

#' Fixed additive scoring model for a trait: marker names and effects, summed over
#' every additive layer (rep 1). Names, not indices, so it maps onto any family
#' sharing the marker map.
#' @keywords internal
#' @noRd
.additive_model <- function(sim, trait) {
  trait <- as.integer(trait)
  add <- Filter(function(l) identical(l$type, "additive"), sim$layers)
  if (length(add) == 0L) {
    stop("cross_usefulness(): `sim` has no additive layer to score progeny with.",
         call. = FALSE)
  }
  nt <- sim$n_traits
  snp <- character(0)
  eff <- numeric(0)
  for (ly in add) {
    # rep-1 fixed loci/effects, matching how the template realized the phenotype.
    if (!is.null(ly$qtn_reps)) {
      idx <- ly$qtn_reps[[1L]][[trait]]
      e   <- ly$effect_reps[[1L]][[trait]]
    } else {
      idx <- ly$qtn[[trait]]
      e   <- ly$effect[[trait]]
    }
    if (is.null(idx) || length(idx) == 0L) next
    # Scale raw effects to the template's REALIZED additive scale -- the same
    # transform .genetic_matrix() applies: comp / sd(comp) * sqrt(prop), with sd
    # taken over the template individuals. Concatenating raw effects would weight
    # layers of different `prop` unequally and can flip cross rankings (they enter
    # mu and sigma on the wrong scale).
    prop_t <- .expand_prop(ly$prop, nt)[trait]
    comp <- .component_raw(ly, sim, trait, 1L)
    s <- stats::sd(comp)
    if (!(is.finite(s) && s > 0 && prop_t > 0)) next   # realizes to zero variance
    snp <- c(snp, sim$map$snp[idx])
    eff <- c(eff, as.numeric(e) * (sqrt(prop_t) / s))
  }
  if (length(snp) == 0L) {
    stop("cross_usefulness(): the additive layer has no usable QTN for trait ",
         trait, " (no polymorphic effect variance).", call. = FALSE)
  }
  list(snp = snp, eff = eff)
}

#' Additive genetic value of a family on the template's fixed effect scale
#' @keywords internal
#' @noRd
.additive_gv <- function(family, model) {
  d <- dosages(family)                       # markers x individuals, -1/0/1
  hit <- match(model$snp, rownames(d))
  if (anyNA(hit)) {
    stop("cross_usefulness(): scoring markers are absent from a progeny family ",
         "(map mismatch).", call. = FALSE)
  }
  as.numeric(crossprod(d[hit, , drop = FALSE], model$eff))   # ind-length vector
}

#' Normalize `pairs` to a 2-column integer matrix of positions into `pop`
#' @keywords internal
#' @noRd
.resolve_pairs <- function(pairs, pop) {
  ids <- pop$ids
  if (is.null(pairs)) {
    if (length(ids) < 2L) {
      stop("cross_usefulness(): need at least two candidate parents.",
           call. = FALSE)
    }
    return(t(utils::combn(length(ids), 2L)))
  }
  pairs <- as.matrix(pairs)
  if (ncol(pairs) != 2L) {
    stop("`pairs` must have two columns (one parent per column).", call. = FALSE)
  }
  idx <- if (is.character(pairs)) {
    matrix(match(pairs, ids), ncol = 2L)
  } else {
    matrix(as.integer(pairs), ncol = 2L)
  }
  if (anyNA(idx) || any(idx < 1L | idx > length(ids))) {
    stop("`pairs` refers to parents not in `pop`.", call. = FALSE)
  }
  idx
}

#' Simulate one progeny family from parents a, b under the chosen scheme
#' @keywords internal
#' @noRd
.make_family <- function(a, b, scheme, n_progeny, generations) {
  if (scheme == "cross") {
    return(cross(a, b, n = n_progeny, seed = NULL))
  }
  f1 <- cross(a, b, n = 1L, seed = NULL)
  if (scheme == "dh") {
    return(double_haploid(f1, n = n_progeny, seed = NULL))
  }
  fam <- selfcross(f1, n = n_progeny, seed = NULL)          # F2
  generations <- .validate_count(generations, "generations", minimum = 1L)
  if (generations > 1L) {
    fam <- single_seed_descent(fam, generations = generations - 1L, seed = NULL)
  }
  fam
}
