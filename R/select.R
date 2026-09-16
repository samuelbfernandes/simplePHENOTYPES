# Selection on a realized phenotype_sim. Methods follow the standard texts
# (Bernardo; Falconer & Mackay; Lynch & Walsh). Selection ranking is deterministic
# and stays in R; the returned Population feeds straight back into the crossing
# core (cross/selfcross/double_haploid) for the next generation.

#' Select individuals from a simulated population
#'
#' Ranks the individuals of a realized [simulate_phenotype()] result on a chosen
#' criterion and returns the selected ones as a [Population][as_population()],
#' ready to cross forward. This is the truncation-selection primitive the
#' breeding-scheme wrappers ([single_seed_descent()], [bulk()], [pedigree()],
#' [recurrent_selection()]) build on.
#'
#' @section Criterion (`on`):
#' \describe{
#'   \item{`"pheno"`}{the observed phenotype -- realistic mass selection. For a
#'     purely additive model this is the classic \eqn{R = i\,h^2\,\sigma_P}
#'     (narrow-sense \eqn{h^2}; Falconer & Mackay). When dominance/epistasis layers
#'     are present the simulated components are not Fisher-orthogonal (SPEC S2), so
#'     only the additive fraction is transmitted and the realized breeding-value
#'     response is governed by the additive variance, not by broad-sense \eqn{H^2}.}
#'   \item{`"gv"`}{the true *total* genetic value (broad-sense: additive plus any
#'     dominance/epistasis). Idealized selection on genetic merit -- an upper bound
#'     on selectable genetic value, but not on breeding-value response, since the
#'     non-additive part is not transmitted to progeny.}
#'   \item{`"bv"`}{the true *breeding* value -- the classical transmissible merit,
#'     \eqn{A_i = \sum_j \alpha_j (x_{ij} - 2p_j)}, summing each causal locus's
#'     average effect of substitution \eqn{\alpha_j = a_j + d_j(q_j - p_j)}
#'     (`.breeding_value_matrix()`). The per-locus effects are reconstructed from
#'     the simulation's own additive/dominance QTN effects (it is a simulation, so
#'     they are known exactly), making this the genetic value transmitted to
#'     random-mated progeny -- robust to both linkage disequilibrium (exact for an
#'     F2) and departures from HWE (after inbreeding/selection). It captures the
#'     additive average effects dominance loci induce away from \eqn{p = 0.5} and
#'     reduces to the additive value for a purely additive model. This is the merit
#'     that governs response to selection in the next generation. Not available
#'     when the model has an epistasis layer or `architecture = "complex"`: an
#'     epistatic term has no per-locus \eqn{a}/\eqn{d}, so its induced additive
#'     average effects cannot be reconstructed and the breeding value would be
#'     incomplete (see `.breeding_value_matrix()`); both cases error rather than
#'     return a partial value. Supply your own predicted values via a
#'     numeric/function criterion there.}
#'   \item{a numeric vector}{one score per individual (named by id or in
#'     population order) -- the hook for **genomic selection, phenomic selection**
#'     and any predicted/estimated value you compute externally.}
#'   \item{a function}{`on(sim)` returning such a vector -- the same hook, resolved
#'     at selection time.}
#' }
#' `on` applies to `method = "mass"`, `"within_family"`, `"among_family"` and
#' `"combined"`. The multi-trait index methods (`"index"`, `"quadratic_index"`)
#' score on all traits' true breeding values and ignore `on` (a single
#' per-individual score cannot supply per-trait breeding values); supplying one
#' with an index method warns. To rank on externally predicted values, compute
#' your index and pass it through `on` with `method = "mass"`.
#'
#' @section Method:
#' `"mass"` truncates on the individual criterion. `"within_family"` keeps the top
#' fraction inside each family, `"among_family"` keeps whole top-ranked families,
#' and `"combined"` ranks on the Lush combined index that optimally weights an
#' individual's own record and its family mean to predict breeding value; the
#' weights are the selection-index solution \eqn{b = V^{-1} c} built from `h2` and
#' the within-family additive relationship `family_relationship` (0.25 half-sibs,
#' 0.5 full-sibs), following Falconer \& Mackay (1996) and Lynch \& Walsh (1998).
#' `"index"` is
#' the Smith--Hazel multi-trait economic index (`weights` = economic weights, one
#' per trait): \eqn{b = P^{-1} G a}, selecting on \eqn{b'y}. `"quadratic_index"` is
#' the nonlinear genomic selection index of Ceron-Rojas et al. (2026),
#' \eqn{\hat I = w'y + y'Wy} (`weights` = linear `w`, `quad_weights` = the
#' symmetric quadratic/cross-product matrix `W`), which captures trait interactions
#' and intermediate optima. `"random"` draws at random (a drift control). Family
#' methods need a `family` grouping; `"combined"` additionally needs `h2`.
#'
#' @param sim a realized `phenotype_sim`, ideally built on a `Population` so the
#'   selected individuals can be crossed on.
#' @param n number of individuals to keep. Exactly one of `n`, `prop`, `intensity`.
#' @param prop proportion of individuals to keep (0-1).
#' @param intensity standardized selection intensity `i`; the number kept is the
#'   count whose expected intensity under normality is closest to `i`.
#' @param on selection criterion (see the Criterion section).
#' @param trait trait index to select on (default 1); ignored for `method =
#'   "index"`, which uses all traits.
#' @param direction `"high"` (default) keeps the largest scores, `"low"` the
#'   smallest.
#' @param method selection method (see the Method section).
#' @param family optional grouping vector (length = individuals) for the family
#'   methods.
#' @param weights economic weights (one per trait) for `method = "index"`, or the
#'   linear weights `w` for `method = "quadratic_index"`.
#' @param quad_weights symmetric `n_traits x n_traits` matrix of quadratic
#'   (diagonal) and cross-product (off-diagonal) weights `W` for
#'   `method = "quadratic_index"` (default: zero, i.e. a linear index).
#' @param h2 narrow-sense heritability of the selection trait, required by
#'   `method = "combined"` to weight family versus individual information.
#' @param family_relationship additive relationship among family members for
#'   `method = "combined"` (default 0.25 = half-sibs; use 0.5 for full-sibs or
#'   selfed families).
#' @param rep replication to select on when several were simulated (default 1).
#' @return the selected individuals as a `Population` (when `sim` is
#'   Population-backed) or their ids, carrying attributes `selected` (ids),
#'   `differential` (selection differential S on the criterion), `intensity`
#'   (realized standardized i), `criterion` and `method`.
#' @references
#' Truncation response and selection intensity: Falconer DS, Mackay TFC (1996)
#'   Introduction to Quantitative Genetics, 4th ed. Longman, Harlow; Lynch M,
#'   Walsh B (1998) Genetics and Analysis of Quantitative Traits. Sinauer
#'   Associates, Sunderland, Massachusetts.
#' Combined (family + individual) selection: Lush JL (1947) Family merit and
#'   individual merit as bases for selection. \emph{The American Naturalist}
#'   81:241--261 (Part I) and 362--379 (Part II). \doi{10.1086/281520}
#' Multi-trait selection index \eqn{b = P^{-1} G a}: Smith HF (1936) A discriminant
#'   function for plant selection. \emph{Annals of Eugenics} 7(3):240--250.
#'   \doi{10.1111/j.1469-1809.1936.tb02143.x}; Hazel LN (1943) The genetic basis
#'   for constructing selection indexes. \emph{Genetics} 28(6):476--490.
#'   \doi{10.1093/genetics/28.6.476}
#' Quadratic (nonlinear) genomic selection index: Ceron-Rojas JJ,
#'   Montesinos-Lopez OA, Montesinos-Lopez A, et al. (2026) Nonlinear genomic
#'   selection index accelerates multi-trait crop improvement. \emph{Nature
#'   Communications} 17:1991. \doi{10.1038/s41467-026-69890-3}
#' Additive average-effect (breeding value) decomposition used as the merit for
#'   the `"index"` / `"quadratic_index"` methods (per-locus average effect of
#'   substitution \eqn{\alpha = a + d(q - p)}): Fisher RA (1918) The correlation between
#'   relatives on the supposition of Mendelian inheritance. \emph{Transactions of
#'   the Royal Society of Edinburgh} 52:399--433. \doi{10.1017/S0080456800012163};
#'   Falconer DS, Mackay TFC (1996) Introduction to Quantitative Genetics, 4th
#'   ed. Longman; Lynch M, Walsh B (1998) Genetics and Analysis of
#'   Quantitative Traits. Sinauer.
#' Breeding schemes: Bernardo R (2020) Breeding for Quantitative Traits in Plants,
#'   3rd ed. Stemma Press.
#' @seealso [single_seed_descent()], [bulk()], [pedigree()],
#'   [recurrent_selection()], [simulate_phenotype()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:60)
#' f1  <- cross(pop[1], pop[2], n = 1, seed = 1)
#' f2  <- selfcross(f1, n = 50, seed = 2)
#' ph  <- simulate_phenotype(f2, h2 = 0.5, seed = 3) |> additive(n_qtn = 20)
#'
#' top <- select_ind(ph, prop = 0.2, on = "pheno")   # keep the best 20%
#' n_individuals(top)
select_ind <- function(sim, n = NULL, prop = NULL, intensity = NULL,
                       on = "pheno", trait = 1L,
                       direction = c("high", "low"),
                       method = c("mass", "within_family", "among_family",
                                  "combined", "index", "quadratic_index",
                                  "random"),
                       family = NULL, weights = NULL, quad_weights = NULL,
                       h2 = NULL, family_relationship = 0.25, rep = 1L) {
  .check_sim(sim)
  # Selecting on a phenotype whose requested h2 was never fully allocated would
  # silently select at the wrong heritability; enforce the same completeness
  # contract the phenotype/QTN accessors do (SPEC 4.1).
  .check_h2_complete(sim)
  if (is.null(sim$pheno)) {
    stop("select_ind(): this phenotype_sim has no realized phenotypes yet. ",
         "Add a genetic layer (e.g. additive()) before selecting.",
         call. = FALSE)
  }
  # Reject an out-of-range replication up front; otherwise a nonexistent rep
  # yields an all-NA criterion and silently "selects" arbitrary individuals.
  rep <- .validate_rep(sim, rep)
  direction <- match.arg(direction)
  method <- match.arg(method)
  ids <- sim$ids
  n_ind <- sim$n_ind

  keep_n <- .resolve_keep(n, prop, intensity, n_ind)

  fam <- if (method %in% c("within_family", "among_family", "combined")) {
    if (is.null(family) || length(family) != n_ind) {
      stop("method = \"", method, "\" needs a `family` grouping vector of ",
           "length ", n_ind, ".", call. = FALSE)
    }
    as.character(family)
  } else NULL

  # The multi-trait index methods score on all traits' breeding values, so a
  # single per-individual `on` cannot feed them; warn rather than silently
  # ignoring an explicitly supplied criterion (which would make an external-GEBV
  # `on` look honoured when it is not).
  if (method %in% c("index", "quadratic_index") && !missing(on) &&
      !identical(on, "pheno")) {
    warning("method = \"", method, "\" ignores `on`: the index scores on all ",
            "traits' breeding values, which a single per-individual `on` cannot ",
            "supply. To rank on your own predictions, pass them via `on` with ",
            "method = \"mass\".", call. = FALSE)
  }
  # The Lush combined index weights own record vs family mean assuming the score
  # is an individual *phenotypic* record: Var(own) = V_P and Cov(A, own) = V_A =
  # h2*V_P. Applying it to a breeding value, genetic value, or arbitrary custom
  # score (for which those identities do not hold) misweights and misranks, so
  # restrict it to on = "pheno".
  if (method == "combined" && !identical(on, "pheno")) {
    stop("method = \"combined\" (Lush index) is defined for phenotypic records ",
         "only: its weights assume Var(own) = V_P and Cov(A, own) = V_A. Use ",
         "on = \"pheno\" (the default), or select on a breeding/genetic value ",
         "with method = \"mass\".", call. = FALSE)
  }

  score <- if (method == "index") {
    .index_score(sim, weights, rep)
  } else if (method == "quadratic_index") {
    .quadratic_index_score(sim, weights, quad_weights, rep)
  } else if (method == "random") {
    stats::setNames(stats::runif(n_ind), ids)
  } else if (method == "combined") {
    .combined_score(.criterion_values(sim, on, trait, rep), fam, h2,
                    family_relationship)
  } else {
    .criterion_values(sim, on, trait, rep)
  }
  if (direction == "low") score <- -score

  sel_idx <- switch(method,
    within_family = .sel_within_family(score, fam, keep_n),
    among_family  = .sel_among_family(score, fam, keep_n),
    .sel_top(score, keep_n))               # mass, combined, index, random

  # Realized selection differential (S) and standardized intensity (i) are computed
  # on the criterion actually used to rank -- the index/quadratic score for those
  # methods, the `on` trait for mass/family selection -- so i = S / sd(criterion) is
  # the realized selection intensity by definition. `score` is already
  # direction-adjusted (selected individuals are its largest values), so S_dir >= 0;
  # i is reported as a positive magnitude and S in the natural (unflipped) sign. A
  # criterion with no spread (e.g. all ties) gives i = 0 rather than NaN.
  S_dir <- mean(score[sel_idx]) - mean(score)
  sd_crit <- stats::sd(score)
  intensity_real <- if (isTRUE(is.finite(sd_crit) && sd_crit > 0)) S_dir / sd_crit else 0

  out <- .selection_result(sim, sel_idx, ids)
  attr(out, "selected") <- ids[sel_idx]
  attr(out, "differential") <- if (direction == "low") -S_dir else S_dir
  attr(out, "intensity") <- intensity_real
  attr(out, "criterion") <- if (is.function(on)) "custom" else
    if (is.numeric(on)) "custom" else on
  attr(out, "method") <- method
  out
}

#' Resolve n / prop / intensity to a count to keep
#' @keywords internal
#' @noRd
.resolve_keep <- function(n, prop, intensity, n_ind) {
  given <- !c(is.null(n), is.null(prop), is.null(intensity))
  if (sum(given) != 1L) {
    stop("Give exactly one of `n`, `prop` or `intensity`.", call. = FALSE)
  }
  keep <- if (!is.null(n)) {
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != floor(n)) {
      stop("`n` must be one whole number.", call. = FALSE)
    }
    as.integer(n)
  } else if (!is.null(prop)) {
    if (!is.numeric(prop) || length(prop) != 1L || !is.finite(prop) ||
        prop <= 0 || prop > 1) {
      stop("`prop` must be one value in (0, 1].", call. = FALSE)
    }
    max(1L, round(prop * n_ind))
  } else {
    # count whose expected standardized intensity under normality is nearest i.
    # For the top fraction p, i = dnorm(qnorm(1 - p)) / p (Falconer & Mackay).
    # A standardized upper-tail intensity is non-negative by definition; a
    # negative value would otherwise silently select the whole population.
    if (!is.numeric(intensity) || length(intensity) != 1L ||
        !is.finite(intensity) || intensity < 0) {
      stop("`intensity` must be one finite, non-negative value (a standardized ",
           "selection intensity).", call. = FALSE)
    }
    p <- seq_len(n_ind) / n_ind
    i_of_p <- stats::dnorm(stats::qnorm(1 - p)) / p
    max(1L, which.min(abs(i_of_p - intensity)))
  }
  if (keep < 1L || keep > n_ind) {
    stop("The number to keep (", keep, ") is outside 1..", n_ind, ".",
         call. = FALSE)
  }
  keep
}

#' Per-individual criterion values (named by id)
#' @keywords internal
#' @noRd
.criterion_values <- function(sim, on, trait, rep) {
  ids <- sim$ids
  if (is.function(on)) {
    v <- on(sim)
  } else if (is.numeric(on)) {
    v <- on
  }
  if (is.function(on) || is.numeric(on)) {
    # A named score is matched to ids by name; an unnamed one is taken in id order.
    # A partially/wrongly named score is an error -- silently treating it as
    # positional would select individuals for whom no score was supplied.
    if (!is.null(names(v))) {
      missing <- setdiff(ids, names(v))
      if (length(missing) > 0L) {
        stop("the named selection criterion is missing scores for: ",
             paste(utils::head(missing, 10L), collapse = ", "),
             if (length(missing) > 10L) ", ..." else "",
             ". Name every individual, or pass an unnamed vector in id order.",
             call. = FALSE)
      }
      v <- v[ids]
    }
  } else if (identical(on, "bv")) {
    v <- .breeding_value_matrix(sim, rep)[, trait]
  } else if (identical(on, "gv")) {
    v <- .genetic_matrix(sim, rep)[, trait]
  } else if (identical(on, "pheno")) {
    ph <- sim$pheno
    ph <- ph[ph$rep == rep & ph$trait == paste0("Trait_", trait), ]
    v <- ph$value[match(ids, ph$id)]
  } else {
    stop("`on` must be \"pheno\", \"gv\", \"bv\", a numeric vector, or a ",
         "function.", call. = FALSE)
  }
  v <- as.numeric(v)
  if (length(v) != length(ids)) {
    stop("The selection criterion has length ", length(v), " but there are ",
         length(ids), " individuals.", call. = FALSE)
  }
  # A non-finite criterion (NA/NaN/Inf) cannot be ranked; ordering it would
  # silently return arbitrary individuals with a NaN differential.
  if (any(!is.finite(v))) {
    stop("The selection criterion has ", sum(!is.finite(v)), " non-finite ",
         "value(s) (NA/NaN/Inf); every individual needs a finite score to be ",
         "ranked.", call. = FALSE)
  }
  stats::setNames(v, ids)
}

#' Quadratic genomic selection index score (Ceron-Rojas et al. 2026)
#'
#' Scores each individual by the nonlinear index
#' \eqn{\hat I = w'\hat\gamma + \hat\gamma' W \hat\gamma}, where
#' \eqn{\hat\gamma} is the individual's vector of (per-trait) breeding values,
#' `w` the linear weights and `W` a symmetric matrix of squared (diagonal) and
#' trait x trait cross-product (off-diagonal) weights. The quadratic term
#' captures nonlinear trait contributions, trait interactions and intermediate
#' optima that the linear Smith-Hazel index cannot.
#'
#' The published QGSI is a *genomic* index: it scores on estimated breeding
#' values (GEBVs), not phenotypes. This package does not fit a genomic-prediction
#' model, so the merit input here is the simulation's *true* breeding value -- the
#' classical transmissible average-effect breeding value (`.breeding_value_matrix()`):
#' \eqn{A_i = \sum_j \alpha_j (x_{ij} - 2p_j)} with per-locus average effects
#' \eqn{\alpha_j = a_j + d_j(q_j - p_j)}. It captures the additive average effects
#' that dominance loci induce (a dominance locus has \eqn{\alpha \ne 0} at
#' \eqn{p \ne 0.5}), and reduces to the additive value for a purely additive model;
#' models with an epistasis layer are refused (see `.breeding_value_matrix()`).
#' With `W = 0` the score reduces to the linear genomic index \eqn{w'\hat\gamma}.
#'
#' This is therefore a simulation of the QGSI's *behaviour* on known-truth merit,
#' not the estimator of Ceron-Rojas et al. (which fits GEBVs from training data):
#' there is no per-trait external-GEBV hook. The `on` argument of [select_ind()]
#' is a single per-individual score and cannot carry the multi-trait breeding
#' values the index needs, so `method = "index"`/`"quadratic_index"` ignore `on`
#' (and warn if one is supplied). To rank on your own predictions, either build
#' the index yourself and pass the result via a numeric/function `on` with
#' `method = "mass"`, or use `method = "mass"` with your GEBVs directly.
#' @keywords internal
#' @noRd
.quadratic_index_score <- function(sim, w, W, rep) {
  .check_sim(sim)
  nt <- sim$n_traits
  if (is.null(w) || length(w) != nt) {
    stop("method = \"quadratic_index\" needs `weights` (linear weights, one per ",
         "trait; ", nt, ").", call. = FALSE)
  }
  if (!is.numeric(w) || any(!is.finite(w))) {
    stop("method = \"quadratic_index\" needs finite numeric `weights`; a ",
         "non-finite (NA/NaN/Inf) linear weight makes the index undefined and ",
         "would silently rank on NA scores.", call. = FALSE)
  }
  if (is.null(W)) {
    W <- matrix(0, nt, nt)
  } else {
    W <- as.matrix(W)
    if (nrow(W) != nt || ncol(W) != nt) {
      stop("`quad_weights` must be an n_traits x n_traits matrix (", nt, "x",
           nt, ").", call. = FALSE)
    }
    if (!is.numeric(W) || any(!is.finite(W))) {
      stop("`quad_weights` must be finite numeric; a non-finite (NA/NaN/Inf) ",
           "quadratic weight makes the index undefined.", call. = FALSE)
    }
    W <- (W + t(W)) / 2                       # use the symmetric part
  }
  # Score on additive breeding values, as the genomic index requires -- NOT on
  # phenotypes, and NOT on the total genetic value (which would let
  # non-transmissible dominance/epistasis drive the ranking). The breeding value
  # is the random-mating transmitting ability (average effect a + d(q - p)), the
  # quantity a GEBV estimates.
  Y <- .breeding_value_matrix(sim, rep)      # ind x trait breeding values
  lin <- as.numeric(Y %*% w)
  quad <- rowSums((Y %*% W) * Y)             # gamma' W gamma per individual
  stats::setNames(lin + quad, sim$ids)
}

#' Smith-Hazel multi-trait index score b = P^{-1} G a
#' @keywords internal
#' @noRd
.index_score <- function(sim, weights, rep) {
  .check_sim(sim)
  nt <- sim$n_traits
  if (is.null(weights) || length(weights) != nt) {
    stop("method = \"index\" needs `weights` (one economic weight per trait; ",
         nt, ").", call. = FALSE)
  }
  if (!is.numeric(weights) || any(!is.finite(weights))) {
    stop("method = \"index\" needs finite numeric `weights`; a non-finite ",
         "(NA/NaN/Inf) economic weight makes the index undefined and would ",
         "silently rank on NA scores.", call. = FALSE)
  }
  # The Smith-Hazel index predicts additive genetic merit, so `G` is the
  # additive breeding-value covariance -- the classical transmissible average
  # effects (which capture the additive effects induced by dominance loci;
  # epistatic models are refused), not the total genotypic covariance.
  G <- .breeding_value_matrix(sim, rep)                # ind x trait breeding values
  wide <- phenotypes_wide(sim)
  wide <- wide[wide$rep == rep, , drop = FALSE]
  P <- as.matrix(wide[, paste0("Trait_", seq_len(nt)), drop = FALSE])
  Pcov <- stats::cov(P)
  Gcov <- stats::cov(G)
  b <- .index_weights(Pcov, Gcov %*% weights)
  stats::setNames(as.numeric(P %*% b), sim$ids)
}

#' Solve the selection-index weights b = P^{-1} r, rank-aware
#'
#' `solve()` errors ("system is exactly singular") when the phenotypic covariance
#' `P` is rank-deficient -- which happens for a perfectly valid model: two traits
#' with correlation 1 (e.g. pleiotropy `cor = 1`, or duplicated/derived traits)
#' realize identical phenotypes and a singular `P`. The Smith-Hazel index is only
#' determined up to the redundant direction there, so any solution ranks
#' individuals identically. Use a Moore-Penrose pseudo-inverse (via SVD, with the
#' standard rank tolerance) so a degenerate but meaningful index still returns the
#' minimum-norm weights instead of crashing; warn once so the user knows to drop
#' the redundant trait. For a full-rank `P` this is the exact inverse, so
#' well-posed indices are unchanged.
#' @keywords internal
#' @noRd
.index_weights <- function(P, r) {
  s <- svd(P)
  tol <- max(dim(P)) * .Machine$double.eps * (if (length(s$d)) s$d[1L] else 0)
  pos <- s$d > tol
  if (!all(pos)) {
    warning("method = \"index\": the phenotypic covariance is singular (rank ",
            sum(pos), " of ", ncol(P), ") -- traits are collinear (e.g. two ",
            "perfectly correlated traits). Using a pseudo-inverse; the index is ",
            "determined only up to the redundant trait(s). Consider dropping ",
            "redundant traits.", call. = FALSE)
  }
  Pinv <- s$v[, pos, drop = FALSE] %*%
    ((1 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
  Pinv %*% r
}

#' Lush combined-selection index score
#'
#' Ranks each individual on the selection-index prediction of its breeding value
#' from two sources of information -- its own record and its family mean -- with
#' weights b = V^{-1} c derived from h2 and the within-family relationship r
#' (Falconer & Mackay; Lynch & Walsh). For a family of size n, with phenotypic
#' variance vP, additive variance vA = h2 vP and phenotypic intraclass
#' correlation t = r h2:
#'   Var(own) = vP;  Var(fam mean) = Cov(own, fam mean) = vP (1+(n-1)t)/n
#'   Cov(A, own) = vA;  Cov(A, fam mean) = vA (1+(n-1)r)/n
#' Both weights are non-negative; the family mean drops out as h2 -> 1.
#' @keywords internal
#' @noRd
.combined_score <- function(values, fam, h2, r) {
  if (is.null(h2) || !is.numeric(h2) || length(h2) != 1L || !is.finite(h2) ||
      h2 <= 0 || h2 > 1) {
    stop("method = \"combined\" needs `h2` in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0 || r > 1) {
    stop("`family_relationship` must be a single value in [0, 1].",
         call. = FALSE)
  }
  vP <- stats::var(values)
  vA <- h2 * vP
  t  <- r * h2
  fam_mean <- tapply(values, fam, mean)
  fam_n    <- tapply(values, fam, length)
  score <- numeric(length(values))
  for (i in seq_along(values)) {
    f <- fam[i]
    n <- fam_n[[f]]
    if (n == 1L) {                       # singleton family: own record only
      score[i] <- values[i]
      next
    }
    Vx2 <- vP * (1 + (n - 1) * t) / n    # = Cov(own, fam mean) too
    c1 <- vA
    c2 <- vA * (1 + (n - 1) * r) / n
    denom <- vP - Vx2                    # = vP (n-1)(1-t)/n
    # denom -> 0 when t = r*h2 -> 1 (the covariance matrix becomes singular) or when
    # vP = 0. There the own record already equals the breeding value and the family
    # mean adds nothing, so rank on the own record rather than dividing by ~0 (NaN).
    if (!is.finite(denom) || abs(denom) <= .Machine$double.eps * max(1, vP) ||
        Vx2 <= 0) {
      score[i] <- values[i]
      next
    }
    b1 <- (c1 - c2) / denom
    b2 <- (vP * c2 - Vx2 * c1) / (Vx2 * denom)
    score[i] <- b1 * values[i] + b2 * fam_mean[[f]]
  }
  stats::setNames(score, names(values))
}

#' Top-keep_n indices of a score vector
#' @keywords internal
#' @noRd
.sel_top <- function(score, keep_n) {
  order(score, decreasing = TRUE)[seq_len(keep_n)]
}

#' Within-family truncation: keep the top fraction inside each family
#' @keywords internal
#' @noRd
.sel_within_family <- function(score, fam, keep_n) {
  groups <- split(seq_along(score), fam)
  sizes  <- lengths(groups)
  ntot   <- length(score)
  # Allocate exactly keep_n across families in proportion to family size
  # (largest-remainder), capped at each family's membership. The old code kept
  # max(1, round(frac * size)) per family, which forces >= 1 from every family
  # and so overshoots keep_n whenever there are more families than the number to
  # keep (e.g. singleton families with keep_n = 1 returned every individual).
  raw   <- keep_n * sizes / ntot
  alloc <- pmin(floor(raw), sizes)
  rem   <- keep_n - sum(alloc)
  if (rem > 0) {
    fracpart <- raw - floor(raw)
    room     <- sizes - alloc
    cand     <- which(room > 0)
    cand     <- cand[order(fracpart[cand], decreasing = TRUE)]
    take     <- utils::head(cand, rem)
    alloc[take] <- alloc[take] + 1L
  }
  out <- integer(0)
  for (g in seq_along(groups)) {
    k <- alloc[[g]]
    if (k > 0L) {
      ix <- groups[[g]]
      out <- c(out, ix[order(score[ix], decreasing = TRUE)[seq_len(k)]])
    }
  }
  out
}

#' Among-family selection: keep every individual of the top-ranked families
#' @keywords internal
#' @noRd
.sel_among_family <- function(score, fam, keep_n) {
  fam_mean <- tapply(score, fam, mean)
  ranked <- names(sort(fam_mean, decreasing = TRUE))
  chosen <- character(0)
  taken <- 0L
  for (f in ranked) {
    members <- which(fam == f)
    chosen <- c(chosen, f)
    taken <- taken + length(members)
    if (taken >= keep_n) break
  }
  which(fam %in% chosen)
}

#' Return the selected individuals as a Population (or ids)
#' @keywords internal
#' @noRd
.selection_result <- function(sim, sel_idx, ids) {
  if (inherits(sim$geno, "Population")) {
    pop_idx <- match(ids[sel_idx], sim$geno$ids)
    return(sim$geno[pop_idx])
  }
  message("select_ind(): the phenotype was not built on a Population, so the ",
          "selected ids are returned rather than a crossable Population. Build ",
          "the population with cross()/selfcross()/double_haploid() to advance ",
          "generations.")
  ids[sel_idx]
}
