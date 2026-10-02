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
#'     (narrow-sense \eqn{h^2}; Falconer & Mackay). In general the expected
#'     response of the mean breeding value is \eqn{R = i\,\mathrm{Cov}(A, P) /
#'     \sigma_P} (the regression of breeding value on the criterion; Smith 1936,
#'     Hazel 1943) **provided \eqn{E[A \mid P]} is linear in \eqn{P}** (as under
#'     joint normality of \eqn{A} and \eqn{P}); it is not an identity for an
#'     arbitrary joint distribution. A discrete counterexample: \eqn{A = (0, 1,
#'     2)} equiprobable with \eqn{P = A^2}, selecting the top third (\eqn{A =
#'     2}), gives an actual response of 1 against a linear-regression value of
#'     about 1.077. Under that linearity assumption the expression equals
#'     \eqn{i\,h^2\,\sigma_P} with \eqn{h^2 = V_A / V_P} only when \eqn{\mathrm{Cov}(A, P - A) = 0}, i.e. when
#'     every non-additive component of the phenotype (dominance, epistasis and the
#'     residual) is uncorrelated with the breeding value. Under random mating /
#'     Hardy-Weinberg proportions with linkage equilibrium this holds for
#'     dominance, but it can fail for epistasis under linkage disequilibrium: the
#'     package's centred epistatic products need not be orthogonal to the additive
#'     values (SPEC S2), so \eqn{\mathrm{Cov}(A, P)} can differ from \eqn{V_A}
#'     (e.g. a two-locus population in Hardy-Weinberg proportions at each locus but
#'     with disequilibrium \eqn{D = 0.05} gives \eqn{\mathrm{Cov}(A, I) \ne 0}).
#'     When dominance/epistasis layers are present the response is therefore
#'     governed by the covariance of the breeding value with the phenotype, not by
#'     broad-sense \eqn{H^2}, and an epistatic term can itself contribute
#'     transmissible (marginal) average effects, which is why `"bv"` is not offered
#'     for epistatic models. Off Hardy-Weinberg (selfed, inbred or selected
#'     populations) the covariance is likewise not \eqn{V_A}.}
#'   \item{`"gv"`}{the true *total* genetic value (broad-sense: additive plus any
#'     dominance/epistasis). Idealized selection on genetic merit -- an upper bound
#'     on selectable genetic value, but not on breeding-value response, since the
#'     non-additive part is not, in general, transmitted to progeny.}
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
#' and `"combined"` ranks on the Lush combined index that weights an
#' individual's own record and its family mean to predict breeding value; the
#' weights are the selection-index solution \eqn{b = V^{-1} c} built from `h2` and
#' `family_relationship`, following Falconer & Mackay (1996) and Lynch & Walsh
#' (1998). `h2` is the candidates' own heritability, so `family_relationship` is
#' the **correlation of breeding values** between two members of a family,
#' \eqn{A_{ij} / \sqrt{A_{ii} A_{jj}}} (\eqn{A} the additive relationship
#' matrix), not \eqn{A_{ij}} alone; the index uses one such value for every
#' family. For families of non-inbred, unrelated parents it is 1/4 for half-sibs,
#' 1/2 for full-sibs, 1/2 for doubled haploids and 2/3 for S1 sibs (for which
#' \eqn{A_{ij} = 1}, \eqn{A_{ii} = 1.5}). Inbred or related parents change
#' these: this package's derivation from the tabular \eqn{A} gives
#' \eqn{2(1 + F) / (3 + F)} for S1 sibs and \eqn{(1 + F) / 2} for doubled
#' haploids of a parent with inbreeding \eqn{F}; for other families compute
#' \eqn{A_{ij} / \sqrt{A_{ii} A_{jj}}} from the pedigree. The index uses one `h2`
#' and one `family_relationship` for every family, so it assumes families of one
#' type from comparable parents; mixing family types (e.g. S1 with full-sib
#' families) is outside this model. The weights are optimal under an **additive**
#' model, in which family members covary only through breeding values
#' (\eqn{t = r h^2}). Dominance, epistasis or a shared family environment add to
#' the covariance of full-sib or selfed family members (e.g. \eqn{V_D / 4} for
#' full sibs), so with those in the phenotype the index is a reasonable but not an
#' optimal predictor. Scores are predictions from deviations of the
#' records from their mean (so families of different sizes are ranked on one
#' scale); a family of one is scored by its own record, `h2` times its deviation.
#' `"index"` is
#' the Smith--Hazel multi-trait economic index (`weights` = economic weights, one
#' per trait): \eqn{b = P^{-1} G a}, selecting on \eqn{b'y}, with \eqn{P} the
#' phenotypic covariance matrix and \eqn{G} the covariance matrix of the true
#' breeding values (`G = Cov(A)`). Hazel's (1943) \eqn{G} is the covariance of each
#' trait's phenotype with each trait's breeding value, \eqn{\mathrm{Cov}(P, A)};
#' the two coincide only when \eqn{\mathrm{Cov}(A, D) = 0} (random mating /
#' Hardy-Weinberg proportions, e.g. an F2 or a random-mated population). In selfed,
#' inbred or previously selected populations with dominance
#' \eqn{\mathrm{Cov}(P, A) \neq \mathrm{Cov}(A)} and the index is then Smith--Hazel
#' under that assumption, not the exact maximiser of the correlation with the
#' aggregate breeding value (DECISION-015). Phenotype rows are matched to
#' individuals by id, so the row order of `sim$pheno` is irrelevant.
#' `"quadratic_index"` is
#' the nonlinear genomic selection index of Ceron-Rojas et al. (2026),
#' \eqn{\hat I = w'y + y'Wy} (`weights` = linear `w`, `quad_weights` = the
#' symmetric quadratic/cross-product matrix `W`), which captures trait interactions
#' and intermediate optima. `"random"` draws at random (a drift control); its
#' `differential` and `intensity` attributes are the *realized* values on the
#' `on` criterion of the individuals drawn (expected 0, either sign), not those of
#' the random draw score, and `on` is validated.
#' `"culling"` is independent culling levels: `culling` gives one proportion per
#' trait (traits in `trait`, default the first `length(culling)`), and an individual
#' is kept iff it is in the top `culling[t]` fraction of every trait (simultaneous;
#' `sequential = TRUE` culls trait 1, then trait 2 among the survivors, and so on).
#' The number kept then follows from the proportions, so `n` / `prop` /
#' `intensity` must be left `NULL`; `direction` may be one value per trait, and `on`
#' may be an individuals x traits matrix of external per-trait predictions. Under
#' Hazel & Lush's (1942) idealized conditions -- `T` uncorrelated traits of equal
#' variance and heritability, equal weights -- this package's derivation gives the
#' per-generation aggregate responses \eqn{i(p)\sqrt{T} h\sigma_A} (index),
#' \eqn{T\,i(p^{1/T}) h\sigma_A} (simultaneous culling at \eqn{p^{1/T}} per trait) and
#' \eqn{i(p) h\sigma_A} (tandem: one trait per generation, via the scheme
#' wrappers' `trait` vector), so index \eqn{\ge} culling \eqn{\ge} tandem; at
#' `T = 2`, `p = 0.1` culling and tandem reach 0.907 and 0.707 of the index. Family
#' methods need a `family` grouping with no missing labels (`NA` is an error);
#' `"combined"` additionally needs `h2`. `"among_family"` keeps **whole** families
#' (highest family mean first) until at least `n` individuals are held, so it
#' returns *at least* `n` -- often more -- and `intensity` describes that whole
#' set; `"within_family"` allocates exactly `n` across families in proportion to
#' family size by largest remainder (remainder ties follow the character sort order
#' of the family labels), or, with `n_per_family`, keeps a stated number from every
#' family (see that argument; this is the form for unequal families). The count
#' from `prop` is `round(prop * N)` (R's
#' half-to-even rounding) with a minimum of 1. The reported `intensity` is
#' `S / sd(criterion)` with the sample (n - 1) standard deviation.
#'
#' @param sim a realized `phenotype_sim`, ideally built on a `Population` so the
#'   selected individuals can be crossed on.
#' @param n number of individuals to keep. Exactly one of `n`, `prop`, `intensity`.
#'   For `method = "among_family"` this is a floor: whole families are kept, so
#'   more than `n` individuals can be returned (see the Method section).
#' @param prop proportion of individuals to keep (0-1); the count is
#'   `round(prop * N)` (half-to-even) and at least 1.
#' @param intensity standardized selection intensity `i`; the number kept is the
#'   count whose asymptotic (infinite-population) truncation-selection intensity
#'   \eqn{i(p) = \phi(z_p)/p} is closest to `i`. The finite-sample expectation of
#'   the realized intensity is somewhat lower than \eqn{i(p)} (small `N`).
#' @param on selection criterion (see the Criterion section).
#' @param trait one trait index in `1..n_traits` to select on (default 1);
#'   not used for `method = "index"`/`"quadratic_index"` (all traits) or for a
#'   numeric/function `on`, but it is still validated in every non-culling method.
#'   A vector is an error here (it is only meaningful for `method = "culling"`, and
#'   as a tandem schedule in [pedigree()] / [recurrent_selection()]).
#' @param direction `"high"` (default) keeps the largest scores, `"low"` the
#'   smallest; one value (a vector is an error, except per trait for
#'   `method = "culling"`).
#' @param method selection method (see the Method section).
#' @param family optional grouping vector (length = individuals) for the family
#'   methods. Every individual needs a non-`NA`, non-empty label: `NA` and empty
#'   (`""`) labels are errors, because they would silently drop the individual or
#'   merge unrelated ones into one pseudo-family.
#' @param weights economic weights (one per trait) for `method = "index"`, or the
#'   linear weights `w` for `method = "quadratic_index"`.
#' @param quad_weights symmetric `n_traits x n_traits` matrix of quadratic
#'   (diagonal) and cross-product (off-diagonal) weights `W` for
#'   `method = "quadratic_index"` (default: zero, i.e. a linear index).
#' @param h2 narrow-sense heritability of the selection trait, required by
#'   `method = "combined"` to weight family versus individual information.
#' @param family_relationship correlation of breeding values among family members
#'   for `method = "combined"`. The default 0.25 is for half-sibs of non-inbred,
#'   unrelated parents; with such parents use 0.5 for full-sibs or doubled
#'   haploids and 2/3 for S1 sibs. Inbred or related parents change these -- see
#'   Details.
#' @param rep replication to select on when several were simulated (default 1).
#' @param culling per-trait kept proportions in (0, 1] for `method = "culling"`.
#' @param sequential for `method = "culling"`: `FALSE` (default) culls every trait
#'   on the whole population at once; `TRUE` culls the traits in order, each among
#'   the survivors of the previous ones.
#' @param n_per_family for `method = "within_family"` only: the number of
#'   individuals to keep **in each family**, in place of the proportional
#'   allocation of a total `n`. One whole number >= 1 keeps that many from every
#'   family; a numeric vector keeps a family-specific number (whole numbers >= 0),
#'   either **named by family label** (every family must be named, none twice or
#'   unknown) or unnamed with one value per family in the order of the sorted
#'   family labels (`sort(unique(as.character(family)))`, the order of
#'   `split()`). `n_per_family` replaces `n`, `prop` and `intensity`: giving
#'   either with it is an error, as is using it with another method. A family
#'   with fewer individuals than requested is an error naming the family (the
#'   count is never capped silently). Within each family the top-scoring
#'   individuals on the criterion are kept; the selection differential S and the
#'   intensity `i = S / sd(criterion)` are the same quantities as for the other
#'   methods (the mean of the kept scores minus the mean of all candidates, over
#'   the sample SD of all candidates). The default `NULL` keeps the proportional
#'   allocation. The scheme wrappers [pedigree()] and [recurrent_selection()] use
#'   mass selection only and do not forward it.
#' @return the selected individuals as a `Population` (when `sim` is
#'   Population-backed) or their ids, carrying attributes `selected` (ids),
#'   `differential` (selection differential S on the criterion), `intensity`
#'   (realized standardized i = S / sample SD of the criterion), `criterion` and
#'   `method` (plus `n_per_family`, the realized named count per family, when that
#'   argument was used). A criterion whose scores, differential or SD overflow to a
#'   non-finite value is an error.
#' @references
#' Truncation response and selection intensity: Falconer DS, Mackay TFC (1996)
#'   Introduction to Quantitative Genetics, 4th ed. Longman, Harlow; Lynch M,
#'   Walsh B (1998) Genetics and Analysis of Quantitative Traits. Sinauer
#'   Associates, Sunderland, Massachusetts.
#' Combined (family + individual) selection: Lush JL (1947) Family merit and
#'   individual merit as bases for selection. \emph{The American Naturalist}
#'   81:241--261 (Part I, \doi{10.1086/281520}) and 362--379 (Part II,
#'   \doi{10.1086/281532}).
#' Multi-trait selection index \eqn{b = P^{-1} G a}: Smith HF (1936) A discriminant
#'   function for plant selection. \emph{Annals of Eugenics} 7(3):240--250.
#'   \doi{10.1111/j.1469-1809.1936.tb02143.x}; Hazel LN (1943) The genetic basis
#'   for constructing selection indexes. \emph{Genetics} 28(6):476--490.
#'   \doi{10.1093/genetics/28.6.476}
#' Index vs independent culling vs tandem selection (`"culling"`, tandem in
#'   [pedigree()] / [recurrent_selection()]): Hazel LN, Lush JL (1942) The
#'   efficiency of three methods of selection. \emph{Journal of Heredity}
#'   33:393--399. \doi{10.1093/oxfordjournals.jhered.a105102}
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
                                  "random", "culling"),
                       family = NULL, weights = NULL, quad_weights = NULL,
                       h2 = NULL, family_relationship = 0.25, rep = 1L,
                       culling = NULL, sequential = FALSE, n_per_family = NULL) {
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
  method <- match.arg(method)
  if (!is.null(n_per_family)) {
    if (method != "within_family") {
      stop("`n_per_family` applies to method = \"within_family\" only.",
           call. = FALSE)
    }
    if (!is.null(n) || !is.null(prop) || !is.null(intensity)) {
      stop("`n_per_family` replaces `n`, `prop` and `intensity`; give only ",
           "`n_per_family`.", call. = FALSE)
    }
  }
  if (method == "culling") {
    return(.select_culling(sim, n, prop, intensity, on, trait,
                           if (missing(direction)) "high" else direction,
                           culling, sequential, rep))
  }
  if (!is.null(culling) || !identical(sequential, FALSE)) {
    stop("`culling` and `sequential` apply to method = \"culling\" only.",
         call. = FALSE)
  }
  # match.arg() silently returns "high" for the whole default vector, so an
  # explicit direction vector must be rejected before it is reduced.
  if (!missing(direction)) .check_direction(direction)
  direction <- match.arg(direction)
  # `trait` is one validated scalar index in every non-culling branch, including
  # those that ignore it (index/quadratic_index score on all traits; a numeric or
  # function `on` supplies its own values): a malformed `trait` is a caller error
  # that must not pass silently just because it happens to be unused.
  .check_trait_index(trait, sim$n_traits)
  ids <- sim$ids
  n_ind <- sim$n_ind

  keep_n <- if (is.null(n_per_family)) {
    .resolve_keep(n, prop, intensity, n_ind)
  } else NA_integer_                      # fixed by the per-family counts below

  fam <- if (method %in% c("within_family", "among_family", "combined")) {
    if (is.null(family) || length(family) != n_ind) {
      stop("method = \"", method, "\" needs a `family` grouping vector of ",
           "length ", n_ind, ".", call. = FALSE)
    }
    # split()/tapply() silently drop NA labels, so an individual with an NA family
    # could never be selected (and the combined index would crash): refuse.
    if (anyNA(family)) {
      stop("`family` has ", sum(is.na(family)), " NA label(s); every individual ",
           "needs a family (NA would silently exclude it from selection).",
           call. = FALSE)
    }
    # an empty (or blank) label is almost always a missing value written as "":
    # refuse it with the fix instead of treating it as a family
    fam_chr <- as.character(family)
    if (any(!nzchar(trimws(fam_chr)))) {
      stop("`family` has ", sum(!nzchar(trimws(fam_chr))), " empty label(s) (\"\"); ",
           "every individual needs a non-empty family label. Replace the empty ",
           "labels with the real family names (for example family[family == \"\"] <- ",
           "\"unknown\", or give each such individual its own label), or remove ",
           "those individuals before selecting.", call. = FALSE)
    }
    fam_chr
  } else NULL
  fam_alloc <- if (!is.null(n_per_family)) .resolve_n_per_family(n_per_family, fam)

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

  crit_real <- NULL     # criterion the statistics are reported on, when the
                        # ranking score is not it (method = "random")
  score <- if (method == "index") {
    .index_score(sim, weights, rep)
  } else if (method == "quadratic_index") {
    .quadratic_index_score(sim, weights, quad_weights, rep)
  } else if (method == "random") {
    # `on` is validated (and evaluated) so S and i can be reported honestly on the
    # phenotype/criterion of the individuals actually drawn; the draw itself uses
    # only the uniform score below.
    crit_real <- .criterion_values(sim, on, trait, rep)
    stats::setNames(stats::runif(n_ind), ids)
  } else if (method == "combined") {
    .combined_score(.criterion_values(sim, on, trait, rep), fam, h2,
                    family_relationship)
  } else {
    .criterion_values(sim, on, trait, rep)
  }
  if (!all(is.finite(score))) {
    stop("The selection score has ", sum(!is.finite(score)), " non-finite ",
         "value(s) (overflow or NA in the index/criterion inputs -- e.g. weights ",
         "of extreme magnitude); rescale the weights so every score is finite.",
         call. = FALSE)
  }
  if (direction == "low") score <- -score

  sel_idx <- switch(method,
    within_family = .sel_within_family(score, fam, keep_n, alloc = fam_alloc),
    among_family  = .sel_among_family(score, fam, keep_n),
    .sel_top(score, keep_n))               # mass, combined, index, random

  # Realized selection differential (S) and standardized intensity (i) are computed
  # on the criterion actually used to rank -- the index/quadratic score for those
  # methods, the `on` trait for mass/family selection -- so i = S / sd(criterion) is
  # the realized selection intensity by definition. `score` is already
  # direction-adjusted (selected individuals are its largest values), so S_dir >= 0;
  # i is reported as a positive magnitude and S in the natural (unflipped) sign. A
  # criterion with no spread (e.g. all ties) gives i = 0 rather than NaN.
  #
  # method = "random": the differential/intensity are the realized values of the
  # individuals drawn on the `on` criterion (natural sign; expected 0), never those
  # of the uniform draw score that ranked them.
  if (method == "random") {
    S_out <- mean(crit_real[sel_idx]) - mean(crit_real)
    sd_crit <- stats::sd(crit_real)
    S_dir <- S_out
  } else {
    S_dir <- mean(score[sel_idx]) - mean(score)
    sd_crit <- stats::sd(score)
    S_out <- if (direction == "low") -S_dir else S_dir
  }
  if (!is.finite(S_dir) || !is.finite(sd_crit)) {
    stop("The selection differential or the criterion's standard deviation ",
         "overflowed (non-finite); the criterion is on an extreme scale -- ",
         "rescale it before selecting.", call. = FALSE)
  }
  intensity_real <- if (sd_crit > 0) S_dir / sd_crit else 0

  out <- .selection_result(sim, sel_idx, ids)
  attr(out, "selected") <- ids[sel_idx]
  attr(out, "differential") <- S_out
  attr(out, "intensity") <- intensity_real
  attr(out, "criterion") <- if (is.function(on)) "custom" else
    if (is.numeric(on)) "custom" else on
  attr(out, "method") <- method
  if (!is.null(fam_alloc)) attr(out, "n_per_family") <- fam_alloc
  out
}

#' Resolve `n_per_family` to one kept count per family
#'
#' Returns a named integer vector in the order `split()` uses for `fam` (sorted
#' family labels). Errors on a malformed value, on names that do not match the
#' family labels, and on a family smaller than its requested count.
#' @keywords internal
#' @noRd
.resolve_n_per_family <- function(n_per_family, fam) {
  groups <- split(seq_along(fam), fam)
  lv <- names(groups)
  sizes <- lengths(groups)
  v <- n_per_family
  if (!is.numeric(v) || !length(v) || anyNA(v) || any(!is.finite(v)) ||
      any(v != floor(v)) || any(v < 0)) {
    stop("`n_per_family` must be a whole number >= 0 (or a vector of them), ",
         "with no NA.", call. = FALSE)
  }
  nm <- names(v)
  if (length(v) == 1L && is.null(nm)) {
    v <- stats::setNames(rep(as.numeric(v), length(lv)), lv)
  } else if (!is.null(nm)) {
    if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
      stop("`n_per_family` names must be non-empty and unique family labels.",
           call. = FALSE)
    }
    unknown <- setdiff(nm, lv)
    missing_f <- setdiff(lv, nm)
    if (length(unknown) || length(missing_f)) {
      stop("`n_per_family` names must match the family labels exactly",
           if (length(unknown)) paste0("; unknown: ", paste(unknown, collapse = ", ")),
           if (length(missing_f)) paste0("; missing: ",
                                         paste(missing_f, collapse = ", ")),
           ".", call. = FALSE)
    }
    v <- v[lv]
  } else {
    if (length(v) != length(lv)) {
      stop("An unnamed `n_per_family` needs one value (same count everywhere) ",
           "or one value per family (", length(lv), ", in sorted family-label ",
           "order); got ", length(v), ". Name the vector by family label to ",
           "avoid depending on the order.", call. = FALSE)
    }
    names(v) <- lv
  }
  v <- v[lv]
  if (sum(v) < 1) {
    stop("`n_per_family` keeps no individual at all; ask for at least one.",
         call. = FALSE)
  }
  short <- which(v > sizes)
  if (length(short)) {
    stop("`n_per_family` asks for more individuals than the family holds in: ",
         paste0(lv[short], " (wants ", v[short], ", has ", sizes[short], ")",
                collapse = "; "),
         ". Lower the count for these families; it is not capped silently.",
         call. = FALSE)
  }
  stats::setNames(as.integer(v), lv)
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

#' Validate the single trait index a criterion is evaluated on
#'
#' A vector would be silently recycled inside `ph$trait == paste0("Trait_",
#' trait)` (an alternating-trait criterion), and an out-of-range index would only
#' surface later as a misleading non-finite-criterion error.
#' @keywords internal
#' @noRd
.check_trait_index <- function(trait, n_traits) {
  if (!is.numeric(trait) || length(trait) != 1L || !is.finite(trait) ||
      trait != floor(trait) || trait < 1 || trait > n_traits) {
    stop("`trait` must be one trait index in 1..", n_traits, "; got ",
         paste(format(trait), collapse = ", "), ". select_ind() ranks on a ",
         "single trait (a vector `trait` is only defined for method = ",
         "\"culling\", or as a tandem schedule in pedigree() / ",
         "recurrent_selection()).", call. = FALSE)
  }
  invisible(trait)
}

#' Validate a single selection direction
#' @keywords internal
#' @noRd
.check_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
      is.na(direction) || !direction %in% c("high", "low")) {
    stop("`direction` must be one of \"high\" or \"low\" (a single value; got ",
         paste(format(direction), collapse = ", "), ").", call. = FALSE)
  }
  invisible(direction)
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
    .check_trait_index(trait, sim$n_traits)
    v <- .breeding_value_matrix(sim, rep)[, trait]
  } else if (identical(on, "gv")) {
    .check_trait_index(trait, sim$n_traits)
    v <- .genetic_matrix(sim, rep)[, trait]
  } else if (identical(on, "pheno")) {
    .check_trait_index(trait, sim$n_traits)
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
  # epistatic models are refused), not the total genotypic covariance. This is
  # Hazel's g_ij = Cov(phenotype_i, breeding value_j) only when Cov(A, D) = 0
  # (random mating / Hardy-Weinberg): in selfed or inbred populations with
  # dominance Cov(P, A) != Cov(A). Kept as Cov(A) by decision (DECISION-015).
  G <- .breeding_value_matrix(sim, rep)                # ind x trait, sim$ids order
  # Phenotypes are joined to individuals by NAME (as .criterion_values() does), not
  # by row position: sim$pheno is a user-visible data frame whose row order carries
  # no guarantee, and G above is in sim$ids order.
  P <- do.call(cbind, lapply(seq_len(nt), function(k) {
    .criterion_values(sim, "pheno", k, rep)
  }))
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
  if (any(!is.finite(P)) || any(!is.finite(r))) {
    stop("method = \"index\": the phenotypic covariance or the genetic ",
         "covariance-weight product has non-finite (NA/NaN/Inf) entries, so the ",
         "index weights are undefined (constant or missing phenotypes?).",
         call. = FALSE)
  }
  # Scale-invariant solve: b = P^{-1} r is unchanged when P and r are divided by
  # any positive constants and the quotient of the constants is restored at the
  # end, so work with P / mP and r / mr (both O(1)) and rescale last. Solving
  # directly overflowed 1 / d (weights NaN) for a finite covariance at denormal
  # scale; the weights are an error only if they are genuinely not representable.
  # Powers of two, so dividing and restoring are exact (results for ordinary
  # scales are bit-identical to the unscaled solve).
  p2 <- function(x) if (x > 0) 2^floor(log2(x)) else 0
  mP <- p2(max(abs(P)))
  mr <- p2(max(abs(r)))
  Pn <- if (mP > 0) P / mP else P
  rn <- if (mr > 0) r / mr else r
  s <- svd(Pn)
  tol <- max(dim(Pn)) * .Machine$double.eps * (if (length(s$d)) s$d[1L] else 0)
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
  b <- Pinv %*% rn
  if (mP > 0 && mr > 0) {
    ratio <- mr / mP
    b <- if (is.finite(ratio) && ratio > 0) b * ratio else (b / mP) * mr
  }
  if (any(!is.finite(b))) {
    stop("method = \"index\": the index weights overflow (non-finite): the ",
         "phenotypic covariance is on an extreme (denormal-scale) magnitude ",
         "relative to the genetic covariance. Rescale the traits.",
         call. = FALSE)
  }
  b
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
#' Both weights are non-negative; the family mean drops out as h2 -> 1 when r < 1
#' (at r = 1 the family mean is the breeding value itself and keeps its weight).
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
  t <- r * h2
  # The index predicts A - E[A] from deviations of the observations from their
  # means (Hazel 1943). With unequal family sizes the weights differ between
  # families, so raw values would make the ranking depend on the arbitrary
  # phenotype origin; centre on the candidate mean.
  mu <- mean(values)
  dev <- values - mu
  # Group by integer code, not by name (a robust grouping; `select_ind()` already
  # rejects NA and empty labels before this point).
  fam_id   <- match(fam, unique(fam))
  fam_mean <- as.numeric(tapply(dev, fam_id, mean))
  fam_n    <- tabulate(fam_id)
  score <- numeric(length(values))
  for (i in seq_along(values)) {
    f <- fam_id[i]
    n <- fam_n[f]
    # 1 - t written as (1 - r) + r (1 - h2): exact near t = 1 (both differences
    # are exact there), so the weights below stay accurate and bounded
    # (b1 <= h2, and b2 finite) arbitrarily close to it. Only t = 1 exactly
    # (r = h2 = 1: the own record and the family mean both equal the breeding
    # value) is singular; there, and for a family of one, use the own record,
    # b = Cov(A, P) / Var(P) = h2.
    one_minus_t <- (1 - r) + r * (1 - h2)
    if (n == 1L || one_minus_t <= 0) {
      score[i] <- h2 * dev[i]
      next
    }
    # b = V^{-1} c for (own record, family mean) with Var(P) = vP, t = r*h2,
    # Cov(own, family mean) = Var(family mean) = vP (1 + (n-1) t) / n,
    # Cov(A, own) = h2 vP and Cov(A, family mean) = h2 vP (1 + (n-1) r) / n.
    # vP cancels, leaving dimensionless weights -- computed directly so that
    # rescaling the records only rescales the scores (products of vP^2 would
    # under/overflow at extreme scales):
    b1 <- h2 * (1 - r) / one_minus_t
    b2 <- h2 * n * r * (1 - h2) / ((1 + (n - 1) * t) * one_minus_t)
    score[i] <- b1 * dev[i] + b2 * fam_mean[f]
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
.sel_within_family <- function(score, fam, keep_n, alloc = NULL) {
  groups <- split(seq_along(score), fam)
  sizes  <- lengths(groups)
  ntot   <- length(score)
  if (!is.null(alloc)) {
    # explicit per-family counts (n_per_family), already validated against the
    # family sizes and given in split() order
    out <- integer(0)
    for (g in seq_along(groups)) {
      k <- alloc[[g]]
      if (k > 0L) {
        ix <- groups[[g]]
        out <- c(out, ix[order(score[ix], decreasing = TRUE)[seq_len(k)]])
      }
    }
    return(out)
  }
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
#'
#' Whole families are taken (highest mean first) until at least `keep_n`
#' individuals are held, so `keep_n` is a floor: more can be returned.
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

#' Independent culling levels (method = "culling")
#' @keywords internal
#' @noRd
.select_culling <- function(sim, n, prop, intensity, on, trait, direction,
                            culling, sequential, rep) {
  if (!is.null(n) || !is.null(prop) || !is.null(intensity)) {
    stop("method = \"culling\": the number kept follows from `culling`; leave ",
         "`n`, `prop` and `intensity` NULL.", call. = FALSE)
  }
  if (is.null(culling) || !is.numeric(culling) || !length(culling) ||
      any(!is.finite(culling)) || any(culling <= 0 | culling > 1)) {
    stop("method = \"culling\" needs `culling`: one kept proportion in (0, 1] per ",
         "trait.", call. = FALSE)
  }
  .validate_flag(sequential, "sequential")
  nt <- length(culling)
  traits <- if (length(trait) == nt) {
    trait
  } else if (length(trait) == 1L && identical(as.numeric(trait), 1)) {
    seq_len(nt)
  } else {
    stop("`trait` must list one trait per `culling` proportion (", nt, ").",
         call. = FALSE)
  }
  if (!is.numeric(traits) || any(!is.finite(traits)) ||
      any(traits != floor(traits)) || any(traits < 1) ||
      any(traits > sim$n_traits) || anyDuplicated(traits)) {
    stop("`trait` must be distinct trait indices in 1..", sim$n_traits, ".",
         call. = FALSE)
  }
  if (!is.character(direction) || !length(direction) ||
      !all(direction %in% c("high", "low")) ||
      !(length(direction) %in% c(1L, nt))) {
    stop("`direction` must be \"high\" or \"low\", one value or one per trait.",
         call. = FALSE)
  }
  direction <- rep_len(direction, nt)
  ids <- sim$ids
  n_ind <- sim$n_ind
  crit <- if (is.matrix(on)) {
    if (!is.numeric(on) || nrow(on) != n_ind || ncol(on) != nt) {
      stop("A matrix `on` for culling must be numeric, individuals x traits (",
           n_ind, " x ", nt, ").", call. = FALSE)
    }
    if (!is.null(rownames(on))) {
      if (!all(ids %in% rownames(on))) {
        stop("The row names of the `on` matrix must include every individual.",
             call. = FALSE)
      }
      on <- on[ids, , drop = FALSE]
    }
    if (any(!is.finite(on))) {
      stop("The `on` matrix has non-finite values.", call. = FALSE)
    }
    lapply(seq_len(nt), function(k) as.numeric(on[, k]))
  } else {
    if (nt > 1L && (is.numeric(on) || is.function(on))) {
      stop("With several culling traits, pass external predictions as an ",
           "individuals x traits matrix `on`.", call. = FALSE)
    }
    lapply(seq_len(nt), function(k) .criterion_values(sim, on, traits[k], rep))
  }
  score <- lapply(seq_len(nt), function(k) {
    if (direction[k] == "low") -crit[[k]] else crit[[k]]
  })
  kept_at <- integer(nt)
  if (sequential) {
    alive <- seq_len(n_ind)
    for (k in seq_len(nt)) {
      keep <- max(1L, round(culling[k] * length(alive)))
      alive <- alive[order(score[[k]][alive], decreasing = TRUE)[seq_len(keep)]]
      kept_at[k] <- length(alive)
    }
    sel_idx <- sort(alive)
  } else {
    tops <- lapply(seq_len(nt), function(k) {
      .sel_top(score[[k]], max(1L, round(culling[k] * n_ind)))
    })
    sel_idx <- sort(Reduce(intersect, tops))
    kept_at <- vapply(tops, length, integer(1))
  }
  if (!length(sel_idx)) {
    stop("method = \"culling\": no individual passes every culling level; ",
         "raise the `culling` proportions.", call. = FALSE)
  }
  diff <- vapply(seq_len(nt), function(k) {
    mean(crit[[k]][sel_idx]) - mean(crit[[k]])
  }, numeric(1))
  sds <- vapply(crit, stats::sd, numeric(1))
  intens <- ifelse(is.finite(sds) & sds > 0, abs(diff) / sds, 0)
  names(diff) <- names(intens) <- paste0("Trait_", traits)
  out <- .selection_result(sim, sel_idx, ids)
  attr(out, "selected") <- ids[sel_idx]
  attr(out, "differential") <- diff
  attr(out, "intensity") <- intens
  attr(out, "criterion") <- if (is.matrix(on) || is.numeric(on) ||
                                is.function(on)) "custom" else on
  attr(out, "method") <- "culling"
  attr(out, "culling") <- list(proportion = culling, traits = traits,
                               sequential = sequential, kept = kept_at)
  out
}
