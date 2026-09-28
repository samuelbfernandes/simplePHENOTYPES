# Marker-based selection (DECISION-029): MAS foreground filtering, staged gene
# pyramiding, and ranking on a marker index (MARS) or any other score.

#' Marker-assisted selection and gene pyramiding
#'
#' Keeps the candidates carrying the favourable allele at target markers --
#' marker-assisted selection (MAS) of major genes -- and ranks the feasible ones.
#' The favourable allele is given directly (`favorable`), so no founders are
#' needed; for introgression from a donor with recovery of the recurrent genome,
#' use [mabc_select()], which this generalizes for the foreground step.
#'
#' A candidate meets the requirement at a marker when it is a `"carrier"` (at
#' least one favourable allele) or a `"homozygote"` for it. It is feasible when it
#' meets the requirement at `min_markers` or more of the targets: the default,
#' all of them, is a full pyramid; fewer is a **staged** pyramid (e.g. 3 of 4
#' genes this cycle), reported as such. Feasible candidates are ranked on
#' `rank_on` (by default the number of targets met), then by a seeded random
#' tie-break.
#'
#' `rank_on` may be any score: a phenotype, predicted breeding values, or a marker
#' index. [additive_value()] with marker names and weights *is* the marker index
#' \eqn{\sum_j w_j x_{ij}} of marker-assisted recurrent selection (MARS). With the
#' simulation's own causal loci weighted by their **average effects** it is an
#' **oracle** index of transmissible merit -- the upper bound of an index on
#' estimated additive marker effects (Lande & Thompson 1990), not a realistic
#' one. For an additive architecture those weights are the additive effects `a`
#' of [template_effects()]; with dominance they are
#' \eqn{\alpha = a + d(q - p)} at the candidates' allele frequencies, not the raw
#' `a` (which ignores dominance's contribution to the breeding value; \eqn{\alpha}
#' is the average effect of Falconer & Mackay 1996). Ranked on such an oracle index,
#' `marker_select()` keeps the individuals `select_ind(on = "bv")` keeps, except
#' among tied scores at the cut-off: `marker_select()` breaks ties at random (see
#' `seed`), `select_ind()` by input order.
#'
#' The favourable allele is identified by the `pop`'s -1/0/1 coding; as for
#' [mabc_select()], that coding cannot be verified from dosages, so `favorable`
#' must refer to the same coding as `pop`.
#'
#' @param pop the candidate `Population`.
#' @param markers target markers, by name or map index.
#' @param favorable which homozygote carries the favourable allele in `pop`'s
#'   coding: `1` or `-1`, one value or one per marker.
#' @param requirement `"carrier"` or `"homozygote"`, one value or one per marker.
#' @param min_markers feasible candidates meet the requirement at this many
#'   targets or more (default all; `0` makes every candidate feasible, e.g. for a
#'   pure index ranking).
#' @param n,prop optional number, or proportion of the candidates, to keep;
#'   neither keeps every feasible candidate. An error if fewer are feasible.
#' @param rank_on optional ranking score among feasible candidates: a numeric
#'   vector (named by id, or in id order) or a function of `pop` returning one.
#' @param direction rank `rank_on` from `"high"` (default) or `"low"`.
#' @param seed optional seed for the tie-break (the ambient RNG state is
#'   restored).
#' @return The selected individuals as a `Population`, with attribute
#'   `"marker_select"`: a data frame over all candidates with `id`, `state`
#'   (per-target genotype, `F` / `H` / `U` for favourable homozygote /
#'   heterozygote / unfavourable homozygote, `;`-separated), `n_met`, `feasible`,
#'   `score`, `tiebreak`, `rank` (among feasible) and `selected`; and attribute
#'   `"marker_select_info"` with the settings, including whether the pyramid was
#'   staged.
#' @references
#' Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- the average effect \eqn{\alpha = a + d(q - p)}.
#'
#' Lande R, Thompson R (1990) Efficiency of marker-assisted selection in the
#'   improvement of quantitative traits. \emph{Genetics} 124:743--756.
#'   \doi{10.1093/genetics/124.3.743}
#' @seealso [mabc_select()], [additive_value()], [select_ind()]
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:4), allele = "A/G", chr = 1:4, pos = 1,
#'                 cm = 0, P1 = 1L, P2 = -1L)
#' pop <- as_population(g)
#' f2 <- selfcross(cross(pop[1], pop[2], seed = 1), n = 200, seed = 2)
#' pyr <- marker_select(f2, markers = c("m1", "m2"), requirement = "homozygote")
#' nrow(attr(pyr, "marker_select")[attr(pyr, "marker_select")$feasible, ])  # ~ 200/16
marker_select <- function(pop, markers, favorable = 1L,
                          requirement = c("carrier", "homozygote"),
                          min_markers = NULL, n = NULL, prop = NULL,
                          rank_on = NULL, direction = c("high", "low"),
                          seed = NULL) {
  .check_population(pop)
  direction <- match.arg(direction)
  idx <- .mabc_markers(pop$map, markers, "markers")
  k <- length(idx)
  if (!k) {
    stop("marker_select(): `markers` must name at least one marker.",
         call. = FALSE)
  }
  if (!is.numeric(favorable) || !(length(favorable) %in% c(1L, k)) ||
      !all(favorable %in% c(-1, 1))) {
    stop("marker_select(): `favorable` must be 1 or -1, one value or one per ",
         "marker.", call. = FALSE)
  }
  favorable <- rep_len(as.integer(favorable), k)
  if (missing(requirement)) requirement <- "carrier"
  if (!is.character(requirement) || !(length(requirement) %in% c(1L, k)) ||
      !all(requirement %in% c("carrier", "homozygote"))) {
    stop("marker_select(): `requirement` must be \"carrier\" or \"homozygote\", ",
         "one value or one per marker.", call. = FALSE)
  }
  requirement <- rep_len(requirement, k)
  if (is.null(min_markers)) min_markers <- k
  if (!is.numeric(min_markers) || length(min_markers) != 1L ||
      min_markers != floor(min_markers) || min_markers < 0 || min_markers > k) {
    stop("marker_select(): `min_markers` must be a whole number in 0..", k, ".",
         call. = FALSE)
  }
  if (!is.null(n) && !is.null(prop)) {
    stop("marker_select(): give at most one of `n` / `prop`.", call. = FALSE)
  }

  ids <- pop$ids
  N <- length(ids)
  z <- dosages(pop)[idx, , drop = FALSE] * favorable     # +1 favourable homozygote
  met <- ifelse(requirement == "homozygote", 1, 0)
  ok <- z >= met                                          # carrier: z >= 0
  n_met <- colSums(ok)
  feasible <- n_met >= min_markers
  state <- apply(matrix(c("U", "H", "F")[z + 2L], nrow = k), 2L, paste,
                 collapse = ";")
  score <- if (is.null(rank_on)) {
    as.numeric(n_met)
  } else {
    v <- if (is.function(rank_on)) rank_on(pop) else rank_on
    if (!is.numeric(v)) {
      stop("marker_select(): `rank_on` must be numeric (or a function returning ",
           "a numeric score).", call. = FALSE)
    }
    if (!is.null(names(v))) {
      if (!all(ids %in% names(v))) {
        stop("marker_select(): a named `rank_on` must score every candidate.",
             call. = FALSE)
      }
      v <- v[ids]
    }
    if (length(v) != N || any(!is.finite(v))) {
      stop("marker_select(): `rank_on` must give one finite score per ",
           "candidate (", N, ").", call. = FALSE)
    }
    as.numeric(v)
  }
  key <- if (direction == "low") -score else score
  n_feasible <- sum(feasible)
  keep <- if (!is.null(n)) {
    .validate_count(n, "n", minimum = 1L)
  } else if (!is.null(prop)) {
    if (!is.numeric(prop) || length(prop) != 1L || !is.finite(prop) ||
        prop <= 0 || prop > 1) {
      stop("marker_select(): `prop` must be one value in (0, 1].", call. = FALSE)
    }
    max(1L, round(prop * N))
  } else {
    n_feasible
  }
  if (n_feasible < keep || n_feasible == 0L) {
    stop("marker_select(): only ", n_feasible, " of ", N, " candidates meet the ",
         "requirement at ", min_markers, " of ", k, " target(s); cannot select ",
         max(keep, 1L), ".", call. = FALSE)
  }
  old <- .Random.seed_safe()
  if (!is.null(seed)) {
    set.seed(.validate_seed(seed))
    on.exit(.restore_seed(old))
  }
  tiebreak <- sample.int(N)
  ord <- order(!feasible, -key, tiebreak)
  rank <- rep(NA_integer_, N)
  rank[ord[seq_len(n_feasible)]] <- seq_len(n_feasible)
  chosen <- ord[seq_len(keep)]
  out <- pop[chosen]
  attr(out, "marker_select") <- data.frame(
    id = ids, state = state, n_met = as.integer(n_met), feasible = feasible,
    score = score, tiebreak = tiebreak, rank = rank,
    selected = seq_len(N) %in% chosen, stringsAsFactors = FALSE,
    row.names = NULL)
  attr(out, "marker_select_info") <- list(
    markers = pop$map$snp[idx], favorable = favorable, requirement = requirement,
    min_markers = as.integer(min_markers), staged = min_markers < k,
    ranked_on = if (is.null(rank_on)) "n_met" else "rank_on")
  out
}
