#' Add an additive variance component
#'
#' @param sim a `phenotype_sim`.
#' @param prop proportion of total phenotypic variance (scalar or length
#'   `n_traits`). May be omitted when `h2` was set in [simulate_phenotype()],
#'   in which case the layer takes whatever of that budget is unspent -- so a
#'   single `additive()` layer gets `prop = h2`. When `prop` is given and `h2`
#'   is set, the layer `prop` values must sum to `h2`.
#' @param n_qtn number of additive QTNs; overrides the baseline `n_qtn` from
#'   [simulate_phenotype()] (with a warning when both are given).
#' @param qtn optional user-supplied QTNs for this layer, so you can fix the
#'   causal loci for one effect type while the others are drawn at random. Give
#'   marker names (matched against the map) or column indices; a vector is used
#'   for every trait, or a length-`n_traits` list sets each trait's loci
#'   separately. `n_qtn` is then taken from what you supply. Fixed loci are not
#'   redrawn by `vary_qtn`.
#' @param effect optional geometric base (scalar) or explicit effect series
#'   (length `n_qtn`).
#' @param orthogonal use the orthogonal genotypic model instead of the
#'   variance-partition coding (default `FALSE`). When `TRUE`, give per-locus
#'   additive effects `a` and dominance deviations `d`; the additive and
#'   dominance components are orthogonal under random mating (see Details).
#' @param a additive effect per locus for `orthogonal = TRUE`: a scalar geometric
#'   base or a length-`n_qtn` vector (the counterpart of `effect`). The locus
#'   value is `-a` / `+d` / `+a` for gene content 0 / 1 / 2.
#' @param d dominance deviation per locus for `orthogonal = TRUE`: a single value
#'   (same at every locus) or a length-`n_qtn` vector. `d = 0` (default) is a
#'   purely additive locus. The heterozygote value is `d`, and the degree of
#'   dominance is `d / abs(a)` (`abs` because `a` may be negative, directly or
#'   through `phase = "repulsion"`): magnitude 0 = additive, 1 = complete
#'   dominance, > 1 = overdominance. Dominance is toward the `+a` homozygote when
#'   `d` and `a` share a sign and toward the `-a` homozygote when they differ --
#'   it is the ratio `d/a`, not the sign of `d` alone, that sets the direction.
#' @param phase QTN effect-sign pattern, `"coupling"` (default) or
#'   `"repulsion"`. Under `"repulsion"` the effect **signs alternate across the
#'   layer's QTNs in draw order** (`+, -, +, -, ...`); coupling leaves the drawn
#'   signs untouched. The realized genetic variance is still `prop` (variance is
#'   pinned by the partition, not by phase); what changes is the sign structure
#'   of the per-QTN effects. This is a convenient stand-in for a repulsion
#'   architecture, **not** a phase derived from the haplotypes: it alternates by
#'   position in the draw, and does not inspect the sign of LD between the loci.
#'   It therefore induces genuine repulsion-style cancellation only when
#'   consecutively drawn QTNs happen to be in positive LD (same chromosome,
#'   coupling-phase alleles); for effectively unlinked QTNs it merely relabels
#'   which allele is "increasing" and leaves the cross-locus covariance
#'   unchanged. `phase` is architecture-independent -- it acts on whichever QTNs
#'   the additive layer draws.
#' @param dist within-layer effect distribution (default "geometric").
#' @return the updated `phenotype_sim`.
#' @details
#' The additive value is the -1/0/1 dosage weighted by the QTN effects, centered,
#' and scaled to `prop` of the phenotypic variance. This is a simulation
#' convention, not Fisher's average-effect decomposition; for an additive-only
#' model `prop` equals the narrow-sense h2 under Hardy-Weinberg.
#'
#' @section Orthogonal genotypic model (`orthogonal = TRUE`):
#' Instead of the dosage coding above, the layer builds each locus's genotypic
#' value from an additive effect `a` and a dominance deviation `d`
#' (value `-a` / `+d` / `+a` for gene content 0 / 1 / 2), and the whole genotypic
#' value is scaled to `prop` of the phenotypic variance. Its additive and
#' dominance **variances then emerge** from `a`, `d` and the allele frequencies
#' rather than from separate layer proportions: the additive (breeding-value)
#' component uses the average effect of substitution
#' \eqn{\alpha_j = a_j + d_j(1 - 2 p_j)}, and the dominance deviation is the
#' realized residual \eqn{D = g - A}. The additive and dominance components are
#' orthogonal (\eqn{Cov(A, D) = 0}) in expectation under **random mating**, which
#' puts each locus in Hardy-Weinberg proportions; \eqn{\alpha} above is the
#' one-generation transmitting-ability form. This does *not* require linkage
#' equilibrium -- between-locus LD is compatible with orthogonality in this
#' no-epistasis model -- but per-locus HWE alone is not sufficient under arbitrary
#' nonrandom multilocus genotype association (e.g. selection or population
#' structure), which can correlate one locus's gene content with another's
#' heterozygosity. A finite or structured sample therefore carries a (usually
#' small) \eqn{Cov(A, D)} that breeders customarily ignore; this implementation
#' does not, and the variance budget reports the **realized** shares
#' \eqn{Var(A)/Var(g)} and \eqn{Var(D)/Var(g)} on separate `additive` and
#' `dominance` rows, plus an `add_dom_cov` row for \eqn{2\,Cov(A, D)/Var(g)}, so
#' the three sum to the layer `prop` (the `add_dom_cov` row is \eqn{\approx 0} for
#' a large random-mating or F2 sample and non-zero otherwise). The realized
#' narrow-sense heritability follows from the
#' degree of dominance `d / abs(a)` (magnitude 1 = complete dominance, > 1 =
#' overdominance; `abs` because `a` may be negative -- which the variance-partition
#' grammar cannot express, because per-component scaling washes a global degree
#' out), and the breeding value used by [select_ind()] / [optimum_contribution()]
#' picks up the dominance-induced average effect automatically. Because the
#' effects are fixed, this mode does not draw a separate `dominance()` layer, is
#' incompatible with `vary_qtn`, and (like `qtn =`) cannot be used with the
#' correlation-controlling `"pleiotropy"` (multi-trait) or `"ld"` architectures.
#' @references
#'   Fisher RA (1918) The correlation between relatives on the supposition of
#'   Mendelian inheritance. \emph{Trans. R. Soc. Edinb.} 52:399-433;
#'   Falconer DS (1985) A note on Fisher's 'average effect' and 'average excess'.
#'   \emph{Genet. Res.} 46(3):337-347 (the average effect of substitution
#'   \eqn{\alpha = a + d(q - p)} and its dependence on random-mating /
#'   Hardy-Weinberg proportions); Falconer DS, Mackay TFC (1996)
#'   \emph{Introduction to Quantitative Genetics}, 4th ed. Longman, Harlow;
#'   Lynch M, Walsh B (1998) \emph{Genetics and Analysis of Quantitative Traits}.
#'   Sinauer, Sunderland, MA (additive / dominance decomposition).
#' @seealso [dominance()], [epistasis()], [vqtl()] for the other layers.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # Additive trait: 3 QTNs explaining half the phenotypic variance,
#' # so h2 = 0.5 and the remaining 0.5 is residual.
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3)
#' ph
#'
#' head(phenotypes_long(ph))
#'
#' # Effect sizes follow a geometric series by default. `effect` sets its base,
#' # so QTN effects here are 0.5, 0.25, 0.125, ...
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3, effect = 0.5)
#'
#' # Or give the effects explicitly, one per QTN.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.5, n_qtn = 3, effect = c(0.6, 0.3, 0.1))
#'
#' # Fix the additive QTNs yourself by marker name.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.5, qtn = c("ss196442916", "ss196439337", "ss196480535"))
#'
#' # Orthogonal genotypic model: give per-locus a and d; the additive and
#' # dominance variances emerge (see the var_budget). Dominance needs
#' # heterozygotes, so simulate on a segregating F2.
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' f2  <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 40, seed = 2)
#' og  <- simulate_phenotype(f2, h2 = 0.6, seed = 3) |>
#'   additive(orthogonal = TRUE, a = 0.5, d = 0.5, n_qtn = 8)
#' og$var_budget
additive <- function(sim, prop = NULL, n_qtn = NULL, qtn = NULL, effect = NULL,
                     phase = c("coupling", "repulsion"), dist = "geometric",
                     orthogonal = FALSE, a = NULL, d = NULL) {
  .check_sim(sim)
  .require_markers(sim, "additive")
  .validate_flag(orthogonal, "orthogonal")
  phase <- match.arg(phase)
  if (orthogonal) {
    if (!is.null(effect)) {
      stop("additive(orthogonal = TRUE): give additive effects via `a` (and ",
           "dominance deviations via `d`), not `effect`.", call. = FALSE)
    }
    if ((sim$architecture == "pleiotropy" && sim$n_traits > 1) ||
        sim$architecture == "ld") {
      stop("additive(orthogonal = TRUE) sets per-locus effects, which ",
           "architecture = \"", sim$architecture, "\" cannot honour (it draws ",
           "its own effects to control the cross-trait correlation). Use ",
           "architecture = \"independent\".", call. = FALSE)
    }
    if (isTRUE(sim$vary_qtn)) {
      stop("additive(orthogonal = TRUE) is not supported with vary_qtn: the ",
           "genotypic model fixes per-locus a and d across replications.",
           call. = FALSE)
    }
    if (!is.null(.last_layer_of_type(sim, "dominance"))) {
      stop("additive(orthogonal = TRUE) already models dominance through `d`; ",
           "remove the separate dominance() layer.", call. = FALSE)
    }
  } else if (!is.null(a) || !is.null(d)) {
    stop("additive(): `a` and `d` configure the orthogonal genotypic model; set ",
         "orthogonal = TRUE to use them (or use `effect` for the standard layer).",
         call. = FALSE)
  }
  prop <- .resolve_prop(sim, prop, "additive")
  user_qtn <- .resolve_qtn_arg(sim, qtn, "additive")
  # Fixing the additive loci is incompatible with the architectures that draw
  # their own loci to build a controlled cross-trait correlation. Under
  # "pleiotropy" the shared/specific partition and the multivariate (PleioArch)
  # effect draw -- with its PSD feasibility check and MAF scaling -- would be
  # bypassed, giving identical per-trait effects and a realized correlation of
  # +1 regardless of `cor`; under "ld" the same locus would become causal for
  # both traits instead of a distinct linked pair. Reject rather than silently
  # mis-simulate (SPEC 4.1; DECISION-007/013/014).
  if (!is.null(user_qtn) &&
      ((sim$architecture == "pleiotropy" && sim$n_traits > 1) ||
       sim$architecture == "ld")) {
    stop("additive(qtn=) fixes the causal loci, which architecture = \"",
         sim$architecture, "\" cannot honour: it draws its own loci (shared + ",
         "trait-specific for \"pleiotropy\"; distinct linked pairs for \"ld\") ",
         "and controls the cross-trait correlation through them. Omit `qtn`, or ",
         "use architecture = \"independent\" to fix loci.", call. = FALSE)
  }
  nq <- if (!is.null(user_qtn)) length(user_qtn[[1]]) else
    .resolve_n_qtn(sim, n_qtn, "additive")
  occ <- .type_occurrence(sim, "additive")

  # In the orthogonal genotypic model the additive effects come from `a`; the
  # standard variance-partition layer uses `effect`.
  eff_arg <- if (orthogonal) a else effect
  d_series <- if (orthogonal) .orthogonal_d_series(d, nq) else NULL

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn)) {
      q <- user_qtn
      e <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, eff_arg))
    } else if (!orthogonal && sim$architecture == "pleiotropy" &&
               sim$n_traits > 1) {
      if (!is.null(effect) || !identical(dist, "geometric")) {
        stop("additive(): `effect` and non-default `dist` cannot be used under ",
             "architecture = \"pleiotropy\" because effects come from the ",
             "multivariate draw that controls `cor`.", call. = FALSE)
      }
      pd <- .pleio_draw(sim, nq, .expand_prop(prop, sim$n_traits), rep_seed)
      q <- pd$qtn
      e <- pd$effect
    } else {
      q <- .draw_qtn(sim, nq, rep_seed)
      e <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, eff_arg))
    }
    list(qtn = q, effect = .apply_phase(e, phase))
  }

  drawn <- .draw_layer(sim, "additive", occ, build, fixed = !is.null(user_qtn))
  layer <- list(type = "additive", prop = prop, n_qtn = nq, dist = dist,
                phase = phase, qtn = drawn$qtn, effect = drawn$effect)
  if (!is.null(drawn$qtn_reps)) {
    layer$qtn_reps <- drawn$qtn_reps
    layer$effect_reps <- drawn$effect_reps
  }
  layer$ld <- attr(drawn$qtn, "ld")
  if (orthogonal) {
    layer$orthogonal <- TRUE
    layer$d_effect <- lapply(seq_len(sim$n_traits), function(t) d_series)
    # A nonzero dominance deviation at a locus with no heterozygotes is silently
    # inert (its het indicator is all zero), so require a heterozygote at *each*
    # locus whose d != 0 -- checked per locus, not collectively over the set.
    if (.orthogonal_hetless_d(sim, layer$qtn, layer$d_effect)) {
      stop("additive(orthogonal = TRUE): a locus with a non-zero dominance ",
           "deviation `d` has no heterozygous individuals, so that `d` cannot be ",
           "simulated -- the genotype is (near-)inbred at the locus. Choose ",
           "het-bearing loci (qtn =), pre-filter with filter_geno(hets = ",
           "\"include\"), set that d = 0, or use an outbred / F2 population.",
           call. = FALSE)
    }
  }
  .add_layer(sim, layer)
}

#' Resolve the orthogonal dominance-deviation series `d` to length n_qtn
#' @keywords internal
#' @noRd
.orthogonal_d_series <- function(d, nq) {
  if (is.null(d)) {
    return(rep(0, nq))
  }
  if (!is.numeric(d) || any(!is.finite(d))) {
    stop("additive(orthogonal = TRUE): `d` must be finite numeric.",
         call. = FALSE)
  }
  if (length(d) == 1L) {
    return(rep(as.numeric(d), nq))
  }
  if (length(d) == nq) {
    return(as.numeric(d))
  }
  stop("additive(orthogonal = TRUE): `d` must be a single value or a vector of ",
       "length n_qtn (", nq, ").", call. = FALSE)
}

#' Add a dominance variance component
#'
#' @inheritParams additive
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`).
#' @return the updated `phenotype_sim`.
#' @details
#' Dominance is modelled as a deviation applied to heterozygotes (the het
#' indicator), with its share of phenotypic variance set by `prop`. In this
#' variance-partition grammar there is no separate "degree of dominance"
#' argument: the ratio of dominance to additive variance is
#' `prop_dominance / prop_additive`, set through the layer proportions. (A single
#' degree-of-dominance scalar would be washed out by the per-component variance
#' scaling and is therefore not offered.)
#'
#' Dominance needs heterozygotes to be identifiable: it is identically zero at a
#' locus with no heterozygous individuals. If the selected loci carry none
#' (common on a near-inbred panel such as the bundled maize lines, ~0.4%
#' heterozygous), `dominance()` **errors** rather than silently substituting
#' loci. Pre-filter to heterozygous markers with `filter_geno(hets = "include")`,
#' fix het-bearing loci with `qtn =`, or simulate dominance on an outbred or
#' F2-type population instead.
#' @seealso [additive()], [epistasis()], [vqtl()], [filter_geno()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # Dominance reuses the additive QTNs by default, so the trait has
#' # h2 = 0.4 + 0.1 = 0.5 spread over one set of loci.
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.4, n_qtn = 5) |>
#'   dominance(prop = 0.1)
#' ph
#'
#' # same_as_add = FALSE draws separate dominance-only loci instead.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 2) |>
#'   additive(prop = 0.4, n_qtn = 5) |>
#'   dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 3)
dominance <- function(sim, prop = NULL, same_as_add = TRUE, n_qtn = NULL,
                      qtn = NULL, dist = "geometric") {
  .check_sim(sim)
  .require_markers(sim, "dominance")
  .validate_flag(same_as_add, "same_as_add")
  if (any(vapply(sim$layers, function(l) isTRUE(l$orthogonal), logical(1)))) {
    stop("dominance() cannot be combined with additive(orthogonal = TRUE): the ",
         "orthogonal genotypic model already carries the dominance deviation ",
         "(via `d`). Specify all dominance there, or use the standard ",
         "additive() + dominance() layers instead.", call. = FALSE)
  }
  prop <- .resolve_prop(sim, prop, "dominance")
  occ <- .type_occurrence(sim, "dominance")
  user_qtn <- .resolve_qtn_arg(sim, qtn, "dominance")

  add_layer <- .last_layer_of_type(sim, "additive")
  if (isTRUE(same_as_add) && is.null(user_qtn) && is.null(add_layer)) {
    stop("dominance(same_as_add = TRUE) requires a prior additive() layer.",
         call. = FALSE)
  }

  if (!is.null(user_qtn)) {
    nq <- length(user_qtn[[1]])
  } else if (isTRUE(same_as_add)) {
    nq <- add_layer$n_qtn
  } else {
    nq <- .resolve_n_qtn(sim, n_qtn, "dominance")
  }

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn)) {
      q <- user_qtn
    } else if (isTRUE(same_as_add)) {
      q <- .rep_qtn(add_layer, rep)      # follow additive's per-rep loci
    } else {
      q <- .draw_qtn(sim, nq, rep_seed)
    }
    e <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(nq, dist))
    list(qtn = q, effect = e)
  }

  fixed <- !is.null(user_qtn) ||
    (isTRUE(same_as_add) && is.null(add_layer$qtn_reps))
  drawn <- .draw_layer(sim, "dominance", occ, build, fixed = fixed)
  # A dominance deviation is a heterozygote effect, so it is identically zero at
  # a locus with no heterozygotes. Rather than silently substitute loci, fail
  # with a clear message so the user knows the (near-inbred) data has none.
  if (.dom_hetless(sim, drawn$qtn)) {
    stop("dominance(): the selected loci have no heterozygous individuals, so a ",
         "dominance deviation (which acts on heterozygotes) cannot be ",
         "simulated -- the genotype is (near-)inbred at those loci. Pre-filter ",
         "to heterozygous markers with filter_geno(hets = \"include\"), choose ",
         "loci with heterozygotes via qtn =, or use an outbred / F2 population.",
         call. = FALSE)
  } else if (.dom_partial_hetless(sim, drawn$qtn)) {
    warning("dominance(): some (but not all) selected loci have no heterozygous ",
            "individuals, so those loci contribute nothing and the remaining ",
            "het-bearing loci absorb the layer's `prop` -- the realized ",
            "dominance rests on fewer loci than requested. Choose het-bearing ",
            "loci (qtn =) or pre-filter with filter_geno(hets = \"include\") if ",
            "that is not intended.", call. = FALSE)
  }
  layer <- list(type = "dominance", prop = prop, n_qtn = nq, dist = dist,
                same_as_add = same_as_add,
                qtn = drawn$qtn, effect = drawn$effect)
  if (!is.null(drawn$qtn_reps)) {
    layer$qtn_reps <- drawn$qtn_reps
    layer$effect_reps <- drawn$effect_reps
  }
  layer$ld <- if (isTRUE(layer$same_as_add)) add_layer$ld else
    attr(drawn$qtn, "ld")
  .add_layer(sim, layer)
}

#' Add an epistatic variance component
#'
#' @inheritParams additive
#' @param n_pairs number of interacting QTN sets.
#' @param interaction number of markers per epistatic QTN (default 2, pairwise).
#' @param interaction_type how each marker in an interacting set contributes:
#'   `"a"` (additive -- the centered dosage) or `"d"` (dominance -- the centered
#'   heterozygote indicator). A single value applies to every position, or a
#'   length-`interaction` vector sets each position: `c("a", "a")` is
#'   additive-by-additive (default), `c("a", "d")` additive-by-dominance,
#'   `c("d", "d")` dominance-by-dominance.
#' @return the updated `phenotype_sim`.
#' @details
#' The epistatic value is the product of the interacting markers' design terms
#' times a per-set effect, then centered and scaled to `prop`. Each design term
#' is centered per locus before multiplying, which removes each locus's mean so
#' the product is not dominated by a constant offset. This does **not**
#' orthogonalize the interaction against the additive/dominance main effects:
#' the centered product can still be correlated with the constituent dosages --
#' strongly so when the interacting loci are in LD or away from allele frequency
#' 0.5 (with two perfectly linked loci it can be collinear with each additive
#' term). It is therefore not a Fisher/NOIA-orthogonal variance component (that
#' is a planned genotypic-model mode), so the "epistatic proportion" is a
#' simulation quantity, not an orthogonal share of variance.
#'
#' `interaction_type` values of `"d"` depend on heterozygotes and are
#' near-degenerate on a fully inbred panel (see [dominance()]); use `"a"` (the
#' default) for inbred data, or an outbred/F2-type population for `"d"`.
#' @seealso [additive()], [dominance()], [vqtl()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # Two pairs of interacting loci contributing 10% of the variance,
#' # on top of an additive background.
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.4, n_qtn = 5) |>
#'   epistasis(prop = 0.1, n_pairs = 2)
#' ph
#'
#' # interaction = 3 makes each epistatic term a three-way interaction
#' # rather than the default pairwise one.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 3) |>
#'   additive(prop = 0.3, n_qtn = 4) |>
#'   epistasis(prop = 0.2, n_pairs = 2, interaction = 3)
#'
#' # Additive-by-dominance and dominance-by-dominance interactions depend on
#' # heterozygotes, so run them on a segregating F2 (the maize282 panel is
#' # inbred and has none).
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' f2 <- selfcross(cross(pop[1], pop[2], n = 1, seed = 1), n = 40, seed = 2)
#' simulate_phenotype(f2, seed = 4) |>
#'   additive(prop = 0.3, n_qtn = 4) |>
#'   epistasis(prop = 0.1, n_pairs = 2, interaction_type = c("a", "d")) |>
#'   epistasis(prop = 0.1, n_pairs = 2, interaction_type = c("d", "d"))
epistasis <- function(sim, prop = NULL, n_pairs = NULL, interaction = 2,
                      interaction_type = "a", qtn = NULL, effect = NULL,
                      dist = "geometric") {
  .check_sim(sim)
  .require_markers(sim, "epistasis")
  interaction <- .validate_count(interaction, "interaction", minimum = 2L)
  prop <- .resolve_prop(sim, prop, "epistasis")
  user_pairs <- .resolve_epi_qtn(sim, qtn, interaction)
  if (!is.null(user_pairs)) interaction <- ncol(user_pairs[[1]])
  itype <- .resolve_interaction_type(interaction_type, interaction)
  np <- if (!is.null(user_pairs)) nrow(user_pairs[[1]]) else
    .resolve_n_qtn(sim, n_pairs, "epistasis", arg = "n_pairs")
  occ <- .type_occurrence(sim, "epistasis")

  build <- function(rep_seed, rep = 0L) {
    q <- if (!is.null(user_pairs)) user_pairs else
      .draw_qtn_pairs(sim, np, interaction, rep_seed)
    e <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(np, dist, effect))
    list(qtn = q, effect = e)
  }

  drawn <- .draw_layer(sim, "epistasis", occ, build,
                       fixed = !is.null(user_pairs))
  # A "d" interaction position acts on heterozygotes, so a hetless locus there
  # makes the pair identically zero and it silently drops out while the rest
  # absorb `prop`. Warn (as dominance() does for its partial case) rather than
  # realize fewer effective pairs than requested without notice.
  if (.epi_hetless_d(sim, drawn$qtn, itype)) {
    warning("epistasis(): a \"d\" interaction position sits on a locus with no ",
            "heterozygous individuals, so that interaction term is identically ",
            "zero and its pair contributes nothing while the remaining pairs ",
            "absorb `prop` -- the genotype is (near-)inbred at that locus. Choose ",
            "het-bearing loci (qtn =), pre-filter with filter_geno(hets = ",
            "\"include\"), or set that position's interaction_type to \"a\" if ",
            "that is not intended.", call. = FALSE)
  }
  layer <- list(type = "epistasis", prop = prop, n_pairs = np,
                interaction = interaction, interaction_type = itype,
                dist = dist,
                qtn = drawn$qtn, effect = drawn$effect)
  if (!is.null(drawn$qtn_reps)) {
    layer$qtn_reps <- drawn$qtn_reps
    layer$effect_reps <- drawn$effect_reps
  }
  .add_layer(sim, layer)
}

#' Add a variance-QTL (vQTL) component
#'
#' A variance QTL affects the *spread* of a trait rather than its mean:
#' genotypes at a vQTL differ in how variable they are. Here `prop` is the
#' marginal phenotypic-variance share assigned to a heterogeneous residual
#' component; it is not genetic variance and is not included in broad-sense
#' heritability. An explicit `prop` is therefore always required, even when the
#' simulation has an `h2` budget.
#'
#' Conditional residual variance uses a log link,
#' `log Var(E_v | genotype) = constant + loading`, where `loading` is the sum of
#' standardized vQTL genotype scores. This guarantees positive conditional
#' variances. The resulting heterogeneous residual is scaled to sample variance
#' `prop`; finite-sample covariance among components means total realized
#' phenotypic variance need not be exactly one.
#'
#' The log-linear variance link is this package's simulation choice (a
#' standard double-generalized-linear-model parameterization of the residual
#' variance); it is **not** a reproduction of the specific phenotype generator or
#' fitted variance model of any one reference. The references below are for the
#' vQTL concept and its detection in experimental crosses and plants -- the
#' setting this layer creates data for -- rather than the source of this exact
#' equation. In particular, Murphy et al. simulate variance heterogeneity with a
#' different generator and fit a DGLM to detect it; only the modeling idea is
#' shared, not the formula used here.
#'
#' @inheritParams additive
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`).
#' @return the updated `phenotype_sim`.
#' @references
#' Ronnegard, L. and Valdar, W. (2011). Detecting major genetic loci
#' controlling phenotypic variability in experimental crosses. \emph{Genetics}
#' 188, 435--447. \doi{10.1534/genetics.111.127068} (vQTL concept / detection.)
#'
#' Murphy, M.D., Fernandes, S.B., Morota, G. et al. (2022). Assessment of two
#' statistical approaches for variance genome-wide association studies in plants.
#' \emph{Heredity} 129, 93--102. \doi{10.1038/s41437-022-00541-1} (plant vGWAS
#' context; its simulation and DGLM differ from the log link used here.)
#' @seealso [additive()], [dominance()], [epistasis()].
#' @rdname vqtl-layer
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#'
#' # A vQTL on the same loci that carry the additive effects (the default).
#' ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
#'   additive(prop = 0.3, n_qtn = 4) |>
#'   vqtl(prop = 0.2)
#' ph
#'
#' # Independent variance-only loci: set same_as_add = FALSE and say how many.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 2) |>
#'   additive(prop = 0.3, n_qtn = 4) |>
#'   vqtl(prop = 0.2, same_as_add = FALSE, n_qtn = 2)
vqtl <- function(sim, prop = NULL, same_as_add = TRUE, n_qtn = NULL,
                 qtn = NULL, dist = "geometric") {
  .check_sim(sim)
  .require_markers(sim, "vqtl")
  .validate_flag(same_as_add, "same_as_add")
  .cite_vqtl()
  prop <- .resolve_prop(sim, prop, "vqtl")
  occ <- .type_occurrence(sim, "vqtl")
  user_qtn <- .resolve_qtn_arg(sim, qtn, "vqtl")

  add_layer <- .last_layer_of_type(sim, "additive")
  if (isTRUE(same_as_add) && is.null(user_qtn) && is.null(add_layer)) {
    stop("vqtl(same_as_add = TRUE) requires a prior additive() layer.",
         call. = FALSE)
  }
  if (!is.null(user_qtn)) {
    nq <- length(user_qtn[[1]])
  } else if (isTRUE(same_as_add)) {
    nq <- add_layer$n_qtn
  } else {
    nq <- .resolve_n_qtn(sim, n_qtn, "vqtl")
  }

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn)) {
      q <- user_qtn
    } else if (isTRUE(same_as_add)) {
      q <- .rep_qtn(add_layer, rep)
    } else {
      q <- .draw_qtn(sim, nq, rep_seed)
    }
    e <- lapply(seq_len(sim$n_traits), function(t) .effect_series(nq, dist))
    list(qtn = q, effect = e)
  }

  fixed <- !is.null(user_qtn) ||
    (isTRUE(same_as_add) && is.null(add_layer$qtn_reps))
  drawn <- .draw_layer(sim, "vqtl", occ, build, fixed = fixed)
  layer <- list(type = "vqtl", prop = prop, n_qtn = nq, dist = dist,
                same_as_add = same_as_add, qtn = drawn$qtn,
                effect = drawn$effect)
  if (!is.null(drawn$qtn_reps)) {
    layer$qtn_reps <- drawn$qtn_reps
    layer$effect_reps <- drawn$effect_reps
  }
  layer$ld <- if (isTRUE(same_as_add)) add_layer$ld else attr(drawn$qtn, "ld")
  .add_layer(sim, layer)
}

# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

#' Validate a phenotype_sim argument
#' @keywords internal
#' @noRd
.check_sim <- function(sim) {
  if (!inherits(sim, "phenotype_sim")) {
    stop("Expected a `phenotype_sim` (the output of simulate_phenotype()).",
         call. = FALSE)
  }
}

#' Resolve a layer's QTN count against the baseline, warning on override
#' @keywords internal
#' @noRd
.resolve_n_qtn <- function(sim, n_qtn, type, arg = "n_qtn") {
  baseline <- sim$n_qtn
  if (is.null(n_qtn)) {
    if (is.null(baseline) || baseline <= 0) {
      stop(type, "() requires ", arg, " (no positive baseline n_qtn was set ",
           "in simulate_phenotype()).", call. = FALSE)
    }
    return(.validate_count(baseline, arg, minimum = 1L))
  }
  n_qtn <- .validate_count(n_qtn, arg, minimum = 1L)
  # Only the additive/dominance/vqtl layers share the n_qtn baseline; epistasis
  # is parameterized by n_pairs and must not warn against it.
  if (arg == "n_qtn" && !is.null(baseline) && baseline > 0 &&
      baseline != n_qtn) {
    warning("Per-layer ", arg, " = ", n_qtn, " overrides the baseline n_qtn = ",
            baseline, " for the ", type, "() layer.", call. = FALSE)
  }
  n_qtn
}

#' Most recent layer of a given type, or NULL
#' @keywords internal
#' @noRd
.last_layer_of_type <- function(sim, type) {
  hits <- Filter(function(l) identical(l$type, type), sim$layers)
  if (length(hits) == 0) NULL else hits[[length(hits)]]
}

#' Explain an exhausted budget when it came from the one-call form
#' @keywords internal
#' @noRd
.one_call_hint <- function(sim) {
  if (!isTRUE(sim$one_call)) {
    return("")
  }
  paste0("\nsimulate_phenotype() was given both `h2` and `n_qtn`, so it ",
         "already built a complete model that uses all of h2. To add more ",
         "layers, drop `n_qtn` from simulate_phenotype() and build the model ",
         "with the layers instead.")
}

#' Resolve a layer's `prop` against an `h2` budget set in simulate_phenotype()
#'
#' With `h2` set, `prop` may be omitted and takes whatever of the budget is
#' still unspent -- so a lone `additive()` layer gets `prop = h2` -- and any
#' `prop` that is supplied must keep the running total within `h2`.
#' @keywords internal
#' @noRd
.resolve_prop <- function(sim, prop, type) {
  nt <- sim$n_traits
  h2 <- sim$h2
  spent <- .total_genetic_prop(sim)

  if (is.null(prop)) {
    if (identical(type, "vqtl")) {
      stop("vqtl() requires an explicit `prop`: residual heterogeneity is not ",
           "part of the broad-sense `h2` budget.", call. = FALSE)
    }
    if (is.null(h2)) {
      stop(type, "() requires `prop` (no `h2` was set in ",
           "simulate_phenotype() to draw the remaining budget from).",
           call. = FALSE)
    }
    remaining <- .expand_prop(h2, nt) - spent
    if (any(remaining <= 1e-10)) {
      stop("No heritability budget is left for the ", type, "() layer: h2 = ",
           paste(sprintf("%.3f", .expand_prop(h2, nt)), collapse = ", "),
           " is already fully allocated. Give the layers explicit `prop` ",
           "values that sum to h2.", .one_call_hint(sim), call. = FALSE)
    }
    return(remaining)
  }

  prop <- .validate_proportion(prop, "prop", nt)

  if (!is.null(h2) && !identical(type, "vqtl")) {
    total <- spent + .expand_prop(prop, nt)
    target <- .expand_prop(h2, nt)
    if (any(total > target + 1e-8)) {
      bad <- which(total > target + 1e-8)
      stop("Genetic proportions exceed the h2 set in simulate_phenotype() ",
           "for trait(s) ", paste(bad, collapse = ", "), ": running total = ",
           paste(sprintf("%.3f", total[bad]), collapse = ", "),
           " against h2 = ", paste(sprintf("%.3f", target[bad]),
                                   collapse = ", "),
           ". Layer `prop` values must sum to h2.",
           .one_call_hint(sim), call. = FALSE)
    }
  }
  prop
}

.add_layer <- function(sim, layer) {
  nt <- sim$n_traits
  layer$prop <- .validate_proportion(layer$prop, "prop", nt)
  prospective <- .total_variance_prop(sim) + .expand_prop(layer$prop, nt)
  if (any(prospective > 1 + 1e-8)) {
    bad <- which(prospective > 1 + 1e-8)
    stop("Variance proportions for the ", layer$type, "() layer push trait(s) ",
         paste(bad, collapse = ", "), " above 1 (running total = ",
         paste(sprintf("%.3f", prospective[bad]), collapse = ", "), ").",
         call. = FALSE)
  }
  sim$layers <- c(sim$layers, list(layer))
  .realize_phenotype(sim)
}

#' Draw a layer once, or once per replication when vary_qtn is TRUE
#'
#' `build(rep_seed, rep)` returns `list(qtn, effect)` for one draw. This calls
#' it for the canonical set (rep 1) and, when the simulation was created with
#' `vary_qtn = TRUE`, once more per replication with a rep-specific sub-seed so
#' every replication gets independent QTNs. `fixed = TRUE` (user-supplied loci)
#' short-circuits: fixed loci are never redrawn.
#' @keywords internal
#' @noRd
.draw_layer <- function(sim, type, occ, build, fixed = FALSE) {
  canonical <- build(.layer_seed(sim$seed, type, occ), 1L)
  if (!isTRUE(sim$vary_qtn) || sim$n_reps <= 1 || fixed) {
    return(list(qtn = canonical$qtn, effect = canonical$effect))
  }
  reps <- vector("list", sim$n_reps)
  reps[[1L]] <- canonical
  for (r in seq_len(sim$n_reps)[-1L]) {
    reps[[r]] <- build(.layer_seed(sim$seed, paste0(type, "_rep", r), occ), r)
  }
  list(
    qtn        = canonical$qtn,
    effect     = canonical$effect,
    qtn_reps   = lapply(reps, function(x) x$qtn),
    effect_reps = lapply(reps, function(x) x$effect)
  )
}

#' TRUE when a reused QTN set cannot support any dominance deviation
#'
#' A dominance layer is a heterozygote effect, so a trait whose loci are *all*
#' homozygous realizes exactly zero variance and cannot fill `prop` at all --
#' that is an error. A set with *some* hetless loci can still realize `prop` from
#' the rest, but the dead loci are silently inert; see [.dom_partial_hetless()],
#' which warns for that case.
#' @keywords internal
#' @noRd
.dom_hetless <- function(sim, q) {
  any(vapply(q, function(idx) {
    if (is.null(idx) || length(idx) == 0L) {
      return(FALSE)
    }
    sum(.geno_cols(sim, idx) == 0, na.rm = TRUE) == 0
  }, logical(1)))
}

#' TRUE when a QTN set has some (but not all) hetless loci
#'
#' The layer still realizes `prop` from the het-bearing loci, but the hetless
#' ones contribute nothing, so the dominance rests on fewer loci than requested.
#' Checked per locus; the all-hetless case is [.dom_hetless()]'s error, not this.
#' @keywords internal
#' @noRd
.dom_partial_hetless <- function(sim, q) {
  any(vapply(q, function(idx) {
    if (is.null(idx) || length(idx) < 2L) {
      return(FALSE)
    }
    hl <- vapply(idx, function(j) {
      sum(.geno_cols(sim, j) == 0, na.rm = TRUE) == 0
    }, logical(1))
    any(hl) && !all(hl)
  }, logical(1)))
}

#' TRUE when any orthogonal locus carrying a non-zero `d` has no heterozygotes
#'
#' Unlike [.dom_hetless()] (which asks whether a whole reused set is homozygous),
#' the orthogonal model attaches a dominance deviation to specific loci, so the
#' check is per locus: a `d != 0` locus with no heterozygote would silently
#' contribute nothing. `qtn` and `d_effect` are per-trait lists in matching order.
#' @keywords internal
#' @noRd
.orthogonal_hetless_d <- function(sim, qtn, d_effect) {
  any(vapply(seq_along(qtn), function(t) {
    idx <- qtn[[t]]
    dv  <- d_effect[[t]]
    if (is.null(idx) || length(idx) == 0L) {
      return(FALSE)
    }
    nz <- which(dv != 0)
    if (length(nz) == 0L) {
      return(FALSE)
    }
    any(vapply(idx[nz], function(j) {
      sum(.geno_cols(sim, j) == 0, na.rm = TRUE) == 0
    }, logical(1)))
  }, logical(1)))
}

#' TRUE when any epistasis pair carries a hetless locus at a `"d"` position
#'
#' A `"d"` interaction position contributes a centered heterozygote indicator, so
#' a locus with no heterozygotes makes that design column identically zero and
#' silently inerts the whole pair -- the remaining pairs then absorb `prop`. The
#' check is per locus at every `"d"` position, mirroring [.orthogonal_hetless_d()]
#' and [.dom_hetless()]. `qtn` is a per-trait list of `n_pairs x interaction`
#' index matrices; `itype` is the length-`interaction` "a"/"d" vector.
#' @keywords internal
#' @noRd
.epi_hetless_d <- function(sim, qtn, itype) {
  d_pos <- which(itype == "d")
  if (length(d_pos) == 0L) {
    return(FALSE)
  }
  any(vapply(qtn, function(mat) {
    if (is.null(mat) || length(mat) == 0L) {
      return(FALSE)
    }
    loci <- unique(as.integer(mat[, d_pos, drop = FALSE]))
    any(vapply(loci, function(j) {
      sum(.geno_cols(sim, j) == 0, na.rm = TRUE) == 0
    }, logical(1)))
  }, logical(1)))
}

#' A prior layer's QTNs for a given replication (per-rep if it varied)
#' @keywords internal
#' @noRd
.rep_qtn <- function(layer, rep) {
  if (rep >= 1L && !is.null(layer$qtn_reps)) {
    return(layer$qtn_reps[[rep]])
  }
  layer$qtn
}

#' Resolve a user-supplied `qtn` argument to a per-trait list of indices
#'
#' Accepts marker names (matched against the map) or 1-based column indices,
#' as a single vector applied to every trait or a length-`n_traits` list. NULL
#' passes through as "draw at random".
#' @keywords internal
#' @noRd
.resolve_qtn_arg <- function(sim, qtn, type) {
  if (is.null(qtn)) {
    return(NULL)
  }
  nt <- sim$n_traits
  as_idx <- function(v) {
    if (is.character(v)) {
      idx <- match(v, sim$map$snp)
      if (anyNA(idx)) {
        bad <- v[is.na(idx)]
        stop(type, "(qtn=): marker(s) not found in the genotype map: ",
             paste(utils::head(bad, 5), collapse = ", "),
             if (length(bad) > 5) ", ..." else "", ".", call. = FALSE)
      }
      if (!length(idx) || anyDuplicated(idx)) {
        stop(type, "(qtn=): loci must be non-empty and must not be ",
             "duplicated within a trait.", call. = FALSE)
      }
      return(as.integer(idx))
    }
    if (!is.numeric(v) || any(!is.finite(v)) || any(v != floor(v))) {
      stop(type, "(qtn=): numeric indices must be finite whole numbers.",
           call. = FALSE)
    }
    v <- as.integer(v)
    if (!length(v) || any(v < 1L | v > sim$n_markers)) {
      stop(type, "(qtn=): index out of range (1..", sim$n_markers, ").",
           call. = FALSE)
    }
    if (anyDuplicated(v)) {
      stop(type, "(qtn=): loci must not be duplicated within a trait.",
           call. = FALSE)
    }
    v
  }
  if (is.list(qtn)) {
    if (length(qtn) != nt) {
      stop(type, "(qtn=): a list must have one element per trait (", nt,
           "); got ", length(qtn), ".", call. = FALSE)
    }
    out <- lapply(qtn, as_idx)
    lens <- vapply(out, length, 0L)
    if (length(unique(lens)) != 1L) {
      stop(type, "(qtn=): every trait must get the same number of QTNs; got ",
           paste(lens, collapse = ", "), ".", call. = FALSE)
    }
    return(out)
  }
  idx <- as_idx(qtn)
  rep(list(idx), nt)
}

#' Resolve a user-supplied epistatic `qtn` argument to an n_pairs x interaction
#' matrix of indices, shared across traits
#' @keywords internal
#' @noRd
.resolve_epi_qtn <- function(sim, qtn, interaction) {
  if (is.null(qtn)) {
    return(NULL)
  }
  to_idx <- function(v) {
    if (is.character(v)) {
      idx <- match(v, sim$map$snp)
      if (anyNA(idx)) {
        stop("epistasis(qtn=): marker(s) not found in the map: ",
             paste(utils::head(v[is.na(idx)], 5), collapse = ", "), ".",
             call. = FALSE)
      }
      return(as.integer(idx))
    }
    if (!is.numeric(v) || any(!is.finite(v)) || any(v != floor(v))) {
      stop("epistasis(qtn=): numeric indices must be finite whole numbers.",
           call. = FALSE)
    }
    as.integer(v)
  }
  m <- if (is.list(qtn)) {
    do.call(rbind, lapply(qtn, to_idx))
  } else if (is.matrix(qtn)) {
    matrix(to_idx(as.vector(qtn)), nrow = nrow(qtn))
  } else {
    if (length(qtn) %% interaction != 0L) {
      stop("epistasis(qtn=): the number of loci must be a multiple of ",
           "`interaction` (", interaction, ").", call. = FALSE)
    }
    matrix(to_idx(qtn), ncol = interaction, byrow = TRUE)
  }
  if (ncol(m) != interaction) {
    stop("epistasis(qtn=): each set must have `interaction` = ", interaction,
         " markers; got ", ncol(m), ".", call. = FALSE)
  }
  if (!length(m) || any(m < 1L | m > sim$n_markers)) {
    stop("epistasis(qtn=): index out of range (1..", sim$n_markers, ").",
         call. = FALSE)
  }
  if (any(apply(m, 1, anyDuplicated) > 0L)) {
    stop("epistasis(qtn=): a locus cannot appear twice in one interaction ",
         "set.", call. = FALSE)
  }
  rep(list(m), sim$n_traits)   # same interacting sets across traits
}

#' Resolve an epistasis `interaction_type` to a length-`interaction` vector
#'
#' `"a"` = additive (centered dosage), `"d"` = dominance (centered heterozygote
#' indicator). A scalar is recycled to every position; a vector must match the
#' number of markers per set.
#' @keywords internal
#' @noRd
.resolve_interaction_type <- function(interaction_type, interaction) {
  it <- as.character(interaction_type)
  if (!all(it %in% c("a", "d"))) {
    stop("epistasis(interaction_type=): each entry must be \"a\" (additive) or ",
         "\"d\" (dominance); got ",
         paste(setdiff(it, c("a", "d")), collapse = ", "), ".", call. = FALSE)
  }
  if (length(it) == 1L) {
    return(rep(it, interaction))
  }
  if (length(it) != interaction) {
    stop("epistasis(interaction_type=): length (", length(it),
         ") must be 1 or `interaction` (", interaction, ").", call. = FALSE)
  }
  it
}

#' Citation notice for the vQTL simulation, once per session
#'
#' The variance-QTL simulation is described in Murphy et al. (2022); credit it
#' the first time a vqtl() layer is added in a session. An [rlang::inform()]
#' message, silenceable with `suppressMessages()`.
#' @keywords internal
#' @noRd
.cite_vqtl <- function() {
  rlang::inform(
    .cite_main("Murphy et al. (2022), Heredity 129:93-102, ",
               "doi:10.1038/s41437-022-00541-1, for the plant variance-QTL ",
               "(vGWAS) context (the log-variance link used here is this ",
               "package's own simulation choice, not Murphy et al.'s ",
               "generator)."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_vqtl_citation"
  )
}

#' Impose coupling / repulsion effect-sign pattern on a per-trait effect list
#'
#' Repulsion alternates the sign of successive QTN effects within each trait, in
#' draw order (`+, -, +, -, ...`). This is a positional sign pattern, not a phase
#' derived from the haplotypes; it produces repulsion-style cancellation only
#' when the alternately-signed loci are actually in positive LD (see the `phase`
#' argument of additive()). Coupling returns the effects unchanged.
#' @keywords internal
#' @noRd
.apply_phase <- function(eff_list, phase) {
  if (identical(phase, "coupling")) {
    return(eff_list)
  }
  lapply(eff_list, function(e) {
    signs <- rep(c(1, -1), length.out = length(e))
    e * signs
  })
}

#' Error if a marker-based layer is added to a genotype-free phenotype
#'
#' A phenotype built from expression alone (`simulate_phenotype(expression = ...)`
#' with no `geno`) has no markers, so additive/dominance/epistasis/vqtl layers
#' cannot be scored. Only `transcriptome()` layers are valid there.
#' @keywords internal
#' @noRd
.require_markers <- function(sim, fn) {
  if (is.null(sim$n_markers) || sim$n_markers < 1L) {
    stop(fn, "() needs genotypes, but this phenotype was built from expression ",
         "alone (no `geno`). Use transcriptome() layers, or rebuild with `geno`.",
         call. = FALSE)
  }
}
