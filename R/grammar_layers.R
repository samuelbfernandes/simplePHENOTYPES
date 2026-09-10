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
#' @param phase QTN linkage phase, `"coupling"` (default) or `"repulsion"`.
#'   Under repulsion the effect signs alternate across the layer's QTNs, so the
#'   increasing allele at one locus is paired with the decreasing allele at the
#'   next -- the classic repulsion architecture, in which a GWAS sees effects
#'   that partially cancel and linked causal loci mask one another. The realized
#'   genetic variance is still `prop` (variance is pinned by the partition, not
#'   by phase); what changes is the sign structure of the per-QTN effects and
#'   the covariance among loci. Coupling leaves the drawn signs untouched.
#'   `phase` is **architecture-independent** -- it acts on whichever QTNs the
#'   additive layer draws, so it is meaningful under `"independent"`,
#'   `"pleiotropy"` and `"ld"` alike. Its effect on the cross-locus covariance is
#'   only material when a trait's additive QTNs are physically linked (on the
#'   same chromosome); for QTNs that are effectively unlinked, alternating the
#'   signs just relabels which allele is "increasing" and leaves the covariance
#'   structure unchanged.
#' @param dist within-layer effect distribution (default "geometric").
#' @return the updated `phenotype_sim`.
#' @details
#' The additive value is the -1/0/1 dosage weighted by the QTN effects, centered,
#' and scaled to `prop` of the phenotypic variance. This is a simulation
#' convention, not Fisher's average-effect decomposition; for an additive-only
#' model `prop` equals the narrow-sense h2 under Hardy-Weinberg.
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
additive <- function(sim, prop = NULL, n_qtn = NULL, qtn = NULL, effect = NULL,
                     phase = c("coupling", "repulsion"), dist = "geometric") {
  .check_sim(sim)
  phase <- match.arg(phase)
  prop <- .resolve_prop(sim, prop, "additive")
  user_qtn <- .resolve_qtn_arg(sim, qtn, "additive")
  nq <- if (!is.null(user_qtn)) length(user_qtn[[1]]) else
    .resolve_n_qtn(sim, n_qtn, "additive")
  occ <- .type_occurrence(sim, "additive")

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn)) {
      q <- user_qtn
      e <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, effect))
    } else if (sim$architecture == "pleiotropy" && sim$n_traits > 1) {
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
                  function(t) .effect_series(nq, dist, effect))
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
  .add_layer(sim, layer)
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
  .validate_flag(same_as_add, "same_as_add")
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
#' is **centered per locus** before multiplying, which subtracts out the
#' lower-order (additive/dominance main-effect) marginals, so the component is a
#' cleaner a x a / a x d / d x d interaction. It is still not a full
#' Fisher-orthogonal variance component (that is a planned genotypic-model mode),
#' so the "epistatic proportion" remains a simulation quantity.
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
#' # Additive-by-dominance and dominance-by-dominance interactions.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 4) |>
#'   additive(prop = 0.3, n_qtn = 4) |>
#'   epistasis(prop = 0.1, n_pairs = 2, interaction_type = c("a", "d")) |>
#'   epistasis(prop = 0.1, n_pairs = 2, interaction_type = c("d", "d"))
epistasis <- function(sim, prop = NULL, n_pairs = NULL, interaction = 2,
                      interaction_type = "a", qtn = NULL, effect = NULL,
                      dist = "geometric") {
  .check_sim(sim)
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
#' @inheritParams additive
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`).
#' @return the updated `phenotype_sim`.
#' @references
#' Ronnegard, L. and Valdar, W. (2011). Detecting major genetic loci
#' controlling phenotypic variability in experimental crosses. \emph{Genetics}
#' 188, 435--447. \doi{10.1534/genetics.111.127068}
#'
#' Murphy, M.D., Fernandes, S.B., Morota, G. et al. (2022). Assessment of two
#' statistical approaches for variance genome-wide association studies in plants.
#' \emph{Heredity} 129, 93--102. \doi{10.1038/s41437-022-00541-1}
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

#' TRUE when a reused QTN set cannot support a dominance deviation
#'
#' A dominance layer is a heterozygote effect, so a trait whose reused loci carry
#' no heterozygotes would contribute exactly zero variance and fail to realize
#' `prop`. Returns TRUE if any trait's loci are entirely homozygous.
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
               "doi:10.1038/s41437-022-00541-1, for plant vQTL simulation."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_vqtl_citation"
  )
}

#' Impose coupling / repulsion phase on a per-trait effect list
#'
#' Repulsion alternates the sign of successive QTN effects within each trait, so
#' the increasing alleles oppose one another and the genetic variance is partly
#' cancelled. Coupling returns the effects unchanged.
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
