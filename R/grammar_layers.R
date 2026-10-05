#' Add an additive variance component
#'
#' @param sim a `phenotype_sim`.
#' @param prop proportion of total phenotypic variance (scalar or length
#'   `n_traits`). May be omitted when `h2` was set in [simulate_phenotype()],
#'   in which case the layer takes whatever of that budget is unspent -- so a
#'   single `additive()` layer gets `prop = h2`. When `prop` is given and `h2`
#'   is set, the layer `prop` values must sum to `h2`.
#' @param n_qtn number of additive QTNs; overrides the baseline `n_qtn` from
#'   [simulate_phenotype()] (with a warning when both are given). Ignored, with a
#'   warning, when `qtn` fixes the loci.
#' @param qtn optional user-supplied QTNs for this layer, so you can fix the
#'   causal loci for one effect type while the others are drawn at random. Give
#'   marker names (matched against the map) or column indices; a vector is used
#'   for every trait, or a length-`n_traits` list sets each trait's loci
#'   separately (every element the same length). `n_qtn` is then taken from what
#'   you supply. Fixed loci are not redrawn by `vary_qtn`. (For [epistasis()],
#'   `qtn` describes interacting sets; see there.) Every architecture accepts
#'   `qtn`, and each keeps its construction: under `"independent"` the loci
#'   carry the effect series; under `"pleiotropy"` (multi-trait) every locus
#'   affects every trait, so give a single vector (or the same loci for each
#'   trait in a list); loci that affect only some traits are partial pleiotropy,
#'   built with [complex_phenotypes()]. If a correlation is controlled (`cor`,
#'   `pi`, ...) the correlated draw sets the effects and the genetic correlation
#'   is its target; otherwise the correlation is just an outcome of the shared
#'   loci, and an explicit `effect` sets the effects (without it the default
#'   draw does, with implicit `cor = 0`). `pi` only matters for controlling a
#'   correlation, and `pi < 1` cannot be used with fixed shared loci. Under
#'   `"ld"` (`ld_type = "direct"` only) give
#'   `qtn = list(trait1_loci, trait2_loci)` of equal length -- the i-th locus of
#'   each is a linked pair on one chromosome, no marker may be causal for both traits, and each
#'   pair's r2 is reported (by [qtn_table()]) and warned about when outside
#'   `[r2_min, r2_max]`. A marker that is monomorphic (or heterozygous in every
#'   individual) carries no variance: it is accepted with a warning, except
#'   where it makes the architecture undefined (an `"ld"` pair, or a pleiotropy
#'   layout with no varying locus), which is an error.
#' @param effect optional geometric base (scalar) or explicit effect series
#'   (length `n_qtn`), used for every trait; or a length-`n_traits` list of these,
#'   one per trait -- e.g. to re-score per-trait effects frozen from an earlier
#'   simulation on fixed `qtn` in a single layer. Under
#'   `architecture = "pleiotropy"` (multi-trait) it is accepted only when no
#'   correlation is controlled (no `cor`, `pi`, ...): it then sets the effects of
#'   the shared loci (the same series for every trait unless a per-trait list is
#'   given) and the genetic correlation is an outcome. Without `effect` (or with
#'   `cor` / `pi`) the pleiotropy draw sets the effects. Not available with
#'   `orthogonal = TRUE` (use `a`).
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
#'   the additive layer draws. For a **genuine haplotype-derived phase** between
#'   two linked traits use `architecture = "ld"` with its `ld_phase` argument
#'   (see [simulate_phenotype()]): it reads the sign of the dosage correlation
#'   of each linked pair and sets trait 2's effect sign so that every pair's
#'   linkage-induced covariance is positive (`"coupling"`) or negative
#'   (`"repulsion"`). The two may be combined: `phase` is applied first, then
#'   `ld_phase` re-signs trait 2's effects pair by pair, so under `"ld"` with a
#'   non-default `ld_phase` the positional alternation is only guaranteed to
#'   survive on trait 1.
#' @param dist within-layer effect distribution (default "geometric").
#' @return the updated `phenotype_sim`.
#' @details
#' The additive value is the -1/0/1 dosage weighted by the QTN effects, centered,
#' and scaled to `prop` of the phenotypic variance. This is a simulation
#' convention, not Fisher's average-effect decomposition; for an additive-only
#' model `prop` equals the narrow-sense h2 under Hardy-Weinberg.
#'
#' Combining this layer with [dominance()] on the same loci (the default) makes
#' the realized genetic variance differ from `prop_A + prop_D` by
#' `2Cov(c_A, c_D)` (with one additive and one dominance layer; several layers of
#' one type add their mutual covariance to Var(c_A) or Var(c_D)), a structural
#' term that depends on the allele frequencies and on which allele is coded +1
#' (see [simulate_phenotype()]); the object reports the realized Var(A), Var(D)
#' and 2Cov(A,D), and Var(c_A), Var(c_D) and 2Cov(c_A,c_D), in `$ad_report` and
#' in `print()`.
#' Use `orthogonal = TRUE` for a Fisher-orthogonal additive/dominance model.
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
#' incompatible with `vary_qtn`, and cannot be used with the
#' correlation-controlling `"pleiotropy"` (multi-trait) or `"ld"` architectures
#' (their effects come from the correlated draw; `qtn =` alone is accepted).
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
#' # A major QTN with a set share of the variance: stack two additive layers,
#' # one locus explaining 20% of V_P and 50 equal-effect loci another 30%
#' g <- SNP55K_maize282_maf04
#' major <- simulate_phenotype(g, seed = 1) |>
#'   additive(prop = 0.2, qtn = g$snp[100]) |>
#'   additive(prop = 0.3, n_qtn = 50, effect = rep(1, 50))
#' tb <- qtn_table(major)
#' head(tb[order(-tb$var_explained), c("snp", "var_explained")], 3)
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
  if (isTRUE(sim$frozen)) return(.frozen_layer_ignored(sim, "additive"))
  .require_markers(sim, "additive")
  .validate_flag(orthogonal, "orthogonal")
  phase <- match.arg(phase)
  if (orthogonal) {
    if (!is.null(effect)) {
      stop("additive(orthogonal = TRUE): give additive effects via `a` (and ",
           "dominance deviations via `d`), not `effect`.", call. = FALSE)
    }
    if (is.list(a)) {
      stop("additive(orthogonal = TRUE): `a` must be a scalar geometric base ",
           "or a length-n_qtn vector (one series for every trait).",
           call. = FALSE)
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
  # A list gives one effect specification per trait (each a geometric base or an
  # explicit series, checked by .effect_series() when the layer is built).
  if (is.list(effect)) {
    if (length(effect) != sim$n_traits) {
      stop("additive(effect=): a list must have one element per trait (",
           sim$n_traits, "); got ", length(effect), ".", call. = FALSE)
    }
    # NULL would silently fall back to the default series for that trait
    bad <- which(!vapply(effect, function(x) is.numeric(x) && length(x) > 0L,
                         logical(1)))
    if (length(bad)) {
      stop("additive(effect=): every list element must be a numeric geometric ",
           "base or effect series; element(s) ", paste(bad, collapse = ", "),
           " are not.", call. = FALSE)
    }
  }
  pleio <- identical(sim$architecture, "pleiotropy") && sim$n_traits > 1
  user_qtn <- .resolve_qtn_arg(sim, qtn, "additive", equal = !pleio)
  if (!is.null(user_qtn) && !is.null(n_qtn)) {
    warning("additive(): `n_qtn` is ignored because `qtn` fixes the loci (",
            max(lengths(user_qtn)), " QTNs).", call. = FALSE)
  }
  # Fixed loci keep each architecture's construction: under "pleiotropy" every
  # locus affects every trait (a locus for only some traits is partial
  # pleiotropy: complex_phenotypes()), and the effects are set as usual -- by
  # `effect` / `dist` when no correlation is being controlled, else by the
  # multivariate (PleioArch) draw that targets `cor` -- so the supplied loci only
  # replace the random choice; "ld" keeps distinct linked pairs: trait 1's i-th
  # locus is linked to trait 2's i-th (SPEC 4.1; DECISION-007/013/014/043).
  series_path <- pleio && !.pleio_controlled(sim) &&
    (!is.null(effect) || !identical(dist, "geometric"))
  pleio_layout <- NULL
  if (!is.null(user_qtn) && pleio) {
    pleio_layout <- .pleio_user_layout(user_qtn, "additive")
    user_qtn <- pleio_layout$q
    if (series_path || isTRUE(orthogonal)) pleio_layout <- NULL
  }
  if (!is.null(user_qtn) && identical(sim$architecture, "ld")) {
    .ld_check_cross_layer(sim, user_qtn, "additive")   # a vqtl()/dominance() may precede
    user_qtn <- .ld_user_pairs(sim, user_qtn, "additive")
  }
  # "ld" pairs each trait's causal loci as distinct linked markers within one
  # draw; a second additive draw knows nothing of the first and could make a
  # marker causal for both traits across layers (DECISION-014/023).
  if (sim$architecture == "ld" && !is.null(.last_layer_of_type(sim, "additive"))) {
    stop("additive(): architecture = \"ld\" supports a single additive layer. ",
         "Each additive() draw pairs distinct linked loci for the two traits, but ",
         "a second draw knows nothing of the first and can make one marker causal ",
         "for both traits across layers. Put all additive QTNs in one layer ",
         "(raise n_qtn).", call. = FALSE)
  }
  nq <- if (!is.null(user_qtn)) max(lengths(user_qtn)) else
    .resolve_n_qtn(sim, n_qtn, "additive")
  occ <- .type_occurrence(sim, "additive")

  # In the orthogonal genotypic model the additive effects come from `a`; the
  # standard variance-partition layer uses `effect`.
  eff_arg <- if (orthogonal) a else effect
  eff_for <- function(t) if (is.list(eff_arg)) eff_arg[[t]] else eff_arg
  eff_name <- if (orthogonal) "a" else "effect"
  d_series <- if (orthogonal) .orthogonal_d_series(d, nq) else NULL

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn) && is.null(pleio_layout)) {
      q <- user_qtn
      e <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, eff_for(t),
                                             arg = eff_name))
    } else if (!orthogonal && !series_path &&
               sim$architecture == "pleiotropy" && sim$n_traits > 1) {
      if (!is.null(effect) || !identical(dist, "geometric")) {
        stop("additive(): `effect` and non-default `dist` cannot be used under ",
             "architecture = \"pleiotropy\" while a correlation is controlled ",
             "(`cor`, `pi`, ...): the effects come from the multivariate draw ",
             "that controls `cor`.", call. = FALSE)
      }
      pd <- .pleio_draw(sim, nq, .expand_prop(prop, sim$n_traits), rep_seed,
                        fixed = pleio_layout)
      q <- pd$qtn
      e <- pd$effect
    } else {
      q <- .draw_qtn(sim, nq, rep_seed)
      e <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, eff_for(t),
                                             arg = eff_name))
    }
    e <- .apply_phase(e, phase)
    # The haplotype-derived phase of the "ld" pairs is imposed last so it holds
    # on the final effects whatever the positional alternation did; it only
    # re-signs trait 2 (see .apply_ld_phase()).
    if (sim$architecture == "ld") {
      e <- .apply_ld_phase(e, q, sim$arch_args[["ld_phase"]])
    }
    list(qtn = q, effect = e)
  }

  drawn <- .draw_layer(sim, "additive", occ, build, fixed = !is.null(user_qtn))
  if (!orthogonal && sim$architecture == "pleiotropy" && sim$n_traits > 1) {
    .pleio_total_cor_check(sim, .expand_prop(prop, sim$n_traits))
  }
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
    active <- .expand_prop(prop, sim$n_traits) > 0     # prop = 0 traits are inert
    if (.orthogonal_hetless_d(sim, layer$qtn[active], layer$d_effect[active])) {
      stop("additive(orthogonal = TRUE): a locus with a non-zero dominance ",
           "deviation `d` has no heterozygous individuals, or is heterozygous ",
           "in every individual, so its heterozygote indicator is constant and ",
           "that `d` cannot be simulated (no heterozygotes: the genotype is ",
           "(near-)inbred at the locus; all heterozygous: an F1-like locus whose ",
           "constant value is absorbed into the mean). Choose loci whose ",
           "heterozygote status varies (qtn =), pre-filter with filter_geno(hets ",
           "= \"include\"), set that d = 0, or use an outbred / F2 population.",
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
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`). Supplying
#'   `qtn` fixes the loci instead, so `same_as_add` is then treated as `FALSE`
#'   (the layer records and prints it as such).
#' @param n_qtn number of dominance QTNs for a fresh draw
#'   (`same_as_add = FALSE`); ignored, with a warning, when the loci are reused
#'   (`same_as_add = TRUE`) or fixed by `qtn`.
#' @param effect optional geometric base (scalar) or explicit effect series
#'   (length `n_qtn`, the number of dominance loci -- the additive layer's
#'   `n_qtn` when `same_as_add = TRUE`), used for every trait; or a
#'   length-`n_traits` list of these, one per trait (the counterpart of v1
#'   `dom_effect`). Under `architecture = "pleiotropy"` (multi-trait) accepted
#'   only when no correlation is controlled (no `cor`, `pi`, ...), as in
#'   [additive()].
#' @return the updated `phenotype_sim`.
#' @details
#' Dominance is modelled as a deviation applied to heterozygotes (the het
#' indicator), with its share of phenotypic variance set by `prop`. In this
#' variance-partition grammar there is no separate "degree of dominance"
#' argument (there is no `degree =`): the ratio of the *scaled components'*
#' variances is `prop_dominance / prop_additive`, set through the layer
#' proportions. (A single degree-of-dominance scalar would be washed out by the
#' per-component variance scaling and is therefore not offered; for a per-locus
#' degree of dominance `d / abs(a)` use `additive(orthogonal = TRUE, a =, d =)`.)
#'
#' **Additive + dominance on the same loci is not orthogonal.** The dosage and
#' the heterozygote indicator covary by \eqn{-(2p - 1) 2pq} at each locus, so the
#' realized genetic variance is `prop_A + prop_D + 2Cov(c_A, c_D)` for one
#' additive and one dominance layer (several layers of a type add their mutual
#' covariance to `Var(c_A)` or `Var(c_D)`) -- a structural,
#' allele-coding-dependent bias, not a finite-sample effect; its per-locus sign is
#' that of `1 - 2p` for the counted-allele frequency `p`, so it is one-signed
#' across loci only when those frequencies are on one side of 0.5 -- and the
#' classical Va:Vd of the block differs from `prop_A:prop_D`. The object reports
#' the realized Var(A), Var(D) and 2Cov(A,D), and Var(c_A), Var(c_D) and
#' 2Cov(c_A,c_D) (`$ad_report`, and the note printed by `print()`); use `additive(orthogonal = TRUE, a =, d =)` for a
#' Fisher-orthogonal partition.
#'
#' Dominance needs heterozygotes to be identifiable: it is identically zero at a
#' locus with no heterozygous individuals. If the selected loci carry none
#' (common on a near-inbred panel such as the bundled maize lines, ~0.4%
#' heterozygous), `dominance()` **errors** rather than silently substituting
#' loci; if only some of them carry none, it warns that the rest absorb `prop`.
#' Pre-filter to heterozygous markers with `filter_geno(hets = "include")`,
#' fix het-bearing loci with `qtn =`, or simulate dominance on an outbred or
#' F2-type population instead.
#'
#' Under `architecture = "pleiotropy"` the dominance effects are drawn with the
#' same covariance as the additive layer (DECISION-023): shared loci carry
#' correlated effects and trait-specific loci independent ones, split by `pi`,
#' so the dominance component targets `cor` -- its realized correlation
#' converges to `cor` as the loci and individuals grow, for loci in approximate
#' linkage equilibrium; see [simulate_phenotype()] `cor` --
#' rather than every trait getting one identical effect series. With `same_as_add = TRUE` the
#' additive layer's shared and trait-specific loci are reused. Each locus's
#' effect is scaled by the realized standard deviation of its heterozygote
#' indicator. Fixing loci with `qtn =` works under every architecture: under
#' "pleiotropy" every locus affects every trait (the same loci for each trait),
#' with an explicit `effect` setting the effects when no correlation is controlled;
#' under "ld" dominance reuses the additive layer's linked loci
#' (`same_as_add = TRUE`) or takes `qtn = list(trait1_loci, trait2_loci)` of
#' disjoint linked pairs (no marker may already be causal for the other trait in
#' an earlier layer), since a fresh random draw could make a marker causal for
#' both traits across layers.
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
#'
#' # Give the dominance effects yourself: a geometric base, or one value per QTN.
#' simulate_phenotype(SNP55K_maize282_maf04, seed = 2) |>
#'   additive(prop = 0.4, n_qtn = 5) |>
#'   dominance(prop = 0.1, same_as_add = FALSE, n_qtn = 3,
#'             effect = c(0.6, 0.3, 0.1))
dominance <- function(sim, prop = NULL, same_as_add = TRUE, n_qtn = NULL,
                      qtn = NULL, dist = "geometric", effect = NULL) {
  .check_sim(sim)
  if (isTRUE(sim$frozen)) return(.frozen_layer_ignored(sim, "dominance"))
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
  pleio <- identical(sim$architecture, "pleiotropy") && sim$n_traits > 1
  user_qtn <- .resolve_qtn_arg(sim, qtn, "dominance", equal = !pleio)
  # Fixed loci keep each architecture's construction (DECISION-043): under
  # "pleiotropy" every locus affects every trait (partial pleiotropy is
  # complex_phenotypes()), with the effects set by `effect` / `dist` when no
  # correlation is controlled and by the correlated draw (DECISION-023) otherwise;
  # under "ld" the loci must be disjoint linked pairs, as in additive().
  series_path <- pleio && !.pleio_controlled(sim) &&
    (!is.null(effect) || !identical(dist, "geometric"))
  if (!is.null(user_qtn) && pleio) {
    .pleio_check_user_units(sim, .expand_prop(prop, sim$n_traits), user_qtn,
                            "dominance", FALSE)
    user_qtn <- .pleio_user_layout(user_qtn, "dominance")$q
  }
  if (!is.null(user_qtn) && identical(sim$architecture, "ld")) {
    .ld_check_cross_layer(sim, user_qtn, "dominance")
    user_qtn <- .ld_user_pairs(sim, user_qtn, "dominance")
  }
  if (pleio && !series_path && (!is.null(effect) || !identical(dist, "geometric"))) {
    stop("dominance(): `effect` and non-default `dist` cannot be used under ",
         "architecture = \"pleiotropy\" while a correlation is controlled ",
         "(`cor`, `pi`, ...): the effects come from the multivariate draw that ",
         "controls `cor`.", call. = FALSE)
  }
  # A list gives one effect specification per trait (each a geometric base or an
  # explicit series, checked by .effect_series() when the layer is built).
  if (is.list(effect)) {
    if (length(effect) != sim$n_traits) {
      stop("dominance(effect=): a list must have one element per trait (",
           sim$n_traits, "); got ", length(effect), ".", call. = FALSE)
    }
    # NULL would silently fall back to the default series for that trait
    bad <- which(!vapply(effect, function(x) is.numeric(x) && length(x) > 0L,
                         logical(1)))
    if (length(bad)) {
      stop("dominance(effect=): every list element must be a numeric geometric ",
           "base or effect series; element(s) ", paste(bad, collapse = ", "),
           " are not.", call. = FALSE)
    }
  }
  eff_for <- function(t) if (is.list(effect)) effect[[t]] else effect
  # Under "ld" a fresh dominance draw knows nothing of the additive layer's loci,
  # so it could make a marker causal for both traits across layers; only reusing
  # the additive linked loci keeps every causal locus trait-specific.
  if (identical(sim$architecture, "ld") && !isTRUE(same_as_add) &&
      is.null(user_qtn)) {
    stop("dominance(same_as_add = FALSE) is not supported under architecture = ",
         "\"ld\" without `qtn`: a fresh draw could make a marker causal for ",
         "both traits across layers. Under \"ld\" dominance reuses the additive ",
         "layer's linked loci (same_as_add = TRUE), or takes disjoint linked ",
         "pairs through `qtn = list(trait1_loci, trait2_loci)`.", call. = FALSE)
  }

  add_layer <- .last_layer_of_type(sim, "additive")
  # User-supplied loci win over reusing the additive layer's: record that.
  if (!is.null(user_qtn)) {
    if (!is.null(n_qtn)) {
      warning("dominance(): `n_qtn` is ignored because `qtn` fixes the loci (",
              max(lengths(user_qtn)), " QTNs).", call. = FALSE)
    }
    same_as_add <- FALSE
  } else if (isTRUE(same_as_add) && !is.null(n_qtn)) {
    warning("dominance(): `n_qtn` is ignored because same_as_add = TRUE reuses ",
            "the additive layer's QTNs; set same_as_add = FALSE to draw `n_qtn` ",
            "fresh dominance loci.", call. = FALSE)
  }
  if (isTRUE(same_as_add) && is.null(user_qtn) && is.null(add_layer)) {
    stop("dominance(same_as_add = TRUE) requires a prior additive() layer.",
         call. = FALSE)
  }

  if (!is.null(user_qtn)) {
    nq <- max(lengths(user_qtn))
  } else if (isTRUE(same_as_add)) {
    nq <- add_layer$n_qtn
  } else {
    nq <- .resolve_n_qtn(sim, n_qtn, "dominance")
  }

  # With same_as_add the locus count is the additive layer's; say so in the
  # length error of an explicit `effect` series.
  count_lab <- if (isTRUE(same_as_add) && is.null(user_qtn))
    "n_qtn (set by the additive layer)" else "n_qtn"

  build <- function(rep_seed, rep = 0L) {
    if (!is.null(user_qtn)) {
      q <- user_qtn
    } else if (isTRUE(same_as_add)) {
      q <- .rep_qtn(add_layer, rep)      # follow additive's per-rep loci
    } else if (pleio) {
      q <- NULL                          # drawn with the correlated effects
    } else {
      q <- .draw_qtn(sim, nq, rep_seed)
    }
    if (pleio && !series_path) {
      # DECISION-023: dominance effects target `cor` across traits (shared
      # units correlated, trait-specific independent), instead of one identical
      # series per trait, which forced a dominance correlation of ~ +1.
      return(.pleio_nonadditive_draw(sim, .expand_prop(prop, sim$n_traits),
                                     rep_seed, "dominance", q = q,
                                     n_units = nq,
                                     shared = attr(q, "pleio_shared")))
    }
    e <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(nq, dist, eff_for(t),
                                           count = count_lab))
    list(qtn = q, effect = e)
  }

  fixed <- !is.null(user_qtn) ||
    (isTRUE(same_as_add) && is.null(add_layer$qtn_reps))
  drawn <- .draw_layer(sim, "dominance", occ, build, fixed = fixed)
  # A dominance deviation is a heterozygote effect, so it is identically zero at
  # a locus with no heterozygotes. Rather than silently substitute loci, fail
  # with a clear message so the user knows the (near-inbred) data has none.
  # every vary_qtn replication's loci, not just the canonical draw; traits with
  # prop = 0 contribute nothing, so their loci are not judged
  active <- .expand_prop(prop, sim$n_traits) > 0
  qsets <- lapply(c(list(drawn$qtn), drawn$qtn_reps), function(q) q[active])
  if (any(vapply(qsets, function(q) .dom_hetless(sim, q), logical(1)))) {
    stop("dominance(): the selected loci have no heterozygous individuals (or ",
         "are heterozygous in every individual), so the heterozygote indicator ",
         "does not vary and a dominance deviation cannot be ",
         "simulated -- the genotype is (near-)inbred (or an F1) at those loci. Pre-filter ",
         "to heterozygous markers with filter_geno(hets = \"include\"), choose ",
         "loci with heterozygotes via qtn =, or use an outbred / F2 population.",
         call. = FALSE)
  } else if (any(vapply(qsets, function(q) .dom_partial_hetless(sim, q),
                        logical(1)))) {
    warning("dominance(): some (but not all) selected loci have no heterozygous ",
            "individuals (or are all heterozygous), so those loci contribute ",
            "nothing and the remaining ",
            "het-bearing loci absorb the layer's `prop` -- the realized ",
            "dominance rests on fewer loci than requested. Choose het-bearing ",
            "loci (qtn =) or pre-filter with filter_geno(hets = \"include\") if ",
            "that is not intended.", call. = FALSE)
  }
  if (pleio && !series_path) {
    .pleio_total_cor_check(sim, .expand_prop(prop, sim$n_traits),
                           drawn$target_cor, drawn$target_cor_reps)
  }
  layer <- list(type = "dominance", prop = prop, n_qtn = nq, dist = dist,
                same_as_add = same_as_add,
                qtn = drawn$qtn, effect = drawn$effect)
  layer$target_cor <- drawn$target_cor
  layer$target_cor_reps <- drawn$target_cor_reps
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
#' @param qtn optional user-supplied interacting sets (marker names or column
#'   indices): a matrix with one row per set and `interaction` columns, or a
#'   vector whose length is a multiple of `interaction` (filled by row). A matrix
#'   or vector is used for **every** trait. A **list** must have one element per
#'   trait (as in [additive()]), each element that trait's sets (a matrix or
#'   vector as above), and every trait must get the same number of sets -- a list
#'   of length `n_traits` is never read as "one set per element". Under
#'   `architecture = "pleiotropy"` every set affects every trait (the same sets
#'   for each trait; sets match as ordered tuples; sets for only some traits are
#'   partial pleiotropy, [complex_phenotypes()]). Not available under
#'   `architecture = "ld"`.
#' @param effect optional geometric base (scalar) or explicit effect series
#'   (length `n_pairs`, one value per interacting set), used for every trait;
#'   or a length-`n_traits` list of these, one per trait (the counterpart of
#'   v1 `epi_effect`). Under `architecture = "pleiotropy"` (multi-trait) accepted
#'   only when no correlation is controlled (no `cor`, `pi`, ...), as in
#'   [additive()].
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
#'
#' Under `architecture = "pleiotropy"` the interacting sets are split by `pi`
#' into shared sets (common to every trait, with correlated effects) and
#' trait-specific sets (independent effects), drawn with the same covariance as
#' the additive layer, so the epistatic component targets `cor` (its realized
#' correlation converges as the sets and individuals grow, for loci in
#' approximate linkage equilibrium; DECISION-023). Each set's effect is scaled by the realized standard
#' deviation of its centered-product design column; `effect =` is rejected there
#' while a correlation is controlled (without `cor` / `pi` it sets the effects, and
#' fixed `qtn =` sets are accepted). `epistasis()` is not supported under
#' `architecture = "ld"`: that architecture makes each trait's causal loci
#' distinct markers in LD within the r2 window, and epistatic sets have no such
#' linked construction -- they would be unlinked and could reuse a marker
#' another layer made causal for the other trait.
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
  if (isTRUE(sim$frozen)) return(.frozen_layer_ignored(sim, "epistasis"))
  .require_markers(sim, "epistasis")
  interaction <- .validate_count(interaction, "interaction", minimum = 2L)
  prop <- .resolve_prop(sim, prop, "epistasis")
  pleio <- identical(sim$architecture, "pleiotropy") && sim$n_traits > 1
  user_pairs <- .resolve_epi_qtn(sim, qtn, interaction, equal = !pleio)
  if (!is.null(user_pairs)) interaction <- ncol(user_pairs[[1]])
  itype <- .resolve_interaction_type(interaction_type, interaction)
  if (!is.null(user_pairs) && !is.null(n_pairs)) {
    warning("epistasis(): `n_pairs` is ignored because `qtn` fixes the sets (",
            max(vapply(user_pairs, nrow, 0L)), " sets).", call. = FALSE)
  }
  np <- if (!is.null(user_pairs)) max(vapply(user_pairs, nrow, 0L)) else
    .resolve_n_qtn(sim, n_pairs, "epistasis", arg = "n_pairs")
  occ <- .type_occurrence(sim, "epistasis")
  # "ld" places each trait's causal loci as DISTINCT markers in LD within the r2
  # window (DECISION-014). An epistatic set has no linked-distinct analog here:
  # sets drawn elsewhere are unlinked (outside the r2 window) and can reuse a
  # marker that another layer made causal for the other trait -- pleiotropy
  # inside an "ld" model. SPEC-nonadditive-correlation 5.3: restrict it.
  if (identical(sim$architecture, "ld")) {
    stop("epistasis() is not supported under architecture = \"ld\": that ",
         "architecture makes each trait's causal loci distinct markers in LD ",
         "within the r2 window (DECISION-014), and epistatic sets have no such ",
         "linked construction -- they would be unlinked and could reuse another ",
         "layer's causal marker for the other trait. Use architecture = ",
         "\"independent\" or \"pleiotropy\" for epistasis (DECISION-023).",
         call. = FALSE)
  }
  # Fixed sets under "pleiotropy" affect every trait (a set for only some traits
  # is partial pleiotropy: complex_phenotypes(), DECISION-043); effects are set by
  # `effect` / `dist` when no correlation is controlled, else by the correlated
  # draw (DECISION-023). A single matrix is used for every trait.
  series_path <- pleio && !.pleio_controlled(sim) &&
    (!is.null(effect) || !identical(dist, "geometric"))
  if (!is.null(user_pairs) && pleio) {
    .pleio_check_user_units(sim, .expand_prop(prop, sim$n_traits), user_pairs,
                            "epistasis", TRUE)
  }
  if (pleio && !series_path && (!is.null(effect) || !identical(dist, "geometric"))) {
    stop("epistasis(): `effect` and non-default `dist` cannot be used under ",
         "architecture = \"pleiotropy\" while a correlation is controlled ",
         "(`cor`, `pi`, ...): the effects come from the multivariate draw that ",
         "controls `cor`.", call. = FALSE)
  }
  # A list gives one effect specification per trait (each a geometric base or an
  # explicit series, checked by .effect_series() when the layer is built), as
  # in additive() and dominance().
  if (is.list(effect)) {
    if (length(effect) != sim$n_traits) {
      stop("epistasis(effect=): a list must have one element per trait (",
           sim$n_traits, "); got ", length(effect), ".", call. = FALSE)
    }
    # NULL would silently fall back to the default series for that trait
    bad <- which(!vapply(effect, function(x) is.numeric(x) && length(x) > 0L,
                         logical(1)))
    if (length(bad)) {
      stop("epistasis(effect=): every list element must be a numeric geometric ",
           "base or effect series; element(s) ", paste(bad, collapse = ", "),
           " are not.", call. = FALSE)
    }
  }
  eff_for <- function(t) if (is.list(effect)) effect[[t]] else effect

  build <- function(rep_seed, rep = 0L) {
    if (pleio && !series_path) {
      # DECISION-023: epistatic effects target `cor` across traits (shared
      # sets correlated, trait-specific sets independent), instead of one
      # identical series on shared sets, which forced a correlation of ~ +1.
      return(.pleio_nonadditive_draw(sim, .expand_prop(prop, sim$n_traits),
                                     rep_seed, "epistasis", q = user_pairs,
                                     n_units = np,
                                     interaction = interaction, itype = itype,
                                     arg = "n_pairs"))
    }
    q <- if (!is.null(user_pairs)) user_pairs else
      .draw_qtn_pairs(sim, np, interaction, rep_seed)
    e <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(np, dist, eff_for(t),
                                           count = "n_pairs"))
    list(qtn = q, effect = e)
  }

  drawn <- .draw_layer(sim, "epistasis", occ, build,
                       fixed = !is.null(user_pairs))
  # A "d" interaction position acts on heterozygotes, so a hetless locus there
  # makes the pair identically zero and it silently drops out while the rest
  # absorb `prop`. Warn (as dominance() does for its partial case) rather than
  # realize fewer effective pairs than requested without notice.
  # (traits with prop = 0 contribute nothing, so their sets are not judged)
  active <- .expand_prop(prop, sim$n_traits) > 0
  dead <- vapply(c(list(drawn$qtn), drawn$qtn_reps),
                 function(q) .epi_dead_status(sim, q[active], itype), 0L)
  if (any(dead == 2L)) {
    stop("epistasis(): every interacting set has a \"d\" position on a locus with ",
         "no heterozygous individuals (or heterozygous in every individual), so ",
         "every interaction term is constant and the layer cannot realize its ",
         "`prop` -- the genotype is (near-)inbred (or an F1) at those loci. ",
         "Choose het-bearing loci (qtn =), pre-filter with filter_geno(hets = ",
         "\"include\"), set the interaction_type to \"a\", or use an outbred / ",
         "F2 population.", call. = FALSE)
  } else if (any(dead == 1L)) {
    warning("epistasis(): a \"d\" interaction position sits on a locus with no ",
            "heterozygous individuals (or all heterozygous), so that interaction ",
            "term is identically zero and its pair contributes nothing while the ",
            "remaining pairs absorb `prop` -- the genotype is (near-)inbred at ",
            "that locus. Choose ",
            "het-bearing loci (qtn =), pre-filter with filter_geno(hets = ",
            "\"include\"), or set that position's interaction_type to \"a\" if ",
            "that is not intended.", call. = FALSE)
  }
  if (pleio && !series_path) {
    .pleio_total_cor_check(sim, .expand_prop(prop, sim$n_traits),
                           drawn$target_cor, drawn$target_cor_reps)
  }
  layer <- list(type = "epistasis", prop = prop, n_pairs = np,
                interaction = interaction, interaction_type = itype,
                dist = dist,
                qtn = drawn$qtn, effect = drawn$effect)
  layer$target_cor <- drawn$target_cor
  layer$target_cor_reps <- drawn$target_cor_reps
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
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`). Supplying
#'   `qtn` fixes the loci instead, so `same_as_add` is then treated as `FALSE`.
#' @param n_qtn number of vQTL loci for a fresh draw (`same_as_add = FALSE`);
#'   ignored, with a warning, when the loci are reused or fixed by `qtn`.
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
  if (isTRUE(sim$frozen)) return(.frozen_layer_ignored(sim, "vqtl"))
  .require_markers(sim, "vqtl")
  .validate_flag(same_as_add, "same_as_add")
  .cite_vqtl()
  prop <- .resolve_prop(sim, prop, "vqtl")
  occ <- .type_occurrence(sim, "vqtl")
  pleio <- identical(sim$architecture, "pleiotropy") && sim$n_traits > 1
  user_qtn <- .resolve_qtn_arg(sim, qtn, "vqtl", equal = !pleio)
  # passed vQTL loci follow the architecture's rules too (DECISION-043):
  # pleiotropy = the same loci for every trait; ld = disjoint linked pairs
  if (!is.null(user_qtn) && pleio) {
    user_qtn <- .pleio_user_layout(user_qtn, "vqtl")$q
  }
  if (!is.null(user_qtn) && identical(sim$architecture, "ld")) {
    .ld_check_cross_layer(sim, user_qtn, "vqtl")
    user_qtn <- .ld_user_pairs(sim, user_qtn, "vqtl")
  }

  add_layer <- .last_layer_of_type(sim, "additive")
  if (!is.null(user_qtn)) {
    if (!is.null(n_qtn)) {
      warning("vqtl(): `n_qtn` is ignored because `qtn` fixes the loci (",
              max(lengths(user_qtn)), " QTNs).", call. = FALSE)
    }
    same_as_add <- FALSE
  } else if (isTRUE(same_as_add) && !is.null(n_qtn)) {
    warning("vqtl(): `n_qtn` is ignored because same_as_add = TRUE reuses the ",
            "additive layer's QTNs; set same_as_add = FALSE to draw `n_qtn` ",
            "fresh vQTL loci.", call. = FALSE)
  }
  if (isTRUE(same_as_add) && is.null(user_qtn) && is.null(add_layer)) {
    stop("vqtl(same_as_add = TRUE) requires a prior additive() layer.",
         call. = FALSE)
  }
  if (!is.null(user_qtn)) {
    nq <- max(lengths(user_qtn))
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
    # One effect per RETAINED locus of each trait. Under architecture =
    # "pleiotropy" the additive layer drops loci that carry exactly zero effect
    # for a trait (pi_t = 0 or 1), so the per-trait locus counts can differ and
    # a common `nq` would leave the effects longer than the loci.
    e <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(length(q[[t]]), dist))
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
    # A trait whose budget is already spent gets prop = 0; refuse only when no
    # trait has anything left (prop = c(0.3, 0) is legal).
    if (all(remaining <= 1e-10)) {
      stop("No heritability budget is left for the ", type, "() layer: h2 = ",
           paste(sprintf("%.3f", .expand_prop(h2, nt)), collapse = ", "),
           " is already fully allocated. Give the layers explicit `prop` ",
           "values that sum to h2.", .one_call_hint(sim), call. = FALSE)
    }
    remaining[remaining <= 1e-10] <- 0
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
  if (identical(sim$architecture, "complex")) {
    stop(.complex_terminal_msg(layer$type), call. = FALSE)
  }
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
    out <- list(qtn = canonical$qtn, effect = canonical$effect)
    out$target_cor <- canonical$target_cor      # pleiotropic non-additive only
    return(out)
  }
  reps <- vector("list", sim$n_reps)
  reps[[1L]] <- canonical
  for (r in seq_len(sim$n_reps)[-1L]) {
    reps[[r]] <- build(.layer_seed(sim$seed, paste0(type, "_rep", r), occ), r)
  }
  out <- list(
    qtn        = canonical$qtn,
    effect     = canonical$effect,
    qtn_reps   = lapply(reps, function(x) x$qtn),
    effect_reps = lapply(reps, function(x) x$effect)
  )
  out$target_cor <- canonical$target_cor
  if (!is.null(canonical$target_cor)) {
    out$target_cor_reps <- lapply(reps, function(x) x$target_cor)
  }
  out
}

#' Per-locus flag: does the heterozygote indicator vary across individuals?
#'
#' FALSE for a locus with no heterozygous individual **and** for one that is
#' heterozygous in every individual (constant indicator): in both cases the
#' centered heterozygote indicator is identically zero, so the locus can carry no
#' dominance deviation.
#' @keywords internal
#' @noRd
.het_varies <- function(sim, idx) {
  blk <- .geno_cols(sim, idx)
  nh <- colSums(blk == 0, na.rm = TRUE)
  nh > 0 & nh < nrow(blk)
}

#' TRUE when a reused QTN set cannot support any dominance deviation
#'
#' A dominance layer is a heterozygote effect, so a trait whose loci *all* have a
#' constant heterozygote indicator (no heterozygotes, or all heterozygous)
#' realizes exactly zero variance and cannot fill `prop` at all -- that is an
#' error. A set with *some* such loci can still realize `prop` from the rest, but
#' the dead loci are silently inert; see [.dom_partial_hetless()], which warns for
#' that case.
#' @keywords internal
#' @noRd
.dom_hetless <- function(sim, q) {
  any(vapply(q, function(idx) {
    if (is.null(idx) || length(idx) == 0L) {
      return(FALSE)
    }
    !any(.het_varies(sim, idx))
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
    hl <- !.het_varies(sim, idx)
    any(hl) && !all(hl)
  }, logical(1)))
}

#' TRUE when any orthogonal locus carrying a non-zero `d` has a constant het indicator
#'
#' (no heterozygotes, or heterozygous in every individual).
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
    any(!.het_varies(sim, idx[nz]))
  }, logical(1)))
}

#' Status of the "d"-position loci of a set of epistatic sets
#'
#' A `"d"` interaction position contributes a centered heterozygote indicator, so
#' a locus whose indicator is constant (no heterozygotes, or all heterozygous)
#' makes that design column identically zero and inerts the whole pair. Returns
#' 0 when no pair is affected, 2 when *every* pair (of some trait) is dead --
#' the layer cannot realize `prop` -- and 1 when only some are, the remaining
#' pairs then absorbing `prop`. `qtn` is a per-trait list of `n_pairs x
#' interaction` index matrices; `itype` the length-`interaction` "a"/"d" vector.
#' @keywords internal
#' @noRd
.epi_dead_status <- function(sim, qtn, itype) {
  d_pos <- which(itype == "d")
  if (length(d_pos) == 0L) {
    return(0L)
  }
  status <- vapply(qtn, function(mat) {
    if (is.null(mat) || length(mat) == 0L) {
      return(0L)
    }
    loci <- unique(as.integer(mat[, d_pos, drop = FALSE]))
    ok <- stats::setNames(.het_varies(sim, loci), loci)
    dead_pair <- apply(mat[, d_pos, drop = FALSE], 1L,
                       function(r) any(!ok[as.character(r)]))
    if (all(dead_pair)) 2L else if (any(dead_pair)) 1L else 0L
  }, 0L)
  if (length(status)) max(status) else 0L
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
.resolve_qtn_arg <- function(sim, qtn, type, equal = TRUE) {
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
    if (equal && length(unique(lens)) != 1L) {
      stop(type, "(qtn=): every trait must get the same number of QTNs; got ",
           paste(lens, collapse = ", "), ".", call. = FALSE)
    }
    return(.warn_monomorphic(sim, out, type))
  }
  idx <- as_idx(qtn)
  .warn_monomorphic(sim, rep(list(idx), nt), type)
}

#' Warn about user-chosen QTNs that cannot carry an effect
#'
#' Randomly drawn QTNs come from `.candidate_markers()` (polymorphic, not
#' heterozygous in every individual); a locus passed through `qtn =` is not
#' screened, but a monomorphic one has a constant design column and so no
#' variance to explain. Say so instead of letting it silently contribute nothing.
#' @keywords internal
#' @noRd
.warn_monomorphic <- function(sim, idx_list, type) {
  idx <- unique(unlist(idx_list, use.names = FALSE))
  bad <- idx[!(idx %in% .candidate_markers(sim))]
  if (length(bad)) {
    warning(type, "(qtn=): ", length(bad), " marker(s) are monomorphic or ",
            "heterozygous in every individual (",
            paste(utils::head(sim$map$snp[bad], 5), collapse = ", "),
            if (length(bad) > 5) ", ..." else "",
            "), so their genotype column is constant and they carry no ",
            "variance; a randomly drawn QTN is never one of these.",
            call. = FALSE)
  }
  idx_list
}

#' Resolve a user-supplied epistatic `qtn` argument to per-trait set matrices
#'
#' A matrix (one row per set, `interaction` columns) or a vector (length a
#' multiple of `interaction`, filled by row) is shared by every trait. A **list**
#' must have one element per trait, each element that trait's sets (matrix or
#' vector), and every trait must get the same number of sets -- exactly as for
#' `additive(qtn = list(...))`. (Earlier versions read a list as "one set per
#' element" and replicated the resulting matrix to every trait, which silently
#' gave identical epistasis and a genetic correlation of 1.)
#' @keywords internal
#' @noRd
.resolve_epi_qtn <- function(sim, qtn, interaction, equal = TRUE) {
  if (is.null(qtn)) {
    return(NULL)
  }
  nt <- sim$n_traits
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
  to_mat <- function(v) {
    if (is.matrix(v)) {
      matrix(to_idx(as.vector(v)), nrow = nrow(v))
    } else {
      if (length(v) %% interaction != 0L) {
        stop("epistasis(qtn=): the number of loci must be a multiple of ",
             "`interaction` (", interaction, ").", call. = FALSE)
      }
      matrix(to_idx(v), ncol = interaction, byrow = TRUE)
    }
  }
  if (is.list(qtn)) {
    if (length(qtn) != nt) {
      stop("epistasis(qtn=): a list must have one element per trait (", nt,
           "); got ", length(qtn), ". Each element is that trait's interacting ",
           "sets (a matrix with one row per set, or a vector); to give several ",
           "sets to every trait pass one matrix instead of a list of sets.",
           call. = FALSE)
    }
    mats <- lapply(qtn, to_mat)
    nr <- vapply(mats, nrow, 0L)
    if (equal && length(unique(nr)) != 1L) {
      stop("epistasis(qtn=): every trait must get the same number of sets; got ",
           paste(nr, collapse = ", "), ".", call. = FALSE)
    }
  } else {
    mats <- rep(list(to_mat(qtn)), nt)      # one matrix, shared by every trait
  }
  for (m in mats) {
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
  }
  .warn_monomorphic(sim, lapply(mats, as.vector), "epistasis")
  mats
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
#' argument of additive()). Coupling returns the effects unchanged. The
#' haplotype-derived phase of the `"ld"` architecture is `.apply_ld_phase()`,
#' applied after this one.
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
  if (identical(sim$architecture, "complex")) {
    stop(.complex_terminal_msg(fn), call. = FALSE)
  }
  if (is.null(sim$n_markers) || sim$n_markers < 1L) {
    stop(fn, "() needs genotypes, but this phenotype was built from expression ",
         "alone (no `geno`). Use transcriptome() layers, or rebuild with `geno`.",
         call. = FALSE)
  }
}

#' Message for a layer added to a complex_phenotypes() result
#' @keywords internal
#' @noRd
.complex_terminal_msg <- function(fn) {
  paste0(fn, "(): a complex_phenotypes() result is terminal -- its genetic ",
         "value is the rescaled sum of its inputs and it has no layers of its ",
         "own, so a layer added afterwards would be silently ignored. Add the ",
         "layer to one of the input models and combine again.")
}
