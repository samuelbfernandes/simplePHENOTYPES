#' PleioArch pleiotropic QTN draw and effects (target genetic correlation)
#'
#' Implements the PleioArch algorithm for the
#' `architecture = "pleiotropy"` additive layer, generalized to any number of
#' traits. The same covariance construction is extended to the dominance and
#' epistasis layers by [.pleio_nonadditive_draw()] (DECISION-023). Pleiotropic
#' effects are drawn from a multivariate normal whose covariance targets the
#' genetic correlations `cor`; trait-specific effects are drawn from univariate
#' normals. Effects are scaled to the per-genotype scale by
#' `1 / sqrt(2 * MAF * (1 - MAF))` (the reference `scaleQTNEffects()` step).
#'
#' The covariance among pleiotropic effects is
#' `Sigma[i, i] = pi_i * V_i` and `Sigma[i, j] = cor_ij * sqrt(V_i * V_j)`.
#' Since trait-specific loci contribute `(1 - pi_i) * V_i` and are independent
#' across traits, the effect draw gives each trait genetic variance `V_i` and a
#' cross-trait covariance of `cor_ij * sqrt(V_i * V_j)` in expectation, all of it
#' from the shared loci. The layer is then rescaled to its `prop`, so what is
#' realized is a correlation `r` (and a covariance `r * sqrt(V_i * V_j)`). `r`
#' converges to `cor_ij` as the shared loci **and the individuals** grow
#' **provided their dosage columns are close to uncorrelated** (linkage
#' equilibrium among the causal loci) **and no locus keeps a non-vanishing share
#' of the variance** (major QTNs set by `n_pleio_major` / `prop_var_major` do, so
#' `r` then does not converge). `r` is a sample correlation over the n
#' individuals, so with n fixed its scatter levels off at that sampling spread
#' however many loci there are. With
#' few shared loci `r` is attenuated toward 0 on average by an amount that depends
#' on `cor`, `pi` and the designs (e.g. about 0.82 x `cor` with two shared loci at
#' `cor = 0.5`, `pi = 1`, uncorrelated designs). Strong LD among them can prevent
#' convergence altogether: under complete LD every realized `r` is +/-1 however
#' many loci, averaging `2 asin(cor) / pi` over replicates. One shared locus gives
#' +/-1, or a single noisy draw alongside trait-specific loci (DECISION-023). For two traits this
#' reduces exactly to the bivariate reference implementation.
#'
#' RNG stays in R. Returns per-trait QTN indices and effects so
#' the additive layer can carry trait-specific loci alongside the shared
#' pleiotropic ones.
#'
#' Provenance: the associated manuscript (Prado et al.) is in preparation and has
#' no public reference, so the authoritative definition of this algorithm is the
#' reference implementation in `context/PleioArch-main/Functions/`
#' (`simulateEffects.R`, `scaleQTNEffects.R`), against which this port was
#' checked (that folder is a development-only resource: it is not in every
#' checkout nor in the package tarball, so parity cannot be re-verified from the
#' package alone). The parity is of the *algorithm*: the port draws the shared
#' effects through a symmetric eigen square root of the covariance where the
#' reference uses `chol()`, so the covariance agrees but the realized values are
#' not bit-identical to the literal reference (no decision requires that).
#' The covariance construction and the `1/sqrt(2*MAF*(1-MAF))` scaling
#' documented here are self-contained; do not cite them to a published paper
#' until one exists.
#'
#' @param sim a `phenotype_sim` with `architecture = "pleiotropy"`.
#' @param n_qtn total QTNs per trait.
#' @param prop_vec per-trait additive variance proportions (length `n_traits`).
#' @param sub_seed deterministic sub-seed.
#' @param fixed optional user-supplied layout from [.pleio_user_layout()]
#'   (`list(shared, spec)`): the shared loci (causal for every trait) and each
#'   trait's own loci. The loci are then not sampled and the shared/specific
#'   partition is the layout's, so `n_qtn` is not used; the variance split
#'   between the two classes is still `pi` and the covariance construction,
#'   scaling and feasibility checks are unchanged. Counts may differ by trait.
#' @return list(qtn = <list len n_traits>, effect = <list len n_traits>).
#' @keywords internal
#' @noRd
.pleio_draw <- function(sim, n_qtn, prop_vec, sub_seed, fixed = NULL) {
  a  <- sim$arch_args
  nt <- sim$n_traits
  # Only when the user actually asked for correlation control: `cor` absent
  # means the engine runs at cor = 0 and no correlation is being controlled.
  if (!is.null(a[["cor"]])) {
    .cite_pleioarch()
  }
  R  <- .pleio_cor_matrix(sim)
  pi_vec <- .pleio_pi_vector(sim)
  # Default: no major QTN. Concentrating variance in one locus makes the
  # realized correlation hinge on that single effect draw, so it no longer
  # converges on `cor` as n_qtn grows -- which defeats the reason for using
  # this engine. A major locus is available, but it must be asked for.
  n_major  <- if (is.null(a[["n_pleio_major"]])) 0 else a[["n_pleio_major"]]
  propMaj  <- if (is.null(a[["prop_var_major"]])) 0 else a[["prop_var_major"]]
  n_major <- .validate_count(n_major, "n_pleio_major", minimum = 0L)
  propMaj <- .validate_proportion(propMaj, "prop_var_major", 1L)

  vg <- prop_vec
  .pleio_check_zero_var(R, vg, "additive")

  # Variance budget (phenotypic variance assumed 1; vg are the target props).
  # Off-diagonals carry the full genetic covariance because trait-specific loci
  # are independent; diagonals carry only the pleiotropic share.
  sigma <- outer(sqrt(vg), sqrt(vg)) * R
  diag(sigma) <- pi_vec * vg

  .check_pleio_feasible(sigma, R, pi_vec)

  # QTN-count partition: pleiotropic (shared) vs trait-specific. A user layout
  # fixes both, and the loci; the guards then apply to the supplied counts.
  if (is.null(fixed)) {
    part <- .pleio_partition(n_qtn, pi_vec, R, vg = vg)
    pleio_n <- part$pleio_n
    spec_n  <- rep(part$spec_n, nt)
  } else {
    # only informative loci (polymorphic in the simulated individuals) carry
    # the covariance: constant columns drop out of the realized value, so the
    # guards count them, while the effect draws below use every supplied locus
    pleio_n <- length(fixed$shared)
    n_inf   <- sum(fixed$shared %in% .candidate_markers(sim))
    spec_n  <- vapply(fixed$spec, length, 0L)
    if (n_inf < 1L && any(vg > 0)) {
      stop("architecture = \"pleiotropy\": none of the loci in `qtn` varies in ",
           "the simulated individuals (monomorphic or heterozygous in every ",
           "one), so no shared locus can carry any variance or correlation.",
           call. = FALSE)
    }
    .pleio_check_layout(R, pi_vec, vg, n_inf, spec_n)
  }
  if (n_major > pleio_n) {
    stop("`n_pleio_major` (", n_major, ") cannot exceed the number of shared ",
         "pleiotropic QTNs ", if (is.null(fixed)) "implied by `pi` (" else
           "supplied in `qtn` (", pleio_n, ").", call. = FALSE)
  }
  n_minor <- pleio_n - n_major
  # With a major/minor split, prop_var_major < 1 allocates (1 - prop_var_major)
  # of the pleiotropic variance to minor loci. If there are none, that share
  # would silently vanish (attenuating the realized correlation); reject rather
  # than lose it.
  if (n_major > 0L && n_minor == 0L && propMaj < 1) {
    stop("`prop_var_major` (", propMaj, ") < 1 assigns ", round(1 - propMaj, 3),
         " of the pleiotropic variance to minor loci, but none remain ",
         "(n_pleio_major equals the shared-QTN count: ",
         pleio_n, "). Increase `n_qtn`, lower `n_pleio_major`, or set ",
         "`prop_var_major = 1`.", call. = FALSE)
  }

  # Keep the major/minor split coherent: variance allocated to a major class
  # with no loci in it would simply vanish from the budget.
  if (n_major == 0 && propMaj > 0) {
    stop("`prop_var_major` is positive but `n_pleio_major` is zero.",
         call. = FALSE)
  }
  if (n_major > 0 && propMaj <= 0) {
    stop("`n_pleio_major` is positive but `prop_var_major` is zero.",
         call. = FALSE)
  }

  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }
  if (is.null(fixed)) {
    cand <- .candidate_markers(sim)
    need <- pleio_n + nt * spec_n[1L]
    if (need > length(cand)) {
      stop("architecture = \"pleiotropy\": the shared + trait-specific QTNs need ",
           need, " distinct markers but only ", length(cand), " candidates exist. ",
           "Lower n_qtn, or raise pi (shared units are reused by every trait, so ",
           "a larger pi needs fewer distinct markers).", call. = FALSE)
    }
    drawn <- sample(cand, need, replace = FALSE)
    pleio_idx <- drawn[seq_len(pleio_n)]
    spec_idx <- lapply(seq_len(nt), function(t) {
      if (spec_n[t] <= 0) {
        return(integer(0))
      }
      drawn[pleio_n + (t - 1L) * spec_n[t] + seq_len(spec_n[t])]
    })
  } else {
    pleio_idx <- fixed$shared
    spec_idx <- fixed$spec
  }

  eff_major <- .draw_mvnorm(n_major, sigma * propMaj)
  eff_minor <- .draw_mvnorm(n_minor, sigma * (1 - propMaj))
  eff_spec <- lapply(seq_len(nt),
                     function(t) .draw_univariate(spec_n[t],
                                                  (1 - pi_vec[t]) * vg[t]))

  maf <- sim$maf
  scale_idx <- function(idx) {
    m <- maf[idx]
    s <- sqrt(2 * m * (1 - m))
    s[!is.finite(s) | s == 0] <- 1
    1 / s
  }
  pleio_scale <- if (pleio_n > 0) scale_idx(pleio_idx) else numeric(0)

  qtn <- vector("list", nt)
  effect <- vector("list", nt)
  for (t in seq_len(nt)) {
    eff_pleio <- c(eff_major[, t], eff_minor[, t]) * pleio_scale
    eff_t <- if (spec_n[t] > 0) eff_spec[[t]] * scale_idx(spec_idx[[t]]) else
      numeric(0)
    qtn_t <- c(pleio_idx, spec_idx[[t]])
    effect_t <- c(eff_pleio, eff_t)
    # A trait with no trait-specific variance (pi_t = 1) or no shared variance
    # (pi_t = 0) is assigned loci that carry exactly zero effect. They were
    # drawn (RNG order is unchanged) but are not causal for that trait, so they
    # are not reported as its QTNs.
    keep <- effect_t != 0
    qtn[[t]] <- qtn_t[keep]
    effect[[t]] <- effect_t[keep]
  }
  # The per-trait retained sets can differ in length and, when a trait has no
  # shared variance (pi_t = 0), no longer intersect in the shared loci. Record
  # the shared loci so a layer that REUSES these loci (dominance) still knows
  # which are shared; a vQTL reuse sizes its effects per trait from the
  # retained sets.
  attr(qtn, "pleio_shared") <- as.integer(pleio_idx)

  list(qtn = qtn, effect = effect)
}

#' Reject or warn on a correlation against a zero-variance trait
#'
#' A correlation with a zero-variance trait is undefined: cov = cor*sqrt(Vi*Vj)
#' is 0 and the realized correlation is 0/0 = NA, so cor cannot be attained.
#' Reject a nonzero requested off-diagonal against any trait whose component
#' proportion is zero, rather than silently returning NA correlations; warn for
#' cor = 0, which is undefined too (reported as NA, not 0).
#' @param component layer named in the messages ("additive", "dominance", ...).
#' @keywords internal
#' @noRd
.pleio_check_zero_var <- function(R, vg, component = "additive") {
  if (!any(vg <= 0)) {
    return(invisible())
  }
  nt_ <- length(vg)
  for (i in seq_len(nt_ - 1L)) {
    for (j in seq(i + 1L, nt_)) {
      if (vg[i] <= 0 || vg[j] <= 0) {
        zero_t <- if (vg[i] <= 0) i else j
        other  <- if (zero_t == i) j else i
        if (abs(R[i, j]) > 0) {
          stop("architecture = \"pleiotropy\": a nonzero `cor` (", signif(R[i, j], 3),
               ") was requested between trait ", i, " and trait ", j,
               ", but trait ", zero_t, " has zero ", component, " ",
               "variance (prop = 0), so that correlation is undefined and cannot ",
               "be realized. Give every correlated trait a positive ", component,
               " `prop`.", call. = FALSE)
        } else {
          warning("architecture = \"pleiotropy\": trait ", zero_t, " has zero ",
                  component, " variance (prop = 0), so its genetic correlation with ",
                  "trait ", other, " is undefined and will be reported as NA (a ",
                  "requested cor = 0 cannot be realized as 0 here). Give every ",
                  "correlated trait a positive ", component, " `prop`.",
                  call. = FALSE)
        }
      }
    }
  }
  invisible()
}

#' Shared (pleiotropic) vs trait-specific unit counts, with their guards
#'
#' `pi` defines the shared/trait-specific variance partition whether or not a
#' correlation is requested, so feasibility is judged against `pi`, not `cor`.
#' Guards: a single shared unit gives exactly +/-1 when there are no
#' trait-specific units and one noisy draw otherwise (warn); a positive `pi`
#' that rounds the shared class to zero, or `pi < 1` with no trait-specific
#' units, would silently drop variance (error). `arg`/`unit`/`units`/`place`/
#' `one` only word the messages (defaults reproduce the additive layer's).
#' @keywords internal
#' @noRd
.pleio_partition <- function(n, pi_vec, R, arg = "n_qtn", unit = "QTN",
                             units = "QTNs", place = "loci", one = "locus",
                             vg = NULL) {
  pleio_n <- max(0L, round(mean(pi_vec) * n))
  pleio_n <- min(pleio_n, n)
  # One shared unit carries the whole cross-trait covariance. A trait with no
  # trait-specific variance (no trait-specific units, or pi_t = 1, which gives
  # them zero variance) is a scalar multiple of that one design column, so a
  # PAIR of such traits correlates exactly +/-1; any other pair hinges on a
  # single effect pair and is one noisy draw around (possibly opposite in sign
  # to) `cor`. Any target strictly inside (-1, 1) is affected -- cor = 0 too.
  if (pleio_n == 1L && any(abs(R[upper.tri(R)]) < 1)) {
    no_spec <- (n - pleio_n == 0L) | pi_vec >= 1
    consequence <- .pleio_single_unit_consequence(R, no_spec, one, vg)
  } else {
    consequence <- ""
  }
  if (nzchar(consequence)) {
    warning("architecture = \"pleiotropy\": only one shared (pleiotropic) ", unit,
            " results from ", arg, " = ", n, ", pi = ",
            paste(round(pi_vec, 3), collapse = ", "), "; ", consequence,
            ". Increase ", arg, " (or pi) so several shared ", place,
            " carry the covariance.", call. = FALSE)
  }
  if (mean(pi_vec) > 0 && n > 0L && pleio_n < 1L) {
    stop("architecture = \"pleiotropy\": ", arg, " = ", n, " and pi = ",
         paste(round(pi_vec, 3), collapse = ", "),
         " round the shared (pleiotropic) class to zero, so the requested ",
         "shared variance", if (any(abs(R[upper.tri(R)]) > 0))
           " and correlation" else "", " cannot be represented. Increase ",
         arg, " (or pi).", call. = FALSE)
  }
  spec_n <- n - pleio_n
  if (spec_n < 1L && any(pi_vec < 1)) {
    stop("architecture = \"pleiotropy\": ", arg, " = ", n, " leaves no ",
         "trait-specific ", units, ", but pi = ",
         paste(round(pi_vec, 3), collapse = ", "),
         " (< 1) requests trait-specific variance. Increase ", arg, " so both ",
         "shared and trait-specific ", place, " fit.", call. = FALSE)
  }
  list(pleio_n = pleio_n, spec_n = spec_n)
}

#' Guard for a user-supplied pleiotropic layout (all loci shared)
#'
#' The counterpart of the guards in [.pleio_partition()] when the loci come from
#' `qtn =`: every supplied locus is shared, so a request for trait-specific
#' variance (`pi_t < 1`) has nothing to live on. `pi` only matters when a
#' correlation is being controlled; the default (`pi = 1`) puts all the variance
#' on the shared loci.
#' @param pleio_n number of shared loci; `spec_n` per-trait specific counts (0).
#' @keywords internal
#' @noRd
.pleio_check_layout <- function(R, pi_vec, vg, pleio_n, spec_n,
                                arg = "qtn", unit = "QTN", units = "QTNs",
                                place = "loci", one = "locus") {
  nt <- length(pi_vec)
  need_spec <- which(pi_vec < 1 & vg > 0 & spec_n < 1L)
  if (length(need_spec)) {
    stop("architecture = \"pleiotropy\": `", arg, "` fixes the ", place,
         ", and each one affects every trait, but pi = ",
         paste(round(pi_vec[need_spec], 3), collapse = ", "),
         " (< 1) requests trait-specific variance. `pi` is only needed to ",
         "control a correlation; use pi = 1 (the default). Trait-specific ",
         units, " are partial pleiotropy: combine models with ",
         "complex_phenotypes().", call. = FALSE)
  }
  if (pleio_n == 1L && nt > 1L && any(abs(R[upper.tri(R)]) < 1)) {
    consequence <- .pleio_single_unit_consequence(R, spec_n < 1L | pi_vec >= 1,
                                                  one, vg)
    if (nzchar(consequence)) {
      warning("architecture = \"pleiotropy\": only one shared (pleiotropic) ",
              unit, " is supplied in `", arg, "`; ", consequence, ". Supply ",
              "several shared ", place, " to carry the covariance.",
              call. = FALSE)
    }
  }
  invisible()
}

#' Validate a user `qtn` under pleiotropy: every locus is shared by every trait
#'
#' Under `architecture = "pleiotropy"` a causal locus affects **every** trait
#' (that is what pleiotropy means here), so a user layout must list the same loci
#' for every trait: a single vector, or a list whose elements hold the same loci
#' (order may differ). A locus that affects only some of the traits is partial
#' pleiotropy, which is built by combining models with [complex_phenotypes()], not
#' inside one pleiotropy model; it is rejected with that pointer instead of being
#' silently treated as independent effects.
#'
#' @param user_qtn per-trait list of marker indices.
#' @param type layer name for messages.
#' @return `list(shared, spec, q)`: the shared loci (trait 1's order), empty
#'   trait-specific sets, and `q`, the per-trait lists carrying the
#'   `"pleio_shared"` attribute the correlated draws expect.
#' @keywords internal
#' @noRd
.pleio_user_layout <- function(user_qtn, type = "additive") {
  nt <- length(user_qtn)
  shared <- user_qtn[[1L]]
  same <- vapply(user_qtn, function(v) setequal(v, shared), logical(1))
  if (!all(same)) {
    stop(type, "(qtn=): under architecture = \"pleiotropy\" every causal locus ",
         "affects every trait, so `qtn` must list the same loci for each trait ",
         "(a single vector does this). Loci that affect only some of the traits ",
         "are partial pleiotropy: build one model per group of traits and ",
         "combine them with complex_phenotypes().", call. = FALSE)
  }
  # each trait keeps its own order: positional effects (effect =) belong to it
  q <- lapply(user_qtn, as.integer)
  attr(q, "pleio_shared") <- as.integer(shared)
  list(shared = as.integer(shared), spec = rep(list(integer(0)), nt), q = q)
}

#' Per-pair consequence of a single shared unit, as warning text
#'
#' `no_spec[t]`: trait t carries no trait-specific variance (no informative
#' trait-specific units, or `pi_t = 1`). A pair of such traits is two scalar
#' multiples of one design column -- exactly +/-1; every other pair is one
#' noisy draw. Pairs with `|cor| = 1` are realizable and skipped, and so are
#' pairs with a zero-variance trait (`vg <= 0`), whose correlation is undefined
#' (.pleio_check_zero_var() reports those). Returns "" when no pair remains.
#' @keywords internal
#' @noRd
.pleio_single_unit_consequence <- function(R, no_spec, one = "locus",
                                           vg = NULL) {
  nt <- nrow(R)
  exact <- character(0)
  noisy <- character(0)
  for (i in seq_len(nt - 1L)) {
    for (j in seq(i + 1L, nt)) {
      if (abs(R[i, j]) >= 1) next
      if (!is.null(vg) && (vg[i] <= 0 || vg[j] <= 0)) next
      lab <- paste0(i, "-", j)
      if (no_spec[i] && no_spec[j]) exact <- c(exact, lab) else noisy <- c(noisy, lab)
    }
  }
  parts <- character(0)
  if (length(exact)) {
    parts <- c(parts, paste0(
      "traits ", paste(exact, collapse = ", "), ": with no trait-specific ",
      "variance the trait values are scalar multiples of one design column, so ",
      "the realized genetic correlation is exactly +/-1 and a target `cor` ",
      "strictly between -1 and 1 (0 included) cannot be realized"))
  }
  if (length(noisy)) {
    parts <- c(parts, paste0(
      "traits ", paste(noisy, collapse = ", "), ": the whole cross-trait ",
      "covariance rests on that single ", one, ", so the realized genetic ",
      "correlation is one noisy draw that can be far from, or opposite in sign ",
      "to, the target `cor`"))
  }
  paste(parts, collapse = "; ")
}

#' Correlated dominance / epistasis draw under the pleiotropy architecture
#'
#' Extends the PleioArch covariance construction from the additive layer to the
#' non-additive mean-effect components (DECISION-023). The same `cor`, `pi` and
#' feasibility rules apply: units (a locus for dominance, an interacting set for
#' epistasis) are partitioned into shared units -- common to every trait, with
#' effects drawn jointly from MVN(0, Sigma / n) -- and trait-specific units with
#' independent effects, where `Sigma[i, i] = pi_i * V_i` and
#' `Sigma[i, j] = cor_ij * sqrt(V_i * V_j)`.
#'
#' Each unit's effect is divided by the realized standard deviation of its design
#' column `z_u` (the heterozygote indicator for dominance; the centered product,
#' [.epi_unit_column()], for epistasis), so every unit contributes equal design
#' variance on the simulated sample. With units uncorrelated with one another
#' (linkage equilibrium, disjoint loci), the raw component `c_t` (before the
#' layer is rescaled to `prop`) has, over effect draws, `E[Cov(c_1, c_2)] =
#' Sigma_12` and `E[Var(c_t)] = V_t` (the shared units contribute
#' `Sigma_tt = pi_t V_t` and the independent trait-specific units the remaining
#' `(1 - pi_t) V_t`; only the covariance is carried by the shared units), so the
#' component targets `cor` in the
#' same sense as the additive layer: after rescaling its realized correlation is a
#' random ratio that converges to `cor` as the units and the individuals grow
#' (given that near-independence), is attenuated toward 0 on average with few
#' units, and
#' under strong LD need not converge (see [.pleio_draw()]). Independently drawn
#' components add no cross-component covariance in expectation, so the TOTAL
#' genetic correlation targets `cor` when the layers' per-trait `prop` profiles
#' are proportional (always so for scalar `prop`) and otherwise targets a value
#' attenuated toward 0 -- see [.pleio_total_cor_check()], which warns with that
#' large-sample target. Constant design columns are left out of the variance allocation
#' ([.pleio_unit_effects()]). The additive layer's `1/sqrt(2 MAF (1 - MAF))` is
#' the HWE expectation of
#' the same normalizer; non-additive design variance is not a clean function of
#' MAF (a near-inbred panel has far fewer heterozygotes than 2pq), so the realized
#' sd is used. This extension is this package's design, not a published method.
#'
#' `q = NULL` draws fresh units; a supplied `q` (e.g. dominance reusing the
#' additive loci) is taken as given, with the shared units identified as those
#' common to every trait -- or, when `shared` is supplied, as exactly those loci
#' (the additive draw drops a trait's zero-effect loci, so a trait with `pi = 0`
#' no longer holds the shared loci and the intersection would be empty).
#' @param shared optional integer vector of the shared loci of a reused `q`.
#' @param component "dominance" or "epistasis".
#' @param n_units units per trait for a fresh draw.
#' @param interaction markers per unit (1 for dominance).
#' @param itype epistasis `interaction_type` vector.
#' @param arg argument named in partition messages.
#' @return list(qtn = <per-trait loci or set matrices>, effect = <per-trait>).
#' @keywords internal
#' @noRd
.pleio_nonadditive_draw <- function(sim, prop_vec, sub_seed, component,
                                    q = NULL, n_units = NULL,
                                    interaction = 1L, itype = NULL,
                                    arg = "n_qtn", shared = NULL) {
  if (!is.null(sim$arch_args[["cor"]])) {
    .cite_pleioarch()
  }
  R <- .pleio_cor_matrix(sim)
  pi_vec <- .pleio_pi_vector(sim)
  vg <- prop_vec
  .pleio_check_zero_var(R, vg, component)
  sigma <- outer(sqrt(vg), sqrt(vg)) * R
  diag(sigma) <- pi_vec * vg
  .check_pleio_feasible(sigma, R, pi_vec, component)

  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }
  fresh <- is.null(q)
  if (fresh) {
    q <- .pleio_units(sim, n_units, interaction, pi_vec, R, arg, vg = vg)
  }
  unit_column <- if (identical(component, "dominance")) {
    function(j) (.geno_cols(sim, j) == 0) * 1
  } else {
    function(loci) .epi_unit_column(sim, loci, itype)
  }
  eff <- .pleio_unit_effects(q, sigma, pi_vec, vg, unit_column, component,
                             R = R, fresh = fresh, shared = shared)
  target <- attr(eff, "target_cor")
  attr(eff, "target_cor") <- NULL
  list(qtn = q, effect = eff, target_cor = target)
}

#' Draw shared + trait-specific units for a pleiotropic non-additive layer
#'
#' Loci are sampled without replacement, so shared and trait-specific units --
#' and the loci within and across interacting sets -- are all disjoint. Uses the
#' caller's RNG state (the caller sets the sub-seed).
#' @keywords internal
#' @noRd
.pleio_units <- function(sim, n_units, interaction, pi_vec, R, arg,
                         vg = NULL) {
  nt <- sim$n_traits
  is_set <- interaction > 1L
  part <- .pleio_partition(n_units, pi_vec, R, arg = arg, vg = vg,
                           unit  = if (is_set) "interacting set" else "QTN",
                           units = if (is_set) "sets" else "QTNs",
                           place = if (is_set) "sets" else "loci",
                           one   = if (is_set) "set" else "locus")
  pleio_n <- part$pleio_n
  spec_n  <- part$spec_n
  cand <- .candidate_markers(sim)
  need <- (pleio_n + nt * spec_n) * interaction
  if (need > length(cand)) {
    stop("architecture = \"pleiotropy\": the shared + trait-specific ",
         if (is_set) "interacting sets" else "QTNs", " need ", need,
         " distinct markers but only ", length(cand), " candidates exist. ",
         "Lower ", arg, ", or raise pi (shared units are reused by every ",
         "trait, so a larger pi needs fewer distinct markers).", call. = FALSE)
  }
  drawn <- sample(cand, need, replace = FALSE)
  spec_rows <- function(t) pleio_n + (t - 1L) * spec_n + seq_len(spec_n)
  if (!is_set) {
    shared <- drawn[seq_len(pleio_n)]
    return(lapply(seq_len(nt), function(t) c(shared, drawn[spec_rows(t)])))
  }
  M <- matrix(drawn, ncol = interaction, byrow = TRUE)
  shared <- M[seq_len(pleio_n), , drop = FALSE]
  lapply(seq_len(nt), function(t) {
    rbind(shared, M[spec_rows(t), , drop = FALSE])
  })
}

#' Per-trait effects: MVN on shared units, univariate on trait-specific ones
#'
#' Shared units are those present for every trait (identified by content, so a
#' reused locus set works). Every effect is divided by its unit's realized design
#' sd. A constant design column (sd 0, e.g. a hetless dominance locus)
#' contributes nothing whatever its effect, so it gets effect 0 and is left OUT
#' of the variance allocation: the shared covariance Sigma is split over the
#' informative shared units only, and each trait's `(1 - pi) V` over its
#' informative specific units only. Counting dead units would shrink one class
#' relative to the other and bias the realized correlation (one constant unit
#' among two shared units takes cor 0.5, pi 0.5 to ~1/3). RNG order: the shared
#' MVN block, then trait 1..n_traits' specific draws.
#' @keywords internal
#' @noRd
.pleio_unit_effects <- function(q, sigma, pi_vec, vg, unit_column,
                                component = "dominance", R = NULL,
                                fresh = TRUE, shared = NULL) {
  nt <- length(q)
  is_set <- is.matrix(q[[1L]])
  as_units <- function(x) {
    if (is_set) lapply(seq_len(nrow(x)), function(i) x[i, ]) else as.list(x)
  }
  units_t <- lapply(q, as_units)
  keys_t  <- lapply(units_t, function(us) {
    vapply(us, paste, character(1), collapse = "-")
  })
  shared <- if (is.null(shared)) Reduce(intersect, keys_t) else
    intersect(as.character(shared), unlist(keys_t, use.names = FALSE))

  flat_keys  <- unlist(keys_t, use.names = FALSE)
  flat_units <- unlist(units_t, recursive = FALSE, use.names = FALSE)
  first <- !duplicated(flat_keys)
  sds <- vapply(flat_units[first], function(u) stats::sd(unit_column(u)),
                numeric(1))
  live <- is.finite(sds) & sds > 0
  inv <- ifelse(live, 1 / sds, 1)
  names(inv) <- flat_keys[first]
  live_keys <- flat_keys[first][live]

  shared_live <- shared[shared %in% live_keys]
  spec_all_any  <- vapply(keys_t, function(kt) any(!(kt %in% shared)),
                          logical(1))
  spec_live_any <- vapply(keys_t, function(kt) {
    any(!(kt %in% shared) & kt %in% live_keys)
  }, logical(1))
  if (length(shared) > 0L && length(shared_live) == 0L) {
    # The lost shared share always matters; the correlation only where a
    # nonzero cor was requested (independent specific effects still target 0).
    consequence <- if (any(sigma[upper.tri(sigma)] != 0)) {
      paste0("the shared variance -- and the requested cross-trait ",
             "correlation -- cannot be carried: this component's correlation ",
             "is no longer controlled by `cor` (its independent trait-specific ",
             "effects target 0, and LD among those loci can still push it ",
             "anywhere in [-1, 1])")
    } else {
      paste0("the shared (pleiotropic) share `pi` of its variance cannot be ",
             "carried and comes from the trait-specific units only (the ",
             "requested cor = 0 is still targeted)")
    }
    warning("architecture = \"pleiotropy\": every shared ", component, " unit ",
            "has a constant design column (e.g. no heterozygotes), so ",
            consequence, ". Use het-bearing loci (filter_geno(hets = ",
            "\"include\")) or an outbred / F2 population.", call. = FALSE)
  }
  # Exactly one informative shared unit: the whole cross-trait covariance rests
  # on one effect pair. .pleio_partition() already warns when a FRESH draw
  # nominally has one shared unit, classifying pairs from the nominal counts;
  # warn here when the count reached one only after dropping constant units,
  # when the units were supplied (reused additive loci), where no partition
  # check ran, or when a trait's specific units all proved constant (which
  # turns that nominal "noisy draw" into an exact +/-1).
  lost_spec <- spec_all_any & !spec_live_any & pi_vec < 1
  if (length(shared_live) == 1L &&
      (!fresh || length(shared) > 1L || any(lost_spec)) &&
      !is.null(R) && any(abs(R[upper.tri(R)]) < 1)) {
    consequence <- .pleio_single_unit_consequence(R, !spec_live_any |
                                                    pi_vec >= 1, "unit", vg)
  } else {
    consequence <- ""
  }
  if (nzchar(consequence)) {
    warning("architecture = \"pleiotropy\": only one informative shared ",
            component, " unit after dropping units with constant design ",
            "columns (e.g. no heterozygotes), so the whole cross-trait ",
            "covariance rests on a single effect pair; ", consequence, ".",
            call. = FALSE)
  }
  # The correlation this component targets once constant units drop out, for
  # .pleio_total_cor_check(): the covariance survives only with a live shared
  # unit, and each variance keeps only its live classes.
  sh_on <- length(shared_live) > 0L
  v_eff <- (if (sh_on) pi_vec * vg else 0) +
    ifelse(spec_live_any, (1 - pi_vec) * vg, 0)
  target <- (if (sh_on) sigma else matrix(0, nt, nt)) / sqrt(outer(v_eff, v_eff))
  target[!is.finite(target)] <- 0
  diag(target) <- 1
  w_shared <- .draw_mvnorm(length(shared_live), sigma)
  out <- lapply(seq_len(nt), function(t) {
    kt <- keys_t[[t]]
    sh <- kt %in% shared_live
    spec_all  <- !(kt %in% shared)
    spec_live <- spec_all & kt %in% live_keys
    # Losing trait t's specific variance leaves its correlations targeting
    # Sigma_tj / sqrt(Sigma_tt Sigma_jj), larger in magnitude than cor_tj -- but
    # only where cor_tj != 0 (a zero target stays 0), so warn only then.
    if (any(spec_all) && !any(spec_live) && pi_vec[t] < 1 &&
        any(sigma[t, -t] != 0)) {
      warning("architecture = \"pleiotropy\": every trait-specific ", component,
              " unit of trait ", t, " has a constant design column, so its ",
              "trait-specific variance is lost and the correlation this ",
              "component targets is inflated in magnitude (further from 0 than ",
              "`cor`; any single draw still scatters around it). Use ",
              "het-bearing loci or an outbred / F2 population.", call. = FALSE)
    }
    e <- numeric(length(kt))                 # constant units stay at 0
    e[sh] <- w_shared[match(kt[sh], shared_live), t]
    e[spec_live] <- .draw_univariate(sum(spec_live), (1 - pi_vec[t]) * vg[t])
    e * unname(inv[kt])
  })
  attr(out, "target_cor") <- target
  out
}

#' Warn when the pleiotropic mean-effect layers cannot give a total of `cor`
#'
#' Each mean-effect layer c targets `cor` on its own (its effect draw's
#' cross-trait covariance is `cor_ij sqrt(V_ci V_cj)` in expectation), and
#' independently drawn layers add no cross-layer
#' covariance in expectation, so the total genetic correlation of traits i, j
#' targets `sum_c cor_ij sqrt(V_ci V_cj) / sqrt(sum_c V_ci * sum_c V_cj)`. By
#' Cauchy-Schwarz this target equals `cor_ij` when the layers' per-trait
#' variance profiles are proportional (always so for scalar `prop`) or
#' `cor_ij = 0`, and is otherwise attenuated toward 0 -- e.g. additive `prop = c(0.49, 0.01)` with
#' dominance `prop = c(0.01, 0.49)` gives 0.28 x `cor` (0.14 at `cor = 0.5`),
#' and no per-layer correlation <= 1 could lift it to `cor`. The layers are added one at a time,
#' so the total is not re-targeted; warn with its target instead (whenever it
#' differs from `cor_ij` by more than 1% of `cor_ij`). That value is
#' a ratio of expected moments -- the large-sample target the realized total
#' converges to as units and individuals grow, under approximate linkage
#' equilibrium -- not `E[Cor(g1, g2)]`:
#' with few shared units the realized total is on average smaller in magnitude
#' (e.g. two shared loci, target 0.40, realized mean ~0.29), as for a single
#' layer.
#' A non-additive layer whose shared or trait-specific units all proved constant
#' no longer targets `cor` itself; its effective target (`layer$target_cor`,
#' from [.pleio_unit_effects()]) replaces `cor` in the sum.
#' @param new_prop_vec per-trait `prop` of the layer being added.
#' @param new_target that layer's effective target correlation matrix (`NULL`:
#'   `cor`).
#' @param new_target_reps per-replication targets under `vary_qtn` (`NULL`: the
#'   layer has one draw). Each replication is checked with every layer's target
#'   for that replication; one that differs from the canonical is reported with
#'   its number.
#' @keywords internal
#' @noRd
.pleio_total_cor_check <- function(sim, new_prop_vec, new_target = NULL,
                                   new_target_reps = NULL) {
  if (is.null(sim$arch_args[["cor"]])) {
    return(invisible())                      # cor = 0: the total is 0 too
  }
  nt <- sim$n_traits
  mean_types <- c("additive", "dominance", "epistasis")
  prior <- Filter(function(l) l$type %in% mean_types && !isTRUE(l$orthogonal),
                  sim$layers)
  if (length(prior) == 0L) {
    return(invisible())
  }
  V <- rbind(do.call(rbind, lapply(prior, function(l) .expand_prop(l$prop, nt))),
             new_prop_vec)
  R <- .pleio_cor_matrix(sim)
  # One target matrix per layer and replication: a layer's own target for that
  # replication (vary_qtn), else its single target, else `cor` (additive).
  pick <- function(single, reps, r) {
    if (!is.null(reps)) reps[[r]] else if (!is.null(single)) single else R
  }
  layer_t <- c(lapply(prior, function(l) list(l$target_cor, l$target_cor_reps)),
               list(list(new_target, new_target_reps)))
  n_rep <- max(1L, vapply(layer_t, function(x) length(x[[2]]), integer(1)))
  # three decimals, but significant digits for tiny values, so a large
  # relative attenuation at a very small `cor` is not rounded to 0.000
  fmt <- function(x) if (abs(x) >= 0.001 || x == 0) sprintf("%.3f", x) else
    sprintf("%.3g", x)
  off <- character(0)
  for (i in seq_len(nt - 1L)) {
    for (j in seq(i + 1L, nt)) {
      den <- sqrt(sum(V[, i]) * sum(V[, j]))
      tot_r <- vapply(seq_len(n_rep), function(r) {
        rho <- vapply(layer_t, function(x) pick(x[[1]], x[[2]], r)[i, j],
                      numeric(1))
        if (den > 0) sum(rho * sqrt(V[, i] * V[, j])) / den else NA_real_
      }, numeric(1))
      for (r in seq_len(n_rep)) {
        tot <- tot_r[r]
        # Relative, not absolute: at small `cor` an absolute 0.01 would hide a
        # large attenuation (cor 0.01 -> 0.0028 is 72%). Proportional profiles
        # give tot == cor up to roundoff, far inside 1%.
        if (!is.finite(tot) || abs(tot - R[i, j]) <= 0.01 * abs(R[i, j])) next
        if (r == 1L) {
          off <- c(off, sprintf("traits %d-%d: %s (cor %s)", i, j, fmt(tot),
                                fmt(R[i, j])))
        } else if (!isTRUE(abs(tot - tot_r[1L]) < 1e-9)) {
          off <- c(off, sprintf("traits %d-%d: %s (cor %s) in replication %d",
                                i, j, fmt(tot), fmt(R[i, j]), r))
        }
      }
    }
  }
  if (length(off)) {
    warning("architecture = \"pleiotropy\": the TOTAL genetic correlation ",
            "targets the variance-weighted combination of the layers' targets, ",
            "not `cor`: ", paste(off, collapse = "; "),
            " (a large-sample target; with few shared QTNs / sets the realized ",
            "total is on average smaller in magnitude, and any single draw ",
            "scatters around it). Each mean-effect layer targets `cor` only ",
            "while its shared and trait-specific units are informative, and the ",
            "total only when the layers' per-trait `prop` profiles are ",
            "proportional (e.g. scalar `prop`).", call. = FALSE)
  }
  invisible()
}

#' Citation notice for the correlation-control engine, once per session
#'
#' The algorithm behind `cor` in the pleiotropy architecture is described in a
#' manuscript still in preparation (Prado et al.), so credit that the first time
#' it is used in a session; the concrete definition is the bundled reference
#' implementation in `context/PleioArch-main/`. Emitted via [rlang::inform()] so
#' it is a message, not a warning, and can be silenced with `suppressMessages()`.
#' @keywords internal
#' @noRd
.cite_pleioarch <- function() {
  rlang::inform(
    .cite_main("Prado et al. (in preparation), when controlling the ",
               "correlation in the pleiotropic architecture."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_pleioarch_citation"
  )
}

#' The primary citation notice, shared by every "please also cite" message
#'
#' Builds the multi-line notice: the lead-in, the main simplePHENOTYPES
#' reference, then the feature-specific reference(s) passed in `...`.
#' @keywords internal
#' @noRd
.cite_main <- function(...) {
  paste0("In addition to citing:\n",
         "Fernandes, S.B., Lipka, A.E. simplePHENOTYPES: SIMulation of ",
         "pleiotropic, linked and epistatic phenotypes. BMC Bioinformatics 21, ",
         "491 (2020). https://doi.org/10.1186/s12859-020-03804-y\n",
         "Please also cite:\n",
         ...)
}

#' Per-trait pleiotropic variance share
#'
#' `pi` sets the share for every trait; `pi_target` / `pi_secondary` are the
#' two-trait spelling kept from the reference implementation.
#' @keywords internal
#' @noRd
.pleio_pi_vector <- function(sim) {
  a <- sim$arch_args
  nt <- sim$n_traits
  if (!is.null(a[["pi"]]) &&
      (!is.null(a[["pi_target"]]) || !is.null(a[["pi_secondary"]]))) {
    stop("Use either `pi` or `pi_target`/`pi_secondary`, not both.",
         call. = FALSE)
  }
  if (nt != 2L &&
      (!is.null(a[["pi_target"]]) || !is.null(a[["pi_secondary"]]))) {
    stop("`pi_target` and `pi_secondary` are the two-trait interface; use ",
         "`pi` when n_traits is not 2.", call. = FALSE)
  }
  if (!is.null(a[["pi"]])) {
    if (!length(a[["pi"]]) %in% c(1L, nt)) {
      stop("`pi` must have length 1 or n_traits (", nt, "); got ",
           length(a[["pi"]]), ".", call. = FALSE)
    }
    p <- rep_len(a[["pi"]], nt)
  } else {
    piT <- if (is.null(a[["pi_target"]])) 1 else a[["pi_target"]]
    p <- rep(piT, nt)
    if (nt >= 2 && !is.null(a[["pi_secondary"]])) {
      p[2] <- a[["pi_secondary"]]
    }
  }
  if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("Pleiotropic shares (`pi`, `pi_target`, `pi_secondary`) must be ",
         "between 0 and 1.", call. = FALSE)
  }
  p
}

#' Check that the requested correlations are attainable
#'
#' The pleiotropic covariance matrix must be positive semi-definite. For two
#' traits this is exactly the reference constraint `cor^2 <= pi_1 * pi_2`; with
#' more traits it additionally rules out mutually inconsistent correlations
#' (for example three traits that are each strongly negatively correlated).
#' @keywords internal
#' @noRd
.check_pleio_feasible <- function(sigma, R, pi_vec, component = NULL) {
  nt <- nrow(sigma)
  # Name the layer when the check runs for a non-additive component, so an
  # infeasible request says which layer it came from (the additive layer keeps
  # its unprefixed message).
  who <- if (is.null(component)) "" else paste0(component, "(): ")
  if (nt == 2) {
    lhs <- R[1, 2]^2
    rhs <- pi_vec[1] * pi_vec[2]
    # Reject at the level of floating-point roundoff only, scaled to the operands
    # actually compared (not an absolute floor like max(1, rhs), which would let
    # cor^2 > 0 slip through when pi1*pi2 = 0, e.g. pi = c(0.25, 0), where any
    # nonzero correlation is impossible). cor^2 and pi1*pi2 are each one product,
    # so a few ULPs of their magnitude is the right slack.
    tol <- 8 * .Machine$double.eps * max(lhs, rhs)
    if (lhs - rhs > tol) {
      stop(who, "Biological constraint violated: cor^2 (", signif(lhs, 6),
           ") cannot exceed pi_target * pi_secondary (", signif(rhs, 6),
           ").", call. = FALSE)
    }
    return(invisible(TRUE))
  }
  # sigma = D^(1/2) M D^(1/2) with D = diag(V), M[i,i] = pi_i, M[i,j] = cor_ij,
  # so for positive variances sigma is PSD exactly when M is (a zero-variance
  # trait has a zero cor row, enforced by .pleio_check_zero_var()). Test M: its
  # entries are O(1) whatever the `prop`s, so the tolerance below cannot be
  # inflated by one large variance and wave through an indefinite block of
  # tiny-variance traits (review round 10: prop = c(1e-16, 1e-16, 0.5)).
  M <- R
  diag(M) <- pi_vec
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  # Reject a negative eigenvalue unless it is at the level of eigen-decomposition
  # roundoff, which is a small multiple of (matrix scale) x (dimension) x eps --
  # NOT sqrt(eps) (~1e-8), which is astronomically larger than real roundoff and
  # would wave a genuinely indefinite matrix (e.g. min eigenvalue -1e-8) straight
  # through to the effect draw. An exactly singular feasible boundary (min
  # eigenvalue ~ 0) still passes and is drawn exactly by .draw_mvnorm().
  scale <- max(abs(ev))
  tol <- scale * nt * 64 * .Machine$double.eps
  if (scale > 0 && min(ev) < -tol) {
    stop(who, "The requested `cor` is not attainable with these pleiotropic ",
         "shares: the implied genetic covariance matrix is not positive ",
         "semi-definite (smallest eigenvalue ", format(min(ev), digits = 3),
         "). Lower the correlations or raise `pi`. For two traits the ",
         "condition is cor^2 <= pi_1 * pi_2.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Multivariate-normal effect draw with per-SNP variance scaling
#'
#' Uses the symmetric eigendecomposition square root rather than a Cholesky of a
#' nudged matrix. Feasibility is already checked upstream, so `sigma` is PSD;
#' the eigen root reproduces it exactly even when it is singular (the feasibility
#' boundary `cor^2 = pi_1 * pi_2`, where `chol()` fails and any nudge would
#' silently perturb a valid request). Only numerical-roundoff negative
#' eigenvalues are floored at zero.
#' @keywords internal
#' @noRd
.draw_mvnorm <- function(n, sigma) {
  nt <- nrow(sigma)
  if (n <= 0 || all(diag(sigma) <= 0)) {
    return(matrix(0, nrow = max(n, 0), ncol = nt))
  }
  sigma_per <- sigma / n
  z <- matrix(stats::rnorm(n * nt), ncol = nt)
  e <- eigen(sigma_per, symmetric = TRUE)
  root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  z %*% root
}

#' Univariate-normal effect draw with per-SNP variance scaling
#' @keywords internal
#' @noRd
.draw_univariate <- function(n, v) {
  if (n <= 0) {
    return(numeric(0))
  }
  stats::rnorm(n, mean = 0, sd = sqrt(max(v, 0) / n))
}

#' Target genetic-correlation matrix for the >2-trait Cholesky fallback
#'
#' Uses a single scalar `cor` for every trait pair, or a supplied n_traits x
#' n_traits matrix passed as `cor`.
#' @keywords internal
#' @noRd
.pleio_cor_matrix <- function(sim) {
  nt <- sim$n_traits
  cor_g <- sim$arch_args[["cor"]]
  if (is.null(cor_g)) cor_g <- 0
  if (is.matrix(cor_g)) {
    if (!all(dim(cor_g) == c(nt, nt))) {
      stop("`cor` matrix must be ", nt, " x ", nt, ".", call. = FALSE)
    }
    if (!is.numeric(cor_g) || any(!is.finite(cor_g)) ||
        any(cor_g < -1 | cor_g > 1)) {
      stop("Every entry of the `cor` matrix must be finite and between -1 ",
           "and 1.", call. = FALSE)
    }
    if (!isTRUE(all.equal(cor_g, t(cor_g), tolerance = 1e-12))) {
      stop("The `cor` matrix must be symmetric.", call. = FALSE)
    }
    if (any(abs(diag(cor_g) - 1) > 1e-12)) {
      stop("The diagonal of the `cor` matrix must equal 1.", call. = FALSE)
    }
    return(cor_g)
  }
  if (!is.numeric(cor_g) || length(cor_g) != 1L || !is.finite(cor_g) ||
      cor_g < -1 || cor_g > 1) {
    stop("`cor` must be one finite value between -1 and 1, or a valid ",
         "correlation matrix.", call. = FALSE)
  }
  R <- matrix(cor_g, nt, nt)
  diag(R) <- 1
  R
}

#' Validate user-supplied dominance loci / epistatic sets under pleiotropy
#'
#' Units (loci for dominance, interacting sets for epistasis) arrive as `q`, a
#' per-trait list of index vectors or set matrices. As for additive loci
#' ([.pleio_user_layout()]) every unit must be listed for every trait: a unit
#' that affects only some traits is partial pleiotropy (`complex_phenotypes()`).
#' Sets match as ordered tuples (`c(3, 9)` and `c(9, 3)` are different sets, as in
#' the correlated draw).
#' @return `q` with every trait given trait 1's units, invisibly.
#' @keywords internal
#' @noRd
.pleio_check_user_units <- function(sim, prop, q, component, is_set) {
  nt <- length(q)
  keys <- lapply(q, function(x) {
    if (is_set) apply(x, 1L, paste, collapse = "-") else as.character(x)
  })
  unit <- if (is_set) "interacting set" else "locus"
  if (any(vapply(keys, anyDuplicated, 0L) > 0L)) {
    stop(component, "(qtn=): a ", unit, " is listed more than once for a ",
         "trait; list each ", unit, " once.", call. = FALSE)
  }
  if (!all(vapply(keys, function(k) setequal(k, keys[[1L]]), logical(1)))) {
    stop(component, "(qtn=): under architecture = \"pleiotropy\" every ", unit,
         " affects every trait, so `qtn` must list the same ",
         if (is_set) "sets" else "loci", " for each trait. A ", unit, " that ",
         "affects only some of the traits is partial pleiotropy: build one ",
         "model per group of traits and combine them with ",
         "complex_phenotypes().", call. = FALSE)
  }
  n_units <- length(keys[[1L]])
  .pleio_check_layout(.pleio_cor_matrix(sim), .pleio_pi_vector(sim),
                      prop, n_units, rep(0L, nt),
                      unit = if (is_set) "interacting set" else "QTN",
                      units = if (is_set) "sets" else "QTNs",
                      place = if (is_set) "sets" else "loci",
                      one = if (is_set) "set" else "locus")
  invisible(q)
}

#' Is a correlation being controlled under pleiotropy?
#'
#' TRUE when the user gave any of `cor`, `pi`, `pi_target`, `pi_secondary`,
#' `n_pleio_major`, `prop_var_major`. Without them the shared loci's effects are
#' just set (by `effect` / `dist`, or the default draw) and the genetic
#' correlation is an outcome of the shared loci, not a target.
#' @keywords internal
#' @noRd
.pleio_controlled <- function(sim) {
  a <- sim$arch_args
  any(!vapply(a[c("cor", "pi", "pi_target", "pi_secondary", "n_pleio_major",
                  "prop_var_major")], is.null, logical(1)))
}
