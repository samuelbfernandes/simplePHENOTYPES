#' Add an additive variance component
#'
#' @param sim a `phenotype_sim`.
#' @param prop proportion of total phenotypic variance (scalar or length
#'   `n_traits`).
#' @param n_qtn number of additive QTNs; overrides the baseline `n_qtn` from
#'   [simulate_phenotype()] (with a warning when both are given).
#' @param effect optional geometric base (scalar) or explicit effect series
#'   (length `n_qtn`).
#' @param dist within-layer effect distribution (default "geometric").
#' @return the updated `phenotype_sim`.
#' @export
additive <- function(sim, prop, n_qtn = NULL, effect = NULL,
                     dist = "geometric") {
  .check_sim(sim)
  nq <- .resolve_n_qtn(sim, n_qtn, "additive")
  occ <- .type_occurrence(sim, "additive")
  sub_seed <- .layer_seed(sim$seed, "additive", occ)

  if (sim$architecture == "pleiotropy" && sim$n_traits == 2) {
    pd <- .pleio_draw(sim, nq, .expand_prop(prop, 2), sub_seed)
    qtn <- pd$qtn
    eff <- pd$effect
  } else if (sim$architecture == "pleiotropy" && sim$n_traits > 2) {
    warning("Pleiotropy with n_traits > 2 uses Cholesky-correlated genetic ",
            "values: the phenotypes are correlated as requested (cor), but ",
            "individual QTN effect sizes are not guaranteed. The exact ",
            "(PleioArch) engine is available for n_traits = 2.", call. = FALSE)
    sim$pleio_cor <- .pleio_cor_matrix(sim)
    # Distinct QTNs per trait so the genetic values differ and Cholesky has
    # something to recorrelate; the target correlation is imposed at realization.
    qtn <- .draw_qtn_per_trait(sim, nq, sub_seed)
    eff <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, effect))
  } else {
    qtn <- .draw_qtn(sim, nq, sub_seed)
    eff <- lapply(seq_len(sim$n_traits),
                  function(t) .effect_series(nq, dist, effect))
  }

  layer <- list(type = "additive", prop = prop, n_qtn = nq, dist = dist,
                qtn = qtn, effect = eff)
  if (sim$architecture == "ld") {
    layer$ld <- .annotate_ld(sim, qtn)
  }
  .add_layer(sim, layer)
}

#' Add a dominance variance component
#'
#' @inheritParams additive
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`).
#' @param degree degree of dominance (optional; modulates the effect series).
#' @return the updated `phenotype_sim`.
#' @export
dominance <- function(sim, prop, same_as_add = TRUE, n_qtn = NULL,
                      degree = NULL, dist = "geometric") {
  .check_sim(sim)
  occ <- .type_occurrence(sim, "dominance")
  sub_seed <- .layer_seed(sim$seed, "dominance", occ)
  if (isTRUE(same_as_add)) {
    add_layer <- .last_layer_of_type(sim, "additive")
    if (is.null(add_layer)) {
      stop("dominance(same_as_add = TRUE) requires a prior additive() layer.",
           call. = FALSE)
    }
    qtn <- add_layer$qtn
    nq  <- add_layer$n_qtn
  } else {
    nq  <- .resolve_n_qtn(sim, n_qtn, "dominance")
    qtn <- .draw_qtn(sim, nq, sub_seed)
  }
  scale <- if (is.null(degree)) 1 else degree
  eff <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(nq, dist) * scale)

  layer <- list(type = "dominance", prop = prop, n_qtn = nq, dist = dist,
                same_as_add = same_as_add, degree = degree,
                qtn = qtn, effect = eff)
  .add_layer(sim, layer)
}

#' Add an epistatic variance component
#'
#' @inheritParams additive
#' @param n_pairs number of interacting QTN sets.
#' @param interaction number of markers per epistatic QTN (default 2, pairwise).
#' @return the updated `phenotype_sim`.
#' @export
epistasis <- function(sim, prop, n_pairs = NULL, interaction = 2,
                      effect = NULL, dist = "geometric") {
  .check_sim(sim)
  np <- .resolve_n_qtn(sim, n_pairs, "epistasis", arg = "n_pairs")
  occ <- .type_occurrence(sim, "epistasis")
  sub_seed <- .layer_seed(sim$seed, "epistasis", occ)
  qtn <- .draw_qtn_pairs(sim, np, interaction, sub_seed)
  eff <- lapply(seq_len(sim$n_traits),
                function(t) .effect_series(np, dist, effect))

  layer <- list(type = "epistasis", prop = prop, n_pairs = np,
                interaction = interaction, dist = dist,
                qtn = qtn, effect = eff)
  .add_layer(sim, layer)
}

#' Add a variance-QTL (vQTL) component
#'
#' @inheritParams additive
#' @param same_as_add reuse the additive layer's QTNs (default `TRUE`).
#' @return the updated `phenotype_sim`.
#' @rdname vqtl-layer
#' @export
vqtl <- function(sim, prop, same_as_add = TRUE, n_qtn = NULL,
                 dist = "geometric") {
  .check_sim(sim)
  occ <- .type_occurrence(sim, "vqtl")
  sub_seed <- .layer_seed(sim$seed, "vqtl", occ)
  if (isTRUE(same_as_add)) {
    add_layer <- .last_layer_of_type(sim, "additive")
    if (is.null(add_layer)) {
      stop("vqtl(same_as_add = TRUE) requires a prior additive() layer.",
           call. = FALSE)
    }
    qtn <- add_layer$qtn
    nq  <- add_layer$n_qtn
  } else {
    nq  <- .resolve_n_qtn(sim, n_qtn, "vqtl")
    qtn <- .draw_qtn(sim, nq, sub_seed)
  }
  eff <- lapply(seq_len(sim$n_traits), function(t) .effect_series(nq, dist))

  layer <- list(type = "vqtl", prop = prop, n_qtn = nq, dist = dist,
                same_as_add = same_as_add, qtn = qtn, effect = eff)
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
    return(baseline)
  }
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

#' Validate the running variance budget and append a layer, then realize
#' @keywords internal
#' @noRd
.add_layer <- function(sim, layer) {
  nt <- sim$n_traits
  .expand_prop(layer$prop, nt)  # validates length
  prospective <- .total_genetic_prop(sim) + .expand_prop(layer$prop, nt)
  if (any(prospective > 1 + 1e-8)) {
    bad <- which(prospective > 1 + 1e-8)
    stop("Genetic proportions for the ", layer$type, "() layer push trait(s) ",
         paste(bad, collapse = ", "), " above 1 (running total = ",
         paste(sprintf("%.3f", prospective[bad]), collapse = ", "), ").",
         call. = FALSE)
  }
  sim$layers <- c(sim$layers, list(layer))
  .realize_phenotype(sim)
}
