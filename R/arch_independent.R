#' Draw QTN marker indices for a layer according to the architecture
#'
#' Parity-critical QTN sampling (RNG stays in R, DECISION-006). Returns a list
#' of length `n_traits`, each element an integer vector of marker-column indices
#' into `sim$G`.
#'
#' - "independent": each trait draws its own distinct QTN set.
#' - "pleiotropy": one shared QTN set is reused across all traits (cross-trait
#'   correlation is produced later by the correlated effect draws).
#' - "ld": causal markers are drawn as for "independent"; LD-companion marker
#'   handling is layered on in [.draw_qtn_ld()] (reporting only — the genetic
#'   value uses the causal markers).
#'
#' @param sim a `phenotype_sim`.
#' @param n_qtn number of QTNs per trait.
#' @param sub_seed deterministic per-layer sub-seed (or NULL).
#' @return list of length `n_traits` of integer index vectors.
#' @keywords internal
#' @noRd
.draw_qtn <- function(sim, n_qtn, sub_seed) {
  cand <- .candidate_markers(sim)
  if (n_qtn > length(cand)) {
    stop("Requested n_qtn (", n_qtn, ") exceeds the number of candidate ",
         "markers (", length(cand), ").", call. = FALSE)
  }
  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }

  nt <- sim$n_traits
  if (sim$architecture == "pleiotropy") {
    shared <- sample(cand, n_qtn, replace = FALSE)
    return(rep(list(shared), nt))
  }
  lapply(seq_len(nt), function(t) sample(cand, n_qtn, replace = FALSE))
}

#' Draw distinct QTN sets per trait regardless of architecture
#'
#' Used by the >2-trait pleiotropy Cholesky fallback, where each trait needs its
#' own genetic values for the recorrelation step.
#' @keywords internal
#' @noRd
.draw_qtn_per_trait <- function(sim, n_qtn, sub_seed) {
  cand <- .candidate_markers(sim)
  if (n_qtn > length(cand)) {
    stop("Requested n_qtn (", n_qtn, ") exceeds candidate markers (",
         length(cand), ").", call. = FALSE)
  }
  lapply(seq_len(sim$n_traits), function(t) {
    s <- .layer_seed(sim$seed, paste0("pleio_qtn_t", t),
                     if (is.null(sub_seed)) 0L else sub_seed %% 1000L)
    old <- .Random.seed_safe()
    if (!is.null(s)) {
      set.seed(s)
      on.exit(.restore_seed(old))
    }
    sample(cand, n_qtn, replace = FALSE)
  })
}

#' Draw epistatic QTN pairs (n_pairs x interaction matrices) per trait
#' @keywords internal
#' @noRd
.draw_qtn_pairs <- function(sim, n_pairs, interaction, sub_seed) {
  cand <- .candidate_markers(sim)
  need <- n_pairs * interaction
  if (need > length(cand)) {
    stop("Requested epistatic markers (", need, ") exceed candidate markers (",
         length(cand), ").", call. = FALSE)
  }
  old <- .Random.seed_safe()
  if (!is.null(sub_seed)) {
    set.seed(sub_seed)
    on.exit(.restore_seed(old))
  }
  nt <- sim$n_traits
  make_one <- function() {
    matrix(sample(cand, need, replace = FALSE), nrow = n_pairs,
           ncol = interaction, byrow = TRUE)
  }
  if (sim$architecture == "pleiotropy") {
    shared <- make_one()
    return(rep(list(shared), nt))
  }
  lapply(seq_len(nt), function(t) make_one())
}

#' Candidate marker columns (optionally constrained); currently all markers
#' @keywords internal
#' @noRd
.candidate_markers <- function(sim) {
  seq_len(ncol(sim$G))
}

#' Count of prior layers of a given type (for sub-seed occurrence index)
#' @keywords internal
#' @noRd
.type_occurrence <- function(sim, type) {
  sum(vapply(sim$layers, function(l) identical(l$type, type), logical(1)))
}
