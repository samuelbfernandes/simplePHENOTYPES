#' Draw QTN marker indices for a layer according to the architecture
#'
#' Parity-critical QTN sampling (RNG stays in R). Returns a list
#' of length `n_traits`, each element an integer vector of marker-column indices
#' into the genotype reference held by `sim`.
#'
#' - "independent": each trait draws its own distinct QTN set.
#' - "pleiotropy": one shared QTN set is reused across all traits (cross-trait
#'   correlation is produced later by the correlated effect draws).
#' - "ld": two traits, each with its own *distinct* causal loci that sit in
#'   linkage disequilibrium, so the traits covary through linkage rather than
#'   pleiotropy (see `.draw_qtn_ld()`). The returned list carries an `"ld"`
#'   attribute describing the linked pairs.
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
  if (sim$architecture == "ld") {
    return(.draw_qtn_ld(sim, n_qtn))
  }
  if (isTRUE(sim$arch_args$distinct_chr) && nt > 1) {
    return(.draw_qtn_distinct_chr(sim, n_qtn, cand))
  }
  lapply(seq_len(nt), function(t) sample(cand, n_qtn, replace = FALSE))
}

#' Draw each trait's QTNs from a disjoint set of chromosomes
#'
#' For `architecture = "independent"` with `distinct_chr = TRUE`: the
#' chromosomes are partitioned among the traits (round-robin) so no chromosome
#' carries QTNs for more than one trait, giving genetically independent traits
#' with no shared linkage. Errors if a trait's assigned chromosomes cannot
#' supply `n_qtn` markers.
#' @keywords internal
#' @noRd
.draw_qtn_distinct_chr <- function(sim, n_qtn, cand) {
  chr <- sim$map$chr[cand]
  levels_chr <- unique(chr)
  nt <- sim$n_traits
  if (length(levels_chr) < nt) {
    stop("distinct_chr = TRUE needs at least n_traits (", nt, ") chromosomes; ",
         "the map has ", length(levels_chr), ".", call. = FALSE)
  }
  assign_to <- rep(seq_len(nt), length.out = length(levels_chr))
  lapply(seq_len(nt), function(t) {
    my_chr <- levels_chr[assign_to == t]
    pool <- cand[chr %in% my_chr]
    if (length(pool) < n_qtn) {
      stop("distinct_chr = TRUE: trait ", t, " has only ", length(pool),
           " markers on its assigned chromosome(s) but needs ", n_qtn, ".",
           call. = FALSE)
    }
    sample(pool, n_qtn, replace = FALSE)
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

#' Candidate polymorphic marker columns
#' @keywords internal
#' @noRd
.candidate_markers <- function(sim) {
  which(is.finite(sim$maf) & sim$maf > 0)
}

#' Count of prior layers of a given type (for sub-seed occurrence index)
#' @keywords internal
#' @noRd
.type_occurrence <- function(sim, type) {
  sum(vapply(sim$layers, function(l) identical(l$type, type), logical(1)))
}
