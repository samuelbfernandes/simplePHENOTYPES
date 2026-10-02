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
    stop("Requested n_qtn (", n_qtn, ") exceeds the number of polymorphic ",
         "markers (", length(cand), ")", .candidate_hint(sim), ".", call. = FALSE)
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
  if (isTRUE(sim$arch_args[["distinct_chr"]]) && nt > 1) {
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
         length(cand), ")", .candidate_hint(sim), ".", call. = FALSE)
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

#' Candidate marker columns for QTN selection
#'
#' Excludes markers with a constant dosage column, which can realize no
#' variance component: monomorphic markers (MAF = 0) and markers that are
#' heterozygous in every individual (MAF = 0.5 but zero dosage variance, e.g. an
#' F1). Polymorphic markers are unaffected. Use [filter_geno()] to constrain the
#' pool further (MAF, LD, heterozygosity) before simulating.
#' @keywords internal
#' @noRd
.candidate_markers <- function(sim) {
  ok <- is.finite(sim$maf) & sim$maf > 0
  if (!is.null(sim$all_het) && length(sim$all_het) == length(ok)) {
    ok <- ok & !sim$all_het
  }
  which(ok)
}

#' Suffix explaining why the candidate pool may be smaller than the marker count
#' @keywords internal
#' @noRd
.candidate_hint <- function(sim) {
  n_het <- sum(sim$all_het, na.rm = TRUE)
  if (n_het > 0L) {
    paste0(" (", n_het, " marker(s) heterozygous in every individual have a ",
           "constant dosage column and are excluded, as are monomorphic ",
           "markers)")
  } else {
    ""
  }
}

#' Count of prior layers of a given type (for sub-seed occurrence index)
#' @keywords internal
#' @noRd
.type_occurrence <- function(sim, type) {
  sum(vapply(sim$layers, function(l) identical(l$type, type), logical(1)))
}
