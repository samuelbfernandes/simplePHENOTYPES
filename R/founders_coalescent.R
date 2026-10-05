# Coalescent founders (breedingDesigner SPEC-0020 item 9; docs/SPEC-coalescent.md,
# DECISION-048). Phase (a): the internal per-chromosome sampler; the exported
# founders_coalescent() with species presets follows in phase (c).

#' Coalescent haplotypes of one chromosome
#'
#' Thin validating wrapper over the Rust SMC' core. Time is in units of 4 N0
#' generations (ms / MaCS); `theta` = 4 N0 mu and `rho` = 4 N0 r per chromosome.
#' `history` is NULL or a data frame (time, size): from each `time` (4 N0 units)
#' the population size is `size` x N0 (ms `-eN`); times in (0, 1e12] and sizes
#' in [1e-9, 1e9], beyond which the rates leave double precision. One 32-bit seed is drawn from
#' R's stream (or taken from `seed`), so set.seed() reproduces the result.
#'
#' A run whose number of mutation + recombination events would pass 2e9 stops
#' with an error (the rates are too large for the sample) instead of hanging.
#'
#' Theory checked by the tests: E[S] = theta sum_{i<n} 1/i (Watterson 1975),
#' E[xi_i] = theta / i (Fu 1995), and the pair TMRCA under a piecewise-constant
#' size, the integral of its survival function (derived here from the pair
#' coalescence rate 2 / lambda(u); no published page is claimed). SMC' as in
#' Marjoram & Wall (2006), SMC as in McVean & Cardin (2005).
#' @references
#' Watterson GA (1975) On the number of segregating sites in genetical models
#' without recombination. Theor Popul Biol 7:256-276 (PMID 1145509).
#' Fu YX (1995) Statistical properties of segregating sites. Theor Popul Biol
#' 48:172-197 (PMID 7482370).
#' Marjoram P, Wall JD (2006) Fast "coalescent" simulation. BMC Genet 7:16
#' (PMID 16539698).
#' McVean GAT, Cardin NJ (2005) Approximating the coalescent with recombination.
#' Philos Trans R Soc Lond B 360:1387-1393 (PMID 16048782).
#' Wilton PR, Carmi S, Hobolth A (2015) The SMC' is a highly accurate
#' approximation to the ancestral recombination graph. Genetics 200:343-355
#' (PMID 25786855; the two-locus correlation gate).
#' @return list(pos, hap = sites x haplotypes 0/1 matrix, n_total, tmrca, length,
#'   tmrca_end = T_MRCA of the tree at the end of the chromosome).
#' @keywords internal
#' @noRd
.coalescent_chromosome <- function(n_hap, theta, rho, history = NULL,
                                   seg_sites = 0L, seed = NULL) {
  n_hap <- .validate_count(n_hap, "n_hap", minimum = 2L)
  seg_sites <- .validate_count(seg_sites, "seg_sites", minimum = 0L)
  for (nm in c("theta", "rho")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 0) {
      stop("`", nm, "` must be a single finite, non-negative number.", call. = FALSE)
    }
  }
  # plain doubles for the Rust core: a classed number (e.g. bit64::integer64)
  # would otherwise pass its raw storage, not its value
  theta <- as.numeric(theta)
  rho <- as.numeric(rho)
  if (!is.finite(theta + rho)) {
    stop("`theta + rho` must be finite.", call. = FALSE)
  }
  if (is.null(history)) {
    ht <- numeric(0)
    hs <- numeric(0)
  } else {
    if (!is.data.frame(history) || !all(c("time", "size") %in% names(history))) {
      stop("`history` must be NULL or a data frame with columns `time` and `size`.",
           call. = FALSE)
    }
    # numeric columns only: as.numeric() on a factor would silently use its codes
    if (!is.numeric(history$time) || !is.numeric(history$size)) {
      stop("`history$time` and `history$size` must be numeric columns.", call. = FALSE)
    }
    ht <- as.numeric(history$time)
    hs <- as.numeric(history$size)
  }
  # The Rust core takes a 32-bit seed; an explicit seed covers the same range as
  # a drawn one, so the recorded `seed` of any run replays it.
  seed <- if (is.null(seed)) {
    floor(stats::runif(1, 0, 2^32))
  } else {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed) || seed < 0 || seed >= 2^32) {
      stop("`seed` must be NULL or one whole number in [0, 2^32).", call. = FALSE)
    }
    as.numeric(seed)
  }
  out <- coalescent_chromosome_core(as.integer(n_hap), theta, rho, ht, hs,
                                    as.integer(seg_sites), as.numeric(seed))
  out$hap <- matrix(out$hap, nrow = length(out$pos), ncol = n_hap)
  out$seed <- seed
  out
}
