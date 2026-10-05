# Coalescent founders (breedingDesigner SPEC-0020 item 9; docs/SPEC-coalescent.md,
# DECISION-050): the internal per-chromosome sampler .coalescent_chromosome() and
# the exported founders_coalescent() with the runMacs() species presets.

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
                                   seg_sites = 0L, seed = NULL,
                                   split_n_first = 0L, split_time = 0) {
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
  split_n_first <- .validate_count(split_n_first, "split_n_first", minimum = 0L)
  if (split_n_first > 0L && (!is.numeric(split_time) || length(split_time) != 1L ||
                             !is.finite(split_time) || split_time <= 0)) {
    stop("`split_time` must be a single positive number.", call. = FALSE)
  }
  out <- coalescent_chromosome_core(as.integer(n_hap), theta, rho, ht, hs,
                                    as.integer(seg_sites), as.numeric(seed),
                                    as.integer(split_n_first),
                                    as.numeric(split_time))
  out$hap <- matrix(out$hap, nrow = length(out$pos), ncol = n_hap)
  out$seed <- seed
  out
}

#' Founder haplotypes from a coalescent with recombination
#'
#' Simulates founders with historical linkage disequilibrium and
#' mutation-drift allele frequencies, as AlphaSimR's `runMacs()` does with MaCS,
#' natively and reproducibly. Each chromosome is walked with a sequential Markov
#' coalescent (SMC', Marjoram & Wall 2006) -- the model AlphaSimR runs, since it
#' calls MaCS with a one-base history window -- under a piecewise-constant
#' population-size history, with infinite-sites mutation. The haplotypes become a
#' [Population][population_from_haplotypes()] ready for [cross()], [selfcross()],
#' [double_haploid()] and the selection engine.
#'
#' @section Model and units:
#' As in ms and MaCS, time is in units of \eqn{4N_0} generations, `theta` =
#' \eqn{4N_0\mu} and `rho` = \eqn{4N_0 r} are per chromosome, and `history`
#' gives the population size relative to \eqn{N_0} from each `time` back (ms
#' `-eN`). A pair of lineages coalesces at rate \eqn{2/\lambda(u)}; along the
#' chromosome, mutations and recombination points arrive at rates
#' \eqn{\theta L} and \eqn{\rho L} for a marginal tree of total length \eqn{L}.
#' The genetic map is linear in physical position, `morgans` long, and starts
#' at 0 cM at the first segregating site (as `runMacs()` rebases it).
#'
#' @section Species presets:
#' `species` reproduces the demographic commands `runMacs()` builds in AlphaSimR
#' 2.1.0 (read from its source; parameter values chosen by AlphaSimR, not
#' estimates made here):
#' \describe{
#'   \item{GENERIC}{\eqn{10^8} bp, \eqn{\theta = 1000}, \eqn{\rho = 400}, 1 Morgan,
#'     \eqn{N_e = 100}, sizes 5 / 15 / 60 / 120 / 1000 from times 0.25 / 2.5 / 25 /
#'     250 / 2500.}
#'   \item{MAIZE}{\eqn{2 \times 10^8} bp, \eqn{\theta = 1000}, \eqn{\rho = 800}, 2 Morgans,
#'     \eqn{N_e = 100}, a 15-step history.}
#'   \item{WHEAT}{\eqn{8 \times 10^8} bp, \eqn{\theta = 320}, \eqn{\rho = 288}, 1.43
#'     Morgans, \eqn{N_e = 50}, a 26-step history.}
#'   \item{CATTLE}{\eqn{2.8 \times 10^9 / 30} bp per chromosome, mutation rate
#'     \eqn{9.4 \times 10^{-9}} and recombination rate \eqn{9.26 \times 10^{-9}} per bp,
#'     \eqn{N_e = 90}, a 12-step history in generations.}
#' }
#' `theta`, `rho`, `history`, `morgans` and `bp` override the preset (the
#' preset's \eqn{N_e} still converts `split` to coalescent time).
#'
#' @param n_ind number of founder individuals (diploid).
#' @param n_chr number of chromosomes.
#' @param seg_sites segregating sites kept per chromosome (one value or one per
#'   chromosome), a uniform random subset of those simulated; `NULL` keeps them
#'   all. Asking for more sites than the coalescent produced is an error (as in
#'   `runMacs()`).
#' @param inbred `TRUE` for fully inbred founders: each individual carries one
#'   simulated haplotype twice (`n_ind` haplotypes are simulated, not `2 n_ind`).
#' @param species one of `"GENERIC"`, `"MAIZE"`, `"WHEAT"`, `"CATTLE"`.
#' @param split optional number of generations ago at which the founders' two
#'   halves split into isolated subpopulations (no migration): individuals
#'   `1..n_ind/2` and the rest. Converted to coalescent time as
#'   `split / (4 Ne) + 1e-6`, as `runMacs()` does. Needs an even `n_ind`, so
#'   that both haplotypes of every individual come from one subpopulation.
#' @param theta,rho optional per-chromosome \eqn{4N_0\mu} and \eqn{4N_0 r},
#'   overriding the preset.
#' @param history optional data frame with numeric columns `time` (4 N0
#'   generations, increasing, in (0, 1e12]) and `size` (relative to N0, in
#'   \[1e-9, 1e9\]); `NULL` keeps the preset's.
#' @param morgans optional genetic length of each chromosome in Morgans (one
#'   value or one per chromosome).
#' @param bp optional physical length of each chromosome in base pairs, used for
#'   the `pos` column of the map.
#' @param pool founder pool label passed to [population_from_haplotypes()].
#' @param seed optional seed (`NULL` or a whole number): with it the founders are
#'   reproducible and the caller's RNG state is restored; without it they follow
#'   the ambient stream (`set.seed()` also reproduces them).
#' @return a `Population` of `n_ind` founders. Allele `1` (the counted allele)
#'   is the derived allele. Its attribute `"coalescent"` records the parameters,
#'   the per-chromosome seeds and the number of segregating sites simulated
#'   before subsetting (`n_total`).
#' @references
#' Marjoram P, Wall JD (2006) Fast "coalescent" simulation. \emph{BMC Genet} 7:16
#' (PMID 16539698).
#' Chen GK, Marjoram P, Wall JD (2009) Fast and flexible simulation of DNA
#' sequence data. \emph{Genome Res} 19:136-142 (PMID 19029539).
#' Hudson RR (2002) Generating samples under a Wright-Fisher neutral model of
#' genetic variation. \emph{Bioinformatics} 18:337-338 (PMID 11847089).
#' Gaynor RC, Gorjanc G, Hickey JM (2021) AlphaSimR: an R package for breeding
#' program simulations. \emph{G3} 11(2):jkaa017. \doi{10.1093/g3journal/jkaa017}
#' (the `runMacs()` presets).
#' @seealso [population_from_haplotypes()], [as_population()], [cross()].
#' @export
#' @examples
#' # 50 founders, 2 chromosomes, 200 segregating sites each (GENERIC preset)
#' f <- founders_coalescent(50, n_chr = 2, seg_sites = 200, seed = 1)
#' f
#' # inbred lines, e.g. for doubled-haploid programs
#' lines <- founders_coalescent(20, seg_sites = 100, inbred = TRUE, seed = 2)
#' all(dosages(lines) != 0)
founders_coalescent <- function(n_ind, n_chr = 1, seg_sites = NULL,
                                inbred = FALSE,
                                species = c("GENERIC", "MAIZE", "WHEAT", "CATTLE"),
                                split = NULL, theta = NULL, rho = NULL,
                                history = NULL, morgans = NULL, bp = NULL,
                                pool = NA_character_, seed = NULL) {
  n_ind <- .validate_count(n_ind, "n_ind", minimum = 1L)
  n_chr <- .validate_count(n_chr, "n_chr", minimum = 1L)
  .validate_flag(inbred, "inbred")
  species <- match.arg(species)
  pre <- .coalescent_preset(species)
  n_hap <- if (inbred) n_ind else 2L * n_ind
  if (n_hap < 2L) {
    stop("founders_coalescent(): at least two haplotypes are needed (n_ind >= 2 ",
         "when inbred = TRUE).", call. = FALSE)
  }
  per_chr <- function(x, arg) {
    if (is.null(x)) return(NULL)
    if (!is.numeric(x) || !length(x) %in% c(1L, n_chr) || any(!is.finite(x)) ||
        any(x <= 0)) {
      stop("founders_coalescent(): `", arg, "` must be positive, one value or one ",
           "per chromosome (", n_chr, ").", call. = FALSE)
    }
    rep_len(as.numeric(x), n_chr)
  }
  # (no %||%: base R has it only from 4.4, the package supports 4.2)
  theta <- if (is.null(theta)) rep(pre$theta, n_chr) else per_chr(theta, "theta")
  rho <- if (is.null(rho)) rep(pre$rho, n_chr) else per_chr(rho, "rho")
  morgans <- if (is.null(morgans)) rep(pre$morgans, n_chr) else per_chr(morgans, "morgans")
  bp <- if (is.null(bp)) rep(pre$bp, n_chr) else per_chr(bp, "bp")
  if (is.null(history)) history <- pre$history
  if (is.null(seg_sites)) {
    keep <- rep(0L, n_chr)
  } else {
    if (!is.numeric(seg_sites) || !length(seg_sites) %in% c(1L, n_chr) ||
        any(!is.finite(seg_sites)) || any(seg_sites < 1) ||
        any(seg_sites != floor(seg_sites)) || any(seg_sites > .Machine$integer.max)) {
      stop("founders_coalescent(): `seg_sites` must be NULL or positive whole ",
           "numbers, one value or one per chromosome (", n_chr, ").", call. = FALSE)
    }
    keep <- rep_len(as.integer(seg_sites), n_chr)
  }
  split_n_first <- 0L
  split_time <- 0
  if (!is.null(split)) {
    if (!is.numeric(split) || length(split) != 1L || !is.finite(split) || split <= 0) {
      stop("founders_coalescent(): `split` must be a single positive number of ",
           "generations.", call. = FALSE)
    }
    if (n_ind %% 2L != 0L) {
      stop("founders_coalescent(): a `split` needs an even `n_ind` (the two ",
           "subpopulations are the halves of the individuals).", call. = FALSE)
    }
    split_n_first <- n_hap %/% 2L
    split_time <- split / (4 * pre$ne) + 1e-6
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != floor(seed) || seed < 0 || seed > .Machine$integer.max) {
      stop("founders_coalescent(): `seed` must be NULL or one non-negative whole ",
           "number no larger than .Machine$integer.max.", call. = FALSE)
    }
    old <- .Random.seed_safe()
    set.seed(seed)
    on.exit(.restore_seed(old), add = TRUE)
  }
  chr_seeds <- floor(stats::runif(n_chr, 0, 2^32))

  haps <- vector("list", n_chr)
  maps <- vector("list", n_chr)
  n_total <- numeric(n_chr)
  for (k in seq_len(n_chr)) {
    x <- .coalescent_chromosome(n_hap, theta[k], rho[k], history = history,
                                seg_sites = keep[k], seed = chr_seeds[k],
                                split_n_first = split_n_first,
                                split_time = split_time)
    n_total[k] <- x$n_total
    s <- length(x$pos)
    if (s == 0L || (keep[k] > 0L && s < keep[k])) {
      stop("founders_coalescent(): chromosome ", k, " has ", x$n_total,
           " segregating site(s), fewer than ",
           if (keep[k] > 0L) paste0("the ", keep[k], " requested in `seg_sites`")
           else "one", "; raise `theta` or lower `seg_sites`.", call. = FALSE)
    }
    haps[[k]] <- x$hap
    maps[[k]] <- data.frame(
      snp = paste0(k, "_", seq_len(s)), chr = k, pos = x$pos * bp[k],
      cm = 100 * morgans[k] * (x$pos - x$pos[1L]), stringsAsFactors = FALSE)
  }
  hap <- do.call(rbind, haps)
  map <- do.call(rbind, maps)
  if (inbred) {
    cis <- trans <- hap
  } else {
    cis <- hap[, seq(1L, n_hap, by = 2L), drop = FALSE]
    trans <- hap[, seq(2L, n_hap, by = 2L), drop = FALSE]
  }
  ids <- as.character(seq_len(n_ind))
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)
  pop <- population_from_haplotypes(cis, trans, map, ids = ids, pool = pool)
  attr(pop, "coalescent") <- list(
    species = species, n_hap = n_hap, inbred = inbred, theta = theta, rho = rho,
    history = history, morgans = morgans, bp = bp, ne = pre$ne, split = split,
    split_time = if (is.null(split)) NULL else split_time, seeds = chr_seeds,
    n_total = n_total)
  pop
}

#' runMacs() species presets (AlphaSimR 2.1.0 source), per chromosome
#'
#' theta and rho are the per-bp ms `-t` / `-r` values times the sequence length;
#' history times are in 4 N0 generations and sizes relative to N0.
#' @keywords internal
#' @noRd
.coalescent_preset <- function(species) {
  switch(species,
    GENERIC = list(
      bp = 1e8, theta = 1e8 * 1e-5, rho = 1e8 * 4e-6, morgans = 1, ne = 100,
      history = data.frame(time = c(0.25, 2.5, 25, 250, 2500),
                           size = c(5, 15, 60, 120, 1000))),
    MAIZE = list(
      bp = 2e8, theta = 2e8 * 5e-6, rho = 2e8 * 4e-6, morgans = 2, ne = 100,
      history = data.frame(
        time = c(0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                 0.50, 2.00, 3.00, 4.00, 5.00),
        size = c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 40, 60, 80, 100))),
    WHEAT = list(
      bp = 8e8, theta = 8e8 * 4e-7, rho = 8e8 * 3.6e-7, morgans = 1.43, ne = 50,
      history = data.frame(
        time = c(0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
                 0.50, 1.00, 2.00, 3.00, 4.00, 5.00, 10.00, 20.00, 30.00, 40.00,
                 50.00, 100.00, 200.00, 300.00, 400.00, 500.00),
        size = c(1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 40, 60, 80, 100, 120,
                 140, 160, 180, 200, 240, 320, 400, 480, 560, 640))),
    CATTLE = {
      chr_bp <- 2.8e9 / 30
      ne <- 90
      rec <- 9.26e-9
      mut <- 9.4e-9
      list(bp = chr_bp, theta = chr_bp * mut * 4 * ne, rho = chr_bp * rec * 4 * ne,
           morgans = rec * chr_bp, ne = ne,
           history = data.frame(
             time = c(3, 6, 12, 18, 24, 154, 454, 654, 1754, 2354, 3354, 33154) / (4 * ne),
             size = c(120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000,
                      17000, 62000) / ne))
    })
}
