# Multi-generation crossing.
#
# R draws every random quantity here, in the order isqg draws them, and the
# Rust core performs only the deterministic remainder (DECISION-012). Keeping
# the draws on R's RNG is what makes set.seed() reproducible and what lets
# tests/testthat/test-isqg-parity.R assert exact agreement with isqg.

#' Draw the randomness for `n_events` whole-genome meiosis events
#'
#' Per event, per chromosome in ascending order:
#'   n_x       ~ rpois(1, L)              L = last map position, in Morgans
#'   chiasmata ~ sort(runif(n_x, 0, L))   not drawn when n_x == 0
#'   flip      ~ rbinom(1, 1, 0.5)        ALWAYS drawn, even when n_x == 0
#'
#' The flip is unconditional: skipping it for crossover-free chromosomes would
#' desynchronise every later draw. Returns the flat (event, chromosome)-ordered
#' vectors the Rust core expects.
#' @keywords internal
#' @noRd
.draw_meiosis <- function(morgans_by_chr, n_events) {
  n_chr <- length(morgans_by_chr)
  counts <- integer(n_events * n_chr)
  flips <- integer(n_events * n_chr)
  chiasmata <- vector("list", n_events * n_chr)

  slot <- 0L
  for (e in seq_len(n_events)) {
    for (c in seq_len(n_chr)) {
      slot <- slot + 1L
      pos <- morgans_by_chr[[c]]
      len <- pos[[length(pos)]]
      k <- stats::rpois(1, len)
      chiasmata[[slot]] <- if (k > 0) sort(stats::runif(k, 0, len)) else numeric(0)
      counts[[slot]] <- k
      flips[[slot]] <- stats::rbinom(1, 1, 0.5)
    }
  }

  list(
    counts    = counts,
    flips     = flips,
    chiasmata = as.numeric(unlist(chiasmata))
  )
}

#' Run one mating design through the Rust core
#' @keywords internal
#' @noRd
.mate <- function(p1, p2, n, design, seed, origin, prefix) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a single positive number of progeny; got ",
         deparse(substitute(n)), ".", call. = FALSE)
  }
  n <- as.integer(n)

  if (!identical(p1$map$snp, p2$map$snp)) {
    stop("The two parents carry different marker maps; they must come from ",
         "the same Population.", call. = FALSE)
  }

  map <- p1$map
  # cm -> Morgans: isqg's Poisson mean is the chromosome length in Morgans.
  by_chr <- split(map$cm / 100, map$chr)
  loci_per_chr <- as.integer(vapply(by_chr, length, integer(1)))
  positions <- as.numeric(unlist(by_chr))

  events_per <- if (design == "dh") 1L else 2L

  if (!is.null(seed)) {
    set.seed(seed)
  }
  draws <- .draw_meiosis(by_chr, n * events_per)

  bits <- function(v) paste(as.integer(v), collapse = "")
  # Ask for haplotypes, not genotypes. A -1/0/1 genotype cannot express the
  # phase of a heterozygote, so deriving progeny strands from one would give an
  # F1 - heterozygous at every locus - a fictitious all-allele-1 / all-allele-2
  # pair, and every later generation would recombine haplotypes that never
  # existed.
  strands <- mate_haplotypes_core(
    loci_per_chr = loci_per_chr,
    positions    = positions,
    p1_cis       = bits(p1$cis[, 1]),
    p1_trans     = bits(p1$trans[, 1]),
    p2_cis       = bits(p2$cis[, 1]),
    p2_trans     = bits(p2$trans[, 1]),
    chiasmata    = draws$chiasmata,
    counts       = draws$counts,
    flips        = draws$flips,
    design       = design,
    n_prog       = n
  )

  # Element 2i-1 is progeny i's first strand, 2i its second.
  unpack <- function(codes) {
    matrix(as.integer(unlist(strsplit(codes, "", fixed = TRUE))),
           nrow = nrow(map), ncol = n)
  }
  ids <- paste0(prefix, seq_len(n))
  cis <- unpack(strands[seq(1, length(strands), by = 2)])
  trans <- unpack(strands[seq(2, length(strands), by = 2)])
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)

  .new_population(map, cis, trans, ids, origin)
}

#' Cross two individuals
#'
#' Produces `n` progeny from a biparental cross. Each progeny receives one
#' recombinant gamete from each parent, so the two parents contribute one
#' homologue apiece.
#'
#' Recombination follows the count-location model: the number of crossovers on
#' a chromosome is Poisson with mean equal to its length in Morgans, and their
#' positions are uniform along it. Chromosomes assort independently. The
#' genetic map is taken from the `cm` column of the population's marker map, so
#' it must not be missing - see [synthetic_map()] if you only have physical
#' positions.
#'
#' @param mother,father single-individual `Population`s (use `[` to select one).
#'   Their roles are symmetric apart from which homologue a progeny inherits
#'   first; there is no sex-specific recombination.
#' @param n number of progeny.
#' @param seed optional RNG seed. All randomness is drawn in R, so `set.seed()`
#'   before the call works equally well.
#' @return A `Population` of `n` progeny.
#' @seealso [selfcross()], [double_haploid()], [as_population()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
#' f1 <- cross(pop[1], pop[2], n = 5, seed = 1)
#' f1
cross <- function(mother, father, n = 1, seed = NULL) {
  .check_single(mother, "mother")
  .check_single(father, "father")
  .mate(mother, father, n, "cross", seed,
        origin = paste0("cross(", mother$ids, " x ", father$ids, ")"),
        prefix = "prog_")
}

#' Self-pollinate an individual
#'
#' Produces `n` progeny by selfing: both gametes come from independent meioses
#' of the same individual. Selfing a heterozygous individual halves
#' heterozygosity each generation, so repeated selfing drives a line toward
#' homozygosity.
#'
#' @inheritParams cross
#' @param parent a single-individual `Population`.
#' @return A `Population` of `n` progeny.
#' @seealso [cross()], [double_haploid()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
#' f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
#' f2 <- selfcross(f1, n = 10, seed = 2)
#' f2
selfcross <- function(parent, n = 1, seed = NULL) {
  .check_single(parent, "parent")
  .mate(parent, parent, n, "selfcross", seed,
        origin = paste0("selfcross(", parent$ids, ")"),
        prefix = "self_")
}

#' Produce doubled haploids from an individual
#'
#' Produces `n` doubled-haploid progeny: a single recombinant gamete is drawn
#' and then doubled, so every individual is completely homozygous at every
#' marker. This is the one-generation route to a fully inbred line, and the
#' result contains no heterozygotes at all.
#'
#' @inheritParams selfcross
#' @return A `Population` of `n` fully homozygous progeny.
#' @seealso [cross()], [selfcross()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
#' f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
#' dh <- double_haploid(f1, n = 10, seed = 3)
#' # No heterozygotes by construction:
#' any(dosages(dh) == 0)
double_haploid <- function(parent, n = 1, seed = NULL) {
  .check_single(parent, "parent")
  .mate(parent, parent, n, "dh", seed,
        origin = paste0("double_haploid(", parent$ids, ")"),
        prefix = "dh_")
}
