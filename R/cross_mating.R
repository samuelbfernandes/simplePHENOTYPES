# Multi-generation crossing.
#
# R draws every random quantity here, in the order isqg draws them, and the
# Rust core performs only the deterministic remainder. Keeping
# the draws on R's RNG is what makes set.seed() reproducible and what lets
# tests/testthat/test-isqg-parity.R assert exact agreement with isqg.

#' Citation notice for the crossing pipeline (isqg algorithms), once per session
#'
#' The meiosis / crossing / double-haploid algorithms are ported from isqg
#' (Toledo et al. 2019), so credit it the first time a crossing function runs in
#' a session. An [rlang::inform()] message, silenceable with
#' `suppressMessages()`.
#' @keywords internal
#' @noRd
.cite_isqg <- function() {
  rlang::inform(
    .cite_main("Toledo, F.H., Perez-Rodriguez, P., Crossa, J. and Burgueno, J. ",
               "(2019). isqg: A Binary Framework for in Silico Quantitative ",
               "Genetics. G3 9(8):2425-2428, doi:10.1534/g3.119.400373, when ",
               "using the crossing pipelines (cross, selfcross, ",
               "double_haploid)."),
    .frequency = "once",
    .frequency_id = "simplePHENOTYPES_isqg_citation"
  )
}

#' Draw the randomness for `n_events` whole-genome meiosis events
#'
#' Per event, per chromosome in ascending order:
#'   n_x       ~ rpois(1, L)              L = LAST map position, in Morgans
#'   chiasmata ~ sort(runif(n_x, 0, L))   not drawn when n_x == 0
#'   flip      ~ rbinom(1, 1, 0.5)        ALWAYS drawn, even when n_x == 0
#'
#' `L` is the chromosome's last map position (`cm / 100`), not its span
#' `last - first`: positions are absolute and are not rebased to a zero origin
#' (the isqg convention, DECISION-012). The two agree when the first marker is at
#' 0; otherwise chiasmata upstream of the first marker exist but only toggle the
#' whole chromosome, which the flip absorbs, so the recombination between markers
#' is still Haldane's.
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

#' Check what is about to be sent to the Rust core
#'
#' The kernel validates everything it receives and returns an error rather than
#' panicking, but a malformed call is cheaper to explain here, where the
#' arguments have names.
#' @keywords internal
#' @noRd
.check_meiosis_call <- function(loci_per_chr, positions, strands, draws, n,
                                events_per) {
  n_loci <- sum(loci_per_chr)
  if (!length(loci_per_chr) || any(loci_per_chr < 1L) ||
      length(positions) != n_loci || any(!is.finite(positions))) {
    stop("Internal error: the chromosome layout does not match the map ",
         "positions.", call. = FALSE)
  }
  for (nm in names(strands)) {
    s <- strands[[nm]]
    if (length(s) != 1L || is.na(s) || nchar(s) != n_loci ||
        !grepl("^[01]*$", s)) {
      stop("The parental strand `", nm, "` must have exactly one 0/1 entry ",
           "per marker (", n_loci, "); the genotypes hold missing or ",
           "out-of-range values.", call. = FALSE)
    }
  }
  n_events <- n * events_per
  if (length(draws$counts) != length(loci_per_chr) * n_events ||
      length(draws$flips) != length(draws$counts) ||
      sum(draws$counts) != length(draws$chiasmata)) {
    stop("Internal error: the drawn meiosis events (", length(draws$counts),
         " counts) do not match ", n, " progeny (", n_events,
         " events x ", length(loci_per_chr), " chromosomes).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Run one mating design through the Rust core
#' @keywords internal
#' @noRd
.mate <- function(p1, p2, n, design, seed, origin, prefix) {
  .cite_isqg()
  n <- .validate_count(n, "n", minimum = 1L)
  seed <- .validate_seed(seed)

  if (!.same_map(p1$map, p2$map)) {
    stop("The two parents carry different marker maps; they must come from ",
         "the same Population.", call. = FALSE)
  }
  .check_orientation(p1$map, p2$map)

  map <- p1$map
  # Match isqg exactly (DECISION-012). isqg sorts the map by (chr, pos) and
  # builds BOTH the map and the parental haplotypes in that order, and takes each
  # chromosome's length as its LAST map position in Morgans -- it does not rebase
  # to a zero origin. So:
  #   * sort markers by (chr, cm) and send the haplotype bits in that same order
  #     (grouping/sorting only the map, as the old code did, mis-assigned
  #     chromosome masks whenever the rows were interleaved or unsorted);
  #   * use absolute positions cm/100 with L = last position (a `(cm - min)/100`
  #     rebase changed both the Poisson mean and the phantom-crossover region
  #     before the first marker, breaking parity on nonzero-origin maps);
  #   * group by run length of the sorted chromosome labels, which ignores unused
  #     factor levels (they otherwise produced empty chromosomes that crashed);
  #   * order the chromosomes by `.chr_rank()` (numeric-aware, locale-free), so
  #     the random draws do not depend on the storage type of `chr` (integer 1, 2,
  #     10 versus text "1", "10", "2") or on the collation locale.
  # The caller's original marker order is restored on the progeny via `inv`.
  ord <- order(.chr_rank(map$chr), map$cm)
  inv <- order(ord)
  chr_s <- as.character(map$chr[ord])
  cm_s  <- map$cm[ord]
  grp <- rle(chr_s)$lengths
  by_chr <- split(cm_s / 100, rep.int(seq_along(grp), grp))
  loci_per_chr <- as.integer(grp)
  positions <- as.numeric(cm_s / 100)

  events_per <- if (design == "dh") 1L else 2L

  if (!is.null(seed)) {
    # the ambient RNG is restored on exit: a seeded call must not disturb the
    # caller's stream
    old_seed <- .Random.seed_safe()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  # The RNG state the meioses are drawn from identifies this mating in the
  # pedigree keys (two matings whose draws happen to coincide, e.g. on a 0 cM
  # map, still get distinct keys); reading it draws nothing.
  rng_state <- .Random.seed_safe()
  draws <- .draw_meiosis(by_chr, n * events_per)
  rng_after <- .Random.seed_safe()

  bits <- function(v) paste(as.integer(v), collapse = "")
  # Ask for haplotypes, not genotypes. A -1/0/1 genotype cannot express the
  # phase of a heterozygote, so deriving progeny strands from one would give an
  # F1 - heterozygous at every locus - a fictitious all-allele-1 / all-allele-2
  # pair, and every later generation would recombine haplotypes that never
  # existed.
  strand_codes <- list(
    p1_cis   = bits(p1$cis[ord, 1]),
    p1_trans = bits(p1$trans[ord, 1]),
    p2_cis   = bits(p2$cis[ord, 1]),
    p2_trans = bits(p2$trans[ord, 1])
  )
  .check_meiosis_call(loci_per_chr, positions, strand_codes, draws, n,
                      events_per)
  strands <- mate_haplotypes_core(
    loci_per_chr = loci_per_chr,
    positions    = positions,
    p1_cis       = strand_codes$p1_cis,
    p1_trans     = strand_codes$p1_trans,
    p2_cis       = strand_codes$p2_cis,
    p2_trans     = strand_codes$p2_trans,
    chiasmata    = draws$chiasmata,
    counts       = draws$counts,
    flips        = draws$flips,
    design       = design,
    n_prog       = n
  )

  # Element 2i-1 is progeny i's first strand, 2i its second. The Rust core
  # returns strands in the sorted marker order; `inv` restores the caller's.
  unpack <- function(codes) {
    m <- matrix(as.integer(unlist(strsplit(codes, "", fixed = TRUE))),
                nrow = nrow(map), ncol = n)
    m[inv, , drop = FALSE]
  }
  ids <- paste0(prefix, seq_len(n))
  cis <- unpack(strands[seq(1, length(strands), by = 2)])
  trans <- unpack(strands[seq(2, length(strands), by = 2)])
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)

  # Pedigree bookkeeping only; it draws nothing, so the RNG stream is unchanged.
  mp <- .mating_pedigree(p1, p2, design,
                         list(rng_state = rng_state, rng_after = rng_after,
                              draws = draws), ids)
  .new_population(map, cis, trans, ids, origin, keys = mp$keys,
                  pedigree = mp$pedigree)
}

#' Cross two individuals
#'
#' Produces `n` progeny from a biparental cross. Each progeny receives one
#' recombinant gamete from each parent, so the two parents contribute one
#' homologue apiece.
#'
#' Recombination follows the count-location model: the number of crossovers on
#' a chromosome is Poisson with mean equal to its length in Morgans, and their
#' positions are uniform along it. That length is the chromosome's **last** map
#' position (`cm / 100`), not its span `max(cm) - min(cm)`: positions are used as
#' given and are not rebased to a zero origin (the isqg convention). The two
#' agree when a chromosome's first marker is at 0; when it is not, the extra
#' crossovers fall upstream of the first marker and only swap the whole
#' chromosome, so the recombination fraction between markers is still
#' Haldane's, \eqn{(1 - e^{-2d})/2}. Chromosomes assort independently. The
#' genetic map is taken from the `cm` column of the population's marker map, so
#' it must not be missing - see [synthetic_map()] if you only have physical
#' positions.
#'
#' Both parents must carry the same marker map (marker names, chromosomes and
#' positions; the `allele` column is not part of it). The dosages are relative
#' to the allele coded `+1`, which [as_numeric()] chooses per data set; crossing
#' populations built from **separately** converted panels can therefore mix up
#' the alleles. Convert the panels together, or with
#' `as_numeric(method = "reference", ref_allele = )`. Where both populations
#' record the `allele` column, a marker whose alleles are listed in opposite
#' order draws a warning and markers with no allele in common an error (see
#' [as_population()]).
#'
#' Chromosomes are processed, and their random draws consumed, in a canonical
#' order that does not depend on the storage type of `chr` or on the locale (see
#' [as_population()]).
#'
#' @param mother,father single-individual `Population`s (use `[` to select one).
#'   Their roles are symmetric apart from which homologue a progeny inherits
#'   first; there is no sex-specific recombination. Crossing an individual with
#'   itself is a self: it draws what [selfcross()] draws and the pedigree records
#'   it as one (`design = "self"`, see [parentage()]).
#' @param n number of progeny.
#' @param seed optional RNG seed. All randomness is drawn in R, so `set.seed()`
#'   before the call works equally well. With a `seed` the caller's RNG state is
#'   restored on exit (a seeded call does not disturb the ambient stream); with
#'   `seed = NULL` the draws consume the ambient stream.
#' @return A `Population` of `n` progeny.
#' @seealso [selfcross()], [double_haploid()], [as_population()]
#' @references
#' Toledo, F.H., Perez-Rodriguez, P., Crossa, J. and Burgueno, J. (2019). isqg:
#' A Binary Framework for in Silico Quantitative Genetics. \emph{G3
#' Genes|Genomes|Genetics} 9(8), 2425--2428. \doi{10.1534/g3.119.400373}
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
#' @references
#' Toledo, F.H., Perez-Rodriguez, P., Crossa, J. and Burgueno, J. (2019). isqg:
#' A Binary Framework for in Silico Quantitative Genetics. \emph{G3
#' Genes|Genomes|Genetics} 9(8), 2425--2428. \doi{10.1534/g3.119.400373}
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
#' @references
#' Toledo, F.H., Perez-Rodriguez, P., Crossa, J. and Burgueno, J. (2019). isqg:
#' A Binary Framework for in Silico Quantitative Genetics. \emph{G3
#' Genes|Genomes|Genetics} 9(8), 2425--2428. \doi{10.1534/g3.119.400373}
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
