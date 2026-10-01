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
#'
#' With `interference = NULL` (the default) this is exactly the isqg stream
#' above, draw for draw. A non-NULL `interference` (a validated
#' `list(nu, p)`, see `.check_interference()`) switches to the two-pathway gamma
#' model of `.draw_meiosis_interference()`, which has its own stream.
#' @keywords internal
#' @noRd
.draw_meiosis <- function(morgans_by_chr, n_events, interference = NULL) {
  if (!is.null(interference)) {
    return(.draw_meiosis_interference(morgans_by_chr, n_events,
                                      interference$nu, interference$p))
  }
  # Bind the three generators once: `stats::rpois` re-resolves the namespace on
  # every call, which was a third of the time of this loop. The sequence of RNG
  # calls (rpois, then runif only when k > 0, then rbinom, per chromosome per
  # event) is unchanged, so the stream is bit-identical to isqg's.
  rpois <- stats::rpois
  runif <- stats::runif
  rbinom <- stats::rbinom
  n_chr <- length(morgans_by_chr)
  last <- vapply(morgans_by_chr, function(pos) pos[[length(pos)]], numeric(1))
  n_slots <- n_events * n_chr
  counts <- integer(n_slots)
  flips <- integer(n_slots)
  chiasmata <- vector("list", n_slots)

  slot <- 0L
  for (e in seq_len(n_events)) {
    for (c in seq_len(n_chr)) {
      slot <- slot + 1L
      len <- last[[c]]
      k <- rpois(1L, len)
      if (k > 0L) {
        u <- runif(k, 0, len)
        # sorted ascending (ties are irrelevant); the common k = 1, 2 cases
        # avoid the generic sort
        if (k > 2L) {
          u <- sort.int(u, method = "radix")
        } else if (k == 2L && u[[1L]] > u[[2L]]) {
          u <- u[2:1]
        }
        chiasmata[[slot]] <- u
        counts[[slot]] <- k
      }
      flips[[slot]] <- rbinom(1L, 1L, 0.5)
    }
  }

  list(
    counts    = counts,
    flips     = flips,
    chiasmata = as.numeric(unlist(chiasmata, use.names = FALSE))
  )
}

#' Validate and normalise the `interference` argument of the crossing functions
#'
#' `NULL` (no interference option: Poisson chiasmata, the isqg stream) or a
#' list / named numeric vector with `nu` (a single number >= 1) and optionally
#' `p` (a single number in `[0, 1]`, default 0). Returns `NULL` or
#' `list(nu =, p =)`.
#' @param x the argument as given.
#' @param fn name of the calling function, for the message.
#' @keywords internal
#' @noRd
.check_interference <- function(x, fn = "cross") {
  if (is.null(x)) return(NULL)
  if (!(is.list(x) || (is.numeric(x) && !is.null(names(x)))) ||
      is.data.frame(x)) {
    stop(fn, "(): `interference` must be NULL or a list(nu =, p =) (a named ",
         "numeric vector also works).", call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm) ||
      !all(nm %in% c("nu", "p")) || !"nu" %in% nm) {
    stop(fn, "(): `interference` must be named with `nu` (required) and ",
         "optionally `p`, e.g. list(nu = 2.6, p = 0).", call. = FALSE)
  }
  num1 <- function(v) {
    is.numeric(v) && length(v) == 1L && is.finite(v)
  }
  nu <- x[["nu"]]
  p <- if ("p" %in% nm) x[["p"]] else 0
  if (!num1(nu) || nu < 1) {
    stop(fn, "(): `interference$nu` must be one finite number >= 1 (the shape ",
         "of the gamma model; 1 is no interference); got ",
         paste(format(nu), collapse = ", "), ".", call. = FALSE)
  }
  if (!num1(p) || p < 0 || p > 1) {
    stop(fn, "(): `interference$p` must be one number in [0, 1] (the share of ",
         "chiasmata from the non-interfering pathway); got ",
         paste(format(p), collapse = ", "), ".", call. = FALSE)
  }
  list(nu = as.numeric(nu), p = as.numeric(p))
}

#' Draw meiosis events under the two-pathway gamma interference model
#'
#' The same return value as the Poisson branch of `.draw_meiosis()` (flat
#' `counts`, `flips`, `chiasmata` in (event, chromosome) order, `chiasmata`
#' sorted within each slot), for a chromosome of last map position `L`
#' (Morgans) and the parameters `nu` (>= 1) and `p` (in [0, 1]).
#'
#' Model (chiasmata on the four-strand bivalent, no chromatid interference):
#'  * the bivalent carries chiasmata at intensity 2 per Morgan, so a gamete,
#'    which takes part in each chiasma independently with probability 1/2,
#'    carries 1 crossover per Morgan on average (the Haldane / isqg convention);
#'  * a fraction `p` of the chiasmata (intensity `2p`) is non-interfering: a
#'    Poisson process;
#'  * the remaining intensity `2(1 - p)` is interfering: a stationary renewal
#'    process whose gaps are Gamma(shape `nu`, rate `2 nu (1 - p)`) (mean gap
#'    `1 / (2 (1 - p))` Morgans);
#'  * the two pathways are independent and their chiasmata superposed.
#'
#' Stationary start: the first interfering chiasma after the chromosome origin
#' is the forward recurrence time of the renewal process, `U * G` with
#' `U ~ Unif(0, 1)` and `G ~ Gamma(nu + 1, rate)` (the interval covering the
#' origin is size-biased, i.e. Gamma(nu + 1), and the origin is uniform in it).
#' Thinning by 1/2 is done per chiasma, and the Poisson pathway is drawn
#' directly at the gamete level (a Poisson process of intensity `2p` thinned by
#' 1/2 is Poisson with intensity `p`). Everything is vectorised over the
#' (event, chromosome) slots.
#' @keywords internal
#' @noRd
.draw_meiosis_interference <- function(morgans_by_chr, n_events, nu, p) {
  n_chr <- length(morgans_by_chr)
  last <- vapply(morgans_by_chr, function(pos) pos[[length(pos)]], numeric(1))
  n_slots <- n_events * n_chr
  len <- rep(last, times = n_events)          # slot = (event - 1) * n_chr + chr
  slot_parts <- list()
  pos_parts <- list()

  # non-interfering pathway: Poisson, gamete-level intensity p per Morgan
  if (p > 0) {
    k <- stats::rpois(n_slots, p * len)
    s <- rep.int(seq_len(n_slots), k)
    slot_parts[[1L]] <- s
    pos_parts[[1L]] <- stats::runif(length(s), 0, len[s])
  }

  # interfering pathway: stationary gamma renewal process on the bivalent,
  # each chiasma kept in the gamete with probability 1/2
  if (p < 1) {
    rate <- 2 * nu * (1 - p)
    x <- stats::runif(n_slots) * stats::rgamma(n_slots, shape = nu + 1,
                                              rate = rate)
    act <- which(x <= len)
    sl <- vector("list", 0L)
    ps <- vector("list", 0L)
    while (length(act)) {
      sl[[length(sl) + 1L]] <- act
      ps[[length(ps) + 1L]] <- x[act]
      x[act] <- x[act] + stats::rgamma(length(act), shape = nu, rate = rate)
      act <- act[x[act] <= len[act]]
    }
    sl <- unlist(sl, use.names = FALSE)
    ps <- unlist(ps, use.names = FALSE)
    keep <- stats::runif(length(sl)) < 0.5
    slot_parts[[length(slot_parts) + 1L]] <- sl[keep]
    pos_parts[[length(pos_parts) + 1L]] <- ps[keep]
  }

  slot_id <- unlist(slot_parts, use.names = FALSE)
  pos <- unlist(pos_parts, use.names = FALSE)
  if (is.null(slot_id)) {
    slot_id <- integer(0)
    pos <- numeric(0)
  }
  o <- order(slot_id, pos)
  list(
    counts    = tabulate(slot_id, nbins = n_slots),
    flips     = stats::rbinom(n_slots, 1L, 0.5),
    chiasmata = as.numeric(pos[o])
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

#' Rust layout of a marker map (isqg order, DECISION-012)
#'
#' isqg sorts the map by (chr, pos) and builds BOTH the map and the parental
#' haplotypes in that order, and takes each chromosome's length as its LAST map
#' position in Morgans -- it does not rebase to a zero origin. So:
#'   * markers are ordered by (chr, cm), and the Rust core gets the permutation
#'     `ord` (the caller's index of the marker of each map rank) instead of
#'     reordered strings, returning progeny in the caller's marker order;
#'   * absolute positions `cm / 100` with L = last position (a `(cm - min)/100`
#'     rebase changed both the Poisson mean and the phantom-crossover region
#'     before the first marker, breaking parity on nonzero-origin maps);
#'   * chromosomes are grouped by run length of the sorted labels, which ignores
#'     unused factor levels (they otherwise produced empty chromosomes that
#'     crashed);
#'   * chromosomes are ordered by `.chr_rank()` (numeric-aware, locale-free), so
#'     the random draws do not depend on the storage type of `chr` (integer 1, 2,
#'     10 versus text "1", "10", "2") or on the collation locale.
#' @return `list(ord, by_chr, loci_per_chr, positions)`.
#' @keywords internal
#' @noRd
.meiosis_layout <- function(map) {
  ord <- order(.chr_rank(map$chr), map$cm)
  chr_s <- as.character(map$chr[ord])
  cm_s  <- map$cm[ord]
  grp <- rle(chr_s)$lengths
  list(
    ord          = ord,
    by_chr       = split(cm_s / 100, rep.int(seq_along(grp), grp)),
    loci_per_chr = as.integer(grp),
    positions    = as.numeric(cm_s / 100)
  )
}

#' Run a batch of matings through the Rust core in one call
#'
#' The vectorised entry point behind every crossing function (SPEC-0020 item 2):
#' `.mate()` (one mating: `cross()`, `selfcross()`, `double_haploid()`) and
#' `mate()` (one row per mating) both come here. The matings are the rows of
#' `mating`/`n`/`design`; their meioses are drawn **in row order, each row's
#' draws exactly those `.draw_meiosis()` would make for it alone** (R's RNG
#' stream is consumed row after row), and the Rust core -- which never draws --
#' turns the pooled events into progeny with a single integer-in, integer-out
#' call. A batch is therefore the same computation, draw for draw and bit for
#' bit, as running the rows one after the other; it only avoids the per-call
#' serialisation of the parents to bit strings, the parsing of the progeny back
#' and the pooling of one Population per row.
#'
#' @param map the shared marker map.
#' @param strands integer matrix, markers (in `map` order) by parental strands,
#'   0/1; may carry column names (used in error messages).
#' @param mating integer matrix, one row per mating, columns the strand indices
#'   (into `strands`) of `p1_cis, p1_trans, p2_cis, p2_trans`.
#' @param n progeny per mating (integer vector).
#' @param design `"cross"`, `"selfcross"` or `"dh"` per mating.
#' @param interference `NULL` or a validated `list(nu, p)`.
#' @return `list(cis, trans, rng_state, rng_after, draws)`: the progeny strands
#'   (integer matrices, markers by progeny, mating after mating, no dimnames)
#'   and, per mating, the RNG state before and after its draws and the draws.
#' @keywords internal
#' @noRd
.mate_many <- function(map, strands, mating, n, design, interference = NULL) {
  n_mating <- length(n)
  lay <- .meiosis_layout(map)
  if (!is.integer(strands)) storage.mode(strands) <- "integer"
  if (anyNA(strands) || any(strands != 0L & strands != 1L)) {
    j <- which(colSums(is.na(strands) | (strands != 0L & strands != 1L)) > 0L)[1L]
    nm <- if (!is.null(colnames(strands))) colnames(strands)[j] else paste0("#", j)
    stop("The parental strand `", nm, "` must have exactly one 0/1 entry ",
         "per marker (", nrow(strands), "); the genotypes hold missing or ",
         "out-of-range values.", call. = FALSE)
  }
  events_per <- ifelse(design == "dh", 1L, 2L)

  # The RNG state the meioses are drawn from identifies each mating in the
  # pedigree keys (two matings whose draws happen to coincide, e.g. on a 0 cM
  # map, still get distinct keys); reading it draws nothing.
  draws <- vector("list", n_mating)
  state <- vector("list", n_mating)
  after <- vector("list", n_mating)
  for (k in seq_len(n_mating)) {
    state[k] <- list(.Random.seed_safe())
    draws[[k]] <- .draw_meiosis(lay$by_chr, n[[k]] * events_per[[k]],
                                interference)
    after[k] <- list(.Random.seed_safe())
  }
  pooled <- if (n_mating == 1L) draws[[1L]] else list(
    counts    = unlist(lapply(draws, `[[`, "counts"), use.names = FALSE),
    flips     = unlist(lapply(draws, `[[`, "flips"), use.names = FALSE),
    chiasmata = unlist(lapply(draws, `[[`, "chiasmata"), use.names = FALSE)
  )
  .check_meiosis_call(lay$loci_per_chr, lay$positions, list(), pooled,
                      sum(n * events_per), 1L)

  res <- mate_many_core(
    loci_per_chr = lay$loci_per_chr,
    positions    = lay$positions,
    order        = lay$ord,
    strands      = strands,
    n_strands    = ncol(strands),
    mating       = as.integer(t(mating)),
    design       = design,
    n_prog       = as.integer(n),
    chiasmata    = pooled$chiasmata,
    counts       = pooled$counts,
    flips        = pooled$flips
  )
  total <- sum(n)
  cis <- res$cis
  trans <- res$trans
  dim(cis) <- dim(trans) <- c(nrow(strands), total)
  list(cis = cis, trans = trans, rng_state = state, rng_after = after,
       draws = draws)
}

#' Run one mating design through the Rust core
#' @keywords internal
#' @noRd
.mate <- function(p1, p2, n, design, seed, origin, prefix, interference = NULL) {
  .cite_isqg()
  n <- .validate_count(n, "n", minimum = 1L)
  seed <- .validate_seed(seed)
  interference <- .check_interference(
    interference, c(cross = "cross", selfcross = "selfcross",
                    dh = "double_haploid")[[design]])

  # a self or doubled haploid has one parent (p2 is p1): nothing to compare
  if (design == "cross") {
    if (!.same_map(p1$map, p2$map)) {
      stop("The two parents carry different marker maps; they must come from ",
           "the same Population.", call. = FALSE)
    }
    .check_orientation(p1$map, p2$map)
  }

  map <- p1$map
  if (!is.null(seed)) {
    # the ambient RNG is restored on exit: a seeded call must not disturb the
    # caller's stream
    old_seed <- .Random.seed_safe()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }

  # Ask for haplotypes, not genotypes. A -1/0/1 genotype cannot express the
  # phase of a heterozygote, so deriving progeny strands from one would give an
  # F1 - heterozygous at every locus - a fictitious all-allele-1 / all-allele-2
  # pair, and every later generation would recombine haplotypes that never
  # existed.
  if (design == "cross") {
    strands <- cbind(p1_cis = p1$cis[, 1], p1_trans = p1$trans[, 1],
                     p2_cis = p2$cis[, 1], p2_trans = p2$trans[, 1])
    idx <- matrix(1:4, nrow = 1L)
  } else {
    # a self or doubled haploid has one parent: p2 is p1
    strands <- cbind(p1_cis = p1$cis[, 1], p1_trans = p1$trans[, 1])
    idx <- matrix(c(1L, 2L, 1L, 2L), nrow = 1L)
  }
  res <- .mate_many(map, strands, idx, n, design, interference)

  ids <- paste0(prefix, seq_len(n))
  cis <- res$cis
  trans <- res$trans
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)

  # Pedigree bookkeeping only; it draws nothing, so the RNG stream is unchanged.
  mp <- .mating_pedigree(p1, p2, design,
                         list(rng_state = res$rng_state[[1L]],
                              rng_after = res$rng_after[[1L]],
                              draws = res$draws[[1L]]), ids)
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
#' @section Crossover interference:
#' By default (`interference = NULL`) crossovers are Poisson, i.e. there is no
#' interference, and the random draws are exactly isqg's. Give
#' `interference = list(nu = , p = )` for a two-pathway gamma model of
#' interference (`nu` plays the role of AlphaSimR's `v`, whose default `2.6`
#' approximates Kosambi's map function):
#'
#' * The four-strand bivalent carries chiasmata at an average of 2 per Morgan.
#'   A gamete takes part in each chiasma independently with probability 1/2
#'   (no chromatid interference), so it carries one crossover per Morgan on
#'   average: **the genetic map (`cm`) keeps its meaning for every `nu` and
#'   `p`**, as under Haldane.
#' * A share `p` of the chiasmata (intensity \eqn{2p} per Morgan) is
#'   non-interfering: a Poisson process.
#' * The rest (intensity \eqn{2(1 - p)}) is interfering: a stationary renewal
#'   process whose gaps are gamma distributed with shape `nu` and rate
#'   \eqn{2 \nu (1 - p)} per Morgan (mean gap \eqn{1 / (2 (1 - p))} Morgans).
#'   `nu >= 1` is the interference strength; `nu = 1` is exponential gaps, i.e.
#'   no interference. The first chiasma is drawn from the stationary
#'   (equilibrium) distribution, so the process does not "start" at the
#'   chromosome end.
#' * The two pathways are independent; the gamete's crossovers are the union.
#'
#' The recombination fraction between two loci \eqn{d} Morgans apart is then
#' \eqn{r(d) = [1 - P_0(d)] / 2}, with \eqn{P_0(d)} the probability of no
#' chiasma in the interval: \eqn{P_0(d) = e^{-2pd} \{1 - F_{\nu+1}(y) -
#' (y / \nu) [1 - F_\nu(y)]\}}, \eqn{y = 2 \nu (1 - p) d}, \eqn{F_a} the
#' Gamma(`a`, 1) distribution function. With `nu = 1`, or `p = 1`, the model is
#' Poisson in distribution (Haldane), although it then uses its own random
#' stream, not isqg's. `p` is in `[0, 1]`; `nu` must be at least 1 (negative
#' interference is not modelled). The gamma model of interference is that of
#' McPeek and Speed (1995, *Genetics*) and the two-pathway extension that of
#' Housworth and Stahl (2003, *American Journal of Human Genetics*) (author,
#' year and journal only; the formulas above were derived and are checked
#' numerically in this package's tests, not copied from either paper or from
#' AlphaSimR). The draws are made in R; the Rust core only applies the drawn
#' crossovers and is unchanged.
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
#' @param interference `NULL` (default: Poisson crossovers, no interference, the
#'   isqg random stream) or `list(nu = , p = )` for the two-pathway gamma model
#'   of crossover interference (see the section "Crossover interference").
#'   `nu >= 1` is the interference strength (1 = none), `p` in `[0, 1]` (default
#'   0 when omitted) the share of chiasmata that do not interfere. The expected
#'   number of crossovers per Morgan is unchanged.
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
cross <- function(mother, father, n = 1, seed = NULL, interference = NULL) {
  .check_single(mother, "mother")
  .check_single(father, "father")
  .mate(mother, father, n, "cross", seed,
        origin = paste0("cross(", mother$ids, " x ", father$ids, ")"),
        prefix = "prog_", interference = interference)
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
selfcross <- function(parent, n = 1, seed = NULL, interference = NULL) {
  .check_single(parent, "parent")
  .mate(parent, parent, n, "selfcross", seed,
        origin = paste0("selfcross(", parent$ids, ")"),
        prefix = "self_", interference = interference)
}

#' Produce doubled haploids from an individual
#'
#' Produces `n` doubled-haploid progeny: a single recombinant gamete is drawn
#' and then doubled, so every individual is completely homozygous at every
#' marker. This is the one-generation route to a fully inbred line, and the
#' result contains no heterozygotes at all.
#'
#' @section Many doubled-haploid families in one call:
#' `double_haploid()` takes one parent. To make doubled haploids from many
#' parents (or many families from one) call [mate()] with a plan whose rows have
#' `design = "dh"`, e.g.
#' `mate(data.frame(mother = ids, father = ids, n = 100, design = "dh"), pop)`:
#' all rows are executed in one pass through the Rust core, with the same random
#' draws, in the same order, and the same progeny and pedigree as calling
#' `double_haploid()` on each parent in turn under one seed, at a fraction of the
#' cost per call (see [mate()]).
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
double_haploid <- function(parent, n = 1, seed = NULL, interference = NULL) {
  .check_single(parent, "parent")
  .mate(parent, parent, n, "dh", seed,
        origin = paste0("double_haploid(", parent$ids, ")"),
        prefix = "dh_", interference = interference)
}
