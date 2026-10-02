# Crossbreeding (DECISION-031): expected breed composition from the pedigree,
# heterosis on a frozen A + D architecture, and the classical crossbreeding
# systems as a scheme wrapper.

#' Expected breed composition from the pedigree
#'
#' Each individual's expected fraction of genes from each founder pool (breed,
#' heterotic group), traced through its recorded pedigree: a founder is 1 for
#' its own pool (the `pool` given to [as_population()]), and a progeny is the
#' mean of its parents (a self or doubled haploid equals its parent). These are
#' **expected** fractions -- the realized fraction of an individual's genome from
#' each breed varies around them and is not tracked (a `Population` carries no
#' per-locus founder origin).
#'
#' @param pop a `Population`.
#' @return A matrix, individuals by pools, rows summing to 1. Founders imported
#'   without a pool label count toward a column `"<unassigned>"`.
#' @seealso [crossbreed()], [heterosis()], [parentage()]
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:10), allele = "A/G", chr = 1,
#'                 pos = 1:10, cm = seq(0, 90, by = 10), P1 = 1L, P2 = -1L)
#' A <- as_population(g[, 1:6], pool = "A"); B <- as_population(g[, -6], pool = "B")
#' f1 <- cross(A, B, n = 2, seed = 1)
#' bc <- cross(f1[1], A, n = 2, seed = 2)
#' breed_composition(bc)                               # 3/4 A, 1/4 B
breed_composition <- function(pop) {
  .check_population(pop)
  pop <- .ensure_pedigree(pop)
  ped <- pop$pedigree
  ped <- ped[order(ped$generation), , drop = FALSE]
  pools <- ifelse(is.na(ped$pool), "<unassigned>", ped$pool)
  lev <- unique(pools[ped$design == "founder"])
  Fm <- matrix(0, nrow(ped), length(lev), dimnames = list(NULL, lev))
  mi <- match(ped$mother, ped$key)
  fi <- match(ped$father, ped$key)
  for (k in seq_len(nrow(ped))) {
    if (ped$design[k] == "founder" || is.na(mi[k])) {
      Fm[k, pools[k]] <- 1
    } else {
      Fm[k, ] <- (Fm[mi[k], ] + Fm[fi[k], ]) / 2
    }
  }
  out <- Fm[match(pop$keys, ped$key), , drop = FALSE]
  rownames(out) <- pop$ids
  out
}

#' Heterosis on a frozen additive + dominance architecture
#'
#' Realized heterosis of a crossbred population: its mean genotypic value minus
#' the mean of the pure breeds weighted by its expected breed composition
#' ([breed_composition()]) -- the mid-parent for an F1. Also returns, for every
#' pair of breeds, the expected F1 heterosis: the expected mean of random
#' matings between the two breed populations minus their mid-parent, computed
#' from the breeds' own genotypes (this package's derivation; it uses only
#' the breeds' gamete frequencies, so linkage does not enter). Per locus, with
#' gamete frequencies \eqn{p_A}, \eqn{p_B} of the `+a` allele, the F1's allele
#' frequency is the breeds' average, so the additive part cancels against the
#' mid-parent and the heterosis is \eqn{d [h_{AB} - (h_A + h_B) / 2]}, where
#' \eqn{h_{AB} = p_A (1 - p_B) + (1 - p_A) p_B} is the F1's heterozygote
#' frequency and \eqn{h_A}, \eqn{h_B} the breeds' own. When each breed is in
#' Hardy-Weinberg proportions (\eqn{h = 2p(1 - p)}) this is
#' \eqn{\sum_j d_j (p_{Aj} - p_{Bj})^2} -- Falconer & Mackay's (1996)
#' \eqn{H_{F1} = \sum d y^2}, with \eqn{y} the difference in allele frequency
#' between the two populations -- zero for a purely additive architecture; for
#' inbred-line breeds (\eqn{h = 0}) it is
#' \eqn{\sum_j d_j (p_{Aj} + p_{Bj} - 2 p_{Aj} p_{Bj})}. The classical fractions of
#' the F1 heterosis retained in later crosses (one half in an F2 or a backcross,
#' two thirds in a two-breed rotation at equilibrium) follow from the probability
#' that an individual's two alleles come from different breeds, a dominance
#' (breed-origin) model, and are exact **only under conditions on the breeds'
#' Hardy-Weinberg deviations** \eqn{h - 2p(1 - p)} at the locus (not, in
#' general, on each breed separately). A single fixed inbred line has no deviation
#' (\eqn{p = 0} or \eqn{1}, so \eqn{h = 0 = 2p(1 - p)}), but a breed made of
#' several fully inbred lines that differ at a locus does (\eqn{h = 0 <
#' 2p(1 - p)}). Exactly, at one locus with non-zero F1 heterosis the retained
#' fraction in a backcross to breed A is one half if and only if
#' \eqn{h_A = 2 p_A (1 - p_A)}, i.e. only the recurrent breed A need be in
#' Hardy-Weinberg proportions (whatever breed B is), and in the F2 if and only if
#' \eqn{[h_A - 2 p_A (1 - p_A)] + [h_B - 2 p_B (1 - p_B)] = 0}, i.e. the two
#' deviations sum to zero (for example \eqn{-0.12} and \eqn{+0.12}, with neither
#' breed in Hardy-Weinberg proportions); for breeds that violate this the F2 and
#' backcross retain a different fraction. (The two-thirds rotation figure is the
#' classical value for breeds in Hardy-Weinberg proportions; its condition for
#' breeds with deviations is not derived here.) For example,
#' an equal mixture of AA and aa lines crossed to an aa line (pure dominance) has
#' F1 heterozygosity 1/2 and backcross heterozygosity 1/2 against baselines of 0,
#' so the retention is 1, not 1/2. `realized` is
#' computed from the actual genotypes, so it does not assume either fraction;
#' epistasis is outside this function.
#'
#' @param pop the crossbred `Population` (pedigree traced to the breed founders).
#' @param breeds a named list of the pure-breed `Population`s, each name the
#'   founder pool its population traces to (`as_population(pool =)`); checked.
#'   Each must hold every founder of its pool that `pop`'s pedigree traces to
#'   (the breed means are taken from these populations only), and share the marker
#'   map of `pop`.
#' @param qtn,a,d frozen loci and effects, as for [genotypic_value()]; `d` may be
#'   a single value.
#' @return A list with `realized` (mean minus composition-weighted breed mean),
#'   `mean`, `breed_means`, `composition` (the population's mean expected breed
#'   fractions) and `expected_f1` (breeds x breeds matrix of the expected F1
#'   heterosis of each pair of distinct breeds; symmetric, `NA` on the diagonal,
#'   since a breed mated with itself is not a cross).
#' @references
#'   Falconer DS, Mackay TFC (1996) \emph{Introduction to Quantitative Genetics},
#'   4th ed. Longman, Harlow -- heterosis of a cross between two populations,
#'   \eqn{H_{F1} = \sum d y^2}.
#' @seealso [breed_composition()], [crossbreed()], [genotypic_value()]
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:20), allele = "A/G",
#'                 chr = rep(1:2, each = 10), pos = rep(1:10, 2),
#'                 cm = rep(seq(0, 90, by = 10), 2))
#' set.seed(1)
#' A <- as_population(cbind(g, matrix(sample(c(-1L, 1L), 200, TRUE), 20,
#'                         dimnames = list(NULL, paste0("A", 1:10)))), pool = "A")
#' B <- as_population(cbind(g, matrix(sample(c(-1L, 1L), 200, TRUE), 20,
#'                         dimnames = list(NULL, paste0("B", 1:10)))), pool = "B")
#' f1 <- crossbreed(list(A = A, B = B), "two_way", n_progeny = 50, seed = 2)
#' heterosis(f1, list(A = A, B = B), qtn = 1:20, a = rep(0.2, 20), d = 0.3)$realized
heterosis <- function(pop, breeds, qtn, a, d = 0) {
  .check_population(pop)
  .check_breeds(breeds, pop)
  nl <- length(.resolve_geno_qtn(pop, qtn, "heterosis")$idx)
  if (length(d) == 1L) d <- rep(d, nl)
  .check_ad(a, d, nl, "heterosis")
  comp <- breed_composition(pop)
  missing <- setdiff(colnames(comp), names(breeds))
  if (length(missing)) {
    stop("heterosis(): the pedigree traces to pool(s) not in `breeds`: ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  mu_b <- vapply(breeds, function(b) mean(genotypic_value(b, qtn, a, d)),
                 numeric(1))
  frac <- colMeans(comp)
  mu <- mean(genotypic_value(pop, qtn, a, d))
  base <- sum(frac * mu_b[names(frac)])
  gam <- lapply(breeds, function(b) {
    rowMeans(dosages(b)[.resolve_geno_qtn(b, qtn, "heterosis")$idx, ,
                        drop = FALSE] + 1) / 2
  })
  nb <- length(breeds)
  # distinct pairs only: a breed mated with itself is not a crossbreed
  H <- matrix(NA_real_, nb, nb, dimnames = list(names(breeds), names(breeds)))
  for (i in seq_len(nb)) for (k in seq_len(nb)) {
    if (i == k) next
    ef1 <- .expected_cross_means(matrix(gam[[i]]), matrix(gam[[k]]), a, d)
    H[i, k] <- as.numeric(ef1) - (mu_b[[i]] + mu_b[[k]]) / 2
  }
  list(realized = mu - base, mean = mu, breed_means = mu_b, composition = frac,
       expected_f1 = H)
}

#' Validate a named list of breed populations
#' @keywords internal
#' @noRd
.check_breeds <- function(breeds, pop = NULL) {
  if (!is.list(breeds) || !length(breeds) || is.null(names(breeds)) ||
      any(!nzchar(names(breeds))) || anyDuplicated(names(breeds)) ||
      !all(vapply(breeds, inherits, logical(1), "Population"))) {
    stop("`breeds` must be a named list of Populations, one per breed.",
         call. = FALSE)
  }
  maps <- lapply(breeds, function(b) b$map)
  if (!all(vapply(maps[-1], .same_map, logical(1), maps[[1]])) ||
      (!is.null(pop) && !.same_map(pop$map, maps[[1]]))) {
    stop("All breeds (and the crossbred population) must share one marker map.",
         call. = FALSE)
  }
  # each name must be the founder pool its population traces to, or
  # compositions and breed means would be attributed to the wrong breed
  for (b in names(breeds)) {
    ped <- .ensure_pedigree(breeds[[b]])$pedigree
    fp <- unique(ped$pool[ped$design == "founder"])
    if (length(fp) != 1L || is.na(fp) || fp != b) {
      stop("`breeds$", b, "` must be a pure breed whose founders all carry the ",
           "pool label \"", b, "\" (as_population(pool = \"", b, "\")); its ",
           "founder pool(s): ", paste(ifelse(is.na(fp), "<none>", fp),
                                      collapse = ", "), ".", call. = FALSE)
    }
  }
  # heterosis() takes each breed mean from `breeds[[b]]` only, so that
  # population must contain the founders of pool b that `pop` descends from;
  # a subset or a different sample would silently change `realized`
  if (!is.null(pop)) {
    ped <- .ensure_pedigree(pop)$pedigree
    for (b in names(breeds)) {
      traced <- ped$key[ped$design == "founder" & !is.na(ped$pool) &
                          ped$pool == b]
      lost <- setdiff(traced, .ensure_pedigree(breeds[[b]])$keys)
      if (length(lost)) {
        stop("`breeds$", b, "` does not contain ", length(lost), " founder(s) of ",
             "pool \"", b, "\" that the crossbred population descends from; ",
             "the breed means must come from the labelled pool the pedigree ",
             "traces to (not a subset or another sample).", call. = FALSE)
      }
    }
  }
  invisible()
}

#' Crossbreeding systems
#'
#' Runs a classical crossbreeding system on pure-breed `Population`s and returns
#' the final crossbred generation, its pedigree traced to the breeds (see
#' [breed_composition()]). Each generation mates `n_progeny` random dam-sire pairs
#' (one progeny each) with [mate()]; sex is not modelled, so "dam" and "sire" are
#' the roles of the two pools in the plan, not sexes.
#'
#' * `"two_way"`: breed 1 dams x breed 2 sires (an F1).
#' * `"backcross"`: F1 dams x breed-1 sires.
#' * `"three_way"`: F1 (breed 1 x breed 2) dams x breed-3 sires.
#' * `"terminal"`: F1 (breed 1 x breed 2) dams x sires of `sire_breed` (a breed
#'   not in the F1); the products are not bred further.
#' * `"rotational"`: F1 dams, then each generation's crossbred dams x sires of
#'   the next breed in rotation (breed 1, 2, ... then back to 1), for
#'   `generations` generations after the F1. The expected composition converges
#'   (for two breeds, toward alternating 2/3 : 1/3).
#'
#' @param breeds a named list of pure-breed `Population`s (founder pools named
#'   as the list).
#' @param system one of the systems above.
#' @param n_progeny progeny per generation.
#' @param generations for `"rotational"`: generations after the F1.
#' @param sire_breed for `"terminal"`: the name of the terminal sire breed
#'   (default the last breed).
#' @param seed optional RNG seed; the caller's RNG state is restored on exit.
#' @param interference `NULL` (default: the option `simplePHENOTYPES.interference`
#'   if set, else Poisson crossovers) or `list(nu = , p = )`,
#'   the crossover interference model of [cross()], applied to every generation.
#' @return The final generation as a `Population`, with attribute `history`: a
#'   data frame with one row per generation (`generation`, `sire_breed` and the
#'   mean expected breed fraction per breed).
#' @seealso [breed_composition()], [heterosis()], [mate()]
#' @export
#' @examples
#' g <- data.frame(snp = paste0("m", 1:10), allele = "A/G", chr = 1,
#'                 pos = 1:10, cm = seq(0, 90, by = 10))
#' mk <- function(pool, v) as_population(cbind(g, matrix(v, 10, 4,
#'   dimnames = list(NULL, paste0(pool, 1:4)))), pool = pool)
#' br <- list(A = mk("A", 1L), B = mk("B", -1L))
#' rot <- crossbreed(br, "rotational", n_progeny = 20, generations = 6, seed = 1)
#' attr(rot, "history")
crossbreed <- function(breeds, system = c("two_way", "backcross", "three_way",
                                          "terminal", "rotational"),
                       n_progeny, generations = 1L, sire_breed = NULL,
                       seed = NULL, interference = NULL) {
  system <- match.arg(system)
  interference <- .check_interference(interference, "crossbreed")
  .check_breeds(breeds)
  n_progeny <- .validate_count(n_progeny, "n_progeny", minimum = 1L)
  nb <- length(breeds); bn <- names(breeds)
  need <- c(two_way = 2L, backcross = 2L, three_way = 3L, terminal = 3L,
            rotational = 2L)[[system]]
  if (nb < need) {
    stop("crossbreed(): system \"", system, "\" needs at least ", need,
         " breeds.", call. = FALSE)
  }
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    # the caller's RNG state is restored on exit
    old_seed <- .Random.seed_safe()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  hist <- list()
  record <- function(gen, pop, sire) {
    comp <- colMeans(breed_composition(pop))
    row <- data.frame(generation = gen, sire_breed = sire,
                      stringsAsFactors = FALSE)
    for (b in bn) row[[b]] <- if (b %in% names(comp)) comp[[b]] else 0
    hist[[length(hist) + 1L]] <<- row
  }
  gen_cross <- function(dams, sire_pop, gen, sire_name) {
    plan <- mating_design(dams, sire_pop, design = "random",
                          n_crosses = n_progeny)
    plan$mother_pool <- "dams"; plan$father_pool <- "sires"
    out <- mate(plan, dams = dams, sires = sire_pop,
                prefix = paste0(system, "_g", gen),
                interference = interference)
    record(gen, out, sire_name)
    out
  }
  cur <- gen_cross(breeds[[1]], breeds[[2]], 1L, bn[2])     # F1
  if (system == "backcross") {
    cur <- gen_cross(cur, breeds[[1]], 2L, bn[1])
  } else if (system == "three_way") {
    cur <- gen_cross(cur, breeds[[3]], 2L, bn[3])
  } else if (system == "terminal") {
    sb <- if (is.null(sire_breed)) bn[nb] else sire_breed
    if (!is.character(sb) || length(sb) != 1L || !sb %in% bn || sb %in% bn[1:2]) {
      stop("crossbreed(): `sire_breed` must name a breed other than the two in ",
           "the F1.", call. = FALSE)
    }
    cur <- gen_cross(cur, breeds[[sb]], 2L, sb)
  } else if (system == "rotational") {
    generations <- .validate_count(generations, "generations", minimum = 1L)
    rot <- bn[seq_len(nb)]
    for (g in seq_len(generations)) {
      sire <- rot[((g - 1L) %% nb) + 1L]
      cur <- gen_cross(cur, breeds[[sire]], g + 1L, sire)
    }
  }
  attr(cur, "history") <- do.call(rbind, hist)
  cur
}
