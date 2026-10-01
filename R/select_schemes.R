# Named breeding-scheme wrappers built on the two primitives: select_ind() and
# the crossing core (cross/selfcross/double_haploid). Each wrapper seeds the RNG
# once and threads the ambient stream through the primitives (which draw on R's
# RNG when seed = NULL), so one `seed` reproduces the whole scheme. The heavy
# meiosis stays in Rust (DECISION-006); these wrappers are thin R orchestration.

#' Pool populations that share a genetic map
#'
#' Column-binds several `Population`s (e.g. the progeny of many crosses) into
#' one. All inputs must share the same marker map: identical `snp` names and
#' chromosomes, and `pos` and `cm` equal to a numerical tolerance (absolute
#' `1e-8` for values below 1, relative `1e-8` above -- so round-off noise is
#' accepted, any real map difference is not). The map of the first input is kept.
#' Individual ids are made unique across the pooled set.
#'
#' @param ... `Population` objects with identical maps.
#' @return a single `Population`.
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:6)
#' a <- selfcross(pop[1], n = 3, seed = 1)
#' b <- selfcross(pop[2], n = 3, seed = 2)
#' c(a, b)
c.Population <- function(...) {
  pops <- list(...)
  pops <- pops[!vapply(pops, is.null, logical(1))]
  if (!length(pops)) stop("c.Population(): nothing to pool.", call. = FALSE)
  if (!all(vapply(pops, inherits, logical(1), "Population"))) {
    stop("c.Population(): all inputs must be `Population` objects.",
         call. = FALSE)
  }
  if (length(pops) == 1L) return(pops[[1]])
  # Compare the WHOLE map, not just SNP names: pooling keeps the first map and
  # discards the rest, so two populations with the same markers but different
  # genetic positions (cm) would silently inherit the wrong recombination
  # distances. The one shared judgement of map identity (.same_map(): snp, chr,
  # pos and cm, element-wise tolerance 1e-8) is used here and by crossing.
  ref <- pops[[1]]$map
  for (p in pops[-1]) {
    if (!.same_map(p$map, ref)) {
      stop("c.Population(): populations have different marker maps (snp, chr, ",
           "pos or cm differ); only populations sharing an identical map can be ",
           "pooled.", call. = FALSE)
    }
    .check_orientation(pops[[1]]$map, p$map)
  }
  cis <- do.call(cbind, lapply(pops, function(p) p$cis))
  trans <- do.call(cbind, lapply(pops, function(p) p$trans))
  ids <- make.unique(unlist(lapply(pops, function(p) p$ids), use.names = FALSE),
                     sep = "_")
  colnames(cis) <- ids
  colnames(trans) <- ids
  # Pedigrees are pooled by key, so renaming colliding display ids above does not
  # break any parent link.
  pops <- lapply(pops, .ensure_pedigree)
  keys <- unlist(lapply(pops, function(p) p$keys), use.names = FALSE)
  ped <- do.call(.pedigree_union, lapply(pops, function(p) p$pedigree))
  .new_population(pops[[1]]$map, cis, trans, ids, "pool", keys = keys,
                  pedigree = .pedigree_relabel(ped, keys, ids))
}

#' Single seed descent
#'
#' Advances every line to near-homozygosity by selfing one seed per line per
#' generation, with no selection during inbreeding (Bernardo 2020). Line count
#' is preserved: `n` lines in, `n` lines out.
#'
#' @param x a `Population`, or a `phenotype_sim` built on one (its genotypes are
#'   used).
#' @param generations number of selfing generations.
#' @param seed optional RNG seed for the whole scheme: `NULL` (default; the
#'   ambient RNG stream is used) or one non-negative whole number. When given, the
#'   caller's RNG state is restored on exit, so a seeded scheme does not disturb
#'   the surrounding random stream.
#' @return a `Population` of inbred lines.
#' @references Bernardo R (2020) \emph{Breeding for Quantitative Traits in
#'   Plants}, 3rd ed. Stemma Press, Woodbury, Minnesota.
#' @seealso [bulk()], [pedigree()], [recurrent_selection()], [select_ind()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:10)
#' f1  <- cross(pop[1], pop[2], n = 20, seed = 1)
#' ril <- single_seed_descent(f1, generations = 5, seed = 2)
#' # Heterozygosity is driven down by repeated selfing:
#' mean(dosages(ril) == 0)
single_seed_descent <- function(x, generations = 5L, seed = NULL) {
  pop <- .as_founder_pop(x)
  generations <- .validate_count(generations, "generations", minimum = 1L)
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  for (g in seq_len(generations)) {
    pop <- .self_each(pop, n_each = 1L, tag = paste0("g", g))
  }
  pop
}

#' Bulk advance
#'
#' Advances a population as an undivided bulk: each generation every plant is
#' selfed and the progeny are pooled, then `n` seeds are drawn at random from the
#' pooled progeny to form the next generation (Bernardo 2020). Every plant is
#' assumed to set ample seed, so the draw is uniform over the current plants: the
#' number of advanced seeds each plant contributes is multinomial(`n`, 1/N)
#' (mean `n / N`, variance about `n / N`), not fixed -- including when `n` is a
#' multiple of the current size N. Line identity is not tracked.
#'
#' @inheritParams single_seed_descent
#' @param n bulk size carried to the next generation (default: keep the current
#'   size); it may be smaller or larger than the current size.
#' @return a `Population` (the final bulk).
#' @references Bernardo R (2020) \emph{Breeding for Quantitative Traits in
#'   Plants}, 3rd ed. Stemma Press, Woodbury, Minnesota.
#' @seealso [single_seed_descent()], [pedigree()], [recurrent_selection()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:10)
#' f1  <- cross(pop[1], pop[2], n = 30, seed = 1)
#' bk  <- bulk(f1, generations = 4, n = 30, seed = 2)
#' n_individuals(bk)
bulk <- function(x, generations = 5L, n = NULL, seed = NULL) {
  pop <- .as_founder_pop(x)
  generations <- .validate_count(generations, "generations", minimum = 1L)
  size <- if (is.null(n)) n_individuals(pop) else
    .validate_count(n, "n", minimum = 1L)
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  for (g in seq_len(generations)) {
    # Bulk advance: draw `size` seeds at random from the pooled selfed progeny of
    # the current plants. Each plant sets ample seed, so every draw picks a parent
    # uniformly: parental contributions are multinomial(size, 1/N) -- random even
    # when size is a multiple of N (the old fixed n/N contribution made that case
    # deterministic, a "multiple-seed descent"). Only the parents that contribute
    # are selfed, once each into their drawn number of seeds.
    counts <- .bulk_counts(n_individuals(pop), size)
    pop <- .self_each(pop, n_each = counts, tag = paste0("bulk_g", g))
    # Line identity is not tracked in a bulk; give anonymous bulk ids rather
    # than carrying the parent lineage embedded in the id string.
    pop <- .relabel(pop, paste0("bulk_g", g, "_", seq_len(n_individuals(pop))))
  }
  pop
}

#' Pedigree selection
#'
#' Selfs and selects each generation: the population is phenotyped, the best
#' fraction is kept, and those are selfed to form the next generation (Bernardo
#' 2020; Falconer & Mackay 1996).
#'
#' @section Scale of the response (fixed-scale caveat):
#' `history$differential` and `history$intensity` are measured on the criterion each
#' generation was ranked on. A [simulate_phenotype()] callback re-standardises its
#' genetic layer to the requested variance budget (`prop`, `h2`) in **every**
#' population it is applied to -- also when the causal loci are fixed with
#' `additive(qtn = )`. The *total* genetic share of the phenotype therefore does
#' not decay across generations. For a purely additive model that share is the
#' additive share, so the accuracy of mass selection, `cor(P, BV) = sqrt(h2)`,
#' stays near its target \eqn{\sqrt{h^2}} (it is not held exactly constant: in
#' an executed example it ranged about 0.695-0.714 against 0.707). With
#' dominance or epistasis layers it declines: the re-standardisation fixes the
#' total genetic variance share, not the part of it that is additive, so
#' `cor(P, BV)` declines as allele frequencies shift and inbreeding rises (in an
#' executed example with additive `prop = 0.3` plus dominance `prop = 0.2` it
#' fell over three pedigree generations, whereas the additive-only model stayed
#' near \eqn{\sqrt{h^2}}). In either case `S` is in each generation's own
#' re-standardised units, so genetic gain does not accumulate on that scale (the
#' loci and their effect ratios are unchanged; only the scale is re-fit). To
#' follow a response on a frozen scale, rank on [phenotype_value()]
#' (fixed effects and a fixed residual variance, DECISION-021) through `on`, e.g.
#' `on = function(s) phenotype_value(s$geno, qtn, effect, var_e = ve)` for a
#' whole-population sim, and measure gain with [additive_value()] on the same
#' `qtn` and `effect` (see the examples). The callback must build its simulation
#' from the population it is given; one that returns a simulation backed by a
#' different population (e.g. a closure over the base population) is an error.
#'
#' @inheritParams single_seed_descent
#' @param phenotype a function mapping a `Population` to a realized
#'   `phenotype_sim` (e.g. `function(p) simulate_phenotype(p, ...) |>
#'   additive(...)`). Called once per generation to score the current
#'   population, which it must build its simulation from. Fix the causal loci with
#'   `additive(qtn = ...)` if the same QTNs should act every generation (see the
#'   fixed-scale caveat above).
#' @param prop,n_select proportion (or count) selected each generation; give one
#'   (supplying both is an error; the default `prop = 0.1` applies only when
#'   `n_select` is `NULL`).
#' @param pop_size number of plants grown each generation (default: the founder
#'   size). Each selected line is selfed into an equal-sized family and the
#'   families are pooled to this size, so the population does not drift down as
#'   it inbreeds.
#' @param on,trait,direction passed to [select_ind()]. `trait` may be a vector
#'   of traits, one per generation (recycled): **tandem selection**, improving one
#'   trait at a time; the history then records the trait selected on. Every
#'   requested trait must exist in the phenotype callback's simulation.
#'   `direction` is one `"high"` or `"low"` (a vector is an error). This scheme
#'   selects on the criterion alone (mass selection), so the family methods and
#'   `n_per_family` of [select_ind()] are not forwarded.
#' @return a `Population`, carrying attribute `history`: a data frame with one
#'   row per generation and columns `generation`, `n_selected`, `differential`
#'   (the selection differential S) and `intensity` (the standardized selection
#'   intensity i), per DECISION-015.
#' @references Bernardo R (2020) \emph{Breeding for Quantitative Traits in
#'   Plants}, 3rd ed. Stemma Press, Woodbury, Minnesota; Falconer DS, Mackay TFC
#'   (1996) \emph{Introduction to Quantitative Genetics}, 4th ed. Longman, Harlow.
#' @seealso [select_ind()], [recurrent_selection()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
#' f1  <- cross(pop[1], pop[2], n = 1, seed = 1)
#' f2  <- selfcross(f1, n = 40, seed = 9)
#' pheno <- function(p) {
#'   simulate_phenotype(p, h2 = 0.5, seed = 7) |> additive(n_qtn = 30)
#' }
#' out <- pedigree(f2, pheno, generations = 3, prop = 0.2, seed = 2)
#' attr(out, "history")
#'
#' # Fixed-scale response: rank on phenotype_value() with frozen effects and a
#' # frozen residual variance, so genetic variance is exhausted as selection
#' # proceeds (a simulate_phenotype() callback would re-standardise it).
#' d   <- dosages(f2)
#' q   <- which(apply(d, 1, stats::sd) > 0)[c(20, 60, 100, 140, 180)]
#' eff <- c(0.5, -0.3, 0.4, 0.2, -0.6)
#' ve  <- attr(phenotype_value(f2, q, eff, h2 = 0.5, seed = 1), "var_e")
#' fixed <- function(s) phenotype_value(s$geno, q, eff, var_e = ve)
#' out2 <- pedigree(f2, pheno, generations = 3, prop = 0.2, on = fixed, seed = 2)
#' attr(out2, "history")
#' var(additive_value(f2, q, eff)); var(additive_value(out2, q, eff))
pedigree <- function(x, phenotype, generations = 5L, prop = 0.1,
                     n_select = NULL, pop_size = NULL,
                     on = "pheno", trait = 1L, direction = "high",
                     seed = NULL) {
  pop <- .as_founder_pop(x)
  .check_phenotyper(phenotype)
  generations <- .validate_count(generations, "generations", minimum = 1L)
  size <- if (is.null(pop_size)) n_individuals(pop) else
    .validate_count(pop_size, "pop_size", minimum = 1L)
  if (!is.null(n_select) && !missing(prop)) {
    stop("pedigree(): give either `prop` or `n_select`, not both (the default ",
         "`prop` applies only when `n_select` is NULL).", call. = FALSE)
  }
  .check_direction(direction)
  .check_tandem(trait)
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  history <- vector("list", generations)
  for (g in seq_len(generations)) {
    sim <- phenotype(pop)
    .check_sim(sim)
    .check_backed(sim, pop)
    .check_tandem_range(trait, sim)
    tr <- trait[((g - 1L) %% length(trait)) + 1L]
    sel <- select_ind(sim, n = n_select,
                      prop = if (is.null(n_select)) prop else NULL,
                      on = on, trait = tr, direction = direction)
    ns <- n_individuals(sel)
    history[[g]] <- data.frame(
      generation = g,
      n_selected = ns,
      differential = attr(sel, "differential"),
      intensity = attr(sel, "intensity")
    )
    if (length(trait) > 1L) history[[g]]$trait <- tr
    # each selected line -> an equal family; pool and trim to the grown size
    prog <- .self_each(sel, n_each = max(1L, ceiling(size / ns)),
                       tag = paste0("ped_g", g))
    if (n_individuals(prog) > size) prog <- prog[sort(sample.int(
      n_individuals(prog), size))]
    pop <- prog
  }
  attr(pop, "history") <- do.call(rbind, history)
  pop
}

#' Recurrent selection
#'
#' Cycles of select-then-intermate for population improvement (Bernardo 2020;
#' Falconer & Mackay 1996): each cycle the population is phenotyped, the best parents
#' are selected, and they are intercrossed to form the next cycle's population.
#' The scale caveat of [pedigree()] applies to `history` here too.
#'
#' @inheritParams pedigree
#' @param cycles number of selection cycles.
#' @param n_parents number of parents selected each cycle. It must be smaller than
#'   the number of plants phenotyped in the cycle (the founder population, then
#'   `n_crosses * progeny_per_cross`); otherwise no selection would take place and
#'   the call is an error.
#' @param n_crosses number of crosses among the selected parents (default:
#'   `n_parents`). Each cross draws an independent random pair of selected
#'   parents (random mating with replacement across crosses), so with few
#'   crosses some selected parents may not be sampled.
#' @param progeny_per_cross progeny produced per cross.
#' @return a `Population` (the final cycle), carrying attribute `history`.
#' @references Bernardo R (2020) \emph{Breeding for Quantitative Traits in
#'   Plants}, 3rd ed. Stemma Press, Woodbury, Minnesota; Falconer DS, Mackay TFC
#'   (1996) \emph{Introduction to Quantitative Genetics}, 4th ed. Longman, Harlow.
#' @seealso [select_ind()], [pedigree()].
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:30)
#' pheno <- function(p) {
#'   simulate_phenotype(p, h2 = 0.5, seed = 7) |> additive(n_qtn = 50)
#' }
#' out <- recurrent_selection(pop, pheno, cycles = 3, n_parents = 10,
#'                            progeny_per_cross = 10, seed = 2)
#' attr(out, "history")
recurrent_selection <- function(x, phenotype, cycles = 3L, n_parents = 10L,
                                n_crosses = NULL, progeny_per_cross = 10L,
                                on = "pheno", trait = 1L, direction = "high",
                                seed = NULL) {
  pop <- .as_founder_pop(x)
  .check_phenotyper(phenotype)
  cycles <- .validate_count(cycles, "cycles", minimum = 1L)
  n_parents <- .validate_count(n_parents, "n_parents", minimum = 2L)
  progeny_per_cross <- .validate_count(progeny_per_cross, "progeny_per_cross",
                                       minimum = 1L)
  n_crosses <- if (is.null(n_crosses)) n_parents else
    .validate_count(n_crosses, "n_crosses", minimum = 1L)
  .check_direction(direction)
  .check_tandem(trait)
  # Selecting n_parents >= N individuals of N is no selection at all (S = i = 0):
  # recurrent *random mating*. Refuse it up front, for the founders and for every
  # later cycle (whose size is n_crosses * progeny_per_cross).
  if (n_parents >= n_individuals(pop)) {
    stop("recurrent_selection(): `n_parents` (", n_parents, ") must be smaller ",
         "than the number of plants phenotyped in a cycle (", n_individuals(pop),
         " in the founder population); otherwise no selection is applied.",
         call. = FALSE)
  }
  if (cycles > 1L && as.numeric(n_crosses) * progeny_per_cross <= n_parents) {
    stop("recurrent_selection(): `n_crosses * progeny_per_cross` (",
         as.numeric(n_crosses) * progeny_per_cross, ") must exceed `n_parents` (",
         n_parents, "); otherwise cycles after the first would apply no selection.",
         call. = FALSE)
  }
  seed <- .validate_seed(seed)
  if (!is.null(seed)) {
    old <- .Random.seed_safe()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  history <- vector("list", cycles)
  for (cy in seq_len(cycles)) {
    sim <- phenotype(pop)
    .check_sim(sim)
    .check_backed(sim, pop)
    .check_tandem_range(trait, sim)
    if (n_parents >= sim$n_ind) {
      stop("recurrent_selection(): cycle ", cy, " phenotyped ", sim$n_ind,
           " plant(s), not more than `n_parents` (", n_parents, "); no selection ",
           "would be applied.", call. = FALSE)
    }
    keep <- n_parents
    tr <- trait[((cy - 1L) %% length(trait)) + 1L]
    parents <- select_ind(sim, n = keep, on = on, trait = tr,
                          direction = direction)
    history[[cy]] <- data.frame(
      cycle = cy,
      n_parents = n_individuals(parents),
      differential = attr(parents, "differential"),
      intensity = attr(parents, "intensity")
    )
    if (length(trait) > 1L) history[[cy]]$trait <- tr
    pop <- .intermate(parents, n_crosses, progeny_per_cross,
                      tag = paste0("cyc", cy))
  }
  attr(pop, "history") <- do.call(rbind, history)
  pop
}

# ---- internal helpers -------------------------------------------------------

#' Validate a (possibly tandem) trait schedule
#' @keywords internal
#' @noRd
.check_tandem <- function(trait) {
  if (!is.numeric(trait) || !length(trait) || any(!is.finite(trait)) ||
      any(trait != floor(trait)) || any(trait < 1)) {
    stop("`trait` must be a trait index, or a vector of them (one per ",
         "generation, recycled) for tandem selection.", call. = FALSE)
  }
  invisible()
}

#' Validate that a tandem trait schedule exists in the callback's simulation
#' @keywords internal
#' @noRd
.check_tandem_range <- function(trait, sim) {
  if (any(trait > sim$n_traits)) {
    stop("`trait` requests trait ", max(trait), " but the `phenotype` function ",
         "returned ", sim$n_traits, " trait(s); every scheduled trait must exist ",
         "in the simulated phenotype (1..", sim$n_traits, ").", call. = FALSE)
  }
  invisible()
}

#' Require the callback's simulation to be backed by the population it was given
#'
#' A callback that closes over another population (typically the base population,
#' a plausible typo) would make every generation re-select from that population
#' and the scheme would never advance, silently. Every individual the simulation
#' holds must be one of the current population's individuals (by pedigree key).
#' @keywords internal
#' @noRd
.check_backed <- function(sim, pop) {
  bad <- function() {
    stop("The `phenotype` function returned a phenotype_sim that is not backed ",
         "by the population it was given, so the scheme would not advance. Build ",
         "the simulation from the function's argument, e.g. function(p) ",
         "simulate_phenotype(p, ...) |> additive(...); a function that closes ",
         "over another population (such as the base population) re-selects from ",
         "that population every generation.", call. = FALSE)
  }
  if (!inherits(sim$geno, "Population")) bad()
  gk <- .ensure_pedigree(sim$geno)$keys
  sk <- if (is.null(sim$ind_idx)) gk else gk[sim$ind_idx]
  if (!all(sk %in% .ensure_pedigree(pop)$keys)) bad()
  invisible()
}

#' Bulk seed counts: how many advanced seeds each of `n_cur` plants contributes
#' when `size` seeds are drawn uniformly (with replacement) from the pooled
#' progeny -- multinomial(size, 1/n_cur).
#' @keywords internal
#' @noRd
.bulk_counts <- function(n_cur, size) {
  tabulate(sample.int(n_cur, size, replace = TRUE), nbins = n_cur)
}

#' Accept a Population or a Population-backed phenotype_sim; return a Population
#' @keywords internal
#' @noRd
.as_founder_pop <- function(x) {
  if (inherits(x, "Population")) return(x)
  if (inherits(x, "phenotype_sim")) {
    if (!inherits(x$geno, "Population")) {
      stop("This phenotype_sim was not built on a Population, so it cannot be ",
           "advanced. Build the founders with as_population()/cross() first.",
           call. = FALSE)
    }
    # Advance only the individuals the sim actually holds: a subset sim
    # (simulate_phenotype(individuals = )) keeps `ind_idx` of the backing
    # population, so the scheme must not silently reintroduce the excluded ones.
    if (is.null(x$ind_idx)) {
      return(x$geno)
    }
    return(x$geno[x$ind_idx])
  }
  stop("Expected a `Population` or a Population-backed `phenotype_sim`.",
       call. = FALSE)
}

#' Require a phenotyping callback
#' @keywords internal
#' @noRd
.check_phenotyper <- function(phenotype) {
  if (!is.function(phenotype)) {
    stop("`phenotype` must be a function mapping a Population to a realized ",
         "phenotype_sim, e.g. function(p) simulate_phenotype(p, ...) |> ",
         "additive(...).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Self every individual once (or n_each times), pooled into one Population.
#' `n_each` is one count for every individual, or one per individual (0 skips it).
#' Progeny ids embed the parent id and a tag so the pedigree stays legible.
#' @keywords internal
#' @noRd
.self_each <- function(pop, n_each = 1L, tag = "self") {
  n <- n_individuals(pop)
  n_each <- rep_len(as.integer(n_each), n)
  kids <- vector("list", n)
  for (j in seq_len(n)) {
    k <- n_each[j]
    if (k < 1L) next
    prog <- selfcross(pop[j], n = k, seed = NULL)
    prog <- .relabel(prog, paste0(pop$ids[j], "_", tag,
                                  if (k > 1L) paste0("_", seq_len(k))
                                  else ""))
    kids[[j]] <- prog
  }
  # drop skipped (zero-count) individuals so c() dispatches on a Population
  do.call(c, Filter(Negate(is.null), kids))
}

#' Intercross selected parents: n_crosses random pairs, progeny_per_cross each.
#' @keywords internal
#' @noRd
.intermate <- function(parents, n_crosses, progeny_per_cross, tag = "cyc") {
  np <- n_individuals(parents)
  if (np < 2L) {
    stop("Recurrent selection needs at least two parents to intercross; ",
         "only ", np, " selected.", call. = FALSE)
  }
  kids <- vector("list", n_crosses)
  for (k in seq_len(n_crosses)) {
    pair <- sample.int(np, 2L)
    prog <- cross(parents[pair[1]], parents[pair[2]],
                  n = progeny_per_cross, seed = NULL)
    prog <- .relabel(prog, paste0(tag, "_x", k, "_",
                                  seq_len(progeny_per_cross)))
    kids[[k]] <- prog
  }
  do.call(c, kids)
}

#' Rebuild a Population with new individual ids
#' @keywords internal
#' @noRd
.relabel <- function(pop, ids) {
  colnames(pop$cis) <- ids
  colnames(pop$trans) <- ids
  .new_population(pop$map, pop$cis, pop$trans, ids, pop$origin, keys = pop$keys,
                  pedigree = .pedigree_relabel(pop$pedigree, pop$keys, ids))
}
