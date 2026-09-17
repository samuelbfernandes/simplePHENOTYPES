# Named breeding-scheme wrappers built on the two primitives: select_ind() and
# the crossing core (cross/selfcross/double_haploid). Each wrapper seeds the RNG
# once and threads the ambient stream through the primitives (which draw on R's
# RNG when seed = NULL), so one `seed` reproduces the whole scheme. The heavy
# meiosis stays in Rust (DECISION-006); these wrappers are thin R orchestration.

#' Pool populations that share a genetic map
#'
#' Column-binds several `Population`s (e.g. the progeny of many crosses) into
#' one. All inputs must share an identical marker map. Individual ids are made
#' unique across the pooled set.
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
  # distances. Require snp, chr, pos and cm to match.
  ref <- pops[[1]]$map
  same_map <- function(m) {
    identical(as.character(m$snp), as.character(ref$snp)) &&
      identical(as.character(m$chr), as.character(ref$chr)) &&
      isTRUE(all.equal(as.numeric(m$pos), as.numeric(ref$pos))) &&
      isTRUE(all.equal(as.numeric(m$cm),  as.numeric(ref$cm)))
  }
  for (p in pops[-1]) {
    if (!same_map(p$map)) {
      stop("c.Population(): populations have different marker maps (snp, chr, ",
           "pos or cm differ); only populations sharing an identical map can be ",
           "pooled.", call. = FALSE)
    }
  }
  cis <- do.call(cbind, lapply(pops, function(p) p$cis))
  trans <- do.call(cbind, lapply(pops, function(p) p$trans))
  ids <- make.unique(unlist(lapply(pops, function(p) p$ids), use.names = FALSE),
                     sep = "_")
  colnames(cis) <- ids
  colnames(trans) <- ids
  .new_population(pops[[1]]$map, cis, trans, ids, "pool")
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
#' @param seed optional RNG seed for the whole scheme.
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
  if (!is.null(seed)) set.seed(seed)
  for (g in seq_len(generations)) {
    pop <- .self_each(pop, n_each = 1L, tag = paste0("g", g))
  }
  pop
}

#' Bulk advance
#'
#' Advances a population as an undivided bulk: each generation every plant is
#' selfed and the progeny are pooled, then `n` are carried forward at random
#' (Bernardo 2020). Line identity is not tracked.
#'
#' @inheritParams single_seed_descent
#' @param n bulk size carried to the next generation (default: keep the current
#'   size).
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
  if (!is.null(seed)) set.seed(seed)
  for (g in seq_len(generations)) {
    # Bulk advance: self every plant into a FAMILY (>= 2 seeds each so the pool
    # genuinely exceeds one-per-plant), pool them all, then draw `size` at
    # random. The >= 2 floor is what makes bulk differ from single_seed_descent:
    # with one seed per plant and n == the current size, no sampling would occur
    # and the two schemes would be genetically identical. Random subsampling
    # from the pool is the drift/selection that defines a bulk.
    n_cur <- n_individuals(pop)
    prog <- .self_each(pop, n_each = max(2L, ceiling(size / n_cur)),
                       tag = paste0("bulk_g", g))
    keep <- if (n_individuals(prog) > size) {
      sort(sample.int(n_individuals(prog), size))
    } else seq_len(n_individuals(prog))
    pop <- prog[keep]
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
#' 2020; Falconer \& Mackay 1996). Response accumulates across generations.
#'
#' @inheritParams single_seed_descent
#' @param phenotype a function mapping a `Population` to a realized
#'   `phenotype_sim` (e.g. `function(p) simulate_phenotype(p, ...) |>
#'   additive(...)`). Called once per generation to score the current
#'   population. Fix the causal loci with `additive(qtn = ...)` if the same QTNs
#'   should act every generation.
#' @param prop,n_select proportion (or count) selected each generation; give one.
#' @param pop_size number of plants grown each generation (default: the founder
#'   size). Each selected line is selfed into an equal-sized family and the
#'   families are pooled to this size, so the population does not drift down as
#'   it inbreeds.
#' @param on,trait,direction passed to [select_ind()].
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
pedigree <- function(x, phenotype, generations = 5L, prop = 0.1,
                     n_select = NULL, pop_size = NULL,
                     on = "pheno", trait = 1L, direction = "high",
                     seed = NULL) {
  pop <- .as_founder_pop(x)
  .check_phenotyper(phenotype)
  generations <- .validate_count(generations, "generations", minimum = 1L)
  size <- if (is.null(pop_size)) n_individuals(pop) else
    .validate_count(pop_size, "pop_size", minimum = 1L)
  if (!is.null(seed)) set.seed(seed)
  history <- vector("list", generations)
  for (g in seq_len(generations)) {
    sim <- phenotype(pop)
    .check_sim(sim)
    sel <- select_ind(sim, n = n_select,
                      prop = if (is.null(n_select)) prop else NULL,
                      on = on, trait = trait, direction = direction)
    ns <- n_individuals(sel)
    history[[g]] <- data.frame(
      generation = g,
      n_selected = ns,
      differential = attr(sel, "differential"),
      intensity = attr(sel, "intensity")
    )
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
#' Falconer \& Mackay 1996): each cycle the population is phenotyped, the best parents
#' are selected, and they are intercrossed to form the next cycle's population.
#'
#' @inheritParams pedigree
#' @param cycles number of selection cycles.
#' @param n_parents number of parents selected each cycle.
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
  if (!is.null(seed)) set.seed(seed)
  history <- vector("list", cycles)
  for (cy in seq_len(cycles)) {
    sim <- phenotype(pop)
    .check_sim(sim)
    keep <- min(n_parents, sim$n_ind)
    parents <- select_ind(sim, n = keep, on = on, trait = trait,
                          direction = direction)
    history[[cy]] <- data.frame(
      cycle = cy,
      n_parents = n_individuals(parents),
      differential = attr(parents, "differential"),
      intensity = attr(parents, "intensity")
    )
    pop <- .intermate(parents, n_crosses, progeny_per_cross,
                      tag = paste0("cyc", cy))
  }
  attr(pop, "history") <- do.call(rbind, history)
  pop
}

# ---- internal helpers -------------------------------------------------------

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
#' Progeny ids embed the parent id and a tag so the pedigree stays legible.
#' @keywords internal
#' @noRd
.self_each <- function(pop, n_each = 1L, tag = "self") {
  n <- n_individuals(pop)
  kids <- vector("list", n)
  for (j in seq_len(n)) {
    prog <- selfcross(pop[j], n = n_each, seed = NULL)
    prog <- .relabel(prog, paste0(pop$ids[j], "_", tag,
                                  if (n_each > 1L) paste0("_", seq_len(n_each))
                                  else ""))
    kids[[j]] <- prog
  }
  do.call(c, kids)
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
  .new_population(pop$map, pop$cis, pop$trans, ids, pop$origin)
}
