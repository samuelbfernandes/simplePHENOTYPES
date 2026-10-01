# Mating plans (DECISION-025): mating_design() writes a plan, mate() executes it.

#' Execute a mating plan
#'
#' Runs every row of a mating plan -- `mother`, `father`, number of progeny `n`
#' -- through [cross()] (or [selfcross()] when a row names the same individual
#' twice, or [double_haploid()] when `design = "dh"`), and pools the progeny into
#' one `Population` with its pedigree recorded (see [parentage()]). Parents may
#' come from several `Population`s sharing one marker map, e.g. two breeds or
#' heterotic pools.
#'
#' The names in `...` only tell `mate()` where to look up each parent; they do
#' not make individuals distinct. Individuals are identified by their pedigree:
#' a founder by its [as_population()] `pool` label, id and haplotypes, a
#' progeny by its mating. Import each breed with its own `pool` label, or two
#' founders that share an id and genotype are one individual, and mating them is
#' a self.
#'
#' With a `seed`, the RNG is seeded once and the rows are run in plan order, so a
#' one-row plan draws the same random numbers and produces the same progeny
#' genotypes as the equivalent call with that seed -- [cross()] for a cross row,
#' [selfcross()] for a self, [double_haploid()] for a `"dh"` row; the returned
#' object differs only in its ids, `origin` and `plan` attribute. The caller's
#' RNG state is restored on exit.
#'
#' The populations must share one marker map and the allele coded `+1` (see
#' [cross()] on allele orientation).
#'
#' @param plan a data frame with columns `mother`, `father` (individual ids) and
#'   `n` (progeny per row, a positive whole number), as produced by
#'   [mating_design()]. With several pools, add `mother_pool` and `father_pool`
#'   naming the pool (an argument name in `...`) of each parent. An optional
#'   `design` column (`"cross"`, `"self"`, `"dh"`) overrides the default, which is
#'   a cross, or a self when mother and father are the same individual; an `NA`
#'   entry takes the default for its row.
#' @param ... the parent `Population`(s). A single population may be unnamed;
#'   several must be named, and the names are the pools the plan refers to.
#' @param seed optional RNG seed (see Details).
#' @param prefix progeny id prefix; progeny are named `<prefix>_1`, `<prefix>_2`,
#'   ... in plan order. Default: the pool name(s) involved, joined by `x`
#'   (e.g. `"A"` or `"AxB"`), or `"prog"` for a single unnamed population.
#' @return A `Population` of all progeny, in plan order, with attribute `plan`
#'   (the plan with the progeny ids of each row in a list column `progeny`).
#' @seealso [mating_design()], [cross()], [parentage()], [families()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:4)
#' plan <- mating_design(pop, design = "half_diallel", progeny_per_cross = 2)
#' prog <- mate(plan, pop, seed = 1)
#' table(families(prog, "full_sib"))
mate <- function(plan, ..., seed = NULL, prefix = NULL) {
  pools <- list(...)
  if (!length(pools) || !all(vapply(pools, inherits, logical(1), "Population"))) {
    stop("mate(): give the parent Population(s) in `...`.", call. = FALSE)
  }
  nm <- names(pools)
  has_pool_cols <- c("mother_pool", "father_pool") %in% names(plan)
  if (length(pools) == 1L && (is.null(nm) || !nzchar(nm))) {
    nm <- ".single"
    # one unnamed population: there is no pool for a plan's pool columns to
    # refer to, and they used to be overwritten (and dropped) without a word
    if (any(has_pool_cols)) {
      stop("mate(): the plan names pools (mother_pool / father_pool) but `...` ",
           "holds one unnamed population; name it (e.g. mate(plan, A = pop)) ",
           "or remove the pool columns.", call. = FALSE)
    }
  } else if (length(pools) == 1L) {
    # one named pool: a plan without pool columns refers to it; a plan with
    # both is checked against it; a plan with only one is malformed, not
    # something to overwrite
    if (sum(has_pool_cols) == 1L) {
      stop("mate(): the plan has ",
           c("mother_pool", "father_pool")[has_pool_cols], " but not ",
           c("mother_pool", "father_pool")[!has_pool_cols],
           "; give both pool columns or neither.", call. = FALSE)
    }
    if (!any(has_pool_cols)) {
      plan$mother_pool <- nm
      plan$father_pool <- nm
    }
  } else if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("mate(): with several populations, name each one (e.g. mate(plan, ",
         "A = popA, B = popB)); the names are the pools the plan refers to.",
         call. = FALSE)
  }
  names(pools) <- nm
  plan <- .check_plan(plan, nm)
  seed <- .validate_seed(seed)
  if (is.null(prefix)) {
    used <- unique(c(plan$mother_pool, plan$father_pool))
    prefix <- if (identical(used, ".single")) "prog" else paste(used, collapse = "x")
  } else if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix)) {
    stop("mate(): `prefix` must be a single string.", call. = FALSE)
  }

  pick <- function(pool, id) {
    p <- pools[[pool]]
    j <- match(id, p$ids)
    if (is.na(j)) {
      stop("mate(): individual \"", id, "\" not found in pool ",
           if (identical(pool, ".single")) "" else paste0("\"", pool, "\" "),
           "(plan row ", which(plan$mother == id | plan$father == id)[1], ").",
           call. = FALSE)
    }
    p[j]
  }
  if (!is.null(seed)) {
    # restore the caller's RNG state on exit
    old_seed <- .Random.seed_safe()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  kids <- vector("list", nrow(plan))
  resolved <- character(nrow(plan))
  for (k in seq_len(nrow(plan))) {
    mo <- pick(plan$mother_pool[k], plan$mother[k])
    fa <- pick(plan$father_pool[k], plan$father[k])
    # self vs cross by pedigree identity, not display id: the same individual
    # under two ids (e.g. after c(pop, pop)) is a self
    same <- identical(.ensure_pedigree(mo)$keys, .ensure_pedigree(fa)$keys)
    if (is.null(plan$design) || is.na(plan$design[k])) {
      design_k <- if (same) "self" else "cross"
    } else {
      design_k <- plan$design[k]
      if (design_k %in% c("self", "dh") && !same) {
        stop("mate(): plan row ", k, " is a \"", design_k, "\" but names two ",
             "different individuals.", call. = FALSE)
      }
      if (design_k == "cross" && same) {
        stop("mate(): plan row ", k, " crosses an individual with itself; use ",
             "design = \"self\".", call. = FALSE)
      }
    }
    resolved[k] <- design_k
    kids[[k]] <- switch(design_k,
      cross = cross(mo, fa, n = plan$n[k]),
      self  = selfcross(mo, n = plan$n[k]),
      dh    = double_haploid(mo, n = plan$n[k])
    )
  }
  out <- do.call(c, kids)
  ids <- paste0(prefix, "_", seq_len(n_individuals(out)))
  out <- .relabel(out, ids)
  out$origin <- paste0("mate(", nrow(plan), " row", if (nrow(plan) > 1L) "s", ")")
  plan$design <- resolved
  ends <- cumsum(plan$n)
  plan$progeny <- lapply(seq_len(nrow(plan)), function(k) {
    ids[seq.int(ends[k] - plan$n[k] + 1L, ends[k])]
  })
  if (identical(nm, ".single")) {
    plan$mother_pool <- NULL
    plan$father_pool <- NULL
  }
  attr(out, "plan") <- plan
  out
}

#' Validate and complete a mating plan
#' @keywords internal
#' @noRd
.check_plan <- function(plan, pools) {
  if (!is.data.frame(plan) || !all(c("mother", "father", "n") %in% names(plan))) {
    stop("mate(): `plan` must be a data frame with columns mother, father and n ",
         "(see mating_design()).", call. = FALSE)
  }
  if (!nrow(plan)) {
    stop("mate(): `plan` has no rows.", call. = FALSE)
  }
  plan$mother <- as.character(plan$mother)
  plan$father <- as.character(plan$father)
  if (anyNA(plan$mother) || anyNA(plan$father)) {
    stop("mate(): `plan` has missing mother or father ids.", call. = FALSE)
  }
  if (!is.numeric(plan$n) || any(!is.finite(plan$n)) || any(plan$n < 1) ||
      any(plan$n != floor(plan$n))) {
    stop("mate(): `plan$n` must be positive whole numbers.", call. = FALSE)
  }
  plan$n <- as.integer(plan$n)
  if (identical(pools, ".single")) {
    plan$mother_pool <- ".single"
    plan$father_pool <- ".single"
  } else {
    if (!all(c("mother_pool", "father_pool") %in% names(plan))) {
      stop("mate(): with several populations the plan needs mother_pool and ",
           "father_pool columns naming each parent's pool.", call. = FALSE)
    }
    bad <- setdiff(unique(c(plan$mother_pool, plan$father_pool)), pools)
    if (length(bad)) {
      stop("mate(): plan refers to pool(s) not given in `...`: ",
           paste(bad, collapse = ", "), ".", call. = FALSE)
    }
    plan$mother_pool <- as.character(plan$mother_pool)
    plan$father_pool <- as.character(plan$father_pool)
  }
  if (!is.null(plan$design)) {
    plan$design <- as.character(plan$design)
    if (!all(is.na(plan$design) | plan$design %in% c("cross", "self", "dh"))) {
      stop("mate(): `plan$design` must be \"cross\", \"self\", \"dh\" or NA.",
           call. = FALSE)
    }
  }
  plan
}

#' Build a mating plan
#'
#' Writes the `mother`, `father`, `n` plan that [mate()] executes, for the
#' standard mating designs.
#'
#' * `"random"`: `n_crosses` pairs drawn at random (a mother from `mothers`, a
#'   father from `fathers`; with one parent set, two distinct individuals).
#'   Without selfs, a mother is drawn among those with at least one father other
#'   than herself, then one of those fathers. The draw is uniform over mothers and
#'   then over each mother's admissible fathers, so it is uniform over the
#'   admissible *pairs* only when every mother has the same number of admissible
#'   fathers; pairs drawn with replacement, a parent can recur.
#' * `"factorial"`: every mother with every father (North Carolina Design II).
#' * `"nested"`: each father with its own set of `mothers_per_father` distinct
#'   mothers, no mother shared between fathers (North Carolina Design I).
#' * `"diallel"`: every ordered pair of distinct parents, reciprocals included
#'   (Griffing's method 3; with `allow_self = TRUE`, method 1).
#' * `"half_diallel"`: every unordered pair of distinct parents, once (Griffing's
#'   method 4; with `allow_self = TRUE`, method 2).
#'
#' Designs pair individuals by id. A pair of the same **individual** is a self,
#' excluded unless `allow_self = TRUE` (`"diallel"` and `"half_diallel"` then add
#' each parent's self). Individuals are compared by their pedigree identity, not
#' their id: two populations drawn from one (e.g. `pop[1:3]` and `pop[2:4]`)
#' share individuals, while two pools that merely reuse ids (e.g. breeds whose
#' animals are both named `P1`) do not, provided they were imported with their
#' own [as_population()] `pool` labels (a founder is identified by that label,
#' its id and its haplotypes). When either set is given as character ids,
#' individuals are compared by id. Each individual enters a set once: the same
#' individual under a second id (e.g. after `c(pop, pop)`) is dropped, keeping
#' its first id, so no pair is listed twice.
#'
#' @param mothers,fathers `Population`s (or character vectors of ids). `fathers`
#'   defaults to `mothers`; the diallels use `mothers` only (`fathers` is
#'   ignored).
#' @param design one of the designs above.
#' @param n_crosses number of random pairs (`"random"` only).
#' @param progeny_per_cross progeny per plan row (the plan's `n`).
#' @param mothers_per_father mothers per father (`"nested"` only).
#' @param allow_self allow a parent to be mated with itself.
#' @param seed optional RNG seed for `"random"`; the caller's RNG state is
#'   restored on exit. The other designs are deterministic, draw nothing and
#'   leave the RNG untouched, so `seed` is ignored for them.
#' @return A data frame with columns `mother`, `father`, `n`, ready for [mate()]
#'   (add `mother_pool` / `father_pool` when mating across populations).
#' @references
#' Griffing, B. (1956). Concept of general and specific combining ability in
#' relation to diallel crossing systems. \emph{Australian Journal of Biological
#' Sciences} 9(4), 463--493. \doi{10.1071/BI9560463}
#' @seealso [mate()]
#' @export
#' @examples
#' data("SNP55K_maize282_maf04")
#' pop <- as_population(SNP55K_maize282_maf04, individuals = 1:4)
#' mating_design(pop, design = "half_diallel")
#' mating_design(pop[1:2], pop[3:4], design = "factorial", progeny_per_cross = 5)
mating_design <- function(mothers, fathers = mothers,
                          design = c("random", "factorial", "nested", "diallel",
                                     "half_diallel"),
                          n_crosses = NULL, progeny_per_cross = 1L,
                          mothers_per_father = NULL, allow_self = FALSE,
                          seed = NULL) {
  design <- match.arg(design)
  # the diallels use `mothers` only; a `fathers` argument must not change them
  if (design %in% c("diallel", "half_diallel")) fathers <- mothers
  ids_of <- function(x, arg) {
    if (inherits(x, "Population")) return(x$ids)
    if (is.character(x) && length(x) && !anyNA(x)) return(x)
    stop("mating_design(): `", arg, "` must be a Population or character ids.",
         call. = FALSE)
  }
  mo <- ids_of(mothers, "mothers")
  fa <- ids_of(fathers, "fathers")
  # identity of each individual: its pedigree key when both sets are
  # Populations, else its id (a character id carries no pedigree)
  both_pop <- inherits(mothers, "Population") && inherits(fathers, "Population")
  key_of <- function(x, ids) {
    if (both_pop) .ensure_pedigree(x)$keys else ids
  }
  km <- key_of(mothers, mo)
  kf <- key_of(fathers, fa)
  # each individual enters a set once (the same individual under two ids, e.g.
  # after c(pop, pop), keeps its first id), so no pair is listed twice
  keep <- !duplicated(km); mo <- mo[keep]; km <- km[keep]
  keep <- !duplicated(kf); fa <- fa[keep]; kf <- kf[keep]
  n <- .validate_count(progeny_per_cross, "progeny_per_cross", minimum = 1L)
  .validate_flag(allow_self, "allow_self")
  no_self <- !allow_self
  seed <- .validate_seed(seed)
  # only the random design draws; the others must not touch the RNG stream
  if (!is.null(seed) && design == "random") {
    # the caller's RNG state is restored on exit
    old_seed <- .Random.seed_safe()
    on.exit(.restore_seed(old_seed), add = TRUE)
    set.seed(seed)
  }
  pairs <- switch(design,
    random = {
      k <- .validate_count(n_crosses, "n_crosses", minimum = 1L)
      # a mother among those with an admissible father, then one of her
      # admissible fathers (any father but herself when selfs are excluded)
      adm <- lapply(seq_along(mo), function(j) {
        if (no_self) which(kf != km[[j]]) else seq_along(fa)
      })
      ok_m <- which(lengths(adm) > 0L)
      if (!length(ok_m)) {
        stop("mating_design(): no admissible (non-self) pair of a mother and a ",
             "father.", call. = FALSE)
      }
      t(vapply(seq_len(k), function(i) {
        j <- ok_m[sample.int(length(ok_m), 1L)]
        f <- adm[[j]]
        c(mo[j], fa[f[sample.int(length(f), 1L)]])
      }, character(2)))
    },
    factorial = {
      g <- expand.grid(i = seq_along(mo), j = seq_along(fa))
      if (no_self) g <- g[km[g$i] != kf[g$j], , drop = FALSE]
      cbind(mo[g$i], fa[g$j])
    },
    nested = {
      m <- .validate_count(mothers_per_father, "mothers_per_father", minimum = 1L)
      # every father gets exactly m mothers, no mother is shared and no father
      # is paired with itself: a bipartite matching of father slots to mothers
      # (augmenting paths, so a feasible assignment is always found)
      slot_f <- rep(seq_along(fa), each = m)
      ok <- function(s, i) !no_self || km[[i]] != kf[[slot_f[s]]]
      owner <- rep(NA_integer_, length(mo))          # mother -> slot
      try_slot <- function(s, seen) {
        for (i in seq_along(mo)) {
          if (!ok(s, i) || seen[i]) next
          seen[i] <- TRUE
          if (is.na(owner[i])) {
            owner[i] <<- s
            return(list(done = TRUE, seen = seen))
          }
          r <- try_slot(owner[i], seen)
          seen <- r$seen
          if (r$done) {
            owner[i] <<- s
            return(list(done = TRUE, seen = seen))
          }
        }
        list(done = FALSE, seen = seen)
      }
      for (s in seq_along(slot_f)) {
        if (!try_slot(s, logical(length(mo)))$done) {
          stop("mating_design(): a nested design with ", m, " mother(s) per ",
               "father (", length(fa), " fathers) has no valid assignment of ",
               "distinct mothers.", call. = FALSE)
        }
      }
      used <- which(!is.na(owner))
      used <- used[order(owner[used])]
      cbind(mo[used], fa[slot_f[owner[used]]])
    },
    diallel = {
      g <- expand.grid(i = seq_along(mo), j = seq_along(mo))
      g <- g[order(g$i, g$j), ]
      # without selfs keep pairs of distinct individuals only (the same
      # individual in two positions is a self)
      if (!allow_self) g <- g[km[g$i] != km[g$j], ]
      cbind(mo[g$i], mo[g$j])
    },
    half_diallel = {
      g <- expand.grid(i = seq_along(mo), j = seq_along(mo))
      g <- g[order(g$i, g$j), ]
      g <- g[if (allow_self) g$i <= g$j else g$i < g$j, ]
      # the same individual in two positions is a self
      if (!allow_self) g <- g[km[g$i] != km[g$j], ]
      cbind(mo[g$i], mo[g$j])
    }
  )
  pairs <- matrix(pairs, ncol = 2L)
  if (!nrow(pairs)) {
    stop("mating_design(): the design produced no matings. Parents are ",
         "compared by pedigree identity: founders imported without distinct ",
         "`as_population(pool =)` labels that share an id and genotype are the ",
         "same individual.", call. = FALSE)
  }
  data.frame(mother = pairs[, 1], father = pairs[, 2], n = n,
             stringsAsFactors = FALSE)
}
