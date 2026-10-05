# test-feat-crossing.R
#
# SPEC-0020 round 5, crossing group:
#   item 2 -- per-call cost: the batched integer-I/O meiosis path (.mate_many(),
#             mate_many_core()) must equal the sequential path draw for draw;
#   item 3 -- crossover interference: the two-pathway gamma model
#             (`interference = list(nu, p)`); default since DECISION-047 (nu = 2.6).

# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

# Interleaved chromosomes (rows sorted within a chromosome only), character
# chromosome labels, a non-zero origin: the layout the isqg conventions are
# hardest on.
# These expectations are for Poisson crossovers (Haldane / the isqg stream),
# the default before DECISION-047; pin it for this file.
withr::local_options(simplePHENOTYPES.interference = "poisson")

.fc_geno <- function(m = 90L, seed = 11) {
  set.seed(seed)
  chr <- rep(c("1", "2", "10"), each = m / 3)
  cm <- stats::ave(rep(0, m), chr, FUN = function(x) sort(stats::runif(length(x), 5, 120)))
  key <- stats::ave(rep(0, m), chr, FUN = function(x) sort(stats::runif(length(x))))
  ii <- order(key)
  data.frame(snp = paste0("s", seq_len(m)), allele = "A/G", chr = chr[ii],
             pos = seq_len(m), cm = cm[ii],
             A = sample(c(-1L, 1L), m, TRUE), B = sample(c(-1L, 1L), m, TRUE),
             C = sample(c(-1L, 0L, 1L), m, TRUE), D = sample(c(-1L, 1L), m, TRUE),
             stringsAsFactors = FALSE)
}
.fc_pop <- function() as_population(.fc_geno())

# HEAD's single-mating implementation, verbatim in substance: parental strands
# serialised to 0/1 strings, mate_haplotypes_core(), progeny parsed back. It is
# the reference the new integer path must reproduce bit for bit.
.legacy_draw <- function(morgans_by_chr, n_events) {
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
  list(counts = counts, flips = flips, chiasmata = as.numeric(unlist(chiasmata)))
}
.legacy_mate <- function(p1, p2, n, design) {
  map <- p1$map
  ord <- order(simplePHENOTYPES:::.chr_rank(map$chr), map$cm)
  inv <- order(ord)
  chr_s <- as.character(map$chr[ord])
  cm_s <- map$cm[ord]
  grp <- rle(chr_s)$lengths
  by_chr <- split(cm_s / 100, rep.int(seq_along(grp), grp))
  events_per <- if (design == "dh") 1L else 2L
  draws <- .legacy_draw(by_chr, n * events_per)
  bits <- function(v) paste(as.integer(v), collapse = "")
  strands <- simplePHENOTYPES:::mate_haplotypes_core(
    loci_per_chr = as.integer(grp), positions = as.numeric(cm_s / 100),
    p1_cis = bits(p1$cis[ord, 1]), p1_trans = bits(p1$trans[ord, 1]),
    p2_cis = bits(p2$cis[ord, 1]), p2_trans = bits(p2$trans[ord, 1]),
    chiasmata = draws$chiasmata, counts = draws$counts, flips = draws$flips,
    design = design, n_prog = n)
  unpack <- function(codes) {
    m <- matrix(as.integer(unlist(strsplit(codes, "", fixed = TRUE))),
                nrow = nrow(map), ncol = n)
    m[inv, , drop = FALSE]
  }
  list(cis = unpack(strands[seq(1, length(strands), by = 2)]),
       trans = unpack(strands[seq(2, length(strands), by = 2)]))
}

# ---------------------------------------------------------------------------
# item 2 -- the default path is bit-identical to HEAD
# ---------------------------------------------------------------------------

test_that(".draw_meiosis(NULL) is the isqg stream, draw for draw (HEAD's loop)", {
  by_chr <- list(c(0.05, 0.4, 1.2), c(0, 0.3), 0.9, c(0.1, 2.5, 3.1, 4.4))
  for (s in 1:5) {
    set.seed(s); a <- simplePHENOTYPES:::.draw_meiosis(by_chr, 300)
    state_a <- .Random.seed
    set.seed(s); b <- .legacy_draw(by_chr, 300)
    expect_identical(a, b)
    expect_identical(state_a, .Random.seed)   # same number of draws consumed
  }
})

test_that("cross / selfcross / double_haploid equal HEAD's string-based path", {
  pop <- .fc_pop()
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 5, seed = 3))
  cases <- list(
    list(f = function(s) cross(pop[3], f1[1], n = 6, seed = s),
         p1 = pop[3], p2 = f1[1], design = "cross"),
    list(f = function(s) selfcross(f1[2], n = 7, seed = s),
         p1 = f1[2], p2 = f1[2], design = "selfcross"),
    list(f = function(s) double_haploid(f1[3], n = 9, seed = s),
         p1 = f1[3], p2 = f1[3], design = "dh"))
  for (cs in cases) {
    for (s in c(1, 2, 3)) {
      got <- cs$f(s)
      set.seed(s)
      ref <- .legacy_mate(cs$p1, cs$p2, n_individuals(got), cs$design)
      expect_identical(unname(got$cis), ref$cis)
      expect_identical(unname(got$trans), ref$trans)
    }
  }
  # the ambient stream (seed = NULL) is consumed identically
  set.seed(40); a <- double_haploid(f1[1], n = 4); after_a <- .Random.seed
  set.seed(40); ref <- .legacy_mate(f1[1], f1[1], 4L, "dh")
  expect_identical(unname(a$cis), ref$cis)
  expect_identical(after_a, .Random.seed)
})

test_that("interference = NULL changes nothing (value, pedigree keys, RNG)", {
  pop <- .fc_pop()
  a <- cross(pop[1], pop[2], n = 5, seed = 6)
  b <- cross(pop[1], pop[2], n = 5, seed = 6, interference = NULL)
  expect_identical(a, b)
})

# ---------------------------------------------------------------------------
# item 2 -- the batched path is the sequential path
# ---------------------------------------------------------------------------

test_that(".mate_many() of several matings equals the matings one after another", {
  pop <- .fc_pop()
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 3, seed = 3))
  map <- pop$map
  S <- cbind(pop$cis[, 1], pop$trans[, 1], f1$cis[, 1], f1$trans[, 1],
             pop$cis[, 3], pop$trans[, 3])
  mating <- rbind(c(1, 2, 3, 4), c(3, 4, 3, 4), c(5, 6, 5, 6), c(1, 2, 5, 6))
  storage.mode(mating) <- "integer"
  n <- c(3L, 4L, 2L, 5L)
  design <- c("cross", "selfcross", "dh", "cross")
  for (interf in list(NULL, list(nu = 3, p = 0.2))) {
    set.seed(77)
    batch <- simplePHENOTYPES:::.mate_many(map, S, mating, n, design, interf)
    set.seed(77)
    rows <- lapply(seq_along(n), function(k) {
      simplePHENOTYPES:::.mate_many(map, S, mating[k, , drop = FALSE], n[k],
                                    design[k], interf)
    })
    state_batch <- NULL
    expect_identical(batch$cis, do.call(cbind, lapply(rows, `[[`, "cis")))
    expect_identical(batch$trans, do.call(cbind, lapply(rows, `[[`, "trans")))
    # each row saw the same RNG state before / after its draws as its solo run
    for (k in seq_along(n)) {
      expect_identical(batch$rng_state[[k]], rows[[k]]$rng_state[[1]])
      expect_identical(batch$rng_after[[k]], rows[[k]]$rng_after[[1]])
      expect_identical(batch$draws[[k]], rows[[k]]$draws[[1]])
    }
    expect_identical(batch$rng_state[[2]], batch$rng_after[[1]])
  }
})

test_that("mate() equals calling cross()/selfcross()/double_haploid() row by row", {
  pop <- .fc_pop()
  both <- c(pop[1:4], suppressMessages(cross(pop[1], pop[2], n = 3, seed = 3)))
  plan <- data.frame(mother = c("A", "C", "prog_1", "prog_2", "D"),
                     father = c("B", "C", "prog_1", "A", "B"),
                     n = c(3, 4, 5, 2, 6),
                     design = c(NA, "dh", "self", NA, "cross"))
  for (interf in list(NULL, list(nu = 2.6))) {
    m <- mate(plan, both, seed = 8, prefix = "k", interference = interf)
    set.seed(8)
    kids <- list(
      cross(both["A"], both["B"], n = 3, interference = interf),
      double_haploid(both["C"], n = 4, interference = interf),
      selfcross(both["prog_1"], n = 5, interference = interf),
      cross(both["prog_2"], both["A"], n = 2, interference = interf),
      cross(both["D"], both["B"], n = 6, interference = interf))
    ref <- do.call(c, kids)
    ref <- simplePHENOTYPES:::.relabel(ref, paste0("k_", seq_len(n_individuals(ref))))
    expect_identical(unname(m$cis), unname(ref$cis))
    expect_identical(unname(m$trans), unname(ref$trans))
    expect_identical(m$keys, ref$keys)
    expect_identical(m$pedigree, ref$pedigree)
    expect_identical(m$ids, ref$ids)
    expect_identical(parentage(m), parentage(ref))
  }
  # one seed, the same ambient stream afterwards is not disturbed
  set.seed(100); before <- .Random.seed
  invisible(mate(plan, both, seed = 8))
  expect_identical(.Random.seed, before)
})

test_that("a many-row doubled-haploid plan is one call and matches the loop", {
  pop <- .fc_pop()
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 6, seed = 5))
  plan <- data.frame(mother = f1$ids, father = f1$ids, n = 4, design = "dh")
  m <- mate(plan, f1, seed = 12)
  set.seed(12)
  loop <- do.call(c, lapply(seq_len(6), function(i) double_haploid(f1[i], n = 4)))
  expect_identical(unname(dosages(m)), unname(dosages(loop)))
  expect_identical(parentage(m)$design, rep("dh", 24))
  expect_identical(parentage(m)$mother, rep(f1$ids, each = 4))
})

test_that("mate() rejects a malformed plan before drawing anything", {
  pop <- .fc_pop()
  set.seed(3); before <- .Random.seed
  expect_error(mate(data.frame(mother = c("A", "B"), father = c("B", "Z"), n = 1),
                    pop), "not found in pool")
  expect_error(mate(data.frame(mother = c("A", "B"), father = c("B", "A"), n = 1,
                               design = c("cross", "dh")), pop),
               "names two different individuals")
  expect_identical(.Random.seed, before)
})

test_that("the Rust batch entry point validates its inputs (errors, not crashes)", {
  mm <- simplePHENOTYPES:::mate_many_core
  ok <- list(loci_per_chr = 3L, positions = c(0, 0.5, 1), order = 1:3,
             strands = c(1L, 0L, 1L, 0L, 1L, 0L), n_strands = 2L,
             mating = c(1L, 2L, 1L, 2L), design = "cross", n_prog = 1L,
             chiasmata = numeric(), counts = c(0L, 0L), flips = c(0L, 1L))
  r <- do.call(mm, ok)
  expect_identical(lengths(r), c(cis = 3L, trans = 3L))
  bad <- function(...) do.call(mm, utils::modifyList(ok, list(...)))
  expect_error(bad(n_prog = 2L), "meiosis events")
  expect_error(bad(design = "bogus"), "design must be")
  expect_error(bad(order = c(1L, 1L, 3L)), "permutation")
  expect_error(bad(order = c(0L, 2L, 3L)), "1-based")
  expect_error(bad(strands = c(1L, 2L, 1L, 0L, 1L, 0L)), "only 0 and 1")
  expect_error(bad(strands = c(1L, NA, 1L, 0L, 1L, 0L)), "only 0 and 1")
  expect_error(bad(n_strands = 3L), "strands has 6 entries")
  expect_error(bad(mating = c(1L, 2L, 1L, 3L)), "strand 3")
  expect_error(bad(mating = c(1L, 2L)), "strand indices")
  expect_error(bad(n_prog = -1L), "non-negative")
  expect_error(bad(chiasmata = 2, counts = c(1L, 0L)), "chiasma")
})

test_that("missing or non-0/1 parental strands are an R error naming the strand", {
  pop <- .fc_pop()
  bad <- pop; bad$cis[3, 1] <- NA_integer_
  expect_error(cross(bad[1], pop[2], n = 1, seed = 1), "p1_cis.*0/1 entry per marker")
  bad <- pop; bad$trans[3, 2] <- 2L
  expect_error(cross(pop[1], bad[2], n = 1, seed = 1), "p2_trans.*0/1 entry per marker")
  expect_error(double_haploid(bad[2], n = 1, seed = 1), "p1_trans.*0/1 entry per marker")
})

test_that("a call at array size is cheap (the per-call string round trip is gone)", {
  skip_on_cran()
  m <- 14000L
  set.seed(1)
  cm <- unlist(lapply(1:20, function(i) sort(stats::runif(m / 20, 0, 150))))
  g <- data.frame(snp = paste0("s", seq_len(m)), allele = "A/G",
                  chr = rep(1:20, each = m / 20), pos = seq_len(m), cm = cm,
                  P1 = sample(c(-1L, 1L), m, TRUE), P2 = sample(c(-1L, 1L), m, TRUE))
  pop <- as_population(g)
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 1, seed = 2))
  invisible(double_haploid(f1, n = 10, seed = 1))
  t <- system.time(for (i in 1:10) double_haploid(f1, n = 100, seed = i))[["elapsed"]]
  # HEAD: ~0.14 s per call on this size; now ~0.02 s. Loose bound against noise.
  expect_lt(t / 10, 0.1)
})

# ---------------------------------------------------------------------------
# item 3 -- crossover interference
# ---------------------------------------------------------------------------

# Analytic P0(d): probability of no chiasma in an interval of d Morgans, from
# the model's equilibrium forward-recurrence distribution; r(d) = (1 - P0)/2.
.fc_P0 <- function(d, nu, p) {
  y <- 2 * nu * (1 - p) * d
  fe <- if (p >= 1) 0 else stats::pgamma(y, nu + 1) + (y / nu) * (1 - stats::pgamma(y, nu))
  exp(-2 * p * d) * (1 - fe)
}
.fc_r <- function(d, nu, p) (1 - .fc_P0(d, nu, p)) / 2
.fc_draw <- function(L, n, nu, p, seed) {
  set.seed(seed)
  simplePHENOTYPES:::.draw_meiosis(list(c(0, L)), n, list(nu = nu, p = p))
}

test_that("interference is validated", {
  chk <- simplePHENOTYPES:::.check_interference
  withr::local_options(simplePHENOTYPES.interference = NULL)
  expect_identical(chk(NULL), list(nu = 2.6, p = 0))
  expect_identical(chk(list(nu = 2.6)), list(nu = 2.6, p = 0))
  expect_identical(chk(c(nu = 2, p = 0.25)), list(nu = 2, p = 0.25))
  expect_identical(chk(list(nu = 1, p = 1)), list(nu = 1, p = 1))
  for (bad in list(list(nu = 0.9), list(nu = -1), list(nu = NA_real_), list(nu = Inf),
                   list(nu = c(2, 3)), list(nu = "2"), list(p = 0.2),
                   list(nu = 2, p = -0.1), list(nu = 2, p = 1.01), list(nu = 2, p = NA),
                   list(nu = 2, q = 0.1), list(2, 0.1), 2.6, "gamma", list(),
                   data.frame(nu = 2))) {
    expect_error(chk(bad, "cross"), "interference", info = paste(deparse(bad), collapse = ""))
  }
  pop <- .fc_pop()
  expect_error(cross(pop[1], pop[2], interference = list(nu = 0.5)), "cross\\(\\).*nu")
  expect_error(selfcross(pop[1], interference = list(nu = 2, p = 2)), "selfcross\\(\\).*`interference\\$p`")
  expect_error(double_haploid(pop[1], interference = 3), "double_haploid\\(\\)")
  expect_error(mate(data.frame(mother = "A", father = "B", n = 1), pop,
                    interference = list(nu = 0)), "mate\\(\\)")
  expect_error(crossbreed(list(A = pop, B = pop), "two_way", n_progeny = 2,
                          interference = list(nu = 0)), "crossbreed\\(\\)")
})

test_that("the draws are well-formed and consumed through the unchanged Rust core", {
  for (cfg in list(c(2.6, 0), c(4, 0.3), c(1, 0), c(3, 1))) {
    dr <- .fc_draw(2.2, 500, cfg[1], cfg[2], 5)
    expect_identical(length(dr$counts), 500L)
    expect_identical(length(dr$flips), 500L)
    expect_true(is.integer(dr$counts) && is.integer(dr$flips))
    expect_identical(sum(dr$counts), length(dr$chiasmata))
    expect_true(all(dr$chiasmata >= 0 & dr$chiasmata <= 2.2))
    blocks <- split(dr$chiasmata, rep(seq_along(dr$counts), dr$counts))
    expect_true(all(vapply(blocks, function(b) !is.unsorted(b), NA)))
    expect_true(all(dr$flips %in% 0:1))
    # the kernel accepts them (range, counts, flips)
    expect_length(simplePHENOTYPES:::gamete_masks_core(
      3L, c(0, 1, 2.2), dr$chiasmata, dr$counts, dr$flips), 500L)
  }
  # several chromosomes, slot order is (event, chromosome)
  set.seed(1)
  dr <- simplePHENOTYPES:::.draw_meiosis(list(c(0, 0.5), c(0, 3)), 400,
                                         list(nu = 3, p = 0.1))
  expect_length(dr$counts, 800)
  by_chr <- rep(1:2, 400)
  expect_lt(mean(dr$counts[by_chr == 1]), mean(dr$counts[by_chr == 2]))
  expect_true(all(dr$chiasmata <= 3))
})

test_that("a zero-length chromosome draws no crossover but still the strand flip", {
  for (cfg in list(c(2, 0), c(2, 0.5), c(2, 1))) {
    set.seed(1)
    dr <- simplePHENOTYPES:::.draw_meiosis(list(c(0, 0), c(0, 0.3)), 50,
                                           list(nu = cfg[1], p = cfg[2]))
    expect_length(dr$flips, 100)
    expect_true(all(dr$counts[seq(1, 100, by = 2)] == 0L))
    expect_true(all(dr$flips %in% 0:1) && length(unique(dr$flips)) == 2L)
  }
  # no crossover anywhere: the draw is still well-formed
  set.seed(1)
  dr <- simplePHENOTYPES:::.draw_meiosis(list(c(0, 0)), 4, list(nu = 2, p = 0))
  expect_identical(dr$counts, integer(4)); expect_identical(dr$chiasmata, numeric(0))
  # and the crossing functions run on a map with no recombination
  pop <- as_population(data.frame(snp = paste0("m", 1:4), allele = "A/G", chr = 1,
                                  pos = 1:4, cm = rep(0, 4), P1 = 1L, P2 = -1L))
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 1, seed = 1))
  expect_equal(n_individuals(double_haploid(f1, n = 3, seed = 2,
                                            interference = list(nu = 3))), 3L)
})

test_that("the expected number of crossovers per Morgan is unchanged for any nu, p", {
  L <- 1.5
  n <- 30000
  for (cfg in list(c(1, 0), c(2.6, 0), c(4, 0.3), c(8, 0), c(8, 0.5), c(20, 0.9),
                   c(3, 1))) {
    dr <- .fc_draw(L, n, cfg[1], cfg[2], 100 + round(cfg[1] * 10 + cfg[2] * 100))
    se <- stats::sd(dr$counts) / sqrt(n)
    expect_lt(abs(mean(dr$counts) - L), 4.5 * se,
              label = sprintf("|mean count - L| for nu = %g, p = %g", cfg[1], cfg[2]))
  }
  # interference makes the count under-dispersed relative to Poisson (Var < mean)
  dr <- .fc_draw(L, n, 8, 0, 31)
  expect_lt(stats::var(dr$counts), 0.85 * L)
})

test_that("nu = 1 and p = 1 are Poisson in distribution (counts, gaps, positions)", {
  L <- 1.5
  n <- 30000
  for (cfg in list(c(1, 0), c(1, 0.4), c(5, 1))) {
    dr <- .fc_draw(L, n, cfg[1], cfg[2], 41)
    # counts ~ Poisson(L): chi-square on the pooled table
    obs <- tabulate(pmin(dr$counts, 7L) + 1L, 8L)
    pr <- c(stats::dpois(0:6, L), stats::ppois(6, L, lower.tail = FALSE))
    pv <- suppressWarnings(stats::chisq.test(obs, p = pr)$p.value)
    expect_gt(pv, 1e-3, label = sprintf("count chi-square, nu = %g, p = %g", cfg[1], cfg[2]))
    # positions are uniform on [0, L]
    expect_gt(suppressWarnings(stats::ks.test(dr$chiasmata / L, "punif")$p.value), 1e-3)
  }
  # gaps between successive crossovers on a very long chromosome are Exp(1)
  for (cfg in list(c(1, 0), c(7, 1))) {
    dr <- .fc_draw(2000, 20, cfg[1], cfg[2], 8)
    ev <- rep(seq_along(dr$counts), dr$counts)
    gaps <- unlist(lapply(split(dr$chiasmata, ev), diff), use.names = FALSE)
    expect_gt(length(gaps), 30000)
    expect_gt(suppressWarnings(stats::ks.test(gaps, "pexp", 1)$p.value), 1e-3,
              label = sprintf("gap KS, nu = %g, p = %g", cfg[1], cfg[2]))
  }
  # ... and are not, for nu > 1: the same KS test rejects decisively
  dr <- .fc_draw(2000, 20, 4, 0, 8)
  ev <- rep(seq_along(dr$counts), dr$counts)
  gaps <- unlist(lapply(split(dr$chiasmata, ev), diff), use.names = FALSE)
  expect_lt(suppressWarnings(stats::ks.test(gaps, "pexp", 1)$p.value), 1e-10)
})

test_that("recombination between two markers follows the model's r(d) (nu = 1: Haldane)", {
  n <- 40000
  check <- function(d, nu, p, seed) {
    dr <- .fc_draw(d, n, nu, p, seed)
    masks <- simplePHENOTYPES:::gamete_masks_core(2L, c(0, d), dr$chiasmata,
                                                  dr$counts, dr$flips)
    r_hat <- mean(substr(masks, 1, 1) != substr(masks, 2, 2))
    r <- .fc_r(d, nu, p)
    se <- sqrt(r * (1 - r) / n)
    expect_lt(abs(r_hat - r), 4.5 * se,
              label = sprintf("|r_hat - r| for nu = %g, p = %g, d = %g M", nu, p, d))
  }
  for (d in c(0.1, 0.3, 1)) {
    check(d, 1, 0, 11)                       # Haldane
    check(d, 2.6, 0, 12)
    check(d, 4, 0.3, 13)
    check(d, 8, 0, 14)
    check(d, 3, 1, 15)                       # p = 1: Haldane
  }
  # the analytic function itself: Haldane at nu = 1 (any p) and at p = 1
  d <- c(0.05, 0.2, 0.7, 2)
  expect_equal(.fc_r(d, 1, 0), (1 - exp(-2 * d)) / 2)
  expect_equal(.fc_r(d, 1, 0.6), (1 - exp(-2 * d)) / 2)
  expect_equal(.fc_r(d, 9, 1), (1 - exp(-2 * d)) / 2)
  # nu = 2.6, p = 0 approximates Kosambi's map function (AlphaSimR's documented
  # default for v): within 0.001 in r over 0-1 M
  dd <- seq(0.02, 1, by = 0.02)
  expect_lt(max(abs(.fc_r(dd, 2.6, 0) - 0.5 * tanh(2 * dd))), 0.001)
  # positive interference never lowers r below Haldane's at the same distance
  expect_true(all(.fc_r(dd, 2.6, 0) >= (1 - exp(-2 * dd)) / 2 - 1e-12))
})

test_that("interference is visible: the coefficient of coincidence is below 1", {
  # Expected number of crossover PAIRS with separation in (x1, x2] on a
  # chromosome of L Morgans: integral of (L - x) c(x) dx, where c is the pair
  # correlation (coincidence) of the model: with a = 2(1 - p), b = 2p,
  # lambda = 2 nu (1 - p), renewal density u(x) = sum_k dgamma(x; k nu, lambda),
  #   c(x) = (a u(x) + b^2 + 2 a b) / 4.            (c = 1 for Poisson)
  target <- function(nu, p, L, x1, x2) {
    a <- 2 * (1 - p); b <- 2 * p; lam <- 2 * nu * (1 - p)
    cfun <- function(x) {
      u <- if (p < 1) vapply(x, function(xx) sum(stats::dgamma(xx, (1:200) * nu, lam)), 0) else 0
      (a * u + b^2 + 2 * a * b) / 4
    }
    stats::integrate(function(x) (L - x) * cfun(x), x1, x2)$value
  }
  observed <- function(nu, p, L, n, x1, x2, seed) {
    dr <- .fc_draw(L, n, nu, p, seed)
    ev <- rep(seq_len(n), dr$counts)
    per <- numeric(n)
    cnt <- vapply(split(dr$chiasmata, ev), function(z) {
      if (length(z) < 2) return(0)
      d <- as.vector(stats::dist(z)); sum(d > x1 & d <= x2)
    }, 0)
    per[as.integer(names(cnt))] <- cnt
    c(mean = mean(per), se = stats::sd(per) / sqrt(n))
  }
  L <- 3
  for (cfg in list(c(4, 0), c(2.6, 0.2))) {
    for (rg in list(c(0, 0.1), c(0.1, 0.3), c(0.3, 0.6))) {
      t <- target(cfg[1], cfg[2], L, rg[1], rg[2])
      o <- observed(cfg[1], cfg[2], L, 20000, rg[1], rg[2], 3)
      expect_lt(abs(o[["mean"]] - t), 4.5 * o[["se"]],
                label = sprintf("pairs in (%g, %g], nu = %g, p = %g", rg[1], rg[2], cfg[1], cfg[2]))
    }
  }
  # decisively below the Poisson value at short separation, near 1 far away
  poisson_short <- L * 0.1 - 0.1^2 / 2
  expect_lt(target(4, 0, L, 0, 0.1), 0.1 * poisson_short)
  expect_lt(target(2.6, 0.2, L, 0, 0.1), 0.6 * poisson_short)
  o <- observed(4, 0, L, 20000, 0, 0.1, 9)
  expect_lt(o[["mean"]] + 4.5 * o[["se"]], 0.1 * poisson_short)
  # Poisson itself matches c = 1
  expect_equal(target(1, 0, L, 0, 0.1), poisson_short, tolerance = 1e-6)
})

test_that("the option reaches cross / selfcross / double_haploid / mate / crossbreed", {
  pop <- .fc_pop()
  itf <- list(nu = 3, p = 0.1)
  f1 <- suppressMessages(cross(pop[1], pop[2], n = 2, seed = 1))
  a <- cross(f1[1], f1[2], n = 6, seed = 2, interference = itf)
  b <- cross(f1[1], f1[2], n = 6, seed = 2, interference = itf)
  d <- cross(f1[1], f1[2], n = 6, seed = 2)
  expect_identical(a, b)                       # reproducible under a seed
  expect_false(identical(unname(a$cis), unname(d$cis)))  # a different stream
  # the seeded call restores the ambient RNG state, interference or not
  set.seed(5); before <- .Random.seed
  invisible(cross(f1[1], f1[2], n = 3, seed = 9, interference = itf))
  expect_identical(.Random.seed, before)
  s <- selfcross(f1[1], n = 5, seed = 4, interference = itf)
  h <- double_haploid(f1[2], n = 5, seed = 4, interference = itf)
  expect_equal(n_individuals(s), 5L)
  expect_false(any(dosages(h) == 0))           # doubled haploids stay homozygous
  m <- mate(data.frame(mother = "A", father = "B", n = 3), pop, seed = 2,
            interference = itf)
  expect_identical(unname(m$cis),
                   unname(cross(pop[1], pop[2], n = 3, seed = 2, interference = itf)$cis))
  mk <- function(pool, v) as_population(
    cbind(data.frame(snp = paste0("m", 1:10), allele = "A/G", chr = 1, pos = 1:10,
                     cm = seq(0, 90, by = 10)),
          matrix(v, 10, 4, dimnames = list(NULL, paste0(pool, 1:4)))), pool = pool)
  br <- list(A = mk("A", 1L), B = mk("B", -1L))
  cb <- crossbreed(br, "backcross", n_progeny = 6, seed = 3, interference = itf)
  expect_equal(n_individuals(cb), 6L)
  # the genetic map keeps its meaning: marker pair 20 cM apart recombines ~ r(0.2)
  set.seed(1)
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = c(0, 20),
                  P1 = c(1L, 1L), P2 = c(-1L, -1L))
  f <- cross(as_population(g[, c(1:5, 6)]), as_population(g[, c(1:5, 7)]), n = 1, seed = 1)
  big <- double_haploid(f, n = 20000, seed = 2, interference = list(nu = 2.6, p = 0))
  r_hat <- mean(big$cis[1, ] != big$cis[2, ])
  r <- .fc_r(0.2, 2.6, 0)
  expect_lt(abs(r_hat - r), 4.5 * sqrt(r * (1 - r) / 20000))
})
