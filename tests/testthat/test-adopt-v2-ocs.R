# test-adopt-v2-ocs.R
#
# Remaining proposals of the independent audit group v2-ocs-usefulness-marker
# (OCS, usefulness, marker selection, MABC) that earlier regression files
# (test-audit-ocs.R, test-ocs.R, test-marker-select.R, test-mabc.R, test-fix*.R)
# do not yet cover: optimizer optimality conditions on random PSD matrices, G
# alignment by name, the slot record of sample_parents(), the scoring model of
# cross_usefulness() on duplicate / orthogonal / replicated / multi-trait layers,
# marker_select() options, exact weighted recovery in recurrent_parent_recovery()
# and mabc_select(), and the documented meaning of merit = "bv".

data("SNP55K_maize282_maf04")

.ad9_f2 <- function(n = 60) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:40)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = 2))
}
.ad9_ph <- function(f2, ...) {
  suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3, ...) |>
                     additive(n_qtn = 30))
}
# a tiny hand-made Population: rows = markers on separate chromosomes, columns =
# individuals (dosages -1/0/1)
.ad9_pop <- function(D, chr = seq_len(nrow(D)), cm = rep(0, nrow(D)),
                     snp = paste0("m", seq_len(nrow(D)))) {
  g <- data.frame(snp = snp, allele = "A/G", chr = chr, pos = seq_len(nrow(D)),
                  cm = cm, stringsAsFactors = FALSE)
  as_population(cbind(g, as.data.frame(D)))
}

# ---- OCS-T4: optimizer optimality on random PSD matrices ---------------------

test_that("Frank-Wolfe returns a feasible, KKT-optimal simplex point on random PSD G (T4)", {
  fw <- simplePHENOTYPES:::.frank_wolfe
  old <- .Random.seed_safe()
  on.exit(.restore_seed(old))
  set.seed(20260930)
  obj <- function(cc, g, G, lam) sum(cc * g) - 0.5 * lam * drop(crossprod(cc, G %*% cc))
  for (rep in 1:4) {
    n <- 8
    Z <- matrix(stats::rnorm(5 * n), 5, n)          # rank 5 < n: singular PSD
    G <- crossprod(Z) / 5
    g <- stats::rnorm(n)
    for (lam in c(0, 0.1, 1, 10, 100)) {
      cc <- fw(g, G, lam, 20000L, 1e-10)
      expect_true(isTRUE(attr(cc, "converged")), info = paste(rep, lam))
      v <- as.numeric(cc)
      # feasibility on the simplex
      expect_gte(min(v), 0)
      expect_lt(abs(sum(v) - 1), 1e-12)
      # first-order optimality: the Frank-Wolfe gap max(grad) - grad'c is ~ 0, so
      # no vertex improves the linearization (concave objective => global optimum)
      grad <- g - lam * drop(G %*% v)
      expect_lt(max(grad) - sum(grad * v), 1e-7 * (1 + max(abs(grad))))
      # and no random simplex point or vertex beats it
      f0 <- obj(v, g, G, lam)
      pts <- matrix(stats::rexp(n * 300), n, 300)
      pts <- sweep(pts, 2, colSums(pts), "/")
      fr <- apply(pts, 2, obj, g = g, G = G, lam = lam)
      fv <- vapply(seq_len(n), function(k) obj(replace(numeric(n), k, 1), g, G, lam),
                   numeric(1))
      expect_gte(f0 + 1e-9, max(fr, fv))
    }
  }
})

# ---- OCS-C5: G identity metadata / alignment by name --------------------------

test_that("a G with permuted column names is rejected; a permuted G is aligned to merit by name (C5)", {
  fnd <- .ad9_pop(matrix(c(1, -1, 1, -1, 0, 0), 3, 2, byrow = TRUE,
                         dimnames = list(NULL, c("x", "y"))))
  Z <- matrix(c(1, 0.5, -0.2, 0.3, 1, 0.1, -0.4, 0.2, 1), 3,
              dimnames = list(NULL, NULL))
  G <- crossprod(Z); dimnames(G) <- list(c("a", "b", "c"), c("a", "b", "c"))
  merit <- c(a = 0.2, b = 1, c = 0.5)
  base <- optimum_contribution(fnd, merit = merit, G = G, lambda = 1.5)
  # same problem, G supplied in a different individual order: contributions per
  # individual are unchanged (alignment is by name, not by position)
  perm <- c("c", "a", "b")
  alt <- optimum_contribution(fnd, merit = merit, G = G[perm, perm], lambda = 1.5)
  expect_identical(names(alt$contributions), perm)
  expect_equal(alt$contributions[names(base$contributions)], base$contributions,
               tolerance = 1e-7)
  expect_equal(alt$merit, base$merit, tolerance = 1e-8)
  expect_equal(alt$coancestry, base$coancestry, tolerance = 1e-8)
  # row and column names present but in a different order: an inconsistent matrix
  Gbad <- G; dimnames(Gbad) <- list(c("a", "b", "c"), c("b", "a", "c"))
  expect_error(optimum_contribution(fnd, merit = merit, G = Gbad, lambda = 1),
               "row and column names")
})

# ---- OCS-C7: the slot record and pedigree of sample_parents() -----------------

test_that("sample_parents(): repeated slots keep the source key, unique ids, valid parentage (C7)", {
  f2 <- .ad9_f2(40)
  ph <- .ad9_ph(f2)
  oc <- optimum_contribution(ph, lambda = 5)
  n <- 30
  for (meth in c("allocate", "multinomial")) {
    m <- sample_parents(oc, f2, n = n, seed = 3, method = meth)
    src <- attr(m, "source")
    expect_equal(nrow(src), n)
    expect_identical(src$slot, seq_len(n))
    expect_identical(src$name, m$ids)                    # display ids = Population ids
    expect_false(anyDuplicated(m$ids) > 0)               # unique names downstream
    expect_identical(src$id, f2$ids[src$index])          # the draw is recorded
    expect_identical(m$keys, f2$keys[src$index])         # genotype identity preserved
    expect_true(anyDuplicated(src$id) > 0)               # (n > support: repeats occur)
    # every slot is a genotype copy of its source individual
    expect_equal(unname(dosages(m)), unname(dosages(f2)[, src$index, drop = FALSE]))
    # pedigree stays valid: one row per slot, copies of one parent share ancestry
    pg <- parentage(m)
    expect_equal(nrow(pg), n)
    expect_identical(pg$id, m$ids)
    expect_identical(pg$key, m$keys)
    for (id in unique(src$id[duplicated(src$id)])) {
      rows <- pg[src$id == id, c("mother_key", "father_key")]
      expect_equal(nrow(unique(rows)), 1L)
    }
  }
})

# ---- cross_usefulness scoring model (C9, C10) ---------------------------------

test_that(".additive_model() sums layers that share markers on the realized scale (C9)", {
  f2 <- .ad9_f2(60)
  base <- suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3))
  # (the F2 has all-heterozygous markers among these loci: the passed-QTN warning)
  ph <- suppressWarnings(suppressMessages(
    base |>
      additive(qtn = c(1, 5, 9), prop = 0.2, effect = c(0.3, 0.2, 0.1)) |>
      additive(qtn = c(5, 9, 20), prop = 0.3, effect = c(0.4, -0.1, 0.2))))
  m <- simplePHENOTYPES:::.additive_model(ph, 1L)
  sc <- function(ly) {
    comp <- simplePHENOTYPES:::.component_raw(ly, ph, 1L, 1L)
    sqrt(ly$prop[1]) / stats::sd(comp)
  }
  s1 <- sc(ph$layers[[1]]); s2 <- sc(ph$layers[[2]])
  coef <- c(
    tapply(m$eff, m$snp, sum)[ph$map$snp[c(1, 5, 9, 20)]])
  expect_equal(unname(coef),
               c(0.3 * s1, 0.2 * s1 + 0.4 * s2, 0.1 * s1 - 0.1 * s2, 0.2 * s2),
               tolerance = 1e-10)
  expect_identical(sum(m$snp == ph$map$snp[5]), 2L)      # kept as separate rows
  # scoring the template with the model reproduces the layer sum (centred)
  gv <- simplePHENOTYPES:::.additive_gv(f2, m)
  tot <- 0
  for (ly in ph$layers) {
    comp <- simplePHENOTYPES:::.component_raw(ly, ph, 1L, 1L)
    tot <- tot + comp * sqrt(ly$prop[1]) / stats::sd(comp)
  }
  expect_equal(gv - mean(gv), tot - mean(tot), tolerance = 1e-8)
})

test_that(".additive_model() scores an orthogonal layer on alpha = a + d(1 - 2p) (C9)", {
  f2 <- .ad9_f2(80)
  og <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.6, seed = 3) |>
      additive(orthogonal = TRUE, a = 0.2, d = 0.6, n_qtn = 12))
  ly <- og$layers[[1]]
  m <- simplePHENOTYPES:::.additive_model(og, 1L)
  x <- dosages(f2)[m$snp, , drop = FALSE] + 1            # gene content 0/1/2
  p <- rowMeans(x) / 2
  comp <- simplePHENOTYPES:::.component_raw(ly, og, 1L, 1L)   # includes d at hets
  a <- ly$effect[[1]]; dd <- ly$d_effect[[1]]
  alpha <- a + dd * (1 - 2 * p)
  expect_equal(m$eff, unname(alpha * sqrt(ly$prop[1]) / stats::sd(comp)),
               tolerance = 1e-10)
  # not the bare a
  expect_gt(max(abs(m$eff - a * sqrt(ly$prop[1]) / stats::sd(comp))), 1e-3)
})

test_that(".additive_model() uses replication 1 of a vary_qtn layer and the requested trait (C9)", {
  f2 <- .ad9_f2(60)
  vr <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, n_reps = 3, vary_qtn = TRUE, seed = 3) |>
      additive(n_qtn = 10))
  ly <- vr$layers[[1]]
  m <- simplePHENOTYPES:::.additive_model(vr, 1L)
  expect_identical(m$snp, vr$map$snp[ly$qtn_reps[[1]][[1]]])
  expect_false(identical(m$snp, vr$map$snp[ly$qtn_reps[[2]][[1]]]))   # reps differ
  # two traits: each trait is scored on its own loci
  two <- suppressMessages(
    simulate_phenotype(f2, n_traits = 2, h2 = c(0.5, 0.5), seed = 3) |>
      additive(n_qtn = 10))
  m1 <- simplePHENOTYPES:::.additive_model(two, 1L)
  m2 <- simplePHENOTYPES:::.additive_model(two, 2L)
  expect_identical(m1$snp, two$map$snp[two$layers[[1]]$qtn[[1]]])
  expect_identical(m2$snp, two$map$snp[two$layers[[1]]$qtn[[2]]])
  expect_length(intersect(m1$snp, m2$snp), 0L)
  u1 <- suppressWarnings(cross_usefulness(two, pairs = rbind(c(1, 2)), n_progeny = 40,
                                          trait = 1, seed = 4))
  u2 <- suppressWarnings(cross_usefulness(two, pairs = rbind(c(1, 2)), n_progeny = 40,
                                          trait = 2, seed = 4))
  expect_false(isTRUE(all.equal(u1$mean, u2$mean)))
})

test_that("cross_usefulness() rejects unknown parents and an unavailable trait (C10)", {
  sim <- suppressMessages(
    simulate_phenotype(as_population(SNP55K_maize282_maf04, individuals = 1:6),
                       h2 = 0.5, seed = 1) |> additive(n_qtn = 20))
  expect_error(cross_usefulness(sim, pairs = rbind(c(1, 9)), n_progeny = 10),
               "not in `pop`")
  expect_error(cross_usefulness(sim, pairs = rbind(c("nope", "also_nope")),
                                n_progeny = 10), "not in `pop`")
  expect_error(cross_usefulness(sim, pairs = rbind(c(0, 2)), n_progeny = 10),
               "not in `pop`")
  # a one-trait simulation has no trait 2 (or 0): it must error (the message
  # should name the argument -- see the known-defect test below)
  expect_error(cross_usefulness(sim, trait = 2, n_progeny = 10))
  expect_error(cross_usefulness(sim, trait = 0, n_progeny = 10))
})

test_that("cross_usefulness(trait = ) is validated with an argument-naming message (C10)", {
  sim <- suppressMessages(
    simulate_phenotype(as_population(SNP55K_maize282_maf04, individuals = 1:6),
                       h2 = 0.5, seed = 1) |> additive(n_qtn = 20))
  expect_error(cross_usefulness(sim, trait = 2, n_progeny = 10), "trait")
  expect_error(cross_usefulness(sim, trait = 0, n_progeny = 10), "trait")
  expect_error(cross_usefulness(sim, trait = 1.5, n_progeny = 10), "trait")
})

# ---- OCS-C14 / A02: what merit = "bv" is, and where it is refused --------------

test_that("the average effect of substitution differs from the least-squares slope off HWE (A02)", {
  # DECISION-019: alpha = a + d(q - p) (transmitting average effect). On the
  # non-HWE locus with gene content (0, 0, 1, 2) and a = 0, d = sqrt(2) the
  # regression slope of the genotypic value on gene content is sqrt(2)/11,
  # i.e. 4/11 of alpha = sqrt(2)/4: the two quantities are not the same, so the
  # documentation must not call "bv" a least-squares projection.
  x <- c(0, 0, 1, 2)
  p <- mean(x) / 2
  d <- sqrt(2)
  alpha <- simplePHENOTYPES:::.avg_effect(0, d, p)
  expect_equal(alpha, sqrt(2) / 4, tolerance = 1e-12)
  gval <- c(-0, -0, d, 0)                       # -a, d, +a at gene content 0, 1, 2 (a = 0)
  slope <- unname(stats::coef(stats::lm(gval ~ x))[2])
  expect_equal(slope, sqrt(2) / 11, tolerance = 1e-12)
  expect_equal(slope / alpha, 4 / 11, tolerance = 1e-12)
  # under Hardy-Weinberg (p = 1/2 here, x = 0, 1, 1, 2) the two coincide for a = 0
  expect_equal(simplePHENOTYPES:::.avg_effect(0, d, 0.5), 0)
})

test_that('merit = "bv" is refused for a model with an epistasis layer (A02)', {
  f2 <- .ad9_f2(40)
  ep <- suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
      additive(n_qtn = 10, prop = 0.3) |> epistasis(n_pairs = 3, prop = 0.2))
  expect_error(optimum_contribution(ep, lambda = 1), "epistasis")
  # an explicit numeric merit bypasses the refusal
  oc <- optimum_contribution(ep, merit = stats::setNames(stats::rnorm(40), f2$ids[1:40]),
                             lambda = 1)
  expect_s3_class(oc, "ocs")
})

test_that("the OCS help text calls bv the transmitting average effect, not a projection (C14)", {
  skip_if_no_source("R", "select_ocs.R")
  src <- readLines(testthat::test_path("..", "..", "R", "select_ocs.R"))
  txt <- gsub("\\s+", " ", paste(sub("^#'\\s?", "", src[grepl("^#'", src)]),
                                 collapse = " "))
  expect_match(txt, "transmitting\\* average effect", fixed = FALSE)
  expect_match(txt, "alpha_j = a_j \\+ d_j\\(q_j - p_j\\)")
  expect_match(txt, "not a least-squares \\(Fisher statistical\\) projection")
  expect_false(grepl("the Fisher average-effect projection", txt, fixed = TRUE))
  # both unsupported cases are named
  expect_match(txt, "epistasis layer and for `architecture = \"complex\"`",
               fixed = FALSE)
})

# ---- marker_select options (C11) ----------------------------------------------

.ad9_ms_pop <- function() {
  # 10 candidates; m1..m3 dosages chosen so that carriers/homozygotes vary
  D <- rbind(m1 = c(1, 1, 0, 0, -1, 1, 0, -1, 1, 0),
             m2 = c(1, 0, 0, 1, -1, 1, 1, 0, -1, 0),
             m3 = c(0, 0, 1, 1, 1, -1, 0, -1, 1, 0))
  colnames(D) <- paste0("I", 1:10)
  .ad9_pop(D)
}

test_that("marker_select(): prop keeps round(prop * N) (at least one) feasible candidates (C11)", {
  pop <- .ad9_ms_pop()
  sc <- stats::setNames(as.numeric(1:10), pop$ids)
  n_sel <- function(p) n_individuals(marker_select(pop, 1, min_markers = 0,
                                                   prop = p, rank_on = sc))
  expect_equal(n_sel(0.3), 3L)
  expect_equal(n_sel(0.04), 1L)                       # max(1, round(0.4))
  expect_equal(n_sel(1), 10L)
  expect_error(marker_select(pop, 1, prop = 0), "prop")
  expect_error(marker_select(pop, 1, prop = 1.2), "prop")
  expect_error(marker_select(pop, 1, n = 2, prop = 0.5), "at most one")
  expect_error(marker_select(pop, 1, min_markers = 0, n = 11), "cannot select")
})

test_that("marker_select(): direction, function / named / unnamed rank_on are applied to the right ids (C11)", {
  pop <- .ad9_ms_pop()
  sc <- stats::setNames(c(5, 9, 1, 7, 3, 10, 2, 8, 4, 6), pop$ids)
  hi <- marker_select(pop, 1, min_markers = 0, n = 3, rank_on = sc)
  lo <- marker_select(pop, 1, min_markers = 0, n = 3, rank_on = sc,
                      direction = "low")
  expect_setequal(hi$ids, names(sc)[sc >= 8])
  expect_setequal(lo$ids, names(sc)[sc <= 3])
  # named, shuffled: aligned by id; unnamed: taken in id order
  shuf <- sc[c(4, 1, 10, 2, 7, 3, 9, 5, 8, 6)]
  hi2 <- marker_select(pop, 1, min_markers = 0, n = 3, rank_on = shuf)
  expect_setequal(hi2$ids, hi$ids)
  hi3 <- marker_select(pop, 1, min_markers = 0, n = 3, rank_on = unname(sc))
  expect_setequal(hi3$ids, hi$ids)
  d <- attr(hi2, "marker_select")
  expect_equal(d$score, unname(sc[d$id]))
  # a function of the population: here the row sum of the three marker dosages
  f <- function(p) colSums(dosages(p))
  fs <- marker_select(pop, 1, min_markers = 0, n = 2, rank_on = f)
  expect_equal(attr(fs, "marker_select")$score, unname(colSums(dosages(pop))))
  expect_true(all(attr(fs, "marker_select")$score[attr(fs, "marker_select")$selected] >=
                    sort(colSums(dosages(pop)), decreasing = TRUE)[2]))
  # incomplete / wrong-length / non-finite scores are rejected
  expect_error(marker_select(pop, 1, min_markers = 0, rank_on = sc[1:5]),
               "every candidate")
  expect_error(marker_select(pop, 1, min_markers = 0, rank_on = 1:4), "one finite")
  expect_error(marker_select(pop, 1, min_markers = 0,
                             rank_on = c(1:9, NA)), "one finite")
  expect_error(marker_select(pop, 1, min_markers = 0, rank_on = letters[1:10]),
               "numeric")
})

test_that("marker_select(): the seed fixes the tie-break and the ambient RNG is restored (C11)", {
  pop <- .ad9_ms_pop()                              # default score n_met: many ties
  a <- marker_select(pop, 1:3, requirement = "carrier", min_markers = 1, n = 4,
                     seed = 7)
  b <- marker_select(pop, 1:3, requirement = "carrier", min_markers = 1, n = 4,
                     seed = 7)
  expect_identical(attr(a, "marker_select"), attr(b, "marker_select"))
  ds <- attr(a, "marker_select")
  # chosen = feasible candidates ordered by (n_met desc, tiebreak asc)
  ord <- order(!ds$feasible, -ds$n_met, ds$tiebreak)
  expect_identical(ds$id[ord[1:4]], a$ids)
  expect_identical(sort(ds$tiebreak), 1:10)
  # seeds change the tie-break
  others <- vapply(1:8, function(s)
    paste(marker_select(pop, 1:3, min_markers = 1, n = 4, seed = s)$ids,
          collapse = ","), character(1))
  expect_gt(length(unique(others)), 1L)
  # ambient RNG untouched by a seeded call
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  invisible(marker_select(pop, 1:3, min_markers = 1, n = 4, seed = 7))
  expect_identical(stats::runif(1), before)
  expect_error(marker_select(pop, 1:3, min_markers = 1, seed = 1.9), "seed")
})

# ---- exact weighted recovery (C12, C13) ----------------------------------------

# Recurrent REC = +1 and donor DON = -1 on four markers of one chromosome; four
# candidates with recurrent-allele fractions s = (x + 1) / 2.
.ad9_rec_fixture <- function(cm = c(0, 10, 20, 100)) {
  D <- cbind(REC = rep(1, 4), DON = rep(-1, 4),
             C1 = c(1, 1, 1, 1), C2 = c(1, 0, -1, 0),
             C3 = c(-1, -1, -1, -1), C4 = c(0, 0, 0, 0))
  rownames(D) <- paste0("m", 1:4)
  pop <- .ad9_pop(D, chr = rep(1, 4), cm = cm)
  list(pop = pop, REC = pop[1], DON = pop[2], cand = pop[3:6],
       S = rbind(C1 = c(1, 1, 1, 1), C2 = c(1, .5, 0, .5),
                 C3 = c(0, 0, 0, 0), C4 = c(.5, .5, .5, .5)))
}

test_that("recurrent_parent_recovery() equals the weighted mean sum(w s) / sum(w) (C12)", {
  fx <- .ad9_rec_fixture()
  r <- function(...) recurrent_parent_recovery(fx$cand, fx$REC, fx$DON, ...)
  expect_equal(as.numeric(r()), unname(rowMeans(fx$S)))                      # equal weights
  w <- c(1, 2, 3, 4)
  expect_equal(as.numeric(r(weights = w)), unname(drop(fx$S %*% w) / sum(w)))
  expect_equal(unname(r(weights = w))[2], 0.4)                   # (1 + 1 + 0 + 2) / 10
  # named weights in any order are matched to the markers by name
  wn <- stats::setNames(w, paste0("m", 1:4))
  expect_equal(r(weights = wn[c(3, 1, 4, 2)]), r(weights = w))
  # a marker subset uses only those weights: (2 * .5 + 3 * 0) / 5
  rs <- r(markers = c("m2", "m3"), weights = w)
  expect_equal(unname(rs)[2], 0.2)
  expect_identical(attr(rs, "n_markers"), 2L)
  # interval weights: cells [0,5], [5,15], [15,60], [60,100] -> 5, 10, 45, 40
  expect_equal(simplePHENOTYPES:::.mabc_weights(fx$cand$map, 1:4, "interval"),
               c(5, 10, 45, 40))
  expect_equal(unname(r(weights = "interval"))[2], (5 * 1 + 10 * .5 + 40 * .5) / 100)
  # invalid weights
  expect_error(r(weights = c(1, 2, 3)), "one value per map marker")
  expect_error(r(weights = c(1, -1, 1, 1)), "non-negative")
  expect_error(r(weights = rep(0, 4)), "sum to zero")
})

test_that("mabc_select(): numeric / named / zero marker weights change the background ranking exactly (C13)", {
  # target m1 (heterozygous in both candidates: a donor carrier); background
  # m2..m4. A is recurrent at m2 only, B at m3 and m4 only.
  D <- cbind(REC = rep(1, 4), DON = rep(-1, 4),
             A = c(0, 1, -1, -1), B = c(0, -1, 1, 1))
  rownames(D) <- paste0("m", 1:4)
  pop <- .ad9_pop(D, chr = rep(1, 4), cm = c(0, 10, 20, 30))
  sel <- function(w, ...) {
    out <- mabc_select(pop[3:4], pop[1], pop[2], target_markers = "m1", n = 1,
                       marker_weights = w, seed = 1, ...)
    list(first = out$ids, d = attr(out, "mabc"))
  }
  eq <- sel(NULL)                         # A = 1/3, B = 2/3
  expect_identical(eq$first, "B")
  expect_equal(eq$d$background_recovery, c(1 / 3, 2 / 3))
  w1 <- sel(c(1, 10, 1, 1))               # weights over all map markers; m1 is the target
  expect_identical(w1$first, "A")
  expect_equal(w1$d$background_recovery, c(10 / 12, 2 / 12))
  wn <- sel(c(m4 = 1, m2 = 10, m3 = 1, m1 = 5))   # named, any order, extra name ok
  expect_identical(wn$first, "A")
  expect_equal(wn$d$background_recovery, c(10 / 12, 2 / 12))
  wz <- sel(c(1, 0, 0, 1))                # zero weight on m2, m3: only m4 counts
  expect_identical(wz$first, "B")
  expect_equal(wz$d$background_recovery, c(0, 1))
  expect_error(sel(c(m2 = 1, m3 = 1)), "missing some scored markers")
  expect_error(sel(c(0, 0, 0, 0)), "sum to zero")
  expect_error(sel(c(1, 1, 1)), "one value per map marker")
})

test_that("mabc_select(): exact ties are broken by the seeded tie-break, not by input order (C13)", {
  # four identical candidates: every ranking key ties
  D <- cbind(REC = rep(1, 3), DON = rep(-1, 3),
             T1 = c(0, 1, 0), T2 = c(0, 1, 0), T3 = c(0, 1, 0), T4 = c(0, 1, 0))
  rownames(D) <- paste0("m", 1:3)
  pop <- .ad9_pop(D, chr = rep(1, 3), cm = c(0, 10, 20))
  pick <- function(s) {
    out <- mabc_select(pop[3:6], pop[1], pop[2], target_markers = "m1", n = 1,
                       seed = s)
    d <- attr(out, "mabc")
    expect_identical(out$ids, d$id[which.min(d$tiebreak)])   # lowest tie-break wins
    out$ids
  }
  first <- vapply(1:30, pick, character(1))
  expect_identical(first, vapply(1:30, pick, character(1)))   # reproducible
  expect_gt(length(unique(first)), 1L)                         # not always T1
})
