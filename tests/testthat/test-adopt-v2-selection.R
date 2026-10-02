# test-adopt-v2-selection.R
#
# Remaining test proposals of the independent v2-selection audit (select_ind()
# and the scheme wrappers) that test-audit-selection.R / test-fix*-selection.R /
# test-culling.R did not adopt. Every test states the CORRECT, documented
# behaviour of the current code; the former known defects are fixed (round 9).

data("SNP55K_maize282_maf04")

.as2 <- function(n = 60, seed = 2, parents = 1:2) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
  f1 <- cross(pop[parents[1]], pop[parents[2]], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = seed))
}
.asp <- function(f2, n_traits = 1, ...) {
  suppressMessages(
    simulate_phenotype(f2, h2 = 0.5, n_traits = n_traits, seed = 3, ...) |>
      additive(n_qtn = 20)
  )
}
.asp1 <- function(p) {
  suppressMessages(simulate_phenotype(p, h2 = 0.5, seed = 7) |>
                     additive(n_qtn = 30))
}

# ---- T4 / C3: intensity -> count mapping ------------------------------------

test_that("intensity maps to the closest expected-normal count, endpoints included (T4/C3)", {
  # N = 1000: i(k/1000) = dnorm(qnorm(1 - p)) / p for k = 1, 100, 500, 1000
  ip <- function(k, N) {
    p <- k / N
    if (k == N) 0 else stats::dnorm(stats::qnorm(1 - p)) / p
  }
  expect_equal(ip(1, 1000), 3.3670901, tolerance = 1e-7)      # analytic anchors
  expect_equal(ip(100, 1000), 1.7549833, tolerance = 1e-7)
  expect_equal(ip(500, 1000), 0.79788456, tolerance = 1e-7)   # sqrt(2/pi)
  for (k in c(1L, 100L, 500L, 1000L)) {
    expect_equal(.resolve_keep(NULL, NULL, ip(k, 1000L), 1000L), k)
  }
  # no selection (i = 0) keeps everyone; an unreachably high intensity keeps one
  expect_equal(.resolve_keep(NULL, NULL, 0, 60L), 60L)
  expect_equal(.resolve_keep(NULL, NULL, 10, 60L), 1L)
  # N = 60, i = 1.4 corresponds to p = 0.2 (i = 1.3998): 12 kept
  expect_equal(.resolve_keep(NULL, NULL, 1.4, 60L), 12L)
  ph <- .asp(.as2(60))
  expect_equal(n_individuals(select_ind(ph, intensity = 1.4)), 12L)
  expect_equal(n_individuals(select_ind(ph, intensity = 0)), 60L)
  expect_equal(n_individuals(select_ind(ph, intensity = 10)), 1L)
})

# ---- C4: criterion values are joined to individuals by name -------------------

test_that(".criterion_values aligns every named criterion with sim$ids, whatever the row order (C4)", {
  f2 <- .as2(60)
  ph2 <- .asp(f2, n_traits = 2)
  ids <- ph2$ids
  # visibly different, exactly known values: trait k, individual j -> 1000 k + j
  tag <- ph2
  tag$pheno$value <- 1000 * as.integer(sub("Trait_", "", tag$pheno$trait)) +
    match(tag$pheno$id, ids)
  set.seed(5)
  tag$pheno <- tag$pheno[sample(nrow(tag$pheno)), ]      # scramble the rows
  expect_equal(unname(.criterion_values(tag, "pheno", 1, 1)), 1000 + seq_along(ids))
  expect_equal(unname(.criterion_values(tag, "pheno", 2, 1)), 2000 + seq_along(ids))
  expect_identical(names(.criterion_values(tag, "pheno", 2, 1)), ids)
  # genetic and breeding values come from the sim itself, in id order
  expect_equal(unname(.criterion_values(ph2, "gv", 2, 1)),
               unname(.genetic_matrix(ph2, 1)[, 2]))
  expect_equal(unname(.criterion_values(ph2, "bv", 1, 1)),
               unname(.breeding_value_matrix(ph2, 1)[, 1]))
  # a named numeric criterion is matched by name; a reversed named vector gives
  # the same ranking as the id-ordered one
  sc <- stats::setNames(as.numeric(seq_along(ids)), ids)
  expect_equal(.criterion_values(ph2, rev(sc), 1, 1), sc)
  for (on in c("pheno", "gv", "bv")) {
    expect_error(.criterion_values(ph2, on, c(1, 2), 1), "one trait index")
  }
})

# ---- T6 / C9: .sel_top tie policy ---------------------------------------------

test_that(".sel_top ranks by score and breaks ties by original position (T6/C9)", {
  expect_identical(.sel_top(c(1, 1, 1, 0), 2L), 1:2)
  expect_identical(.sel_top(c(3, 1, 2), 2L), c(1L, 3L))
  expect_identical(.sel_top(c(0, 2, 2, 2, 1), 3L), 2:4)
  expect_identical(.sel_top(rep(7, 5), 5L), 1:5)
})

# ---- T10: within-family allocation with equal remainders ------------------------

test_that("within_family: equal remainders follow the character order of the labels (T10)", {
  # 6 families x 10, keep 3: raw = 0.5 each, no whole slot, three remainders tie
  fam <- rep(LETTERS[1:6], each = 10)
  score <- c(5, 1, 9, 3, 8, 2, 7, 6, 4, 10,  rep(1, 50))
  score[11:20] <- 21:30; score[21:30] <- 31:40; score[31:60] <- rep(1:10, 3)
  idx <- .sel_within_family(score, fam, 3L)
  expect_length(idx, 3L)
  expect_equal(fam[idx], c("A", "B", "C"))                # one each, label order
  expect_equal(score[idx], c(10, 30, 40))                 # the best of each family
  # labels "10" and "2" sort as character: "10" < "2"
  fam2 <- rep(c("2", "10"), each = 5)
  s2 <- c(1, 2, 3, 4, 5, 10, 9, 8, 7, 6)
  idx2 <- .sel_within_family(s2, fam2, 1L)
  expect_equal(fam2[idx2], "10")
  expect_equal(s2[idx2], 10)
  # a whole-slot allocation does not depend on the remainder rule
  idx3 <- .sel_within_family(s2, fam2, 2L)
  expect_equal(sort(fam2[idx3]), c("10", "2"))
  # the same through select_ind(): exactly n kept
  ph <- .asp(.as2(60))
  out <- select_ind(ph, n = 3, method = "within_family", family = rep(1:6, each = 10))
  expect_equal(n_individuals(out), 3L)
  expect_equal(sort(ceiling(match(attr(out, "selected"), ph$ids) / 10)), 1:3)
})

# ---- T9 / C11: .sel_among_family ------------------------------------------------

test_that(".sel_among_family takes whole families by family mean until n is reached (T9/C11)", {
  fam <- rep(c("a", "b"), c(3, 5))
  # family a (3 plants) has the higher mean: 3 < 4, so family b is added as well
  s1 <- c(10, 10, 10, 1, 1, 1, 1, 1)
  expect_equal(.sel_among_family(s1, fam, 4L), 1:8)
  expect_equal(.sel_among_family(s1, fam, 3L), 1:3)
  # family b has the higher mean and already holds 5 >= 4: only b is returned
  s2 <- c(1, 1, 1, 9, 9, 9, 9, 9)
  expect_equal(.sel_among_family(s2, fam, 4L), 4:8)
  # ranking is by the family MEAN, not by the best member
  s3 <- c(100, 0, 0, 3, 3, 3, 3, 3)                       # means: 33.3 vs 3
  expect_equal(.sel_among_family(s3, fam, 1L), 1:3)
  s4 <- c(100, -90, -90, 3, 3, 3, 3, 3)                   # means: -26.7 vs 3
  expect_equal(.sel_among_family(s4, fam, 1L), 4:8)
})

# ---- T11: Lush combined index over a parameter grid ------------------------------

test_that("Lush weights are non-negative, b1 <= h2, and match V^-1 c over a grid (T11)", {
  # recover (b1, b2) from the scores of designed records, independently of the
  # closed form: with a zero family mean score_i = b1 * dev_i; with family mean
  # equal to dev, score_i = (b1 + b2) * dev_i
  lush_b <- function(n, h2, r) {
    fam <- rep(c("a", "b"), each = n)
    y1 <- c(1, -1, rep(0, 2 * n - 2))
    y2 <- c(rep(1, n), rep(-1, n))
    b1 <- unname(.combined_score(y1, fam, h2, r)[1])
    b12 <- unname(.combined_score(y2, fam, h2, r)[1])
    c(b1 = b1, b2 = b12 - b1)
  }
  for (h2 in c(0.01, 0.3, 0.7, 0.99)) {
    for (r in c(0, 0.25, 0.5, 2 / 3, 1)) {
      for (n in c(2L, 5L, 50L)) {
        b <- lush_b(n, h2, r)
        info <- paste0("h2=", h2, " r=", round(r, 3), " n=", n)
        expect_gte(b[["b1"]], -1e-12, label = paste("b1", info))
        expect_gte(b[["b2"]], -1e-12, label = paste("b2", info))
        expect_lte(b[["b1"]], h2 + 1e-12, label = paste("b1<=h2", info))
        # V b = c with V, c from Falconer & Mackay's family-mean variances
        t <- r * h2
        m <- (1 + (n - 1) * t) / n
        V <- matrix(c(1, m, m, m), 2)
        cc <- h2 * c(1, (1 + (n - 1) * r) / n)
        expect_equal(unname(b), as.numeric(solve(V, cc)), tolerance = 1e-8,
                     info = info)
      }
    }
  }
  # r = 0: no family information, the own record carries all the weight
  b0 <- lush_b(5L, 0.4, 0)
  expect_equal(unname(b0), c(0.4, 0), tolerance = 1e-12)
  # the family mean drops out as h2 -> 1 when r < 1
  expect_lt(lush_b(5L, 0.99, 0.5)[["b2"]], 0.1 * lush_b(5L, 0.3, 0.5)[["b2"]])
})

# ---- C13: culling attributes -----------------------------------------------------

test_that("culling reports per-trait differentials and intensities on each natural scale (C13)", {
  f2 <- .as2(60)
  ph2 <- .asp(f2, n_traits = 2)
  y1 <- .criterion_values(ph2, "pheno", 1, 1)
  y2 <- .criterion_values(ph2, "pheno", 2, 1)
  sel <- suppressMessages(select_ind(ph2, method = "culling", culling = c(0.6, 0.6),
                                     direction = c("high", "low")))
  ix <- match(attr(sel, "selected"), ph2$ids)
  expect_gt(length(ix), 0L)
  d <- attr(sel, "differential")
  expect_equal(names(d), c("Trait_1", "Trait_2"))
  expect_equal(unname(d), c(mean(y1[ix]) - mean(y1), mean(y2[ix]) - mean(y2)),
               tolerance = 1e-12)
  expect_gt(d[["Trait_1"]], 0)             # kept the high tail of trait 1
  expect_lt(d[["Trait_2"]], 0)             # kept the low tail of trait 2
  expect_equal(unname(attr(sel, "intensity")),
               c(abs(d[["Trait_1"]]) / stats::sd(y1),
                 abs(d[["Trait_2"]]) / stats::sd(y2)), tolerance = 1e-12)
  expect_true(all(attr(sel, "intensity") >= 0))
  expect_equal(attr(sel, "culling")$kept, rep(36L, 2))
})

# ---- T20: culling at tiny N -------------------------------------------------------

test_that("culling at N = 5 pins the zero-survivor and sequential contracts (T20)", {
  f2 <- .as2(60)
  sub <- suppressMessages(
    simulate_phenotype(f2, individuals = 1:5, h2 = 0.5, n_traits = 2, seed = 3) |>
      additive(n_qtn = 10))
  expect_equal(sub$n_ind, 5L)
  # round(0.5 * 5) = 2 (half-to-even): the top 2 of each trait
  aligned <- cbind(5:1, 5:1)                  # same two best individuals
  s1 <- suppressMessages(select_ind(sub, method = "culling", culling = c(0.5, 0.5),
                                    on = aligned))
  expect_setequal(attr(s1, "selected"), sub$ids[1:2])
  expect_equal(attr(s1, "culling")$kept, c(2L, 2L))
  # opposed criteria: the two top sets are disjoint -> a clear error, no partial result
  opposed <- cbind(5:1, 1:5)
  expect_error(select_ind(sub, method = "culling", culling = c(0.5, 0.5), on = opposed),
               "no individual passes")
  # sequential culling applies trait 2 among the 2 survivors of trait 1 and
  # always keeps at least one: round(0.5 * 2) = 1 survivor, the better on trait 2
  s2 <- suppressMessages(select_ind(sub, method = "culling", culling = c(0.5, 0.5),
                                    on = opposed, sequential = TRUE))
  expect_equal(attr(s2, "selected"), sub$ids[2])
  expect_equal(attr(s2, "culling")$kept, c(2L, 1L))
})

# ---- C12: empty family label -----------------------------------------------------

test_that("an empty-string family label is an error with a suggested fix, for every family method (C12)", {
  ph <- .asp(.as2(60))
  fam_e <- rep(c("", "a"), each = 30)
  for (m in c("within_family", "among_family", "combined")) {
    expect_error(select_ind(ph, n = 6, method = m, family = fam_e, h2 = 0.5),
                 "empty label.*Replace the empty")
  }
  # blank labels are empty too; real labels (including a replacement for "") work
  expect_error(select_ind(ph, n = 6, method = "within_family",
                          family = rep(c("  ", "a"), each = 30)), "empty label")
  fam_ok <- rep(c("unknown", "a"), each = 30)
  for (m in c("within_family", "among_family", "combined")) {
    expect_no_error(select_ind(ph, n = 6, method = m, family = fam_ok, h2 = 0.5))
  }
})

# ---- C6: .index_weights ----------------------------------------------------------

test_that(".index_weights is the exact solve for a full-rank P and finite for ill-conditioned P (C6)", {
  set.seed(4)
  A <- crossprod(matrix(stats::rnorm(9), 3)) + diag(3)
  r <- c(1, -2, 0.5)
  expect_equal(as.numeric(.index_weights(A, r)), as.numeric(solve(A, r)),
               tolerance = 1e-10)
  # singular P: warn, and the minimum-norm solution reproduces the right-hand side
  P <- matrix(c(1, 1, 1, 1), 2)
  expect_warning(w <- .index_weights(P, c(2, 2)), "singular")
  expect_equal(as.numeric(P %*% w), c(2, 2), tolerance = 1e-12)
  expect_equal(as.numeric(w), c(1, 1), tolerance = 1e-12)
  # a badly scaled but finite covariance gives finite weights
  w2 <- suppressWarnings(.index_weights(diag(c(1e-300, 1)), c(1e10, 1)))
  expect_true(all(is.finite(w2)))
})

test_that(".index_weights never returns non-finite weights for finite denormal-scale input (C6)", {
  # P and r both at denormal scale: b = P^{-1} r = (1, 2) is representable
  w <- .index_weights(diag(c(1e-310, 1e-310)), c(1e-310, 2e-310))
  expect_true(all(is.finite(w)))
  expect_equal(as.numeric(w), c(1, 2), tolerance = 1e-3)  # denormals carry few bits
  # P at denormal scale against an O(1) r: the true weights (1e310) are not
  # representable, so a clear error replaces silent NaN weights
  expect_error(.index_weights(diag(c(1e-310, 1e-310)), c(1, 1)), "overflow")
})

test_that(".index_weights is invariant to a common positive scaling of P and r (C6)", {
  set.seed(11)
  A <- crossprod(matrix(stats::rnorm(9), 3)) + diag(3)
  r <- c(1, -2, 0.5)
  b0 <- as.numeric(.index_weights(A, r))
  for (cc in c(1e-200, 1e-100, 1e100, 1e200)) {
    expect_equal(as.numeric(.index_weights(A * cc, r * cc)), b0, tolerance = 1e-10)
  }
})

test_that("an index on phenotypes at an extreme scale fails loudly, never silently (C2/C6)", {
  ph2 <- .asp(.as2(60), n_traits = 2)
  sm <- ph2
  sm$pheno$value <- sm$pheno$value * 1e-160     # P covariance ~ 1e-320 (denormal)
  expect_error(select_ind(sm, n = 5, method = "index", weights = c(1, 1)),
               "non-finite")
})

# ---- T16 / C16: scheme edges -----------------------------------------------------

test_that("scheme edge sizes: one line, one seed, one generation, pop_size below the selected count (T16)", {
  f2 <- .as2(60)
  one <- single_seed_descent(f2[1], generations = 1, seed = 1)
  expect_equal(n_individuals(one), 1L)
  expect_match(one$ids, "_g1$")
  b1 <- suppressMessages(bulk(f2, generations = 2, n = 1, seed = 1))
  expect_equal(n_individuals(b1), 1L)
  p1 <- suppressMessages(pedigree(f2, .asp1, generations = 1, prop = 0.2, seed = 1))
  expect_equal(nrow(attr(p1, "history")), 1L)
  expect_equal(n_individuals(p1), 60L)
  # more lines selected (30) than plants grown (10): the next generation is
  # trimmed to pop_size, and the history still records the 30 selected
  p2 <- suppressMessages(pedigree(f2, .asp1, generations = 2, prop = 0.5,
                                  pop_size = 10, seed = 1))
  expect_equal(n_individuals(p2), 10L)
  expect_equal(attr(p2, "history")$n_selected[1], 30L)
  expect_equal(attr(p2, "history")$n_selected[2], 5L)
  # the phenotyping callback must be a function
  expect_error(pedigree(f2, "not a function", generations = 1), "must be a function")
  expect_error(recurrent_selection(f2, NULL, cycles = 1, n_parents = 5),
               "must be a function")
})

test_that("bulk() reaches the requested size and keeps a connected, reproducible pedigree (C16)", {
  f2 <- .as2(60)[1:10]
  for (size in c(4L, 10L, 25L)) {                 # smaller, equal, larger than N
    bk <- suppressMessages(bulk(f2, generations = 3, n = size, seed = 8))
    expect_equal(n_individuals(bk), size)
    expect_false(anyDuplicated(bk$ids) > 0)
    pp <- parentage(bk, ancestors = TRUE)
    # every recorded parent resolves to a pedigree row; founders have none
    pk <- c(pp$mother_key, pp$father_key)
    expect_true(all(pk[!is.na(pk)] %in% pp$key))
    # the bulk descends from the 10 starting plants through 3 selfing generations
    expect_equal(max(pp$generation) - min(pp$generation[pp$key %in% parentage(f2)$key]),
                 3L)
  }
  a <- suppressMessages(bulk(f2, generations = 2, n = 25, seed = 8))
  b <- suppressMessages(bulk(f2, generations = 2, n = 25, seed = 8))
  expect_identical(dosages(a), dosages(b))
})

# ---- C22: scheme internals keep keys unique and parent links resolvable -------------

test_that(".self_each / .intermate give distinct keys and parents that resolve (C22)", {
  f2 <- .as2(20)
  kids <- suppressMessages(.self_each(f2[1:3], n_each = 2L, tag = "t"))
  pk <- parentage(kids)
  expect_false(anyDuplicated(pk$key) > 0)
  expect_true(all(pk$mother_key %in% f2$keys[1:3]))
  expect_true(all(pk$father_key == pk$mother_key))          # selfing
  expect_identical(pk$design, rep("self", 6L))
  x <- suppressMessages(.intermate(f2[1:4], n_crosses = 3L, progeny_per_cross = 4L,
                                   tag = "c"))
  px <- parentage(x)
  expect_equal(nrow(px), 12L)
  expect_false(anyDuplicated(px$key) > 0)
  expect_true(all(px$mother_key %in% f2$keys[1:4]))
  expect_true(all(px$father_key %in% f2$keys[1:4]))
  expect_true(all(px$mother_key != px$father_key))          # no selfing in a cross
  # each cross contributes exactly progeny_per_cross siblings
  expect_equal(sort(unique(as.integer(table(paste(px$mother_key, px$father_key,
                                                  sep = "|"))) %% 4L)), 0L)
})

# ---- C20: subset sims advance exactly their own individuals ---------------------------

test_that(".as_founder_pop returns exactly the simulated individuals, in order (C20)", {
  f2 <- .as2(40)
  sub <- suppressMessages(
    simulate_phenotype(f2, individuals = c(7, 3, 12, 20), h2 = 0.5, seed = 3) |>
      additive(n_qtn = 10))
  fp <- .as_founder_pop(sub)
  expect_equal(n_individuals(fp), 4L)
  expect_identical(fp$ids, sub$ids)
  expect_identical(fp$ids, f2$ids[sub$ind_idx])
  expect_identical(dosages(fp), dosages(f2)[, sub$ind_idx, drop = FALSE])
})

# ---- T1: response to selection, tight pooled check ----------------------------------

test_that("pooled over 20 populations, the realized response is R = Cov(A, P) / V_P * S (T1)", {
  skip_on_cran()
  f1 <- cross(as_population(SNP55K_maize282_maf04, individuals = 1:20)[1],
              as_population(SNP55K_maize282_maf04, individuals = 1:20)[2],
              n = 1, seed = 1)
  one <- function(s) {
    f <- suppressMessages(selfcross(f1, n = 300, seed = 100 + s))
    ph <- suppressMessages(simulate_phenotype(f, h2 = 0.5, seed = 200 + s) |>
                             additive(n_qtn = 20))
    sel <- select_ind(ph, prop = 0.2)
    A <- stats::setNames(.breeding_value_matrix(ph, 1)[, 1], ph$ids)
    P <- .criterion_values(ph, "pheno", 1, 1)
    c(obs = mean(A[attr(sel, "selected")]) - mean(A),
      hat = stats::cov(A, P) / stats::var(P) * attr(sel, "differential"))
  }
  res <- vapply(1:20, one, numeric(2))
  expect_equal(sum(res["obs", ]) / sum(res["hat", ]), 1, tolerance = 0.05)
  expect_true(all(res["obs", ] > 0))
})
