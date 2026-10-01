# test-crossbreed.R -- DECISION-031: breed composition, heterosis, systems.

.breed <- function(pool, n, p, seed, label = pool) {
  old <- .Random.seed_safe()
  on.exit(.restore_seed(old), add = TRUE)
  set.seed(seed)
  m <- length(p)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0(pool, seq_len(n))
  as_population(cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                                 chr = rep(1:10, each = m / 10),
                                 pos = rep(seq_len(m / 10), 10),
                                 cm = rep(seq(0, 100, length.out = m / 10), 10),
                                 stringsAsFactors = FALSE), as.data.frame(g)),
                pool = label)
}
# allele frequencies drawn under a local seed: no test file changes the RNG
# state of the session
.freqs <- function() {
  old <- .Random.seed_safe()
  on.exit(.restore_seed(old))
  set.seed(1)
  list(A = stats::runif(200, 0.05, 0.95), B = stats::runif(200, 0.05, 0.95),
       C = stats::runif(200, 0.05, 0.95))
}
FR <- .freqs()
PA <- FR$A; PB <- FR$B; PC <- FR$C
BR <- list(A = .breed("A", 150, PA, 2), B = .breed("B", 150, PB, 3),
           C = .breed("C", 150, PC, 4))

test_that("expected breed composition: F1, backcross, three-way, rotation", {
  f1 <- crossbreed(BR[1:2], "two_way", n_progeny = 20, seed = 1)
  expect_equal(unname(colMeans(breed_composition(f1))), c(0.5, 0.5))
  bc <- crossbreed(BR[1:2], "backcross", n_progeny = 20, seed = 2)
  expect_equal(unname(colMeans(breed_composition(bc))), c(0.75, 0.25))
  tw <- crossbreed(BR, "three_way", n_progeny = 20, seed = 3)
  expect_equal(unname(colMeans(breed_composition(tw))), c(0.25, 0.25, 0.5))
  te <- crossbreed(BR, "terminal", n_progeny = 20, sire_breed = "C", seed = 4)
  expect_equal(unname(colMeans(breed_composition(te))), c(0.25, 0.25, 0.5))
  rot <- crossbreed(BR[1:2], "rotational", n_progeny = 20, generations = 10,
                    seed = 5)
  h <- attr(rot, "history")
  expect_equal(nrow(h), 11L)
  a_frac <- h$A[h$generation >= 10]
  expect_lt(min(abs(a_frac - 2 / 3), abs(a_frac - 1 / 3)), 1e-3)
  expect_equal(unname(rowSums(breed_composition(rot))), rep(1, 20))
})

test_that("realized F1 heterosis tracks its expectation; additive gives zero", {
  set.seed(6); a <- stats::rnorm(200, sd = 0.2); d <- abs(stats::rnorm(200, sd = 0.3))
  r <- vapply(1:12, function(s) {
    f1 <- crossbreed(BR[1:2], "two_way", n_progeny = 150, seed = 10 + s)
    heterosis(f1, BR[1:2], qtn = 1:200, a = a, d = d)$realized
  }, numeric(1))
  h <- heterosis(crossbreed(BR[1:2], "two_way", n_progeny = 10, seed = 1),
                 BR[1:2], qtn = 1:200, a = a, d = d)
  expect_lt(abs(mean(r) - h$expected_f1["A", "B"]),
            3 * stats::sd(r) / sqrt(12) + 1e-8)
  # under within-breed HWE the expectation is close to sum d (pA - pB)^2
  pa <- rowMeans(dosages(BR$A) + 1) / 2; pb <- rowMeans(dosages(BR$B) + 1) / 2
  expect_equal(h$expected_f1["A", "B"], sum(d * (pa - pb)^2), tolerance = 0.1)
  h0 <- heterosis(crossbreed(BR[1:2], "two_way", n_progeny = 10, seed = 1),
                  BR[1:2], qtn = 1:200, a = a, d = 0)
  expect_equal(h0$expected_f1["A", "B"], 0, tolerance = 1e-10)
  expect_true(all(is.na(diag(h$expected_f1))))
  expect_equal(h$expected_f1["A", "B"], h$expected_f1["B", "A"])
  # an inbred-line pair: the expectation is sum d (pA + pB - 2 pA pB)
  g <- data.frame(snp = "q", allele = "A/G", chr = 1, pos = 1, cm = 0, L1 = 1L,
                  L2 = -1L, stringsAsFactors = FALSE)
  l1 <- as_population(g[, 1:6], pool = "A"); l2 <- as_population(g[, -6], pool = "B")
  hl <- heterosis(cross(l1, l2, seed = 1), list(A = l1, B = l2), qtn = 1, a = 1,
                  d = 0.4)
  expect_equal(hl$expected_f1["A", "B"], 0.4)
  expect_equal(hl$realized, 0.4)
})

test_that("F2 keeps half the F1 heterosis; a two-breed rotation about two thirds", {
  set.seed(7); a <- stats::rnorm(200, sd = 0.2); d <- abs(stats::rnorm(200, sd = 0.3))
  hf1 <- heterosis(crossbreed(BR[1:2], "two_way", n_progeny = 10, seed = 1),
                   BR[1:2], qtn = 1:200, a = a, d = d)$expected_f1["A", "B"]
  f2 <- vapply(1:10, function(s) {
    f1 <- crossbreed(BR[1:2], "two_way", n_progeny = 150, seed = 20 + s)
    plan <- mating_design(f1, design = "random", n_crosses = 150, seed = 40 + s)
    heterosis(mate(plan, f1, seed = 60 + s), BR[1:2], qtn = 1:200, a = a,
              d = d)$realized
  }, numeric(1))
  expect_equal(mean(f2) / hf1, 0.5, tolerance = 0.1 / 0.5)
  rot <- vapply(1:6, function(s) {
    r <- crossbreed(BR[1:2], "rotational", n_progeny = 150, generations = 8,
                    seed = 80 + s)
    heterosis(r, BR[1:2], qtn = 1:200, a = a, d = d)$realized
  }, numeric(1))
  expect_equal(mean(rot) / hf1, 2 / 3, tolerance = 0.12 / (2 / 3))
})

test_that("input checks", {
  expect_error(crossbreed(BR[1], "two_way", n_progeny = 5), "at least 2")
  expect_error(crossbreed(BR, "terminal", n_progeny = 5, sire_breed = "A"),
               "other than the two")
  expect_error(heterosis(crossbreed(BR, "three_way", n_progeny = 5, seed = 1),
                         BR[1:2], qtn = 1:5, a = rep(1, 5)), "not in `breeds`")
})

test_that("breed names must be the founder pools (review D)", {
  sw <- list(A = BR$B, B = BR$A)                          # swapped labels
  expect_error(crossbreed(sw, "two_way", n_progeny = 5, seed = 1),
               "pool label \"A\"")
  f1 <- crossbreed(BR[1:2], "two_way", n_progeny = 5, seed = 1)
  expect_error(heterosis(f1, sw, qtn = 1:5, a = rep(1, 5)), "pure breed")
  mixed <- list(A = c(BR$A, BR$B), B = BR$B)
  expect_error(crossbreed(mixed, "two_way", n_progeny = 5), "founder pool\\(s\\): A, B")
  unl <- .breed("A", 20, PA, 2, label = NA_character_)
  expect_error(crossbreed(list(A = unl, B = BR$B), "two_way", n_progeny = 5),
               "<none>")
})

test_that("reserved pool labels cannot collide with unlabelled founders (review D r3)", {
  g <- data.frame(snp = "q", allele = "A/G", chr = 1, pos = 1, cm = 0, L1 = 1L,
                  stringsAsFactors = FALSE)
  expect_error(as_population(g, pool = ""), "non-empty")
  expect_error(as_population(g, pool = "<unassigned>"), "reserved")
  u <- as_population(g)                                  # no label
  b <- as_population(transform(g, L1 = -1L), pool = "B")
  x <- cross(u, b, seed = 1)
  expect_equal(colnames(breed_composition(x)), c("<unassigned>", "B"))
  expect_error(heterosis(x, list(B = b), qtn = 1, a = 1), "not in `breeds`")
})
