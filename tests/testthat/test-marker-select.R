# test-marker-select.R -- DECISION-029: MAS / pyramiding / marker index.

.two_founders <- function(chr, cm) {
  k <- length(chr)
  g <- data.frame(snp = paste0("m", seq_len(k)), allele = "A/G", chr = chr,
                  pos = seq_len(k), cm = cm, P1 = 1L, P2 = -1L,
                  stringsAsFactors = FALSE)
  as_population(g)
}

test_that("unlinked targets: F2 feasibility follows 1/16 and 9/16", {
  fnd <- .two_founders(chr = 1:2, cm = c(0, 0))
  f1 <- cross(fnd[1], fnd[2], seed = 1)
  for (s in 1:5) {
    f2 <- selfcross(f1, n = 2000, seed = 10 + s)
    hom <- attr(marker_select(f2, c("m1", "m2"), requirement = "homozygote"),
                "marker_select")
    car <- attr(marker_select(f2, c("m1", "m2"), requirement = "carrier"),
                "marker_select")
    se <- function(p) sqrt(p * (1 - p) / 2000)
    expect_lt(abs(mean(hom$feasible) - 1 / 16), 3 * se(1 / 16))
    expect_lt(abs(mean(car$feasible) - 9 / 16), 3 * se(9 / 16))
  }
})

# cross() draws Poisson(d Morgans) crossovers (no interference), so two loci d
# Morgans apart recombine when the count is odd: r = (1 - exp(-2 d)) / 2
# (Haldane's map function).
test_that("linked targets follow the Poisson map: P(double favourable homozygote) = ((1-r)/2)^2", {
  for (d in c(5, 10, 50)) {
    fnd <- .two_founders(chr = c(1, 1), cm = c(0, d))
    f1 <- cross(fnd[1], fnd[2], seed = 2)                 # AB/ab, coupling
    r <- 0.5 * (1 - exp(-2 * d / 100))
    target <- ((1 - r) / 2)^2
    ph <- vapply(1:4, function(s) {
      f2 <- selfcross(f1, n = 3000, seed = 100 * d + s)
      mean(attr(marker_select(f2, c("m1", "m2"), requirement = "homozygote"),
                "marker_select")$feasible)
    }, numeric(1))
    expect_lt(abs(mean(ph) - target), 3 * sqrt(target * (1 - target) / 12000))
  }
})

test_that("staged pyramid, n, errors and determinism", {
  fnd <- .two_founders(chr = 1:4, cm = rep(0, 4))
  f2 <- selfcross(cross(fnd[1], fnd[2], seed = 1), n = 400, seed = 3)
  st <- marker_select(f2, paste0("m", 1:4), requirement = "homozygote",
                      min_markers = 3)
  d <- attr(st, "marker_select")
  expect_identical(d$feasible, d$n_met >= 3)
  expect_true(attr(st, "marker_select_info")$staged)
  a <- marker_select(f2, paste0("m", 1:4), min_markers = 3, n = 20, seed = 5)
  b <- marker_select(f2, paste0("m", 1:4), min_markers = 3, n = 20, seed = 5)
  expect_identical(attr(a, "marker_select"), attr(b, "marker_select"))
  expect_equal(n_individuals(a), 20L)
  expect_error(marker_select(f2, paste0("m", 1:4), requirement = "homozygote",
                             n = 400), "cannot select")
  expect_error(marker_select(f2, "m1", favorable = 2), "1 or -1")
  expect_error(marker_select(f2, "m1", min_markers = 2), "0..1")
  # favorable = -1 flips which homozygote counts
  fl <- attr(marker_select(f2, "m1", favorable = -1, requirement = "homozygote"),
             "marker_select")
  expect_identical(fl$feasible, unname(dosages(f2)["m1", ] == -1))
})

test_that("foreground matches mabc_select() with the donor allele favourable", {
  g <- data.frame(snp = paste0("m", 1:40), allele = "A/G", chr = rep(1:2, each = 20),
                  pos = rep(1:20, 2), cm = rep(seq(0, 95, by = 5), 2),
                  P1 = 1L, P2 = -1L)
  fnd <- as_population(g)
  bc1 <- cross(cross(fnd[1], fnd[2], seed = 1), fnd[1], n = 200, seed = 2)
  ms <- attr(marker_select(bc1, "m5", favorable = -1), "marker_select")
  mb <- attr(mabc_select(bc1, fnd[1], fnd[2], target_markers = "m5", n = 1,
                         seed = 1), "mabc")
  expect_identical(ms$feasible, mb$feasible)
  rec <- recurrent_parent_recovery(bc1, fnd[1], fnd[2])
  rk <- attr(marker_select(bc1, "m5", favorable = -1, rank_on = rec, seed = 1),
             "marker_select")
  expect_equal(order(rk$rank)[seq_len(sum(rk$feasible))][1:5],
               order(mb$rank)[1:5])
})

test_that("oracle MARS: an index on all causal effects reproduces on = 'bv'", {
  set.seed(4)
  m <- 200; n <- 300
  p <- stats::runif(m, 0.2, 0.8)
  gm <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(gm) <- paste0("I", seq_len(n))
  G <- cbind(data.frame(snp = paste0("s", seq_len(m)), allele = "A/G",
                        chr = rep(1:4, each = m / 4), pos = rep(seq_len(m / 4), 4),
                        cm = rep(seq(0, 100, length.out = m / 4), 4)), as.data.frame(gm))
  pop <- as_population(G)
  sim <- simulate_phenotype(pop, h2 = 0.5, seed = 5) |> additive(n_qtn = 40)
  te <- template_effects(sim)
  idx <- additive_value(pop, te$qtn, te$a)
  mars <- marker_select(pop, markers = 1, min_markers = 0, n = 30, rank_on = idx)
  bvsel <- select_ind(sim, n = 30, on = "bv")
  # the same individuals up to ties at the cut-off (the two break ties
  # differently): compare the kept scores
  expect_equal(sort(unname(idx[bvsel$ids])), sort(unname(idx[mars$ids])))
  expect_equal(sort(unname(idx[mars$ids])), sort(unname(idx), decreasing = TRUE)[30:1])
})

test_that("with dominance the oracle weights are the average effects, not raw a (review B)", {
  set.seed(6)
  m <- 200; n <- 300
  p0 <- stats::runif(m, 0.2, 0.8)
  gm <- t(vapply(p0, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(gm) <- paste0("I", seq_len(n))
  G <- cbind(data.frame(snp = paste0("s", seq_len(m)), allele = "A/G",
                        chr = rep(1:4, each = m / 4), pos = rep(seq_len(m / 4), 4),
                        cm = rep(seq(0, 100, length.out = m / 4), 4)), as.data.frame(gm))
  pop <- as_population(G)
  sim <- simulate_phenotype(pop, h2 = 0.5, seed = 7) |>
    additive(prop = 0.3, n_qtn = 40) |> dominance(prop = 0.2)
  te <- template_effects(sim)
  x <- dosages(pop)[te$qtn, , drop = FALSE] + 1
  p <- rowMeans(x) / 2
  alpha <- te$a + te$d * (1 - 2 * p)
  idx <- additive_value(pop, te$qtn, alpha)
  oracle <- marker_select(pop, markers = 1, min_markers = 0, n = 30, rank_on = idx)
  bvsel <- select_ind(sim, n = 30, on = "bv")
  expect_equal(sort(unname(idx[bvsel$ids])), sort(unname(idx[oracle$ids])))
  naive <- additive_value(pop, te$qtn, te$a)            # raw a: not the oracle
  expect_false(setequal(names(sort(naive, decreasing = TRUE))[1:30], bvsel$ids))
})

test_that("oracle selection agrees with on = 'bv' up to ties (review B r2)", {
  g <- data.frame(snp = "s1", allele = "A/G", chr = 1, pos = 1, cm = 0,
                  I1 = 1L, I2 = 0L, I3 = 0L, I4 = -1L, stringsAsFactors = FALSE)
  pop <- as_population(g)
  sc <- additive_value(pop, 1, 1)                      # 1, 0, 0, -1
  for (s in 1:4) {
    m <- marker_select(pop, markers = 1, min_markers = 0, n = 2, rank_on = sc,
                       seed = s)
    expect_true("I1" %in% m$ids)
    expect_equal(sort(unname(sc[m$ids])), c(0, 1))       # a tie at the cut-off
  }
})
