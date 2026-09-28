# test-culling.R -- DECISION-028: independent culling levels in select_ind().

.cull_geno <- function(n = 500, m = 600, seed = 1) {
  set.seed(seed)
  p <- stats::runif(m, 0.2, 0.8)
  g <- t(vapply(p, function(pp) stats::rbinom(n, 2, pp) - 1L, integer(n)))
  colnames(g) <- paste0("I", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:10, each = m / 10), pos = rep(seq_len(m / 10), 10),
                   cm = rep(seq(0, 100, length.out = m / 10), 10),
                   stringsAsFactors = FALSE), as.data.frame(g))
}
G <- .cull_geno()
.two_trait <- function(seed) {
  suppressMessages(simulate_phenotype(G, n_traits = 2, h2 = 0.5, seed = seed) |>
                     additive(n_qtn = 150))
}

test_that("kept set is the intersection of per-trait top sets", {
  sim <- .two_trait(1)
  sel <- suppressMessages(select_ind(sim, method = "culling", culling = c(0.3, 0.5)))
  y1 <- sim$pheno$value[sim$pheno$trait == "Trait_1"]
  y2 <- sim$pheno$value[sim$pheno$trait == "Trait_2"]
  top1 <- order(y1, decreasing = TRUE)[1:150]; top2 <- order(y2, decreasing = TRUE)[1:250]
  expect_setequal(sel, sim$ids[intersect(top1, top2)])
  expect_equal(attr(sel, "culling")$kept, c(150L, 250L))
  # sequential: trait 2 among trait-1 survivors
  seq_sel <- suppressMessages(select_ind(sim, method = "culling", culling = c(0.3, 0.5),
                                         sequential = TRUE))
  expect_setequal(seq_sel, sim$ids[top1[order(y2[top1], decreasing = TRUE)[1:75]]])
  # per-trait direction
  lo <- suppressMessages(select_ind(sim, method = "culling", culling = c(0.3, 0.5),
                                    direction = c("high", "low")))
  expect_setequal(lo, sim$ids[intersect(top1, order(y2)[1:250])])
})

test_that("one-trait culling equals mass selection at the same proportion", {
  sim <- .two_trait(2)
  a <- suppressMessages(select_ind(sim, method = "culling", culling = 0.2))
  b <- suppressMessages(select_ind(sim, prop = 0.2))
  expect_setequal(a, b)
})

test_that("an individuals x traits matrix supplies per-trait predictions", {
  sim <- .two_trait(3)
  M <- cbind(seq_len(500), seq_len(500) %% 7)
  sel <- suppressMessages(select_ind(sim, method = "culling", culling = c(0.5, 0.5),
                                     on = M))
  top1 <- order(M[, 1], decreasing = TRUE)[1:250]
  top2 <- order(M[, 2], decreasing = TRUE)[1:250]
  expect_setequal(sel, sim$ids[intersect(top1, top2)])
  expect_equal(attr(sel, "criterion"), "custom")
})

test_that("culling errors clearly", {
  sim <- .two_trait(4)
  expect_error(select_ind(sim, method = "culling", culling = c(0.3, 0.3), n = 10),
               "leave `n`")
  expect_error(select_ind(sim, method = "culling"), "needs `culling`")
  for (bad in list(c(1, NA), c(1, NaN), c(1, Inf))) {
    expect_error(select_ind(sim, method = "culling", culling = c(0.5, 0.5),
                            trait = bad), "distinct trait indices")
  }
  expect_error(select_ind(sim, method = "culling", culling = c(0.3, 1.2)), "\\(0, 1\\]")
  M <- cbind(seq_len(500), rev(seq_len(500)))
  expect_error(select_ind(sim, method = "culling", culling = c(0.4, 0.4), on = M),
               "no individual passes")
  expect_error(select_ind(sim, prop = 0.2, culling = c(0.3, 0.3)), "apply to")
})

test_that("index >= culling >= tandem, near 1 : 0.907 : 0.707 (T = 2, p = 0.1)", {
  gain <- function(sim, sel) {
    bv <- rowSums(.breeding_value_matrix(sim))
    mean(bv[match(sel, sim$ids)]) - mean(bv)
  }
  res <- t(vapply(1:30, function(s) {
    sim <- .two_trait(100 + s)
    idx <- suppressMessages(select_ind(sim, prop = 0.1, method = "index",
                                       weights = c(1, 1)))
    cul <- suppressMessages(select_ind(sim, method = "culling",
                                       culling = rep(sqrt(0.1), 2)))
    tan <- suppressMessages(select_ind(sim, prop = 0.1, trait = 1))
    c(gain(sim, idx), gain(sim, cul), gain(sim, tan))
  }, numeric(3)))
  m <- colMeans(res)
  expect_true(m[1] > m[2] && m[2] > m[3])
  expect_equal(m[2] / m[1], 0.907, tolerance = 0.06 / 0.907)
  expect_equal(m[3] / m[1], 0.707, tolerance = 0.06 / 0.707)
})
