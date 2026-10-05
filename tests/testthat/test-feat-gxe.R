# Fixed-scale G x E: gxe_value() and phenotype_value(gxe =, env =, var_env =)
# (breedingDesigner SPEC-0020 item 4; AlphaSimR addTraitAG / setPheno(p =)).

gxe_pop <- function() {
  data("SNP55K_maize282_maf04", envir = environment())
  as_population(SNP55K_maize282_maf04, individuals = 1:30)
}
q3 <- c(1, 5, 9)
a3 <- c(0.5, -1, 2)
b3 <- c(0.2, -0.1, 0.3)

test_that("gxe_value is intercept + sum of dosage x slope effect, unscaled", {
  pop <- gxe_pop()
  z <- dosages(pop)[q3, , drop = FALSE]
  expect_equal(gxe_value(pop, q3, b3), colSums(z * b3))
  expect_equal(gxe_value(pop, q3, b3, intercept = 0.7), colSums(z * b3) + 0.7)
  # same scale as additive_value(): the slope is an additive value of b
  expect_equal(gxe_value(pop, q3, b3), additive_value(pop, q3, b3))
  expect_named(gxe_value(pop, q3, b3), pop$ids)
})

test_that("gxe_value validates effect and intercept", {
  pop <- gxe_pop()
  expect_error(gxe_value(pop, q3, b3[1:2]), "one value per locus")
  expect_error(gxe_value(pop, q3, c(0.1, NA, 0.2)), "finite numeric")
  expect_error(gxe_value(pop, q3, b3, intercept = c(0, 1)), "single finite")
  expect_error(gxe_value(pop, q3, b3, intercept = Inf), "single finite")
})

test_that("phenotype_value without gxe is unchanged (value and stream)", {
  pop <- gxe_pop()
  y <- phenotype_value(pop, q3, a3, var_e = 0.5, seed = 4)
  set.seed(4)
  e <- rnorm(n_individuals(pop), sd = sqrt(0.5))
  expect_equal(unname(as.numeric(y)), unname(additive_value(pop, q3, a3) + e))
  expect_null(attr(y, "gxe_value"))
  expect_null(attr(y, "env"))
})

test_that("phenotype_value(gxe) follows AlphaSimR calcPheno: y = g + s * qnorm(p, sd) + e", {
  pop <- gxe_pop()
  g <- additive_value(pop, q3, a3)
  n <- n_individuals(pop)
  for (ve_env in c(0, 2.5)) {
    s <- gxe_value(pop, q3, b3, intercept = 1)
    y <- phenotype_value(pop, q3, a3, var_e = 0.5, gxe = b3, gxe_intercept = 1,
                         env = 0.3, var_env = ve_env, seed = 9)
    sd_w <- if (ve_env == 0) 1 else sqrt(ve_env)   # AlphaSimR: envVar = 1 when varEnv = 0
    w <- qnorm(0.3, sd = sd_w)
    set.seed(9)
    e <- rnorm(n, sd = sqrt(0.5))
    expect_equal(unname(as.numeric(y)), unname(g + s * w + e))
    expect_equal(attr(y, "gxe_value"), s)
    expect_equal(attr(y, "env"), 0.3)
    expect_equal(attr(y, "env_value"), w)
    expect_equal(attr(y, "genetic_value"), g)   # G x E is not genetic value
  }
})

test_that("a NULL env is drawn by runif before the residual (setPheno order)", {
  pop <- gxe_pop()
  n <- n_individuals(pop)
  y <- phenotype_value(pop, q3, a3, var_e = 0.5, gxe = b3, seed = 11)
  set.seed(11)
  p <- runif(1)
  e <- rnorm(n, sd = sqrt(0.5))
  expect_equal(attr(y, "env"), p)
  expect_equal(unname(as.numeric(y)),
               unname(additive_value(pop, q3, a3) +
                        gxe_value(pop, q3, b3) * qnorm(p) + e))
  # the caller's RNG state is restored
  set.seed(1); before <- runif(1)
  set.seed(1); phenotype_value(pop, q3, a3, var_e = 0.5, gxe = b3, seed = 11)
  expect_identical(runif(1), before)
})

test_that("the same env and seed give a common covariate across calls (one trial)", {
  pop <- gxe_pop()
  y1 <- phenotype_value(pop[1:10], q3, a3, var_e = 0, gxe = b3, env = 0.8)
  y2 <- phenotype_value(pop, q3, a3, var_e = 0, gxe = b3, env = 0.8)
  expect_equal(y1, y2[1:10], ignore_attr = TRUE)
})

test_that("h2 sets var_e from the genetic variance only (G x E excluded)", {
  pop <- gxe_pop()
  y0 <- phenotype_value(pop, q3, a3, h2 = 0.4, seed = 1)
  y1 <- phenotype_value(pop, q3, a3, h2 = 0.4, gxe = b3, env = 0.9, seed = 1)
  expect_equal(attr(y1, "var_e"), attr(y0, "var_e"))
})

test_that("G x E arguments are validated before any draw", {
  pop <- gxe_pop()
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, env = 0.5),
               "give the per-locus G x E effects")
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, var_env = 2),
               "give the per-locus G x E effects")
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, gxe_intercept = 1),
               "give the per-locus G x E effects")
  for (bad in list(0, 1, -0.1, NA_real_, c(0.2, 0.3), "0.5")) {
    expect_error(phenotype_value(pop, q3, a3, var_e = 1, gxe = b3, env = bad),
                 "single probability")
  }
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, gxe = b3, var_env = -1),
               "non-negative")
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, gxe = b3[1:2]),
               "`gxe` must be a finite numeric vector")
  expect_error(phenotype_value(pop, q3, a3, var_e = 1, gxe = b3,
                               gxe_intercept = NA_real_),
               "`gxe_intercept` must be a single finite")
  set.seed(5); before <- runif(1)
  set.seed(5)
  try(phenotype_value(pop, q3, a3, var_e = 1, gxe = b3, env = 2), silent = TRUE)
  expect_identical(runif(1), before)
})

test_that("across environments the G x E variance is Var(s) * var_env (AlphaSimR scaling)", {
  # AlphaSimR scales slopes so Var(s) = varGxE / varEnv; then over trials
  # (w ~ N(0, varEnv)) the G x E term s_i * w has variance varGxE around the
  # environment main effect when the mean slope is 1. Check the identity
  # Var_w(s_i w - mean(s) w) = Var(s) * var_env on a quadrature over p.
  pop <- gxe_pop()
  s <- gxe_value(pop, q3, b3, intercept = 1)
  p <- (seq_len(4000) - 0.5) / 4000
  w <- qnorm(p, sd = sqrt(2))
  dev <- outer(s - mean(s), w)            # individuals x trials
  expect_equal(mean(colMeans(dev^2)), mean((s - mean(s))^2) * 2, tolerance = 2e-3)
})

test_that("h2 with numeric qtn scores a reordered ref at the same markers (Codex V2)", {
  x <- matrix(c(-1, 0, 1, 1, 1, 1, -1, 1), 2, byrow = TRUE,
              dimnames = list(c("m1", "m2"), paste0("i", 1:4)))
  ref <- x[2:1, ]
  by_index <- phenotype_value(x, 1:2, c(1, 10), h2 = 0.5, ref = ref, seed = 1)
  by_name <- phenotype_value(x, c("m1", "m2"), c(1, 10), h2 = 0.5, ref = ref, seed = 1)
  expect_equal(attr(by_index, "var_e"), attr(by_name, "var_e"))
  expect_equal(attr(by_index, "var_e"), var(additive_value(x, 1:2, c(1, 10))))
  # without marker names, indices are the only identity (unchanged behaviour)
  xu <- x; ru <- ref; rownames(xu) <- rownames(ru) <- NULL
  expect_equal(attr(phenotype_value(xu, 1:2, c(1, 10), h2 = 0.5, ref = ru, seed = 1), "var_e"),
               var(additive_value(ru, 1:2, c(1, 10))))
})

test_that("print.Population ignores unused chromosome factor levels (Codex I1)", {
  m <- data.frame(snp = c("a", "b"), chr = factor(c("1", "1"), levels = c("1", "2")),
                  pos = 1:2, cm = c(0, 100))
  p <- population_from_haplotypes(cbind(i = c(1, 0)), cbind(i = c(0, 1)), m)
  expect_no_warning(out <- capture.output(print(p)))
  expect_match(out, "100 cM total span", all = FALSE)
})

test_that("a named ref may be a marker subset or carry NA at non-causal loci", {
  x <- matrix(c(-1, 0, 1, 1, 0, -1, -1, 1, 0, 0, 1, -1), nrow = 3, byrow = TRUE,
              dimnames = list(c("m1", "m2", "m3"), paste0("i", 1:4)))
  ref <- x[c("m2", "m1", "m3"), ]
  ref["m2", 1] <- NA
  expect_equal(attr(phenotype_value(x, 1, 2, h2 = 0.5, ref = ref, seed = 1), "var_e"),
               var(additive_value(ref, "m1", 2)))
  ref2 <- x["m3", , drop = FALSE]
  expect_equal(attr(phenotype_value(x, 3, 2, h2 = 0.5, ref = ref2, seed = 1), "var_e"),
               var(additive_value(ref2, "m3", 2)))
})

test_that("duplicated marker names fall back to row indices", {
  x <- matrix(c(-1, -1, 1, 1, -1, 0, 0, 1), nrow = 2, byrow = TRUE,
              dimnames = list(c("m", "m"), paste0("i", 1:4)))
  expect_equal(attr(phenotype_value(x, 2, 2, h2 = 0.5, ref = x, seed = 1), "var_e"),
               var(additive_value(x, 2, 2)))
})

test_that("trial means order by env only with a main effect (mean slope 1, var_env > 0)", {
  pop <- gxe_pop()
  args <- list(x = pop, qtn = q3, effect = a3, var_e = 0.5, gxe = b3, seed = 3)
  # interaction only (defaults): the mean shift is mean(s) * (w_hi - w_lo), not a main effect
  lo <- do.call(phenotype_value, c(args, env = 0.1))
  hi <- do.call(phenotype_value, c(args, env = 0.9))
  s0 <- gxe_value(pop, q3, b3)
  expect_equal(mean(hi) - mean(lo), mean(s0) * (qnorm(0.9) - qnorm(0.1)))
  # main effect: mean slope 1 in the base population, as AlphaSimR sets gxeInt
  int <- 1 - mean(s0)
  lo <- do.call(phenotype_value, c(args, env = 0.1, gxe_intercept = int, var_env = 2))
  hi <- do.call(phenotype_value, c(args, env = 0.9, gxe_intercept = int, var_env = 2))
  expect_equal(mean(hi) - mean(lo), qnorm(0.9, sd = sqrt(2)) - qnorm(0.1, sd = sqrt(2)))
  expect_gt(mean(hi), mean(lo))
})
