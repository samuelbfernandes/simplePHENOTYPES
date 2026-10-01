# Round-5 feature tests for the grammar group (SPEC-0020 items 1 and 6).

# --- Item 1: cheap, bounded naming of `geno` -------------------------------

# A Population large enough that deparse() of it would take many seconds,
# built cheaply from the internal constructor (no meiosis, no crossing).
.big_population <- function(n_ind = 1000L, n_mark = 6000L) {
  set.seed(1)
  map <- simplePHENOTYPES:::.make_map(
    snp = paste0("m", seq_len(n_mark)), chr = rep(1L, n_mark),
    pos = seq_len(n_mark), cm = seq_len(n_mark) * 0.1)  # centiMorgans
  cis <- matrix(stats::rbinom(n_ind * n_mark, 1L, 0.4), n_mark, n_ind)
  trans <- matrix(stats::rbinom(n_ind * n_mark, 1L, 0.4), n_mark, n_ind)
  ids <- paste0("i", seq_len(n_ind))
  dimnames(cis) <- dimnames(trans) <- list(map$snp, ids)
  simplePHENOTYPES:::.new_population(map, cis, trans, ids, "founder")
}

test_that("item 1: a large inline geno is named in milliseconds, never deparsed", {
  pop <- .big_population()
  namer <- function(geno) simplePHENOTYPES:::.geno_label(substitute(geno))
  elapsed <- system.time(lab <- do.call(namer, list(geno = pop)))[["elapsed"]]
  expect_lt(elapsed, 0.5)                       # deparse() of this takes seconds
  expect_identical(lab, "<inline Population>")
  mat <- matrix(0, 500, 400)
  elapsed <- system.time(lab <- do.call(namer, list(geno = mat)))[["elapsed"]]
  expect_lt(elapsed, 0.5)
  expect_identical(lab, "<inline matrix 500 x 400>")
  df <- as.data.frame(matrix(0, 50, 2000))
  expect_lt(system.time(lab <- do.call(namer, list(geno = df)))[["elapsed"]], 0.5)
  expect_identical(lab, "<inline data.frame 50 x 2000>")
})

test_that("item 1: the label of ordinary calls is the old deparse() text", {
  lab <- simplePHENOTYPES:::.geno_label
  old <- function(e) deparse(e)
  exprs <- list(quote(geno), quote(SNP55K_maize282_maf04), quote(df[1:3, ]),
                quote(head(x, 3)), quote(sub$geno), quote(m[, c(1, 2, 3)]),
                quote(pkg::obj))
  for (e in exprs) expect_identical(lab(e), old(e))
  expect_identical(lab(NULL), "NULL")
  # a long call keeps only its first deparse line (bounded), still identical to
  # the first element of the old multi-line deparse
  long <- quote(f(aaaaaaaaaaaa, bbbbbbbbbbbbbb, cccccccccccccc, dddddddddddddd,
                  eeeeeeeeeeeeee, ffffffffffffff, gggggggggggggg))
  expect_identical(lab(long), old(long)[1L])
  # a call that embeds a large object is not deparsed either
  big_call <- as.call(list(as.name("head"), matrix(0, 300, 300)))
  expect_identical(lab(big_call), "<inline call>")
  # the time is bounded for a call carrying a huge constant, too
  expect_lt(system.time(lab(big_call))[["elapsed"]], 0.5)
})

test_that("item 1: simulate_phenotype() through do.call() with an inline object works", {
  data("SNP55K_maize282_maf04")
  geno <- SNP55K_maize282_maf04
  direct <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 2)
  inline <- do.call(simulate_phenotype,
                    list(geno = geno, h2 = 0.5, n_qtn = 3, seed = 2))
  expect_identical(direct$geno_name, "geno")
  expect_identical(inline$geno_name,
                   sprintf("<inline data.frame %d x %d>", nrow(geno), ncol(geno)))
  expect_identical(direct$pheno, inline$pheno)
  # a Population keeps its own label and the naming step does not deparse it
  pop <- .big_population(60L, 300L)
  elapsed <- system.time(
    ph <- do.call(simulate_phenotype, list(geno = pop, h2 = 0.5, n_qtn = 3, seed = 1))
  )[["elapsed"]]
  expect_identical(ph$geno_name, "<Population: founder>")
  expect_lt(elapsed, 20)
})

# --- Item 6: entry-mean replication (`reps`) -------------------------------

.fx <- function() {
  data("SNP55K_maize282_maf04", envir = environment())
  get("SNP55K_maize282_maf04", envir = environment())
}

test_that("item 6: reps = 1 (default or explicit) is bit-identical to HEAD", {
  geno <- .fx()
  # reference values recorded from the pre-`reps` code (same seed and call)
  a <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 1)
  expect_equal(a$pheno$value[1:4],
               c(1.0748089703966, -0.158761061241073, 1.15789014330759,
                 -2.22997737383096), tolerance = 1e-13)
  expect_equal(sum(a$pheno$value * seq_along(a$pheno$value)),
               -2690.81707705573, tolerance = 1e-12)
  a1 <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 1, reps = 1)
  expect_identical(a1$pheno, a$pheno)
  expect_identical(a1$var_budget, a$var_budget)

  suppressWarnings({
    b <- simulate_phenotype(geno, n_traits = 2, architecture = "pleiotropy",
                            cor = 0.3, seed = 7)
    b <- additive(b, prop = 0.4, n_qtn = 4)
    b <- vqtl(b, prop = 0.1, n_qtn = 2)
    b1 <- simulate_phenotype(geno, n_traits = 2, architecture = "pleiotropy",
                             cor = 0.3, seed = 7, reps = c(1, 1))
    b1 <- additive(b1, prop = 0.4, n_qtn = 4)
    b1 <- vqtl(b1, prop = 0.1, n_qtn = 2)
  })
  expect_equal(b$pheno$value[1:4],
               c(1.33445703484397, 0.109673804273308, -1.38364684046686,
                 -0.138993794787653), tolerance = 1e-13)
  expect_equal(sum(b$pheno$value * seq_along(b$pheno$value)),
               -3613.65028357235, tolerance = 1e-12)
  expect_identical(b1$pheno, b$pheno)
  # the default print is unchanged: no replication note
  expect_false(any(grepl("Entry means", utils::capture.output(print(a)))))
})

test_that("item 6: residual variance is exactly var_e / reps; draws are unchanged", {
  geno <- .fx()
  for (rp in c(2L, 4L, 10L)) {
    for (s in 1:25) {
      noise <- simulate_phenotype(geno, seed = s, reps = rp)   # pure residual, V_E = 1
      expect_equal(stats::var(noise$pheno$value), 1 / rp, tolerance = 1e-12)
    }
  }
  # same RNG stream: reps = r is the reps = 1 residual times 1 / sqrt(r)
  n1 <- simulate_phenotype(geno, seed = 11)
  n4 <- simulate_phenotype(geno, seed = 11, reps = 4)
  expect_equal(n4$pheno$value, n1$pheno$value / 2, tolerance = 1e-14)
  # with a genetic part: V_E = 1 - prop, so Var(y - g) = (1 - prop) / reps
  ph <- simulate_phenotype(geno, h2 = 0.4, n_qtn = 3, seed = 5, reps = 5)
  g <- genetic_values(ph)[, 1]
  expect_equal(stats::var(g), 0.4, tolerance = 1e-12)
  expect_equal(stats::var(ph$pheno$value - g), 0.6 / 5, tolerance = 1e-12)
  # the pipe carries reps through later layers
  p <- simulate_phenotype(geno, seed = 5, reps = 5)
  p <- additive(p, prop = 0.4, n_qtn = 3)
  expect_equal(stats::var(p$pheno$value - genetic_values(p)[, 1]), 0.6 / 5,
               tolerance = 1e-12)
  expect_identical(p$reps, 5L)
})

test_that("item 6: reps = r matches the mean of r independent records in variance", {
  # Semantics check (AlphaSimR setPheno(varE, reps)): the entry mean of r iid
  # records of variance V_E has variance V_E / r. Average r reps = 1 residuals
  # from independent seeds and compare the Monte Carlo variance of the mean with
  # the single-draw reps = r value, with a standard-error tolerance. The records
  # here are the exact-variance residual draws of the grammar, so the Monte
  # Carlo variance of their mean is 1 / r up to the (small) between-record
  # sample covariance.
  geno <- .fx()
  r <- 4L
  recs <- vapply(1:r, function(s) simulate_phenotype(geno, seed = 100 + s)$pheno$value,
                 numeric(nrow(simulate_phenotype(geno, seed = 1)$pheno)))
  v_mean <- stats::var(rowMeans(recs))
  n <- nrow(recs)
  se <- (1 / r) * sqrt(2 / (n - 1))     # SE of a variance estimate, normal approx
  expect_lt(abs(v_mean - 1 / r), 4 * se)
  v_reps <- stats::var(simulate_phenotype(geno, seed = 100, reps = r)$pheno$value)
  expect_equal(v_reps, 1 / r, tolerance = 1e-12)
})

test_that("item 6: entry-mean H2 = Vg / (Vg + Ve/reps); h2 stays single-record", {
  geno <- .fx()
  ph1 <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 3)
  ph4 <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 3, reps = 4)
  g <- genetic_values(ph4)[, 1]
  y <- ph4$pheno$value
  e <- y - g
  vg <- stats::var(g); ve <- stats::var(e)
  # realized (stored phenotype = entry mean) H2 is Vg / Var(g + e)
  expect_equal(simplePHENOTYPES:::.realized_h2(ph4), vg / stats::var(y),
               tolerance = 1e-12)
  # analytic target 0.5 / (0.5 + 0.5 / 4) = 0.8 (cov(g, e) is only sampling noise)
  expect_equal(vg / (vg + ve), 0.5 / (0.5 + 0.5 / 4), tolerance = 1e-12)
  expect_lt(abs(simplePHENOTYPES:::.realized_h2(ph4) - 0.8), 0.1)
  # the single-record reconstruction recovers exactly the reps = 1 realized H2
  expect_equal(simplePHENOTYPES:::.realized_h2(ph4, scale = "record"),
               simplePHENOTYPES:::.realized_h2(ph1), tolerance = 1e-12)
  expect_equal(genetic_values(ph4), genetic_values(ph1))
  expect_gt(simplePHENOTYPES:::.realized_h2(ph4),
            simplePHENOTYPES:::.realized_h2(ph4, scale = "record"))
  # h2 itself is the single-record request and still validated on [0, 1]
  expect_equal(ph4$h2, 0.5)
  # reps = 1: both scales coincide
  expect_equal(simplePHENOTYPES:::.realized_h2(ph1, scale = "record"),
               simplePHENOTYPES:::.realized_h2(ph1))
})

test_that("item 6: per-trait reps and vqtl heterogeneity are rescaled", {
  geno <- .fx()
  suppressWarnings({
    mk <- function(reps) {
      s <- simulate_phenotype(geno, n_traits = 2, seed = 9, reps = reps)
      s <- additive(s, prop = 0.4, n_qtn = 3)
      vqtl(s, prop = 0.2, n_qtn = 2)
    }
    base <- mk(1); mixed <- mk(c(1, 4))
  })
  y1 <- split(base$pheno$value, base$pheno$trait)
  ym <- split(mixed$pheno$value, mixed$pheno$trait)
  g <- genetic_values(base)
  # trait 1 untouched; trait 2: genetic part kept, residual (incl. vqtl) / sqrt(4)
  expect_identical(ym$Trait_1, y1$Trait_1)
  expect_equal(unname(ym$Trait_2), unname(g[, 2] + (y1$Trait_2 - g[, 2]) / 2),
               tolerance = 1e-13)
  # residual = exact-variance noise (0.4) + vqtl component (0.2): variance is
  # (0.4 + 0.2 + 2 cov) / 4, so ~0.15 up to the sample covariance of the two parts
  expect_equal(stats::var(ym$Trait_2 - g[, 2]), (1 - 0.4) / 4, tolerance = 0.05)
  expect_identical(mixed$reps, c(1L, 4L))
  # the printed note names both scales
  out <- utils::capture.output(print(mixed))
  expect_true(any(grepl("Entry means of reps = \\[1, 4\\]", out)))
  expect_true(any(grepl("entry-mean", out)))
  expect_true(any(grepl("Single-record realized", out)))
})

test_that("item 6: print() states which scale is shown", {
  geno <- .fx()
  ph <- simulate_phenotype(geno, h2 = 0.5, n_qtn = 3, seed = 3, reps = 4)
  out <- utils::capture.output(print(ph))
  expect_true(any(grepl("Entry means of reps = 4 records", out)))
  expect_true(any(grepl("entry-mean H", out)))
  expect_true(any(grepl("Single-record realized", out)))
})

test_that("item 6: invalid reps is rejected", {
  geno <- .fx()
  bad <- list(0, -1, 1.5, NA, NA_real_, Inf, "a", TRUE, c(1, 2, 3), numeric(0),
              NULL, 2^31)
  for (b in bad) {
    expect_error(simulate_phenotype(geno, reps = b), "`reps` must be",
                 info = deparse(b))
  }
  expect_error(simulate_phenotype(geno, n_traits = 2, reps = c(2, 3, 4)),
               "`reps` must be")
  expect_error(simulate_phenotype(geno, n_traits = 2, reps = c(2, 0)),
               "`reps` must be")
  # valid: scalar, whole-valued double, per-trait vector
  expect_identical(simulate_phenotype(geno, reps = 3)$reps, 3L)
  expect_identical(simulate_phenotype(geno, n_traits = 3, reps = 2)$reps, c(2L, 2L, 2L))
})

test_that("item 6: complex_phenotypes(reps =) scales the common residual", {
  geno <- .fx()
  suppressWarnings({
    s1 <- additive(simulate_phenotype(geno, n_traits = 2, seed = 21),
                   prop = 0.3, n_qtn = 3)
    s2 <- additive(simulate_phenotype(geno, n_traits = 2, seed = 22),
                   prop = 0.3, n_qtn = 3)
  })
  suppressWarnings({    # inputs built with different seeds: first one is used
    c1 <- complex_phenotypes(s1, s2, h2 = 0.5)
    c1b <- complex_phenotypes(s1, s2, h2 = 0.5, reps = 1)
    c4 <- complex_phenotypes(s1, s2, h2 = 0.5, reps = 4)
  })
  expect_identical(c1b$pheno, c1$pheno)
  g <- c4$complex_genetic[, , 1]
  y <- split(c4$pheno$value, c4$pheno$trait)
  for (t in 1:2) {
    expect_equal(stats::var(y[[t]] - g[, t]), 0.5 / 4, tolerance = 1e-12)
  }
  y1 <- split(c1$pheno$value, c1$pheno$trait)
  expect_equal(y[[1]], g[, 1] + (y1[[1]] - g[, 1]) / 2, tolerance = 1e-13)
  expect_equal(simplePHENOTYPES:::.realized_h2(c4, scale = "record"),
               simplePHENOTYPES:::.realized_h2(c1), tolerance = 1e-12)
  expect_true(any(grepl("Entry means of reps = 4",
                        utils::capture.output(print(c4)))))
  expect_error(suppressWarnings(complex_phenotypes(s1, s2, h2 = 0.5, reps = 0)),
               "`reps` must be")
})
