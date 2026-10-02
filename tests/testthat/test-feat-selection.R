# test-feat-selection.R
#
# Round-5 feature requests (breedingDesigner SPEC-0020 items 7 and 8):
#   * sample_parents() records the source of every drawn parent (attribute
#     "source"); and
#   * select_ind(method = "within_family") accepts a per-family count,
#     `n_per_family`, for unequal families.

data("SNP55K_maize282_maf04")

.fs_f2 <- function(n = 40) {
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:20)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  suppressMessages(selfcross(f1, n = n, seed = 2))
}
.fs_ph <- function(f2) {
  suppressMessages(simulate_phenotype(f2, h2 = 0.5, seed = 3) |>
                     additive(n_qtn = 20))
}

# ---- item 7: sample_parents() source attribute ------------------------------

test_that("sample_parents() records the source id/slot of every drawn parent", {
  f2 <- .fs_f2(40)
  oc <- optimum_contribution(.fs_ph(f2), merit = "bv", lambda = 5)
  for (meth in c("allocate", "multinomial")) {
    m <- sample_parents(oc, f2, n = 25, seed = 4, method = meth)
    src <- attr(m, "source")
    expect_s3_class(src, "data.frame")
    expect_identical(names(src), c("slot", "id", "index", "name"))
    expect_identical(src$slot, 1:25)
    expect_identical(nrow(src), n_individuals(m))
    expect_identical(src$name, m$ids)                      # unique ids in the result
    expect_identical(src$id, f2$ids[src$index])            # source id <-> position
    expect_identical(src$index, match(src$id, f2$ids))
    # the slot carries exactly that individual's genotype
    d <- dosages(m); d0 <- dosages(f2)
    for (j in seq_len(nrow(src))) {
      expect_identical(unname(d[, j]), unname(d0[, src$index[j]]))
    }
    # repeated parents share the source id but have distinct names
    expect_true(anyDuplicated(src$id) > 0L)
    expect_false(anyDuplicated(src$name) > 0L)
  }
})

test_that("sample_parents() draws are unchanged by recording the source", {
  f2 <- .fs_f2(40)
  oc <- optimum_contribution(.fs_ph(f2), merit = "bv", lambda = 5)
  for (meth in c("allocate", "multinomial")) {
    set.seed(123); before <- runif(1)
    set.seed(123)
    m <- sample_parents(oc, f2, n = 18, seed = 7, method = meth)
    expect_identical(runif(1), before)                     # caller's RNG untouched
    # same seed -> same parents; source agrees with the key-based record
    m2 <- sample_parents(oc, f2, n = 18, seed = 7, method = meth)
    expect_identical(m$ids, m2$ids)
    expect_identical(attr(m, "source"), attr(m2, "source"))
    expect_identical(m$keys, f2$keys[attr(m, "source")$index])
  }
})

test_that("sample_parents() source with one holder of all the contribution", {
  f2 <- .fs_f2(10)
  contr <- stats::setNames(rep(0, 10), f2$ids); contr[3] <- 1
  oc <- structure(list(contributions = contr, parents = f2$ids[3], n_parents = 1L,
                       merit = 0, coancestry = 0, lambda = 0), class = "ocs")
  m <- sample_parents(oc, f2, n = 3, seed = 1)
  src <- attr(m, "source")
  expect_identical(src$index, c(3L, 3L, 3L))
  expect_identical(src$id, rep(f2$ids[3], 3))
  expect_identical(src$name, m$ids)
})

# ---- item 8: select_ind(within_family, n_per_family) ------------------------

.fs_fam <- function(sizes) {
  labs <- rep(names(sizes), sizes)
  n <- sum(sizes)
  list(f2 = .fs_f2(n), fam = labs)
}

test_that("n_per_family: one count keeps that many from every family", {
  x <- .fs_fam(c(A = 10, B = 14, C = 16))
  ph <- .fs_ph(x$f2)
  s <- select_ind(ph, method = "within_family", family = x$fam, n_per_family = 3)
  expect_equal(n_individuals(s), 9L)
  sel <- attr(s, "selected")
  fam_of <- stats::setNames(x$fam, ph$ids)
  expect_equal(unname(table(fam_of[sel])), unname(table(rep(c("A", "B", "C"), 3))),
               ignore_attr = TRUE)
  expect_identical(attr(s, "n_per_family"), c(A = 3L, B = 3L, C = 3L))
  # the kept ones are the family-wise top scorers on the criterion
  crit <- simplePHENOTYPES:::.criterion_values(ph, "pheno", 1L, 1L)
  for (f in c("A", "B", "C")) {
    ix <- which(x$fam == f)
    top <- ph$ids[ix][order(crit[ix], decreasing = TRUE)[1:3]]
    expect_setequal(intersect(sel, ph$ids[ix]), top)
  }
})

test_that("n_per_family: vector named by family, unequal families", {
  x <- .fs_fam(c(A = 5, B = 12, C = 23))
  ph <- .fs_ph(x$f2)
  npf <- c(C = 6L, A = 1L, B = 4L)                   # deliberately not in label order
  s <- select_ind(ph, method = "within_family", family = x$fam, n_per_family = npf)
  expect_equal(n_individuals(s), 11L)
  fam_of <- stats::setNames(x$fam, ph$ids)
  got <- table(fam_of[attr(s, "selected")])
  expect_equal(as.integer(got[c("A", "B", "C")]), c(1L, 4L, 6L))
  expect_identical(attr(s, "n_per_family"), c(A = 1L, B = 4L, C = 6L))
  # unnamed vector = order of sorted family labels
  s2 <- select_ind(ph, method = "within_family", family = x$fam,
                   n_per_family = c(1, 4, 6))
  expect_setequal(attr(s2, "selected"), attr(s, "selected"))
  # zero for a family drops it
  s3 <- select_ind(ph, method = "within_family", family = x$fam,
                   n_per_family = c(A = 0, B = 2, C = 2))
  expect_equal(n_individuals(s3), 4L)
  expect_false(any(fam_of[attr(s3, "selected")] == "A"))
  # direction = "low" picks the lowest per family
  sl <- select_ind(ph, method = "within_family", family = x$fam,
                   n_per_family = 2, direction = "low")
  crit <- simplePHENOTYPES:::.criterion_values(ph, "pheno", 1L, 1L)
  for (f in c("A", "B", "C")) {
    ix <- which(x$fam == f)
    expect_setequal(intersect(attr(sl, "selected"), ph$ids[ix]),
                    ph$ids[ix][order(crit[ix])[1:2]])
  }
})

test_that("n_per_family: S and i follow the standard definitions", {
  x <- .fs_fam(c(A = 8, B = 12, C = 20))
  ph <- .fs_ph(x$f2)
  npf <- c(A = 2, B = 3, C = 5)
  s <- select_ind(ph, method = "within_family", family = x$fam, n_per_family = npf)
  crit <- simplePHENOTYPES:::.criterion_values(ph, "pheno", 1L, 1L)
  sel <- attr(s, "selected")
  S <- mean(crit[match(sel, ph$ids)]) - mean(crit)
  expect_equal(attr(s, "differential"), unname(S))
  expect_equal(attr(s, "intensity"), unname(S / stats::sd(crit)))
  expect_identical(attr(s, "method"), "within_family")
  expect_identical(attr(s, "criterion"), "pheno")
  # all individuals asked for in every family -> S = 0
  all_in <- select_ind(ph, method = "within_family", family = x$fam,
                       n_per_family = c(A = 8, B = 12, C = 20))
  expect_equal(n_individuals(all_in), 40L)
  expect_equal(attr(all_in, "differential"), 0)
})

test_that("n_per_family: invalid inputs are errors", {
  x <- .fs_fam(c(A = 4, B = 8, C = 12))
  ph <- .fs_ph(x$f2)
  sel <- function(...) select_ind(ph, method = "within_family",
                                  family = x$fam, ...)
  # together with n / prop / intensity
  expect_error(sel(n = 6, n_per_family = 2), "replaces")
  expect_error(sel(prop = 0.5, n_per_family = 2), "replaces")
  expect_error(sel(intensity = 1, n_per_family = 2), "replaces")
  # only within_family
  expect_error(select_ind(ph, method = "mass", n_per_family = 2),
               "within_family")
  expect_error(select_ind(ph, method = "culling", culling = 0.5,
                          n_per_family = 2), "within_family")
  expect_error(select_ind(ph, method = "among_family", family = x$fam,
                          n_per_family = 2), "within_family")
  # a family smaller than requested: error naming it, never capped
  expect_error(sel(n_per_family = 5), "A \\(wants 5, has 4\\)")
  expect_error(sel(n_per_family = c(A = 2, B = 9, C = 1)), "B \\(wants 9, has 8\\)")
  # malformed values
  expect_error(sel(n_per_family = 1.5), "whole number")
  expect_error(sel(n_per_family = -1), "whole number")
  expect_error(sel(n_per_family = NA_real_), "whole number")
  expect_error(sel(n_per_family = "2"), "whole number")
  expect_error(sel(n_per_family = 0), "no individual")
  # names not matching the labels / wrong length
  expect_error(sel(n_per_family = c(A = 1, B = 1)), "missing: C")
  expect_error(sel(n_per_family = c(A = 1, B = 1, C = 1, D = 1)), "unknown: D")
  expect_error(sel(n_per_family = c(A = 1, A = 1, B = 1)), "unique")
  expect_error(sel(n_per_family = c(1, 2)), "one value per family")
  # family still required
  expect_error(select_ind(ph, method = "within_family", n_per_family = 2),
               "family")
})

test_that("select_ind() defaults are unchanged without n_per_family", {
  x <- .fs_fam(c(A = 10, B = 14, C = 16))
  ph <- .fs_ph(x$f2)
  s <- select_ind(ph, n = 12, method = "within_family", family = x$fam)
  expect_equal(n_individuals(s), 12L)
  expect_null(attr(s, "n_per_family"))
  # proportional allocation 3 / 4 / 5 (largest remainder), as before
  fam_of <- stats::setNames(x$fam, ph$ids)
  expect_equal(as.integer(table(fam_of[attr(s, "selected")])), c(3L, 4L, 5L))
  # the argument is appended: positional order of the old signature is intact
  nm <- names(formals(select_ind))
  expect_identical(utils::tail(nm, 3), c("n_per_family", "lambda", "min_gain"))
  expect_identical(nm[1:3], c("sim", "n", "prop"))
  # mass selection still needs one of n / prop / intensity
  expect_error(select_ind(ph, method = "within_family", family = x$fam),
               "exactly one")
})

test_that("pedigree()/recurrent_selection() do not take n_per_family", {
  expect_false("n_per_family" %in% names(formals(pedigree)))
  expect_false("n_per_family" %in% names(formals(recurrent_selection)))
})
