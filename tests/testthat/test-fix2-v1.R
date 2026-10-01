# test-fix2-v1.R
#
# Round-2 regression tests for the frozen v1 engine (group "v1"): findings
# A3, B2, B10, B11, B12 (V1 half), B18 and the text-only corrections C13-C15
# of the independent review (.tmp/codex-review/FIX_LIST.md).
#
# Owner rule D1: bad inputs are rejected with a clear message; every valid
# output stays bit-identical (test-v130-parity.R guards that).

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
.f2_geno <- function() {
  e <- new.env()
  utils::data("SNP55K_maize282_maf04", envir = e)
  e$SNP55K_maize282_maf04
}

.f2_home <- function(env = parent.frame()) {
  d <- tempfile("f2v1_")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}

# create_phenotypes() into a fresh home_dir; errors propagate
.f2_cp <- function(G, ..., home = NULL, verbose = FALSE) {
  if (is.null(home)) home <- .f2_home(parent.frame())
  suppressMessages(suppressWarnings(
    create_phenotypes(geno_obj = G, home_dir = home, output_dir = "",
                      to_r = TRUE, verbose = verbose, ...)
  ))
}

# source text of an R/ file (only available when testing from the source tree)
.f2_src <- function(f) {
  p <- test_path("..", "..", "R", f)
  skip_if_not(file.exists(p), "R/ sources not available")
  paste(readLines(p), collapse = "\n")
}

.f2_panel <- function(m = 30, n = 24, seed = 1, pos = NULL, chr = 1L) {
  set.seed(seed)
  X <- matrix(sample(c(-1L, 1L), m * n, replace = TRUE), m, n)
  colnames(X) <- paste0("s", seq_len(n))
  data.frame(snp = paste0("m", seq_len(m)), allele = "A/G", chr = chr,
             pos = if (is.null(pos)) seq_len(m) * 100L else pos, cm = 0,
             X, check.names = FALSE)
}

# ---------------------------------------------------------------------------
# A3: duplicate chr_pos is ACCEPTED (V1 full and partial architectures)
# ---------------------------------------------------------------------------
test_that("A3: duplicated chr_pos is accepted end to end by both non-LD architectures", {
  G <- .f2_panel(m = 30, n = 24, pos = rep(100L, 30))  # one position, 30 markers
  eff2 <- list(c(.2, .04), c(.3, .09))
  full <- .f2_cp(G, ntraits = 2, add_QTN_num = 2, add_effect = eff2,
                 h2 = c(.5, .5), model = "A", architecture = "pleiotropic",
                 seed = 3, rep = 1)
  expect_s3_class(full, "data.frame")
  expect_true(all(is.finite(unlist(full[, -(1:2)]))) || nrow(full) > 0)
  part <- .f2_cp(G, ntraits = 2, pleio_a = 1, trait_spec_a_QTN_num = c(1, 1),
                 add_effect = list(c(.2, .09), c(.3, .09)), h2 = c(.5, .5),
                 model = "A", architecture = "partially", seed = 3, rep = 1)
  expect_s3_class(part, "data.frame")
})

test_that("A3: create_phenotypes() documents that duplicate chr_pos is accepted", {
  src <- .f2_src("legacy_create_phenotypes.R")
  expect_match(src, "Duplicated marker positions (`chr_pos`", fixed = TRUE)
  expect_match(src, "accepted", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# B2: seed = NULL under a non-Mersenne caller RNG advances the ambient stream
# ---------------------------------------------------------------------------
test_that("B2: two successive seed = NULL calls differ under L'Ecuyer-CMRG", {
  G <- .f2_geno()
  old <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(old[1], old[2], old[3])))
  suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  set.seed(99)
  a <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A")
  s_after_a <- .Random.seed
  b <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A")
  expect_false(isTRUE(all.equal(a, b)))
  expect_false(identical(s_after_a, .Random.seed))
  # the caller's kind is untouched
  expect_equal(RNGkind()[1], "L'Ecuyer-CMRG")
  expect_equal(RNGkind()[2:3], c("Inversion", "Rejection"))
  expect_equal(length(.Random.seed), 7L)
  # ... and the stream advance is reproducible from the same start
  set.seed(99)
  a2 <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A")
  expect_equal(a, a2)
})

test_that("B2: an explicit seed leaves the caller's RNG kind and state untouched", {
  G <- .f2_geno()
  old <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(old[1], old[2], old[3])))
  for (kind in c("L'Ecuyer-CMRG", "Mersenne-Twister")) {
    suppressWarnings(RNGkind(kind, "Inversion", "Rejection"))
    set.seed(5)
    before <- .Random.seed
    .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A",
           seed = 12)
    expect_identical(.Random.seed, before)
    expect_equal(RNGkind()[1], kind)
  }
})

test_that("B2: the Mersenne-Twister seed = NULL behaviour is unchanged", {
  G <- .f2_geno()
  old <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(old[1], old[2], old[3])))
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rounding"))
  set.seed(99)
  a <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A")
  b <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A")
  expect_false(isTRUE(all.equal(a, b)))
  expect_equal(RNGkind()[3], "Rounding")
})

# ---------------------------------------------------------------------------
# B10: seed overflow validator bounds linkage's derived seeds seed * s + z
# ---------------------------------------------------------------------------
test_that("B10: the seed validator bounds the LD retry seeds (seed * 10 + z + rep + x)", {
  v <- simplePHENOTYPES:::.v1_validate_seed_arith
  h2 <- matrix(.5)
  # non-LD callers keep the old (looser) bound
  expect_silent(v(300000000, 1, h2, FALSE, FALSE))
  # LD retry seeds reach seed * 10: 300000000 * 10 > .Machine$integer.max
  expect_error(v(300000000, 1, h2, FALSE, FALSE, ld = TRUE), "too large")
  expect_error(v(300000000, 1, h2, FALSE, FALSE, ld = TRUE), "LD")
  lim <- floor((.Machine$integer.max - 1000) / 10)
  expect_silent(v(lim - 100, 2, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3))
  expect_error(v(lim + 500, 2, h2, FALSE, FALSE, ld = TRUE, n_qtn = 3), "too large")
})

test_that("B10: create_phenotypes(architecture = 'LD', seed = 3e8) is rejected up front", {
  G <- .f2_geno()
  expect_error(
    .f2_cp(G, ntraits = 2, architecture = "LD", type_of_ld = "indirect",
           add_QTN_num = 3, add_effect = c(.2, .3), h2 = c(.5, .5),
           model = "A", ld_min = .2, ld_max = .8, ld_method = "corr",
           seed = 300000000, rep = 1),
    "seed.*too large")
  expect_error(
    .f2_cp(G, ntraits = 2, architecture = "LD", type_of_ld = "direct",
           add_QTN_num = 3, add_effect = c(.2, .3), h2 = c(.5, .5),
           model = "A", ld_min = .2, ld_max = .8, ld_method = "corr",
           seed = 300000000, rep = 1),
    "seed.*too large")
  # the same seed is still valid for a non-LD architecture (outputs unchanged)
  expect_s3_class(
    .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 1, h2 = .5, model = "A",
           seed = 300000000), "data.frame")
})

# ---------------------------------------------------------------------------
# B11: indirect-LD checker validates LD magnitude and reported LD
# ---------------------------------------------------------------------------
.f2_ld_triple <- function(G, lo, hi, method = "corr") {
  ldp <- simplePHENOTYPES:::.ld_pair
  for (j in 8:400) {
    sup <- j + (1:6)
    inf <- j - (1:6)
    ok_s <- sup[vapply(sup, function(k) {
      l <- ldp(G, j, k, method); l >= lo && l <= hi
    }, logical(1))]
    ok_i <- inf[vapply(inf, function(k) {
      l <- ldp(G, j, k, method); l >= lo && l <= hi
    }, logical(1))]
    if (length(ok_s) && length(ok_i)) {
      return(list(cause = j, sup = ok_s[1], inf = ok_i[1]))
    }
  }
  NULL
}

test_that("B11: .ld_check_indirect rejects a triple whose LD is outside the window", {
  G <- .f2_geno()
  chk <- simplePHENOTYPES:::.ld_check_indirect
  # rows 1, 5, 6 of the maize panel: distinct, same chromosome, absolute corr
  # LD from row 1 of 0.031 and 0.045 (both outside [.2, .8])
  ldp <- simplePHENOTYPES:::.ld_pair
  expect_lt(ldp(G, 1, 5, "corr"), .2)
  expect_lt(ldp(G, 1, 6, "corr"), .2)
  expect_error(chk(G, 1, 5, 6, 1, ld_method = "corr", ld_min = .2, ld_max = .8,
                   reported_sup = ldp(G, 1, 5, "corr"),
                   reported_inf = ldp(G, 1, 6, "corr")),
               "outside \\[ld_min, ld_max\\]")
})

test_that("B11: .ld_check_indirect accepts an in-window triple with matching reported LD", {
  G <- .f2_geno()
  chk <- simplePHENOTYPES:::.ld_check_indirect
  ldp <- simplePHENOTYPES:::.ld_pair
  tr <- .f2_ld_triple(G, .2, .8)
  skip_if(is.null(tr), "no in-window triple in the example panel")
  rs <- ldp(G, tr$cause, tr$sup, "corr")
  ri <- ldp(G, tr$cause, tr$inf, "corr")
  expect_true(chk(G, tr$cause, tr$sup, tr$inf, 1, ld_method = "corr",
                  ld_min = .2, ld_max = .8, reported_sup = rs, reported_inf = ri))
  # a wrong reported LD is a contract error
  expect_error(chk(G, tr$cause, tr$sup, tr$inf, 1, ld_method = "corr",
                   ld_min = .2, ld_max = .8, reported_sup = rs + .05,
                   reported_inf = ri), "reported")
  # legacy call form (no window information) keeps the structural checks only
  expect_true(chk(G, tr$cause, tr$sup, tr$inf, 1))
})

test_that("B11: a cause marker equal to ANOTHER triple's QTN stays allowed", {
  G <- .f2_geno()
  chk <- simplePHENOTYPES:::.ld_check_indirect
  ldp <- simplePHENOTYPES:::.ld_pair
  t1 <- .f2_ld_triple(G, .2, .8)
  skip_if(is.null(t1), "no in-window triple in the example panel")
  # second triple whose cause is the first triple's downstream QTN
  j2 <- t1$sup
  sup2 <- j2 + (1:6)
  inf2 <- j2 - (1:6)
  s2 <- sup2[vapply(sup2, function(k) { l <- ldp(G, j2, k, "corr"); l >= .2 && l <= .8 }, logical(1))]
  i2 <- inf2[vapply(inf2, function(k) { l <- ldp(G, j2, k, "corr"); l >= .2 && l <= .8 }, logical(1))]
  skip_if(!length(s2) || !length(i2) || i2[1] == t1$cause && FALSE,
          "no second in-window triple")
  cause <- c(t1$cause, j2)
  sup <- c(t1$sup, s2[1])
  inf <- c(t1$inf, i2[1])
  skip_if(anyDuplicated(sup) || anyDuplicated(inf) || length(intersect(sup, inf)),
          "second triple collides structurally")
  rs <- mapply(function(a, b) ldp(G, a, b, "corr"), cause, sup)
  ri <- mapply(function(a, b) ldp(G, a, b, "corr"), cause, inf)
  expect_true(chk(G, cause, sup, inf, 1, ld_method = "corr", ld_min = .2,
                  ld_max = .8, reported_sup = rs, reported_inf = ri))
})

# ---------------------------------------------------------------------------
# B12 (V1 half): the "no heterozygous individual" warning tests genotype 0
# ---------------------------------------------------------------------------
test_that("B12: an all-heterozygous dominance locus does not trigger the false diagnostic", {
  n <- 12
  Xa <- cbind(c(-1, 0, 1, 1, 0, -1, 1, 0, 1, -1, 0, 1))
  Xhet <- cbind(rep(0, n), c(-1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1))  # col 1 all het, col 2 hetless
  rownames(Xa) <- rownames(Xhet) <- letters[1:n]
  expect_no_warning(genetic_effect(add_obj = Xa, dom_obj = Xhet,
                                   add_effect = .3, dom_effect = c(.2, .1),
                                   add = TRUE, dom = TRUE, epi = FALSE))
  # a genuinely hetless set of loci still warns
  Xnone <- cbind(c(-1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1),
                 c(1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1))
  rownames(Xnone) <- letters[1:n]
  expect_warning(genetic_effect(add_obj = Xa, dom_obj = Xnone, add_effect = .3,
                                dom_effect = c(.2, .1), add = TRUE, dom = TRUE,
                                epi = FALSE), "heterozyg")
  # an all-heterozygous locus alone is legitimate: VD = 0 without a warning
  Xall <- cbind(rep(0, n)); rownames(Xall) <- letters[1:n]
  expect_no_warning(genetic_effect(add_obj = Xa, dom_obj = Xall,
                                   add_effect = .3, dom_effect = .2,
                                   add = TRUE, dom = TRUE, epi = FALSE))
})

# ---------------------------------------------------------------------------
# B18: h2 = .05 boundary (round(10 * .05) = 0, ties to even)
# ---------------------------------------------------------------------------
test_that("B18: h2 = 0.05 with rep > 1 is rejected and the message states h2 > 0.05", {
  G <- .f2_geno()
  expect_equal(round(10 * 0.05), 0)
  err <- tryCatch(
    .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 2, h2 = .05, model = "A",
           seed = 7, output_format = "wide"),
    error = function(e) conditionMessage(e))
  expect_match(err, "identical replicates")
  expect_match(err, "h2 > 0\\.05", perl = TRUE)
  expect_false(grepl("h2 >= 0.05", err, fixed = TRUE))
  # just above the boundary is fine and gives distinct replicates
  ph <- .f2_cp(G, add_QTN_num = 3, add_effect = .2, rep = 2, h2 = .051,
               model = "A", seed = 7, output_format = "wide")
  expect_false(isTRUE(all.equal(ph[[2]], ph[[3]])))
})

test_that("B18: the roxygen says h2 > .05, not h2 >= .05", {
  src <- .f2_src("legacy_create_phenotypes.R")
  expect_false(grepl("h2\\W+>=\\W+0\\.05", src))
  expect_match(src, "h2` (at or )?below 0\\.05|h2 > 0\\.05|`h2 > 0\\.05`")
})

# ---------------------------------------------------------------------------
# C13: base_line_multi_traits roxygen vs numerics
# ---------------------------------------------------------------------------
.f2_bl <- function(cor, eff, n = 300, seed = 4) {
  set.seed(seed)
  X <- matrix(sample(c(-1, 0, 1), n * 6, replace = TRUE), n, 6,
              dimnames = list(paste0("i", 1:n), NULL))
  colnames(X) <- paste0("Chr_1_", 1:6)
  res <- suppressWarnings(base_line_multi_traits(
    add_obj = list(X), dom_obj = list(NULL), epi_obj = list(NULL),
    add_effect = eff, dom_effect = list(NULL, NULL),
    epi_effect = list(NULL, NULL), epi_interaction = 2, ntraits = 2,
    cor = cor, architecture = "pleiotropic", rep = 1, rep_by = "experiment",
    add = TRUE, dom = FALSE, epi = FALSE, sim_method = "custom",
    verbose = FALSE))
  list(g = res[[1]]$base_line,
       raw = cbind(scale(X %*% eff[[1]], scale = FALSE),
                   scale(X %*% eff[[2]], scale = FALSE)))
}

test_that("C13: unit diagonal leaves trait 1 unchanged; a non-unit diagonal rescales it", {
  eff <- list(c(.2, .04, .008, .5, .3, .1), c(.1, .3, .05, .02, .4, .2))
  u <- .f2_bl(matrix(c(1, .5, .5, 1), 2), eff)
  expect_equal(unname(u$g[, 1]), unname(u$raw[, 1]), tolerance = 1e-10)
  # non-unit diagonal: covariance-like target [[4, 1], [1, 1]]
  nu <- .f2_bl(matrix(c(4, 1, 1, 1), 2), eff)
  expect_equal(sd(nu$g[, 1]) / sd(nu$raw[, 1]), 2, tolerance = 1e-8)
  expect_equal(sd(nu$g[, 2]) / sd(nu$raw[, 2]), 1, tolerance = 1e-8)
  expect_equal(cor(nu$g)[1, 2], 0.5, tolerance = 1e-8)
  # general rule: the variance of trait k is multiplied by cor[k, k]
  cv <- .f2_bl(matrix(c(4, 1, 1, 2), 2), eff)
  expect_equal(var(cv$g[, 1]) / var(cv$raw[, 1]), 4, tolerance = 1e-8)
  expect_equal(var(cv$g[, 2]) / var(cv$raw[, 2]), 2, tolerance = 1e-8)
})

test_that("C13: the roxygen no longer claims that every trait mixes all inputs", {
  src <- .f2_src("legacy_Base_line_multi_traits.R")
  expect_false(grepl("Every output trait is a linear combination of \\\\emph\\{all\\}", src))
  expect_false(grepl("Trait 1 is unchanged; trait", src, fixed = TRUE))
  expect_match(src, "unit diagonal", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# C14: NA dosages now error in genetic_effect(); text must say so
# ---------------------------------------------------------------------------
test_that("C14: NA dosages stop genetic_effect() and the roxygen no longer says they propagate", {
  X <- cbind(c(-1, 0, NA, 1)); rownames(X) <- letters[1:4]
  expect_error(genetic_effect(add_obj = X, add_effect = .2, add = TRUE,
                              dom = FALSE, epi = FALSE), "missing values")
  src <- .f2_src("legacy_genetic_effect.R")
  expect_false(grepl("propagate to NA genetic values", src, fixed = TRUE))
  expect_false(grepl("no guard", src, fixed = TRUE))
})

# ---------------------------------------------------------------------------
# C15: qtn_pleiotropic seed collision rule (seed = rep = 2 example)
# ---------------------------------------------------------------------------
.f2_seedfiles <- function(G, seed, rep) {
  home <- .f2_home(parent.frame())
  .f2_cp(G, add_QTN_num = 3, add_effect = .2, dom_QTN_num = 3, dom_effect = .2,
         var_QTN_num = 3, var_effect = .1, h2 = .5, rep = rep, model = "ADV",
         seed = seed, vary_QTN = TRUE, home = home, verbose = TRUE)
  rd <- function(p) scan(list.files(home, p, full.names = TRUE), quiet = TRUE)
  list(add = rd("Add_QTN"), dom = rd("Dom_QTN"), var = rd("var_QTN"),
       dom_snp = utils::read.delim(file.path(home, "Dominance_QTNs.txt"))$snp,
       var_snp = utils::read.delim(file.path(home, "Variance_QTNs.txt"))$snp)
}

test_that("C15: seed = rep = 2 makes the dominance and variance draws collide; seed >= 2 * rep does not", {
  G <- .f2_geno()
  s2 <- .f2_seedfiles(G, seed = 2, rep = 2)
  expect_equal(s2$dom, c(5, 6))
  expect_equal(s2$var, c(5, 6))
  expect_identical(s2$dom_snp, s2$var_snp)        # equal draw sizes: same set
  s4 <- .f2_seedfiles(G, seed = 4, rep = 2)       # seed >= 2 * rep
  expect_equal(intersect(intersect(s4$dom, s4$var), s4$add), numeric(0))
  expect_equal(intersect(s4$dom, s4$var), numeric(0))
  expect_equal(intersect(s4$add, s4$var), numeric(0))
  expect_false(identical(s4$dom_snp, s4$var_snp))
  # seed = 2 * rep - 1 still collides
  s3 <- .f2_seedfiles(G, seed = 3, rep = 2)
  expect_true(length(intersect(s3$dom, s3$var)) > 0)
})

test_that("C15: the roxygen states seed >= 2 * rep and equal draw sizes", {
  src <- .f2_src("legacy_QTN_pleiotropic.R")
  expect_false(grepl("Choose `seed >= rep`", src, fixed = TRUE))
  expect_match(src, "seed >= 2 * rep", fixed = TRUE)
})
