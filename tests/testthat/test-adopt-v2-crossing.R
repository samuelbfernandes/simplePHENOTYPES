# test-adopt-v2-crossing.R
#
# Test-gap proposals of the 2026-09 independent audit (groups v2-crossing and
# v2-rust-core) that no earlier test file adopted. Each block names the audit
# proposal it adopts (CROSS-T* / CROSS-C* = Fable / Codex rows of v2-crossing,
# RUST-T* / RUST-C* = Fable / Codex rows of v2-rust-core). All assertions
# describe the behaviour of the current code; none touches the Rust sources.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# run `expr` under a local seed and leave the caller's RNG state untouched
.adopt_with_seed <- function(seed, expr) {
  old <- .Random.seed_safe()
  on.exit(.restore_seed(old), add = TRUE)
  set.seed(seed)
  force(expr)
}

.adopt_geno <- function(n_mk = 30, chr = 1L, cm = seq(0, 100, length.out = n_mk),
                        allele = "A/G") {
  data.frame(snp = paste0("m", seq_len(n_mk)), allele = allele, chr = chr,
             pos = seq_len(n_mk) * 1e6, cm = cm,
             P1 = rep(1L, n_mk), P2 = rep(-1L, n_mk), stringsAsFactors = FALSE)
}

# one-locus breeds of fully inbred lines: breed A fixed for +1, B for -1, C for +1
.adopt_lines <- function(pool, value, n = 12) {
  g <- matrix(value, 1L, n, dimnames = list(NULL, paste0(pool, seq_len(n))))
  as_population(cbind(data.frame(snp = "q", allele = "A/G", chr = 1L, pos = 1L,
                                 cm = 0, stringsAsFactors = FALSE),
                      as.data.frame(g)), pool = pool)
}

# ---------------------------------------------------------------------------
# CROSS-T7: three-way cross heterosis = 1/2 (H_AC + H_BC)
# ---------------------------------------------------------------------------

test_that("CROSS-T7: three-way heterosis is half the sum of the two F1 heterosis values", {
  skip_on_cran()
  m <- 100L
  breed <- function(pool, p, seed) {
    g <- .adopt_with_seed(seed, t(vapply(p, function(pp) stats::rbinom(80, 2, pp) - 1L,
                                         integer(80))))
    colnames(g) <- paste0(pool, seq_len(80))
    as_population(cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                                   chr = rep(1:10, each = m / 10),
                                   pos = rep(seq_len(m / 10), 10),
                                   cm = rep(seq(0, 100, length.out = m / 10), 10),
                                   stringsAsFactors = FALSE),
                        as.data.frame(g)), pool = pool)
  }
  fr <- .adopt_with_seed(1, list(A = stats::runif(m, 0.05, 0.95),
                                 B = stats::runif(m, 0.05, 0.95),
                                 C = stats::runif(m, 0.05, 0.95)))
  br <- list(A = breed("A", fr$A, 2), B = breed("B", fr$B, 3), C = breed("C", fr$C, 4))
  ad <- .adopt_with_seed(5, list(a = stats::rnorm(m, sd = 0.2),
                                 d = abs(stats::rnorm(m, sd = 0.3))))
  # (A x B) dams x C sires: per locus, the progeny heterozygosity is
  # (h_AC + h_BC) / 2 and the composition-weighted mid-parent is
  # h_A / 4 + h_B / 4 + h_C / 2, so H_ABC = (H_AC + H_BC) / 2 for ANY breeds
  expected <- function(h) 0.5 * (h$expected_f1["A", "C"] + h$expected_f1["B", "C"])
  r <- vapply(1:8, function(s) {
    tw <- suppressMessages(crossbreed(br, "three_way", n_progeny = 300, seed = 100 + s))
    h <- heterosis(tw, br, qtn = seq_len(m), a = ad$a, d = ad$d)
    c(realized = h$realized, expected = expected(h))
  }, numeric(2))
  # the expectation is a function of the breeds only (identical in every replicate)
  expect_equal(max(r["expected", ]) - min(r["expected", ]), 0, tolerance = 1e-12)
  se <- stats::sd(r["realized", ]) / sqrt(ncol(r))
  expect_gt(r["expected", 1], 0)
  expect_lt(abs(mean(r["realized", ]) - r["expected", 1]), 4 * se)
  # the composition really is 1/4 A + 1/4 B + 1/2 C
  tw <- suppressMessages(crossbreed(br, "three_way", n_progeny = 20, seed = 1))
  expect_equal(colMeans(breed_composition(tw))[c("A", "B", "C")],
               c(A = 0.25, B = 0.25, C = 0.5))
})

# ---------------------------------------------------------------------------
# CROSS-T8 / CROSS-C10: breed-origin recurrences, exactly
# ---------------------------------------------------------------------------

test_that("CROSS-T8: three-breed rotation obeys C_t = (C_{t-1} + sire)/2 and cycles (4/7, 2/7, 1/7)", {
  br <- list(A = .adopt_lines("A", 1L), B = .adopt_lines("B", -1L),
             C = .adopt_lines("C", 1L))
  rot <- suppressMessages(crossbreed(br, "rotational", n_progeny = 10,
                                     generations = 15, seed = 1))
  h <- attr(rot, "history")
  expect_identical(nrow(h), 16L)
  comp <- as.matrix(h[, c("A", "B", "C")])
  expect_equal(unname(rowSums(comp)), rep(1, 16))
  # exact recurrence: the composition of generation t is the mean of the dams'
  # composition (generation t-1) and the pure sire breed's indicator
  sire <- match(h$sire_breed, c("A", "B", "C"))
  for (t in 2:16) {
    e <- numeric(3); e[sire[t]] <- 1
    expect_equal(unname(comp[t, ]), unname((comp[t - 1, ] + e) / 2), tolerance = 1e-12)
  }
  # the sires cycle through the breeds in order, and the limit cycle is a
  # permutation of (4/7, 2/7, 1/7), reached to within 2^-14
  expect_identical(h$sire_breed[1], "B")                   # the F1 is A x B
  expect_identical(h$sire_breed[2:16], rep(c("A", "B", "C"), 5))
  for (t in 14:16) {
    expect_equal(sort(unname(comp[t, ])), sort(c(4, 2, 1) / 7), tolerance = 1e-3)
  }
})

test_that("CROSS-C10: F1, backcross, F2 and rotation retain 1, 1/2, 1/2 and 2/3 of the F1 heterosis (inbred lines, one locus)", {
  skip_on_cran()
  # fixed inbred lines: h = 0 = 2p(1 - p), so the classical fractions are exact;
  # a = 0, d = 1 makes the genotypic value the heterozygosity indicator and the
  # pure-breed means 0
  A <- .adopt_lines("A", 1L); B <- .adopt_lines("B", -1L)
  br <- list(A = A, B = B)
  hz <- function(pop) heterosis(pop, br, qtn = 1, a = 0, d = 1)
  n <- 1000L
  tol <- 4 * sqrt(0.25 / n)           # 4 binomial SE at p = 1/2

  f1 <- suppressMessages(crossbreed(br, "two_way", n_progeny = n, seed = 1))
  h_f1 <- hz(f1)
  expect_equal(h_f1$expected_f1["A", "B"], 1)
  expect_equal(h_f1$realized, 1)                       # every F1 is heterozygous

  bc <- suppressMessages(crossbreed(br, "backcross", n_progeny = n, seed = 2))
  expect_equal(colMeans(breed_composition(bc))[c("A", "B")], c(A = 0.75, B = 0.25))
  expect_lt(abs(hz(bc)$realized / h_f1$realized - 0.5), tol)

  plan <- mating_design(f1, design = "random", n_crosses = n, seed = 3)
  f2 <- suppressMessages(mate(plan, f1, seed = 4))
  expect_equal(colMeans(breed_composition(f2))[c("A", "B")], c(A = 0.5, B = 0.5))
  expect_lt(abs(hz(f2)$realized / h_f1$realized - 0.5), tol)

  # two-breed rotation: the expected heterozygosity of generation t is the chance
  # that the dam gamete's breed differs from the sire breed, read from the
  # composition history; at equilibrium it is 2/3
  rot <- suppressMessages(crossbreed(br, "rotational", n_progeny = n,
                                     generations = 8, seed = 5))
  hist <- attr(rot, "history")
  last <- nrow(hist)
  a_dam <- hist$A[last - 1L]
  sire_is_a <- identical(hist$sire_breed[last], "A")
  expected_het <- if (sire_is_a) 1 - a_dam else a_dam
  expect_equal(expected_het, 2 / 3, tolerance = 5e-3)
  expect_lt(abs(hz(rot)$realized / h_f1$realized - expected_het), 4 * sqrt(0.25 / n))
})

# ---------------------------------------------------------------------------
# CROSS-T10: families by father, and the pedigree line of print()
# ---------------------------------------------------------------------------

test_that("CROSS-T10: families(by = 'paternal_half_sib') and the pedigree line of print()", {
  pop <- as_population(.adopt_geno(), pool = "A")
  f1 <- cross(pop[1], pop[2], n = 2, seed = 1)
  f2 <- selfcross(f1[1], n = 2, seed = 2)
  pat <- families(c(f1, f2), "paternal_half_sib")
  expect_identical(levels(pat), c("P2", "prog_1"))
  expect_identical(as.character(pat), c("P2", "P2", "prog_1", "prog_1"))
  expect_identical(names(pat), c(f1$ids, f2$ids))
  # founders have no father
  expect_true(all(is.na(families(pop, "paternal_half_sib"))))
  # a selfed individual's father is its parent
  expect_identical(as.character(families(f2, "paternal_half_sib")), c("prog_1", "prog_1"))

  out <- capture.output(print(f2))
  expect_true(any(grepl("Pedigree: 5 recorded individuals; generation 2$", out)))
  out2 <- capture.output(print(c(pop[1], f2)))
  expect_true(any(grepl("Pedigree: [0-9]+ recorded individuals; generation 0-2$", out2)))
})

# ---------------------------------------------------------------------------
# CROSS-T15: synthetic_map, unnamed per-chromosome arguments, first-appearance order
# ---------------------------------------------------------------------------

test_that("CROSS-T15: an unnamed total_cm / centromere follows the order chromosomes first appear", {
  chr <- c(2, 2, 2, 1, 1, 1)
  pos <- rep(c(1, 2, 3), 2) * 1e6
  # chromosome 2 appears first: it gets 30 cM, chromosome 1 gets 60 cM
  expect_equal(synthetic_map(chr, pos, total_cm = c(30, 60), suppression = 0),
               c(0, 15, 30, 0, 30, 60))
  # named values are matched by chromosome label whatever the order
  expect_equal(synthetic_map(chr, pos, total_cm = c(`1` = 60, `2` = 30), suppression = 0),
               c(0, 15, 30, 0, 30, 60))
  # an unnamed centromere is positional too: 2 Mb is the midpoint of chr 2, and
  # chr 1 puts its centromere at 1 Mb, so the two maps differ only on chr 1
  a <- synthetic_map(chr, pos, total_cm = c(30, 60), centromere = c(2e6, 3e6))
  b <- synthetic_map(chr, pos, total_cm = c(30, 60), centromere = c(2e6, 1e6))
  expect_equal(a[1:3], b[1:3])
  expect_false(isTRUE(all.equal(a[4:6], b[4:6])))
  expect_equal(c(a[3], a[6]), c(30, 60))
})

# ---------------------------------------------------------------------------
# CROSS-C8 (supplement): the count draw is Poisson(last position), not Poisson(span)
# ---------------------------------------------------------------------------

test_that("CROSS-C8: a span-based Poisson(0.1) stream is NOT what .draw_meiosis consumes", {
  a <- .adopt_with_seed(2, .draw_meiosis(list(c(0.5, 0.6)), 5)$counts)
  span <- .adopt_with_seed(2, {
    out <- integer(5)
    for (i in 1:5) {
      out[i] <- k <- stats::rpois(1, 0.1)
      if (k > 0) stats::runif(k, 0, 0.1)
      stats::rbinom(1, 1, 0.5)
    }
    out
  })
  expect_false(identical(a, span))
})

# ---------------------------------------------------------------------------
# RUST-C5 (supplement): flip 0 / 1 with no breakpoints give an all-0 / all-1 mask
# ---------------------------------------------------------------------------

test_that("RUST-C5: flips 0 and 1 without breakpoints give the all-zero and all-one masks", {
  pos <- c(0, 0.5, 1, 1.5, 2)
  expect_identical(gamete_masks_core(5L, pos, numeric(), 0L, 0L), "00000")
  expect_identical(gamete_masks_core(5L, pos, numeric(), 0L, 1L), "11111")
  expect_identical(gamete_masks_core(5L, pos, numeric(), c(0L, 0L), c(0L, 1L)),
                   c("00000", "11111"))
})

# ---------------------------------------------------------------------------
# RUST-T7: the MAF = 0.5 tie rule, pinned
# ---------------------------------------------------------------------------

test_that("RUST-T7: compute_flip() swaps the alleles only when allele 2 is strictly more frequent", {
  cf <- function(...) simplePHENOTYPES:::compute_flip(matrix(c(...), nrow = 1L))[[1L]]
  expect_false(cf(0L, 2L))                    # 1 vs 1: tie, allele 1 stays the reference
  expect_false(cf(2L, 0L))
  expect_false(cf(0L, 2L, 1L, 1L, NA))        # hets and missing calls do not break the tie
  expect_false(cf(0L, 0L, 2L, 2L))            # 2 vs 2
  expect_true(cf(0L, 2L, 2L))                 # 1 vs 2: allele 2 strictly more frequent
  expect_false(cf(0L, 0L, 2L))                # 2 vs 1
  # a het-only (or all-missing) marker: no homozygote class, no flip
  expect_false(cf(1L, 1L, 1L))
  expect_false(cf(NA_integer_, NA_integer_))
})

# ---------------------------------------------------------------------------
# RUST-T10: Haldane at the mask level, every marker pair, absolute bound
# ---------------------------------------------------------------------------

test_that("RUST-T10: gamete masks follow Haldane for every distance on a 41-marker, 2 M chromosome", {
  n <- 20000L
  pos <- seq(0, 2, by = 0.05)
  dr <- .adopt_with_seed(11, .draw_meiosis(list(pos), n))
  masks <- gamete_masks_core(length(pos), pos, dr$chiasmata, dr$counts, dr$flips)
  expect_length(masks, n)
  bits <- do.call(rbind, strsplit(masks, "", fixed = TRUE))
  first <- bits[, 1L]
  d <- pos[-1L]
  r <- (1 - exp(-2 * d)) / 2
  r_hat <- vapply(seq_along(d), function(j) mean(bits[, j + 1L] != first), numeric(1))
  se <- sqrt(r * (1 - r) / n)
  expect_lt(max(abs(r_hat - r) / se), 4)
  # the strand an end of the chromosome starts on is a fair coin; the number of
  # exchanges is Poisson with mean L = 2 (the last map position, in Morgans)
  expect_lt(abs(mean(dr$flips) - 0.5), 4 * sqrt(0.25 / n))
  expect_lt(abs(mean(dr$counts) - 2), 4 * sqrt(2 / n))
  expect_lt(abs(stats::var(dr$counts) / 2 - 1), 0.08)   # Poisson: var = mean
})

# ---------------------------------------------------------------------------
# RUST-T11: a zero-length chromosome through the production path
# ---------------------------------------------------------------------------

test_that("RUST-T11: an all-0 cM chromosome is inherited as one block and segregates 1:2:1", {
  g <- .adopt_geno(n_mk = 9, chr = c(rep(1L, 4), rep(2L, 5)),
                   cm = c(rep(0, 4), seq(0, 80, length.out = 5)))
  pop <- as_population(g)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  expect_true(all(dosages(f1) == 0L))
  n <- 2000L
  f2 <- selfcross(f1, n = n, seed = 2)
  d <- dosages(f2)
  # the four markers of chromosome 1 never separate
  expect_true(all(apply(d[1:4, , drop = FALSE], 2L, function(x) length(unique(x)) == 1L)))
  # coupling phase: F1 haplotypes are all-A / all-G, so F2 dosage is -1, 0, 1 as 1:2:1
  se <- function(p) 4 * sqrt(p * (1 - p) / n)
  expect_lt(abs(mean(d[1, ] == -1L) - 0.25), se(0.25))
  expect_lt(abs(mean(d[1, ] == 0L) - 0.50), se(0.50))
  expect_lt(abs(mean(d[1, ] == 1L) - 0.25), se(0.25))
  # the normal chromosome still recombines: its end markers are not always equal
  expect_gt(mean(d[5, ] != d[9, ]), 0.1)
  # doubled haploids of the F1: one fair coin per gamete for the whole chromosome
  dh <- dosages(double_haploid(f1, n = n, seed = 3))
  expect_true(all(dh[1:4, ] %in% c(-1L, 1L)))
  expect_true(all(apply(dh[1:4, , drop = FALSE], 2L, function(x) length(unique(x)) == 1L)))
  expect_lt(abs(mean(dh[1, ] == 1L) - 0.5), 4 * sqrt(0.25 / n))
})

# ---------------------------------------------------------------------------
# RUST-C9: as_numeric() dispatch and the all-missing marker
# ---------------------------------------------------------------------------

test_that("RUST-C9: a character matrix is data, never a path", {
  m <- matrix(c("AA", "AG", "GG", "AA", "AA", "AA", "GG", "GG"), 2L, 4L, byrow = TRUE,
              dimnames = list(c("m1", "m2"), paste0("s", 1:4)))
  # a bare string is a path (a missing file is an error about the file) ...
  expect_error(as_numeric(file.path(tempdir(), "no-such-file.vcf")), "exist|not found|No such")
  # ... a character matrix is data: it is never opened as a file, and with no
  # recognisable format the error is about the format, not a file
  expect_error(suppressMessages(as_numeric(m, to_r = TRUE, verbose = FALSE)),
               "format was not detected")
  expect_false(any(grepl("exist|No such", tryCatch(as_numeric(m), error = conditionMessage))))
})

test_that("RUST-C9: an all-missing HapMap marker follows the documented imputation and coding", {
  skip_if_not_installed("data.table")
  calls <- rbind(c("NN", "NN", "NN", "NN"), c("AA", "AA", "GG", "GG"))
  hm <- data.frame(`rs#` = c("m1", "m2"), alleles = "A/G", chrom = 1L, pos = c(100L, 200L),
                   strand = "+", `assembly#` = NA, center = NA, protLSID = NA,
                   assayLSID = NA, panelLSID = NA, QCcode = NA, check.names = FALSE,
                   stringsAsFactors = FALSE)
  for (j in 1:4) hm[[paste0("s", j)]] <- calls[, j]
  conv <- function(...) suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE, ...))
  g <- function(res) as.matrix(res[, 6:9])
  # the second marker is a MAF = 0.5 tie: allele 1 (A) is the reference
  tie <- c(1, 1, -1, -1)
  for (code_as in c("-101", "012")) {
    v <- if (code_as == "012") c(maj = 2, het = 1, min = 0) else c(maj = 1, het = 0, min = -1)
    for (impute in c("None", "Middle", "Minor", "Major")) {
      res <- conv(impute = impute, code_as = code_as)
      want1 <- switch(impute, None = rep(NA_real_, 4), Middle = rep(v[["het"]], 4),
                      Minor = rep(v[["min"]], 4), Major = rep(v[["maj"]], 4))
      expect_equal(unname(as.numeric(g(res)[1, ])), unname(want1),
                   info = paste(code_as, impute))
      expect_equal(unname(as.numeric(g(res)[2, ])),
                   if (code_as == "012") tie + 1 else tie, info = paste(code_as, impute))
    }
  }
})

# ---------------------------------------------------------------------------
# RUST-C10: backend-contract signature manifest (names and order, all exports)
# ---------------------------------------------------------------------------

test_that("RUST-C10: contract functions keep their formal names and order (new formals may only be appended)", {
  # test-backend-contract.R freezes the crossing surface exactly; this covers the
  # rest of the contract. A manifest is a PREFIX of the current formals: renaming,
  # removing or reordering a formal fails, appending an optional one (the package's
  # backward-compatible way to grow an API) does not.
  manifest <- list(
    population_from_haplotypes = c("cis", "trans", "map", "ids", "pool", "individuals_in_rows"),
    haplotypes = "x", dosages = "x", n_individuals = "x",
    simulate_phenotype = c("geno", "architecture", "n_traits", "n_qtn", "n_reps", "vary_qtn",
                           "seed", "h2", "mean", "individuals", "model", "expression",
                           "transcriptome", "reps", "resid_cor", "refit", "..."),
    additive = c("sim", "prop", "n_qtn", "qtn", "effect", "phase", "dist", "orthogonal", "a", "d"),
    dominance = c("sim", "prop", "same_as_add", "n_qtn", "qtn", "dist"),
    epistasis = c("sim", "prop", "n_pairs", "interaction", "interaction_type", "qtn",
                  "effect", "dist"),
    vqtl = c("sim", "prop", "same_as_add", "n_qtn", "qtn", "dist"),
    complex_phenotypes = c("...", "h2", "reps"),
    genetic_values = c("sim", "rep"), qtn_table = c("sim", "rep"),
    phenotypes_long = "sim", phenotypes_wide = "sim",
    write_phenotypes = c("sim", "file", "format", "sep", "file_type"),
    select_ind = c("sim", "n", "prop", "intensity", "on", "trait", "direction", "method",
                   "family", "weights", "quad_weights", "h2", "family_relationship", "rep",
                   "culling", "sequential", "n_per_family"),
    single_seed_descent = c("x", "generations", "seed", "interference"),
    bulk = c("x", "generations", "n", "seed", "interference"),
    pedigree = c("x", "phenotype", "generations", "prop", "n_select", "pop_size", "on",
                 "trait", "direction", "seed", "interference"),
    recurrent_selection = c("x", "phenotype", "cycles", "n_parents", "n_crosses",
                            "progeny_per_cross", "on", "trait", "direction", "seed",
                            "interference"),
    g_matrix = c("x", "ridge", "base_freq"),
    optimum_contribution = c("x", "merit", "trait", "direction", "lambda", "target_coancestry",
                             "max_coancestry", "G", "rep", "min_contribution", "max_iter", "tol"),
    sample_parents = c("ocs", "pop", "n", "seed", "method"),
    cross_usefulness = c("sim", "pairs", "scheme", "n_progeny", "generations", "select_top",
                         "trait", "direction", "seed", "interference"),
    mabc_select = c("pop", "recurrent", "donor", "target_markers", "target_requirement",
                    "background_markers", "exclude_interval", "flanking_markers", "n",
                    "marker_weights", "seed"),
    recurrent_parent_recovery = c("pop", "recurrent", "donor", "markers", "weights"),
    parentage = c("x", "ancestors"), families = c("x", "by"),
    combining_ability = c("candidates", "testers", "qtn", "a", "d", "design", "method",
                          "n_progeny", "h2", "var_e", "ref", "seed", "interference"),
    template_effects = c("sim", "trait", "rep"),
    progeny_test = c("parents", "mates", "qtn", "a", "d", "n_progeny", "h2", "var_e", "ref",
                     "seed", "interference"),
    marker_select = c("pop", "markers", "favorable", "requirement", "min_markers", "n", "prop",
                      "rank_on", "direction", "seed"),
    predict_ebv = c("x", "pheno", "method", "h2", "var_a", "var_e", "ref", "K", "base_freq",
                    "ridge"),
    a_matrix = c("pop", "ids", "founder_f"), prediction_accuracy = c("ebv", "truth"),
    selection_methods = character(),
    breed_composition = "pop",
    additive_value = c("x", "qtn", "effect"), genotypic_value = c("x", "qtn", "a", "d"),
    phenotype_value = c("x", "qtn", "effect", "h2", "var_e", "ref", "seed", "d"),
    filter_geno = c("geno", "maf_above", "maf_below", "hets", "remove_monomorphic",
                    "indep_pairwise", "indep_pairphase", "indep", "blocks", "block_max_kb",
                    "window_unit", "code_as", "verbose")
  )
  ns <- asNamespace("simplePHENOTYPES")
  for (f in names(manifest)) {
    expect_true(f %in% getNamespaceExports("simplePHENOTYPES"), info = f)
    now <- names(formals(get(f, envir = ns)))
    want <- manifest[[f]]
    if (!length(want)) {
      expect_null(now, info = paste("formals of", f))      # takes no arguments
    } else {
      expect_identical(now[seq_along(want)], want, info = paste("formals of", f))
    }
  }
})

# ---------------------------------------------------------------------------
# RUST-C11: one small cross-format fixture, exact expected matrix, never skipped on CI
# ---------------------------------------------------------------------------

test_that("RUST-C11: HapMap, VCF, GDS, BED and PED of one fixture give the same hand-computed matrix", {
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("data.table")
  # six samples, four markers, no MAF = 0.5 ties and no missing calls, so the
  # major allele is unambiguous at every marker
  calls <- rbind(m1 = c("AA", "AA", "AA", "AG", "GG", "AA"),   # A major
                 m2 = c("GG", "GG", "AG", "AA", "GG", "GG"),   # G major
                 m3 = c("AA", "AG", "GG", "GG", "GG", "AG"),   # G major
                 m4 = c("AA", "AA", "GG", "AA", "AA", "AA"))   # A major
  samples <- paste0("s", 1:6)
  want <- rbind(m1 = c(1, 1, 1, 0, -1, 1),
                m2 = c(1, 1, 0, -1, 1, 1),
                m3 = c(-1, 0, 1, 1, 1, 0),
                m4 = c(1, 1, -1, 1, 1, 1))
  dir <- tempfile("fixture"); dir.create(dir)

  # HapMap file
  hm <- data.frame(`rs#` = rownames(calls), alleles = "A/G", chrom = 1L,
                   pos = c(100L, 200L, 300L, 400L), strand = "+", `assembly#` = NA,
                   center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA,
                   QCcode = NA, check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_along(samples)) hm[[samples[j]]] <- calls[, j]
  hmp <- file.path(dir, "fixture.hmp.txt")
  utils::write.table(hm, hmp, sep = "\t", quote = FALSE, row.names = FALSE)

  # the same genotypes as VCF (REF = A, ALT = G), then GDS, BED and PED from it
  gt <- ifelse(calls == "AA", "0/0", ifelse(calls == "GG", "1/1", "0/1"))
  vcf <- file.path(dir, "fixture.vcf")
  writeLines(c("##fileformat=VCFv4.2",
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                       "FORMAT", samples), collapse = "\t"),
               vapply(seq_len(nrow(calls)), function(i)
                 paste(c("1", 100L * i, rownames(calls)[i], "A", "G", ".", "PASS", ".",
                         "GT", gt[i, ]), collapse = "\t"), character(1))), vcf)
  gds <- file.path(dir, "fixture.gds")
  SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "copy.num.of.ref", verbose = FALSE)
  gf <- SNPRelate::snpgdsOpen(gds)
  bed_base <- file.path(dir, "fixture_bed")
  ped_base <- file.path(dir, "fixture_ped")
  SNPRelate::snpgdsGDS2BED(gf, bed_base, verbose = FALSE)
  SNPRelate::snpgdsGDS2PED(gf, ped_base, verbose = FALSE)
  SNPRelate::snpgdsClose(gf)

  conv <- function(x) suppressMessages(as_numeric(x, to_r = TRUE, verbose = FALSE))
  mat <- function(res) unname(as.matrix(res[, -(1:5), drop = FALSE]))
  for (src in list(hmp, vcf, gds, paste0(bed_base, ".bed"), paste0(ped_base, ".ped"))) {
    res <- conv(src)
    expect_identical(nrow(res), 4L, info = basename(src))
    expect_equal(mat(res), unname(want), ignore_attr = TRUE, info = basename(src))
    # BED/PED (through SNPRelate) carry running numbers, not the marker names
    if (!grepl("\\.(bed|ped)$", src)) {
      expect_identical(res$snp, rownames(calls), info = basename(src))
    }
    expect_identical(as.integer(res$pos), c(100L, 200L, 300L, 400L), info = basename(src))
  }
})
