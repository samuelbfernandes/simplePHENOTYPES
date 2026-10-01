# test-audit-rust-crossing.R
#
# Regression tests for the 2026-09 audit of the Rust core and the crossing
# layer (reconciliation v2-rust-core / v2-crossing).
#
# Anything that used to PANIC in the Rust kernel aborts the whole R process on
# some toolchains (gcc-linked macOS builds), so those calls are made in a
# subprocess: the test then also proves that the process survives (exit status
# 0) and that the failure is an ordinary R error.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

.probe_preamble <- function() {
  if (requireNamespace("pkgload", quietly = TRUE) &&
      getExportedValue("pkgload", "is_dev_package")("simplePHENOTYPES")) {
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE, compile = FALSE))",
            deparse(as.character(getExportedValue("pkgload", "pkg_path")())))
  } else {
    "suppressMessages(library(simplePHENOTYPES))"
  }
}

# Run `code` in a fresh Rscript; returns list(status, out).
.run_probe <- function(code) {
  script <- tempfile(fileext = ".R")
  writeLines(c(.probe_preamble(), code), script)
  on.exit(unlink(script), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--no-save", "--no-restore", shQuote(script)),
    stdout = TRUE, stderr = TRUE,
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
  ))
  list(status = attr(out, "status") %||% 0L, out = as.character(out))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# Evaluate named calls in the package namespace, one line of output each:
#   CASE <name> ERROR <message>   |   CASE <name> OK
.run_cases <- function(cases) {
  body <- vapply(names(cases), function(nm) {
    sprintf(paste0(
      "r <- tryCatch({ eval(quote(%s), asNamespace('simplePHENOTYPES')); 'OK' },",
      " error = function(e) paste('ERROR', conditionMessage(e)))\n",
      "cat('CASE %s', r, '\\n', sep = ' ')"), cases[[nm]], nm)
  }, character(1))
  res <- .run_probe(c(body, "cat('DONE\\n')"))
  lines <- grep("^CASE ", res$out, value = TRUE)
  nm <- sub("^CASE ([^ ]+) .*$", "\\1", lines)
  msg <- sub("^CASE [^ ]+ ", "", lines)
  list(status = res$status, done = any(res$out == "DONE"),
       result = stats::setNames(msg, nm), out = res$out)
}

.geno <- function(n_mk = 30, chr = 1L, cm = seq(0, 100, length.out = n_mk),
                  allele = "A/G") {
  data.frame(snp = paste0("m", seq_len(n_mk)), allele = allele, chr = chr,
             pos = seq_len(n_mk) * 1e6, cm = cm,
             P1 = rep(1L, n_mk), P2 = rep(-1L, n_mk), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# RUST-F1 / RC-01: panics must be R errors, not process aborts
# ---------------------------------------------------------------------------

test_that("the six abort probes are ordinary R errors and the process survives", {
  skip_on_cran()
  cases <- list(
    A1 = 'mate_haplotypes_core(1L,0,"1","1","1","1",numeric(),integer(),integer(),"cross",1L)',
    A2 = 'numericalize_core(0L,1L,2L,FALSE,"-101","Add","None")',
    A3 = 'gamete_masks_core(3L, c(0,.5,1), c(0.1,0.2), 1L, 0L)',
    A4 = 'mate_haplotypes_core(3L, c(0,.1,.2), "111","000","111","000", numeric(0), c(0L,0L), c(0L,0L), "cross", 2L)',
    A5 = 'mate_haplotypes_core(3L, c(0,.1,.2), "11111111","00000000","11111111","00000000", numeric(0), c(0L,0L), c(0L,0L), "cross", 1L)',
    A6 = 'mate_haplotypes_core(3L, c(0,.1,.2), "111","000","111","000", numeric(0), c(0L,0L), c(0L,0L), "cross", -1L)'
  )
  r <- .run_cases(cases)
  expect_identical(r$status, 0L)
  expect_true(r$done)
  expect_length(r$result, 6L)
  expect_true(all(startsWith(r$result, "ERROR")), info = paste(r$result, collapse = "\n"))
})

# ---------------------------------------------------------------------------
# RUST-F2/F3/F6/F7, RC-01..05: the kernel validates what it receives
# ---------------------------------------------------------------------------

test_that("malformed kernel calls raise R errors that name the problem", {
  skip_on_cran()
  mh <- function(strands = c("111", "000", "111", "000"), counts = c(0L, 0L),
                 flips = c(0L, 0L), chi = "numeric(0)", n = 1L, design = '"cross"',
                 pos = "c(0,.1,.2)", loci = "3L") {
    sprintf('mate_haplotypes_core(%s, %s, "%s","%s","%s","%s", %s, %s, %s, %s, %s)',
            loci, pos, strands[1], strands[2], strands[3], strands[4], chi,
            deparse(counts), deparse(flips), design, n)
  }
  cases <- list(
    short_strand   = mh(strands = c("11", "000", "111", "000")),
    long_strand    = mh(strands = c("1111", "000", "111", "000")),
    alphabet_2     = mh(strands = c("121", "000", "111", "000")),
    surplus_events = mh(counts = rep(0L, 4), flips = rep(0L, 4)),
    deficit_events = mh(counts = 0L, flips = 0L),
    dh_deficit     = mh(counts = integer(), flips = integer(), design = '"dh"'),
    negative_n     = mh(n = -1L),
    na_n           = mh(n = "NA_integer_"),
    bad_design     = mh(design = '"selfing"'),
    chi_nan        = mh(counts = c(1L, 0L), chi = "NaN"),
    chi_negative   = mh(counts = c(1L, 0L), chi = "-1"),
    chi_beyond_L   = mh(counts = c(1L, 0L), chi = "5"),
    chi_inf        = mh(counts = c(1L, 0L), chi = "Inf"),
    counts_sum     = 'gamete_masks_core(3L, c(0,.5,1), c(0.1,0.2), 1L, 0L)',
    counts_neg     = 'gamete_masks_core(1L, 0, numeric(), c(-1L, 1L), c(0L, 0L))',
    counts_na      = 'gamete_masks_core(1L, 0, numeric(), c(NA_integer_, 0L), c(0L, 0L))',
    flip_2         = 'gamete_masks_core(1L, 0, numeric(), c(0L,0L), c(2L,-1L))',
    flip_na        = 'gamete_masks_core(1L, 0, numeric(), c(0L,0L), c(NA_integer_,0L))',
    counts_ragged  = 'gamete_masks_core(c(1L,1L), c(0,0), numeric(), c(0L,0L,0L), c(0L,0L,0L))',
    unsorted_pos   = 'gamete_masks_core(3L, c(0,.5,.2), numeric(), 1L, 0L)',
    nonfinite_pos  = 'gamete_masks_core(2L, c(0,NaN), numeric(), 1L, 0L)',
    empty_chr      = 'gamete_masks_core(c(2L,0L), c(0,1), numeric(), c(0L,0L), c(0L,0L))',
    nonpositive_na = 'gamete_masks_core(NA_integer_, numeric(), numeric(), integer(), integer())',
    no_chr         = 'gamete_masks_core(integer(), numeric(), numeric(), integer(), integer())',
    pos_length     = 'gamete_masks_core(3L, c(0,.5), numeric(), 1L, 0L)',
    num_short_flip = 'numericalize_core(c(0L,2L,0L,2L), 2L, 2L, TRUE, "-101","Add","None")',
    num_na_flip    = 'numericalize_core(c(0L,2L), 2L, 1L, c(TRUE,NA), "-101","Add","None")',
    num_raw3       = 'numericalize_core(c(0L,3L), 1L, 2L, FALSE, "-101","Add","None")',
    num_raw_long   = 'numericalize_core(c(0L,1L,2L,0L,1L), 1L, 2L, FALSE, "-101","Add","None")',
    num_raw_short  = 'numericalize_core(0L, 1L, 2L, FALSE, "-101","Add","None")',
    num_code_as    = 'numericalize_core(c(0L,2L), 1L, 2L, FALSE, "abc","Add","None")',
    num_model      = 'numericalize_core(c(0L,2L), 1L, 2L, FALSE, "-101","add","None")',
    num_impute     = 'numericalize_core(c(0L,2L), 1L, 2L, FALSE, "-101","Add","none")',
    num_n_neg      = 'numericalize_core(integer(), -1L, 2L, logical(), "-101","Add","None")'
  )
  r <- .run_cases(cases)
  expect_identical(r$status, 0L)
  expect_true(r$done)
  expect_setequal(names(r$result), names(cases))
  expect_true(all(startsWith(r$result, "ERROR")),
              info = paste(names(r$result), r$result, sep = ": ", collapse = "\n"))
  m <- r$result
  expect_match(m[["short_strand"]], "characters but the layout has 3 loci")
  expect_match(m[["long_strand"]], "characters but the layout has 3 loci")
  expect_match(m[["alphabet_2"]], "only '0' and '1'")
  expect_match(m[["surplus_events"]], "need 2 meiosis events")
  expect_match(m[["deficit_events"]], "need 2 meiosis events")
  expect_match(m[["negative_n"]], "non-negative")
  expect_match(m[["bad_design"]], "design must be")
  expect_match(m[["chi_nan"]], "chiasma")
  expect_match(m[["chi_negative"]], "chiasma")
  expect_match(m[["chi_beyond_L"]], "chiasma")
  expect_match(m[["counts_sum"]], "counts sum to 1 but 2 chiasmata")
  expect_match(m[["flip_2"]], "flips must be 0 or 1")
  expect_match(m[["flip_na"]], "flips must be 0 or 1")
  expect_match(m[["unsorted_pos"]], "non-decreasing")
  expect_match(m[["empty_chr"]], "at least one locus")
  expect_match(m[["num_short_flip"]], "flip has 1 entries")
  expect_match(m[["num_na_flip"]], "flip must not contain NA")
  expect_match(m[["num_raw3"]], "0, 1, 2 or NA")
  expect_match(m[["num_raw_long"]], "raw_dosage has 5 values")
  expect_match(m[["num_code_as"]], "code_as must be")
  expect_match(m[["num_model"]], "model must be")
  expect_match(m[["num_impute"]], "impute must be")
})

test_that("the production path checks its arguments before calling the kernel", {
  expect_error(
    .check_meiosis_call(1L, 0, list(p1_cis = "12"), list(counts = c(0L, 0L),
      flips = c(0L, 0L), chiasmata = numeric()), 1L, 2L),
    "parental strand")
  expect_error(
    .check_meiosis_call(2L, c(0, 1), list(p1_cis = "11"), list(counts = 0L,
      flips = 0L, chiasmata = numeric()), 1L, 2L),
    "do not match 1 progeny")
  expect_true(.check_meiosis_call(2L, c(0, 1), list(p1_cis = "11"),
    list(counts = c(0L, 0L), flips = c(0L, 0L), chiasmata = numeric()), 1L, 2L))
  # missing genotypes become an R error, not a crash
  pop <- as_population(.geno())
  pop$cis[3, 1] <- NA_integer_
  expect_error(cross(pop[1], as_population(.geno())[2], n = 1, seed = 1),
               "0/1 entry per marker")
})

# ---------------------------------------------------------------------------
# valid-input behaviour is unchanged (RC-07 / T-4 grid, boundaries)
# ---------------------------------------------------------------------------

test_that("numericalize_core: all 64 valid combinations match a reference mapper", {
  raw <- c(0L, 1L, 2L, NA)
  ref <- function(r, flip, code_as, model, impute) {
    v <- if (code_as == "012") c(maj = 2L, het = 1L, min = 0L) else c(maj = 1L, het = 0L, min = -1L)
    add <- if (is.na(r)) {
      switch(impute, None = NA_integer_, Middle = v[["het"]], Minor = v[["min"]],
             Major = v[["maj"]])
    } else if (r == 1L) v[["het"]] else if ((r == 0L) != flip) v[["maj"]] else v[["min"]]
    if (is.na(add)) return(NA_integer_)
    switch(model, Add = add,
           Dom = if (add != v[["het"]]) v[["min"]] else v[["het"]],
           Left = if (add == v[["het"]]) v[["min"]] else add,
           Right = if (add == v[["het"]]) v[["maj"]] else add)
  }
  n <- 0L
  for (code_as in c("-101", "012")) for (model in c("Add", "Dom", "Left", "Right")) {
    for (impute in c("None", "Middle", "Minor", "Major")) for (flip in c(FALSE, TRUE)) {
      got <- numericalize_core(raw, 1L, 4L, flip, code_as, model, impute)
      want <- vapply(raw, ref, integer(1), flip = flip, code_as = code_as,
                     model = model, impute = impute)
      expect_identical(got, unname(want), info = paste(code_as, model, impute, flip))
      n <- n + 1L
    }
  }
  expect_identical(n, 64L)
})

test_that("gamete masks: exact-marker boundary, duplicates cancel, order irrelevant", {
  pos <- c(0, 0.5, 1, 1.5, 2)
  mask <- function(chi, flip = 0L) {
    gamete_masks_core(5L, pos, chi, length(chi), flip)
  }
  # a chiasma exactly on a marker leaves that marker upstream (isqg upper_bound)
  expect_identical(mask(0.5), "00111")
  # at the last position (L): no toggle, and it is a legal value
  expect_identical(mask(2), "00000")
  expect_identical(mask(0), "01111")               # a marker AT the chiasma stays upstream
  expect_identical(mask(c(0.75, 0.75)), "00000")   # duplicates cancel
  expect_identical(mask(c(1.75, 0.75)), mask(c(0.75, 1.75)))
  expect_identical(mask(c(0.75, 1.75)), "00110")
  expect_identical(mask(numeric(), flip = 1L), "11111")
})

test_that("known-answer vectors for the FNV-1a-128 pedigree hash (RFC 9923)", {
  expect_identical(stable_hash_core(c("a", "foobar")),
                   c("d228cb696f1a8caf78912b704e4a8964",
                     "343e1662793c64bf6f0d3597ba446f18"))
  # empty input is the 128-bit offset basis
  expect_identical(stable_hash_core(""), "6c62272e07bb014262b821756295c58d")
  expect_identical(stable_hash_core(character()), character())
})

# ---------------------------------------------------------------------------
# Haldane property (T-9/T-10): absolute bound, >= 20 000 gametes
# ---------------------------------------------------------------------------

test_that("recombination between two markers is Haldane's, within 3 SE, 20 000 gametes", {
  n <- 20000L
  check <- function(first, d_m, seed) {
    pos <- c(first, first + d_m)              # non-zero origin allowed
    set.seed(seed)
    dr <- .draw_meiosis(list(pos), n)
    masks <- gamete_masks_core(2L, pos, dr$chiasmata, dr$counts, dr$flips)
    r_hat <- mean(substr(masks, 1, 1) != substr(masks, 2, 2))
    r <- (1 - exp(-2 * d_m)) / 2
    se <- sqrt(r * (1 - r) / n)
    expect_lt(abs(r_hat - r), 3 * se,
              label = sprintf("|r_hat - r| for origin %.2f, d = %.2f M", first, d_m))
  }
  check(0, 0.09, 101)
  check(0, 0.5, 102)
  check(0.37, 1, 103)    # origin != 0: the Poisson mean is the LAST position
})

test_that(".draw_meiosis draws Poisson(last position), not the span", {
  set.seed(1)
  dr <- .draw_meiosis(list(c(0.5, 1.5)), 6000)
  # span 1, last position 1.5
  se <- sqrt(1.5 / 6000)
  expect_lt(abs(mean(dr$counts) - 1.5), 4 * se)
  # identical stream to a manual rpois(1, last) / runif / rbinom loop
  set.seed(2); a <- .draw_meiosis(list(c(0.5, 0.6)), 5)$counts
  set.seed(2)
  b <- integer(5)
  for (i in 1:5) {
    b[i] <- k <- stats::rpois(1, 0.6)
    if (k > 0) stats::runif(k, 0, 0.6)
    stats::rbinom(1, 1, 0.5)
  }
  expect_identical(a, b)
})

test_that("an L = 0 chromosome consumes exactly one Bernoulli per gamete (T2)", {
  pop <- as_population(.geno(n_mk = 4, cm = rep(0, 4)))
  set.seed(5)
  invisible(suppressMessages(cross(pop[1], pop[2], n = 3)))
  after <- .Random.seed
  set.seed(5); stats::rbinom(6, 1, 0.5)
  expect_identical(after, .Random.seed)
})

# ---------------------------------------------------------------------------
# RUST-F8 / CROSS-F4: chromosome order does not depend on the chr type
# ---------------------------------------------------------------------------

test_that(".chr_rank is numeric-aware and locale-free", {
  expect_identical(.chr_rank(c(10L, 2L, 1L)), c(3L, 2L, 1L))
  expect_identical(.chr_rank(c("10", "2", "1")), c(3L, 2L, 1L))
  expect_identical(.chr_rank(factor(c("1", "10", "2"))), c(1L, 3L, 2L))
  expect_identical(.chr_rank(c(1e5, 2)), c(2L, 1L))
  expect_identical(.chr_rank(c("chr10", "chr2", "chrX", "chr1")), c(3L, 2L, 4L, 1L))
  # numbers before other labels; uppercase before lowercase (byte order)
  expect_identical(.chr_rank(c("X", "chr1", "2", "MT")), c(3L, 4L, 1L, 2L))
})

test_that("integer and character chr give the same seeded progeny (RUST-F8)", {
  g <- .geno(n_mk = 30, chr = rep(c(2L, 10L, 1L), each = 10),
             cm = rep(seq(0, 80, length.out = 10), 3))
  gc <- g
  gc$chr <- as.character(gc$chr)
  run <- function(x) {
    pop <- as_population(x)
    f1 <- cross(pop[1], pop[2], n = 1, seed = 7)
    list(f1 = dosages(f1), dh = dosages(double_haploid(f1, n = 20, seed = 8)),
         f2 = dosages(selfcross(f1, n = 10, seed = 9)))
  }
  expect_identical(run(g), run(gc))
  # and a chromosome order that string sorting would get wrong ("1","10","2")
  gf <- g; gf$chr <- factor(gf$chr)
  expect_identical(run(g), run(gf))
})

# ---------------------------------------------------------------------------
# CROSS-F7: seed= restores the ambient RNG
# ---------------------------------------------------------------------------

test_that("seeded crossing functions restore the caller's RNG state", {
  pop <- as_population(.geno())
  f1 <- cross(pop[1], pop[2], n = 2, seed = 1)
  next_draw <- function(f) {
    set.seed(123)
    invisible(f())
    stats::runif(1)
  }
  base <- {set.seed(123); stats::runif(1)}
  expect_equal(next_draw(function() cross(pop[1], pop[2], n = 3, seed = 9)), base)
  expect_equal(next_draw(function() selfcross(f1[1], n = 3, seed = 9)), base)
  expect_equal(next_draw(function() double_haploid(f1[1], n = 3, seed = 9)), base)
  plan <- mating_design(pop, design = "half_diallel", progeny_per_cross = 1)
  expect_equal(next_draw(function() mate(plan, pop, seed = 9)), base)
  expect_equal(next_draw(function() mating_design(pop, design = "random",
                                                  n_crosses = 3, seed = 9)), base)
  A <- as_population(.geno(), pool = "A")
  B <- as_population(.geno(), pool = "B")
  B$ids <- colnames(B$cis) <- colnames(B$trans) <- c("Q1", "Q2")
  expect_equal(next_draw(function() crossbreed(list(A = A, B = B), "two_way",
                                                n_progeny = 3, seed = 9)), base)
  # results are still reproducible from the seed
  expect_identical(dosages(cross(pop[1], pop[2], n = 3, seed = 9)),
                   dosages(cross(pop[1], pop[2], n = 3, seed = 9)))
  # with no prior RNG state the seeded call leaves none behind
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(cross(pop[1], pop[2], n = 1, seed = 3))
  expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
})

# ---------------------------------------------------------------------------
# CROSS-F3 / VC-002: one judgement of map identity
# ---------------------------------------------------------------------------

test_that("integer vs double chr is one map for pooling, crossing and breeds", {
  g <- .geno()
  gd <- g; gd$chr <- as.numeric(gd$chr)
  p1 <- as_population(g)
  p2 <- as_population(gd)
  expect_false(identical(p1$map, p2$map))
  expect_true(.same_map(p1$map, p2$map))
  expect_s3_class(cross(p1[1], p2[1], n = 1, seed = 1), "Population")
  pool <- suppressWarnings(c(p1, p2))
  expect_s3_class(cross(pool[1], p2[1], n = 1, seed = 1), "Population")
  # a genuinely different map is refused by all of them
  g3 <- g; g3$cm <- g3$cm * 1.1
  p3 <- as_population(g3)
  expect_false(.same_map(p1$map, p3$map))
  expect_error(cross(p1[1], p3[1]), "different marker maps")
  expect_error(c(p1, p3), "different marker maps")
})

# ---------------------------------------------------------------------------
# VC-001: cross-pool allele orientation
# ---------------------------------------------------------------------------

test_that("populations that list a marker's alleles in opposite order draw a warning", {
  a <- as_population(.geno(n_mk = 6, cm = seq(0, 50, length.out = 6), allele = "A/G"))
  b <- as_population(.geno(n_mk = 6, cm = seq(0, 50, length.out = 6), allele = "G/A"))
  expect_warning(cross(a[1], b[2], n = 1, seed = 1), "opposite alleles as \\+1")
  # the same allele order: silent
  expect_no_warning(cross(a[1], a[2], n = 1, seed = 1))
  # disjoint alleles: not the same panel
  d <- as_population(.geno(n_mk = 6, cm = seq(0, 50, length.out = 6), allele = "C/T"))
  expect_error(cross(a[1], d[2], n = 1, seed = 1), "different alleles")
  # a population without the column is not checked
  e <- a; e$map$allele <- NULL
  expect_no_warning(cross(a[1], e[2], n = 1, seed = 1))
  # the allele column is not part of the map identity
  expect_true(.same_map(a$map, b$map))
})

test_that("as_numeric(method = 'reference') / joint conversion fix the orientation", {
  skip_if_not_installed("data.table")
  hmp <- function(calls_a, calls_b) {
    data.frame(`rs#` = c("s1", "s2"), alleles = "A/G", chrom = 1L,
               pos = c(100L, 200L), strand = "+", `assembly#` = NA,
               center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA,
               QCcode = NA, A1 = calls_a, B1 = calls_b,
               check.names = FALSE, stringsAsFactors = FALSE)
  }
  conv <- function(x, ...) suppressMessages(as_numeric(x, to_r = TRUE, ...))
  attach_cm <- function(z) { z$cm <- c(0, 50); z }
  pa <- hmp(c("AA", "AA"), c("AA", "AA"))
  pb <- hmp(c("GG", "GG"), c("GG", "GG"))
  ra <- attach_cm(conv(pa, method = "reference", ref_allele = c("A", "A")))
  rb <- attach_cm(conv(pb, method = "reference", ref_allele = c("A", "A")))
  f1 <- cross(as_population(ra, individuals = 1), as_population(rb, individuals = 1),
              n = 1, seed = 1)
  expect_true(all(dosages(f1) == 0))
})

# ---------------------------------------------------------------------------
# CROSS-F12 / VC-004: subsetting
# ---------------------------------------------------------------------------

test_that("duplicated subscripts get unique ids; an empty subscript prints cleanly", {
  pop <- as_population(.geno(), pool = "A")
  d <- pop[c(1, 1, 2)]
  expect_identical(d$ids, c("P1", "P1_1", "P2"))
  expect_identical(colnames(dosages(d)), d$ids)
  expect_identical(d$keys[1], d$keys[2])
  e <- pop[integer(0)]
  expect_identical(n_individuals(e), 0L)
  expect_no_warning(out <- capture.output(print(e)))
  expect_false(any(grepl("Inf", out)))
})

# ---------------------------------------------------------------------------
# CROSS-F6: units guard
# ---------------------------------------------------------------------------

test_that("a genetic map in Morgans draws a warning; a centiMorgan map does not", {
  expect_warning(as_population(.geno(n_mk = 41, cm = seq(0, 2, length.out = 41))),
                 "Morgans")
  expect_no_warning(as_population(.geno(n_mk = 41, cm = seq(0, 200, length.out = 41))))
  expect_no_warning(as_population(.geno(n_mk = 5, cm = seq(0, 2, length.out = 5))))
})

# ---------------------------------------------------------------------------
# VC-003 / VC-005: synthetic_map guards
# ---------------------------------------------------------------------------

test_that("synthetic_map never returns NaN and rejects duplicated chromosome names", {
  res <- tryCatch(synthetic_map(c(1, 1, 1), c(0, 5e-301, 1e-300), total_cm = 1,
                                width = 1e-300),
                  error = function(e) e)
  if (inherits(res, "error")) {
    expect_match(conditionMessage(res), "representable|too small|not finite")
  } else {
    expect_false(anyNA(res))
    expect_false(is.unsorted(res))
  }
  expect_error(synthetic_map(c(1, 1, 2, 2), c(1, 2, 1, 2) * 1e6,
                             total_cm = c(`1` = 10, `1` = 99, `2` = 5)),
               "more than one entry")
  # realistic scales are unchanged
  expect_equal(synthetic_map(c(1, 1, 1), c(0, 5e5, 1e6), total_cm = 1,
                             suppression = 0), c(0, 0.5, 1))
})

# ---------------------------------------------------------------------------
# VC-006 / VC-007: mating_design / mate input handling
# ---------------------------------------------------------------------------

test_that("allow_self is validated before it is used", {
  pop <- as_population(.geno())
  for (bad in list("no", NA, 1L, c(TRUE, FALSE))) {
    expect_error(mating_design(pop, design = "diallel", allow_self = bad),
                 "`allow_self` must be TRUE or FALSE")
  }
})

test_that("mate() does not overwrite or drop pool columns", {
  pop <- as_population(.geno(), pool = "A")
  plan <- data.frame(mother = "P1", father = "P2", n = 1L, mother_pool = "WRONG",
                     stringsAsFactors = FALSE)
  expect_error(mate(plan, A = pop), "give both pool columns or neither")
  plan2 <- data.frame(mother = "P1", father = "P2", n = 1L, mother_pool = "X",
                      father_pool = "Y", stringsAsFactors = FALSE)
  expect_error(mate(plan2, pop), "names pools")
  expect_error(mate(plan2, A = pop), "not given")
  ok <- data.frame(mother = "P1", father = "P2", n = 1L, mother_pool = "A",
                   father_pool = "A", stringsAsFactors = FALSE)
  expect_s3_class(mate(ok, A = pop, seed = 1), "Population")
  expect_s3_class(mate(ok[, 1:3], A = pop, seed = 1), "Population")
})

# ---------------------------------------------------------------------------
# CROSS-F2 / CROSS-F1: heterosis
# ---------------------------------------------------------------------------

.breeds3 <- function() {
  m <- 3
  base <- data.frame(snp = paste0("m", 1:m), allele = "A/G", chr = 1L,
                     pos = 1:m, cm = c(0, 30, 60), stringsAsFactors = FALSE)
  mk <- function(pool, g) {
    colnames(g) <- paste0(pool, seq_len(ncol(g)))
    as_population(cbind(base, as.data.frame(g)), pool = pool)
  }
  # NOT in Hardy-Weinberg proportions within breed
  A <- mk("A", cbind(c(0L, 1L, -1L), c(0L, 1L, -1L), c(1L, 1L, -1L), c(-1L, 0L, 1L)))
  B <- mk("B", cbind(c(1L, 1L, 1L), c(1L, -1L, 0L), c(-1L, -1L, 1L), c(1L, 0L, 1L)))
  list(A = A, B = B)
}

test_that("heterosis() takes breed means from the labelled pool only", {
  br <- .breeds3()
  f1 <- suppressMessages(crossbreed(br, "two_way", n_progeny = 30, seed = 4))
  qtn <- 1:3; a <- c(0.7, -0.3, 1.1); d <- c(0.5, 0.2, -0.4)
  ok <- heterosis(f1, br, qtn, a, d)
  expect_true(is.finite(ok$realized))
  # a subset of breed A silently changed `realized` before
  expect_error(heterosis(f1, list(A = br$A[1:2], B = br$B), qtn, a, d),
               "does not contain")
  # hand value of the expected F1 heterosis on non-HWE breeds
  gam <- lapply(br, function(b) rowMeans(dosages(b) + 1) / 2)
  h <- function(p) 2 * p * (1 - p)
  hA <- rowMeans(dosages(br$A) == 0); hB <- rowMeans(dosages(br$B) == 0)
  hAB <- gam$A * (1 - gam$B) + (1 - gam$A) * gam$B
  expect_equal(ok$expected_f1["A", "B"], sum(d * (hAB - (hA + hB) / 2)))
  expect_false(isTRUE(all.equal(ok$expected_f1["A", "B"],
                                sum(d * (gam$A - gam$B)^2))))
})

test_that("one-locus non-HWE case: expected F1 heterosis is 0, not the HWE shortcut", {
  base <- data.frame(snp = "m1", allele = "A/G", chr = 1L, pos = 1L, cm = 0,
                     stringsAsFactors = FALSE)
  A <- as_population(cbind(base, A1 = 0L, A2 = 0L), pool = "A")   # all heterozygous
  B <- as_population(cbind(base, B1 = 1L, B2 = 1L), pool = "B")   # all +1
  f1 <- suppressMessages(crossbreed(list(A = A, B = B), "two_way", n_progeny = 4,
                                    seed = 1))
  h <- heterosis(f1, list(A = A, B = B), qtn = 1, a = 0, d = 1)
  expect_equal(h$expected_f1["A", "B"], 0)
})

# ---------------------------------------------------------------------------
# small redistributable format check (RC-08): bundled data -> HapMap file -> NUM
# ---------------------------------------------------------------------------

test_that("a HapMap file written from the bundled panel converts to numeric format", {
  skip_if_not_installed("data.table")
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES", envir = environment())
  g <- SNP55K_maize282_maf04[1:40, 1:9]
  hp <- data.frame(`rs#` = g$snp, alleles = "A/G", chrom = g$chr, pos = g$pos,
                   strand = "+", `assembly#` = NA, center = NA, protLSID = NA,
                   assayLSID = NA, panelLSID = NA, QCcode = NA,
                   check.names = FALSE, stringsAsFactors = FALSE)
  calls <- ifelse(g[, 6:9] == 1, "AA", ifelse(g[, 6:9] == -1, "GG", "AG"))
  hp <- cbind(hp, calls)
  f <- tempfile(fileext = ".hmp.txt")
  utils::write.table(hp, f, sep = "\t", quote = FALSE, row.names = FALSE)
  out <- suppressMessages(as_numeric(f, to_r = TRUE))
  expect_identical(names(out)[1:5], c("snp", "allele", "chr", "pos", "cm"))
  expect_identical(nrow(out), 40L)
  expect_identical(ncol(out), 9L)
  expect_true(all(unlist(out[, 6:9]) %in% c(-1L, 0L, 1L, NA)))
})

# ---------------------------------------------------------------------------
# selfing / DH of an inbred parent; interleaved rows (T1, T4)
# ---------------------------------------------------------------------------

test_that("selfing and DH of a fully inbred parent return the parent exactly", {
  pop <- as_population(.geno())
  s <- selfcross(pop[1], n = 5, seed = 1)
  dh <- double_haploid(pop[2], n = 5, seed = 2)
  expect_true(all(dosages(s) == 1L))
  expect_true(all(dosages(dh) == -1L))
})

test_that("interleaved chromosome rows give the same progeny as sorted rows", {
  g <- .geno(n_mk = 20, chr = rep(1:2, each = 10), cm = rep(seq(0, 90, length.out = 10), 2))
  set.seed(3); g$P1 <- sample(c(-1L, 1L), 20, TRUE); g$P2 <- -g$P1
  shuffled <- g[c(rbind(1:10, 11:20)), ]
  a <- dosages(double_haploid(cross(as_population(g)[1], as_population(g)[2], 1, seed = 1),
                              n = 10, seed = 2))
  b <- dosages(double_haploid(cross(as_population(shuffled)[1], as_population(shuffled)[2], 1, seed = 1),
                              n = 10, seed = 2))
  expect_identical(a, b[rownames(a), ])
})

# ---------------------------------------------------------------------------
# edge maps and founder phasing (Fable T-11, T-12, T-16)
# ---------------------------------------------------------------------------

test_that("single-marker chromosomes and tied cM positions work through selfcross()", {
  g <- .geno(n_mk = 5, chr = c(1L, 1L, 1L, 2L, 3L), cm = c(0, 10, 10, 50, 0))
  pop <- as_population(g)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  f2 <- selfcross(f1, n = 200, seed = 2)
  d <- dosages(f2)
  # tied markers (m2, m3, both at 10 cM, coupling in the F1) are always co-inherited
  expect_true(all(d["m2", ] == d["m3", ]))
  # the single-marker chromosomes segregate 1:2:1 (mean het ~ 0.5)
  expect_equal(mean(d["m4", ] == 0), 0.5, tolerance = 0.2)
  expect_equal(mean(d["m5", ] == 0), 0.5, tolerance = 0.2)
})

test_that("founder heterozygotes are phased allele-1 on the first strand", {
  g <- .geno(n_mk = 3, cm = c(0, 10, 20))
  g$P1 <- c(1L, 0L, -1L)
  pop <- as_population(g, individuals = "P1")
  expect_identical(as.integer(pop$cis[, 1]), c(1L, 1L, 0L))
  expect_identical(as.integer(pop$trans[, 1]), c(1L, 0L, 0L))
})
