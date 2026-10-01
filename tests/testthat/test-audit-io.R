# test-audit-io.R
#
# Regression tests for the genotype I/O audit (IO-F1..F16, Codex O1..O10):
# allele orientation, raw contract of every reader, call vocabulary, input
# validation of filter_geno(), and round-trip fidelity. The bundled test panel
# has no tied (MAF = 0.5) marker and no allele-column mismatch, so the fixtures
# here are small synthetic ones written to tempfile()s.

test_dir_audit <- normalizePath(file.path(testthat::test_path(), ".."))

# A one-marker (or n-marker) HapMap data frame. `calls` is a list/vector for a
# single marker, or a character matrix (markers x samples).
mk_hmp <- function(calls, alleles = "A/G", snp = NULL, chrom = 1L,
                   pos = NULL) {
  calls <- if (is.matrix(calls)) calls else matrix(calls, nrow = 1L)
  n <- nrow(calls)
  alleles <- rep_len(alleles, n)
  df <- data.frame(
    `rs#` = if (is.null(snp)) paste0("m", seq_len(n)) else snp,
    alleles = alleles, chrom = chrom,
    pos = if (is.null(pos)) 100L * seq_len(n) else pos,
    strand = "+", `assembly#` = NA, center = NA, protLSID = NA,
    assayLSID = NA, panelLSID = NA, QCcode = NA,
    check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_len(ncol(calls))) df[[paste0("s", j)]] <- calls[, j]
  df
}
hmp_num <- function(df, ...) {
  as_numeric(df, to_r = TRUE, verbose = FALSE, ...)
}
geno_of <- function(res) as.matrix(res[, -(1:5), drop = FALSE])

# ---------------------------------------------------------------------------
# IO-F1 / Codex O2: HapMap `alleles` column validated against the calls
# ---------------------------------------------------------------------------

test_that("a stale HapMap alleles column cannot collapse two homozygote classes", {
  # declared A/G, but the calls are C/T: previously 1 1 1 0 (CC and TT both +1)
  expect_warning(
    res <- hmp_num(mk_hmp(c("CC", "CC", "TT", "CT"), alleles = "A/G")),
    "not in the declared allele pair")
  expect_equal(as.numeric(geno_of(res)), c(1, 1, -1, 0))
  expect_identical(res$allele, "C/T")
})

test_that("a foreign homozygote under a declared pair is not coded as allele 2", {
  # declared A/G, calls AA and TT (Codex O2): previously TT was coded as G/G
  expect_warning(res <- hmp_num(mk_hmp(c("AA", "TT"), alleles = "A/G")),
                 "not in the declared allele pair")
  expect_identical(res$allele, "A/T")
  expect_equal(sort(as.numeric(geno_of(res))), c(-1, 1))
})

test_that("a consistent alleles column gives no warning and keeps allele order", {
  expect_no_warning(res <- hmp_num(mk_hmp(c("AA", "GG", "AG", "AA"))))
  expect_identical(res$allele, "A/G")
  expect_equal(as.numeric(geno_of(res)), c(1, -1, 0, 1))
  # only allele 2 observed: still a valid A/G marker, no warning
  expect_no_warning(res2 <- hmp_num(mk_hmp(c("GG", "GG", "GG"))))
  expect_identical(res2$allele, "A/G")
})

test_that("method = 'reference' is checked against the alleles actually observed", {
  df <- mk_hmp(c("CC", "CC", "TT", "CT"), alleles = "A/G")
  expect_warning(
    res <- hmp_num(df, method = "reference", ref_allele = "T"),
    "declared allele pair")
  expect_equal(as.numeric(geno_of(res)), c(-1, -1, 1, 0))   # T is the reference
  expect_warning(
    expect_error(hmp_num(df, method = "reference", ref_allele = "A"),
                 "ref_allele"),
    "declared allele pair")
})

# ---------------------------------------------------------------------------
# IO-F2: every reader delivers the same raw contract (0 = hom allele 1)
# ---------------------------------------------------------------------------

write_tie_vcf <- function(path, missing_id_site = TRUE) {
  writeLines(c(
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4",
    "1\t100\tsnpA\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/0\t1/1\t1/1",   # tie: 2 REF / 2 ALT
    "1\t200\tsnpB\tC\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1",   # not a tie
    if (missing_id_site)
      "1\t300\t.\tG\tT\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/0"),   # missing ID
    path)
  path
}

test_that("a tied marker gives identical dosages from a VCF path and a data frame", {
  skip_if_not_installed("SNPRelate")
  vcf <- write_tie_vcf(tempfile(fileext = ".vcf"))
  from_file <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE))
  vdf <- read.table(vcf, comment.char = "", skip = 1, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(vdf)[1L] <- "#CHROM"
  from_df <- suppressMessages(
    as_numeric(vdf, from = "vcf", to_r = TRUE, verbose = FALSE))

  expect_identical(geno_of(from_file), geno_of(from_df))
  expect_identical(from_file$allele, from_df$allele)
  expect_identical(from_file$snp, from_df$snp)
  # tie rule: allele 1 (REF) is the reference, so the REF homozygote is +1
  expect_equal(as.numeric(geno_of(from_file)[1L, ]), c(1, 1, -1, -1))
  # non-tied marker: REF is the majority and is +1 as well
  expect_equal(as.numeric(geno_of(from_file)[2L, ]), c(1, 1, 1, -1))
  # a VCF site with no ID gets chr:pos in both paths
  expect_identical(from_file$snp[3L], "1:300")
})

test_that("GDS, BED and PED readers agree with the VCF on tied markers", {
  skip_if_not_installed("SNPRelate")
  vcf <- write_tie_vcf(tempfile(fileext = ".vcf"), missing_id_site = FALSE)
  gds <- tempfile(fileext = ".gds")
  SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "copy.num.of.ref",
                           verbose = FALSE)
  gf <- SNPRelate::snpgdsOpen(gds)
  bed_base <- tempfile("tie_bed")
  ped_base <- tempfile("tie_ped")
  SNPRelate::snpgdsGDS2BED(gf, bed_base, verbose = FALSE)
  SNPRelate::snpgdsGDS2PED(gf, ped_base, verbose = FALSE)
  SNPRelate::snpgdsClose(gf)

  ref <- suppressMessages(as_numeric(vcf, to_r = TRUE, verbose = FALSE))
  for (res in list(
    suppressMessages(as_numeric(gds, to_r = TRUE, verbose = FALSE)),
    suppressMessages(as_numeric(paste0(bed_base, ".bed"), to_r = TRUE,
                                verbose = FALSE)),
    suppressMessages(as_numeric(paste0(ped_base, ".ped"), to_r = TRUE,
                                verbose = FALSE)))) {
    expect_identical(geno_of(res), geno_of(ref))
  }
})

test_that(".read_gds_to_raw() delivers raw 0 = homozygous for allele 1", {
  skip_if_not_installed("SNPRelate")
  vcf <- write_tie_vcf(tempfile(fileext = ".vcf"))
  gds <- tempfile(fileext = ".gds")
  SNPRelate::snpgdsVCF2GDS(vcf, gds, method = "copy.num.of.ref",
                           verbose = FALSE)
  gf <- SNPRelate::snpgdsOpen(gds)
  on.exit(SNPRelate::snpgdsClose(gf), add = TRUE)
  parts <- simplePHENOTYPES:::.read_gds_to_raw(gf)
  expect_identical(parts$allele1[1:2], c("A", "C"))
  expect_equal(unname(parts$raw[1L, ]), c(0L, 0L, 2L, 2L))
  expect_equal(unname(parts$raw[3L, ]), c(0L, 1L, 2L, 0L))
})

test_that("a tied HapMap marker follows the allele order of the alleles column", {
  a <- hmp_num(mk_hmp(c("AA", "AA", "GG", "GG"), alleles = "A/G"))
  b <- hmp_num(mk_hmp(c("AA", "AA", "GG", "GG"), alleles = "G/A"))
  expect_equal(as.numeric(geno_of(a)), c(1, 1, -1, -1))
  expect_equal(as.numeric(geno_of(b)), c(-1, -1, 1, 1))
})

# ---------------------------------------------------------------------------
# Codex O1 / IO-F10: FinalReport
# ---------------------------------------------------------------------------

write_fr <- function(rows, header = "SNP Name\tSample ID\tAllele1 - AB\tAllele2 - AB\tChr\tPosition") {
  path <- tempfile(fileext = ".txt")
  writeLines(c("[Header]", "Num SNPs\t1", "[Data]", header, rows), path)
  path
}

test_that("FinalReport literal NA allele fields stay missing", {
  fr <- write_fr(c("m1\ts1\tA\tA\t1\t100",
                   "m1\ts2\t\t\t1\t100",
                   "m1\ts3\tNA\tNA\t1\t100"))
  res <- suppressMessages(as_numeric(fr, to_r = TRUE, verbose = FALSE))
  expect_equal(as.numeric(geno_of(res)), c(1, NA, NA))
  # one-sided blank or NA is a no-call, not a homozygote
  fr2 <- write_fr(c("m1\ts1\tA\tB\t1\t100",
                    "m1\ts2\tA\t\t1\t100",
                    "m1\ts3\tNA\tB\t1\t100",
                    "m1\ts4\tB\tB\t1\t100"))
  res2 <- suppressMessages(as_numeric(fr2, to_r = TRUE, verbose = FALSE))
  # only the B homozygote is observed among the called homozygotes: it is the
  # majority (+1), the heterozygote 0
  expect_equal(as.numeric(geno_of(res2)), c(0, NA, NA, 1))
  # the allele column is filled from the calls
  expect_identical(res2$allele, "B/A")
})

test_that("FinalReport with missing required columns gives its own diagnostic", {
  bad <- tempfile(fileext = ".txt")
  writeLines(c("[Header]", "Num SNPs\t1", "[Data]", "foo\tbar", "x\ty"), bad)
  expect_error(as_numeric(bad, from = "finalreport", to_r = TRUE,
                          verbose = FALSE),
               "cannot find required columns")
})

# ---------------------------------------------------------------------------
# IO-F3 / IO-F4 / Codex O3: call vocabulary and table metadata
# ---------------------------------------------------------------------------

test_that("format_conversion() default hets is the shared vocabulary", {
  expect_identical(eval(formals(simplePHENOTYPES:::format_conversion)$hets,
                        asNamespace("simplePHENOTYPES")),
                   simplePHENOTYPES:::.HETS)
  tb <- data.frame(s1 = "AA", s2 = "GA", s3 = "GG", s4 = "AG",
                   stringsAsFactors = FALSE)
  res <- suppressMessages(as_numeric(tb, from = "table", to_r = TRUE,
                                     verbose = FALSE))
  expect_equal(as.numeric(geno_of(res)), c(1, 0, -1, 0))   # reversed digraph GA
})

test_that("calls are matched case-insensitively and '-' / '0' are missing", {
  p <- simplePHENOTYPES:::parse_hapmap_chars_to_raw
  raw <- p(matrix(c("aa", "AG", "gg", "r"), 1L))
  expect_equal(as.integer(raw), c(0L, 1L, 2L, 1L))
  expect_identical(attr(raw, "alleles")[1L, ], c("A", "G"))
  for (miss in c("-", "0", ".", "n")) {
    raw <- p(matrix(c("A", "A", miss, "G"), 1L))
    expect_equal(as.integer(raw), c(0L, 0L, NA, 2L), info = miss)
  }
})

test_that("a non-biallelic marker is reported by a counted warning, not a message", {
  m <- rbind(c("AA", "TT", "CC", "AA"), c("AA", "GG", "AG", "AA"),
             c("AA", "TT", "AT", "AA"))
  expect_message(
    expect_warning(raw <- simplePHENOTYPES:::parse_hapmap_chars_to_raw(m),
                   "1 marker\\(s\\) are not biallelic"),
    NA)
  expect_true(all(is.na(raw[1L, ])))
})

test_that("unlisted calls under a custom homo vector are counted in a warning", {
  tb <- data.frame(s1 = "AA", s2 = "X", s3 = "GG", s4 = "AG",
                   stringsAsFactors = FALSE)
  expect_warning(
    res <- suppressMessages(as_numeric(tb, from = "table", to_r = TRUE,
                                       verbose = FALSE)),
    "recognised")
  expect_equal(as.numeric(geno_of(res)), c(1, NA, -1, 0))
})

test_that("table allele metadata is letters, also with hets or one homozygote", {
  a <- suppressMessages(as_numeric(
    data.frame(s1 = "AA", s2 = "AG", s3 = "GG", stringsAsFactors = FALSE),
    from = "table", to_r = TRUE, verbose = FALSE))
  b <- suppressMessages(as_numeric(
    data.frame(s1 = "AA", s2 = "AG", stringsAsFactors = FALSE),
    from = "table", to_r = TRUE, verbose = FALSE))
  expect_identical(a$allele, "A/G")
  expect_identical(b$allele, "A/G")
  expect_equal(as.numeric(geno_of(b)), c(1, 0))
})

# ---------------------------------------------------------------------------
# Exhaustive coding table (Codex proposal): raw dosage x model x code_as x impute
# ---------------------------------------------------------------------------

test_that("the coding step matches the documented transform for every combination", {
  raw_vals <- c(0L, 1L, 2L, NA)
  for (flip in c(FALSE, TRUE)) for (code_as in c("-101", "012"))
    for (model in c("Add", "Dom", "Left", "Right"))
      for (impute in c("None", "Middle", "Minor", "Major")) {
        cv <- if (code_as == "012") c(major = 2L, het = 1L, minor = 0L) else
          c(major = 1L, het = 0L, minor = -1L)
        add <- vapply(raw_vals, function(r) {
          if (is.na(r)) {
            switch(impute, None = NA_integer_, Middle = cv[["het"]],
                   Minor = cv[["minor"]], Major = cv[["major"]])
          } else if (r == 1L) cv[["het"]]
          else if ((r == 0L) != flip) cv[["major"]] else cv[["minor"]]
        }, integer(1))
        expected <- vapply(add, function(a) {
          if (is.na(a)) return(NA_integer_)
          switch(model,
                 Add = a,
                 Dom = if (a != cv[["het"]]) cv[["minor"]] else a,
                 Left = if (a == cv[["het"]]) cv[["minor"]] else a,
                 Right = if (a == cv[["het"]]) cv[["major"]] else a)
        }, integer(1))
        got <- simplePHENOTYPES:::numericalize_core(
          raw_vals, 1L, 4L, flip, code_as, model, impute)
        expect_identical(as.integer(got), as.integer(expected),
                         info = paste(flip, code_as, model, impute))
      }
})

# ---------------------------------------------------------------------------
# IO-F5 / IO-F12: detection
# ---------------------------------------------------------------------------

test_that("a numeric matrix with an all-zero 12th column is not HapMap", {
  m <- matrix(0L, 3, 15)
  m[, 1] <- c(-1L, 0L, 1L)
  expect_false(identical(simplePHENOTYPES:::detect_format(m), "hapmap"))
  # a real headerless nucleotide block is still recognised
  df <- as.data.frame(matrix("A", 3, 13), stringsAsFactors = FALSE)
  df[[12]] <- c("AA", "AG", "GG")
  expect_identical(simplePHENOTYPES:::detect_format(df), "hapmap")
})

test_that("detect_format() is case-insensitive and safe on short .txt headers", {
  hmp <- mk_hmp(c("AA", "GG"))
  lower <- hmp
  names(lower) <- tolower(names(lower))
  f <- tempfile(fileext = ".txt")
  data.table::fwrite(lower, f, sep = "\t")
  expect_identical(simplePHENOTYPES:::detect_format(f), "hapmap")
  expect_identical(simplePHENOTYPES:::detect_format(lower), "hapmap")
  short <- tempfile(fileext = ".txt")
  writeLines(c("a\tb\tc", "1\t2\t3"), short)
  expect_identical(simplePHENOTYPES:::detect_format(short), "unknown")
})

test_that("duplicated marker IDs and sample names are an error", {
  dup_snp <- mk_hmp(rbind(c("AA", "GG"), c("AA", "GG")), snp = c("m1", "m1"))
  expect_error(hmp_num(dup_snp), "Duplicated marker ID")
  dup_smp <- mk_hmp(c("AA", "GG"))
  names(dup_smp)[13] <- "s1"
  expect_error(hmp_num(dup_smp), "Duplicated sample name")
  num <- hmp_num(mk_hmp(rbind(c("AA", "GG"), c("AA", "AG"))))
  num$snp[2] <- num$snp[1]
  expect_error(as_numeric(num, from = "numeric", to_r = TRUE, verbose = FALSE),
               "Duplicated marker ID")
})

# ---------------------------------------------------------------------------
# IO-F11: partial / haploid VCF calls
# ---------------------------------------------------------------------------

test_that("partial or haploid GT calls are not reported as multiallelic", {
  vdf <- data.frame(`#CHROM` = 1L, POS = 100L, ID = "a", REF = "A", ALT = "G",
                    QUAL = ".", FILTER = "PASS", INFO = ".", FORMAT = "GT",
                    s1 = "0/0", s2 = "./1", s3 = "0", s4 = "1/1",
                    check.names = FALSE, stringsAsFactors = FALSE)
  w <- testthat::capture_warnings(
    res <- as_numeric(vdf, from = "vcf", to_r = TRUE, verbose = FALSE))
  expect_true(any(grepl("partially missing", w)))
  expect_false(any(grepl("multiallelic", w)))
  expect_equal(as.numeric(geno_of(res)), c(1, NA, NA, -1))
  vdf$s2 <- "0/2"
  w2 <- testthat::capture_warnings(
    as_numeric(vdf, from = "vcf", to_r = TRUE, verbose = FALSE))
  expect_true(any(grepl("multiallelic", w2)))
})

# ---------------------------------------------------------------------------
# IO-F8 / IO-F9: file names and round-trip types
# ---------------------------------------------------------------------------

test_that("the auto-generated output name only rewrites the file extension", {
  d <- file.path(tempfile("root"), "x.txt.d")
  dir.create(d, recursive = TRUE)
  f <- file.path(d, "geno.hmp.txt")
  data.table::fwrite(mk_hmp(rbind(c("AA", "GG", "AG"), c("CC", "CT", "TT")),
                            alleles = c("A/G", "C/T")), f, sep = "\t")
  suppressMessages(as_numeric(f, verbose = FALSE))
  expect_true(file.exists(file.path(d, "geno_numeric.txt")))
  # no extension at all: the suffix is simply appended
  g <- file.path(d, "genonoext")
  file.copy(f, g)
  suppressMessages(as_numeric(g, from = "hapmap", verbose = FALSE))
  expect_true(file.exists(file.path(d, "genonoext_numeric.txt")))
})

test_that("as_numeric(file written from as_numeric(x)) reproduces as_numeric(x)", {
  hmp <- mk_hmp(rbind(c("AA", "GG", "AG", "AA"), c("CC", "CT", "TT", "CC")),
                alleles = c("A/G", "C/T"), chrom = c(1L, 2L))
  first <- hmp_num(hmp)
  f <- tempfile(fileext = ".txt")
  data.table::fwrite(first, f, sep = "\t", na = NA)
  again <- suppressMessages(as_numeric(f, to_r = TRUE, verbose = FALSE))
  # The counted-allele record (`attr(, "counted_allele")`) lives in the R object
  # only; a text file cannot carry it, so the table itself is what must agree.
  expect_false(is.null(attr(first, "counted_allele")))
  attr(first, "counted_allele") <- NULL
  expect_identical(again, first)
  expect_type(first$chr, "character")
  expect_type(first$cm, "double")
})

test_that("chr is character and cm numeric for every reader", {
  skip_if_not_installed("SNPRelate")
  d <- normalizePath(file.path(testthat::test_path(), ".."))
  res <- lapply(c("test.hmp.txt", "test.vcf", "test.bed", "test.ped",
                  "test.gds"), function(f)
    suppressMessages(as_numeric(file.path(d, f), to_r = TRUE, verbose = FALSE)))
  for (r in res) {
    expect_type(r$chr, "character")
    expect_type(r$cm, "double")
    expect_type(r$snp, "character")
  }
})

test_that("as_numeric() rejects malformed path arguments", {
  expect_error(as_numeric(character(0)), "one non-empty path")
  expect_error(as_numeric(c("a.vcf", "b.vcf")), "one non-empty path")
  expect_error(as_numeric(NA_character_), "one non-empty path")
  expect_error(as_numeric(""), "one non-empty path")
})

test_that(".apply_coding() validates what the Rust kernel would trust", {
  ac <- simplePHENOTYPES:::.apply_coding
  meta <- data.frame(snp = "m1", allele = "A/G", chr = "1", pos = 1L, cm = NA_real_)
  expect_error(ac(matrix(c(0L, 3L), 1L), meta, c("a", "b")), "0, 1, 2 or NA")
  expect_error(ac(matrix(c(0L, 2L), 1L), meta, c("a", "b"), code_as = "x"),
               "should be one of")
  ok <- ac(matrix(c(0L, 2L), 1L), meta, c("a", "b"))
  expect_equal(as.numeric(ok[1, 6:7]), c(1, -1))
})

# ---------------------------------------------------------------------------
# Codex O7: export contract
# ---------------------------------------------------------------------------

test_that("as_numeric() is the public entry; format_conversion() is internal glue", {
  ex <- getNamespaceExports("simplePHENOTYPES")
  expect_true("as_numeric" %in% ex)
  expect_false("format_conversion" %in% ex)
})

# ---------------------------------------------------------------------------
# filter_geno(): argument validation (IO-F6, IO-F7, Codex O4/O5/O6/O10)
# ---------------------------------------------------------------------------

mk_num <- function(mat, chr = 1L, pos = NULL) {
  m <- nrow(mat)
  data.frame(snp = paste0("m", seq_len(m)), allele = "A/G", chr = chr,
             pos = if (is.null(pos)) 1000L * seq_len(m) else pos, cm = 0,
             stats::setNames(as.data.frame(mat), paste0("i", seq_len(ncol(mat)))),
             check.names = FALSE, stringsAsFactors = FALSE)
}
d_maf <- mk_num(rbind(c(-1, -1, 1, 1, 1, 1, 1, 1),     # both MAF 0.25
                      c(-1, -1, 1, 1, 1, 1, 1, 1)))

test_that("maf_above / maf_below must be a single finite number in [0, 0.5]", {
  expect_error(filter_geno(d_maf, maf_above = c(0.3, 0), verbose = FALSE),
               "single finite number")
  expect_error(filter_geno(d_maf, maf_above = NA_real_, verbose = FALSE),
               "single finite number")
  expect_error(filter_geno(d_maf, maf_below = 0.7, verbose = FALSE),
               "single finite number")
  expect_error(filter_geno(d_maf, maf_above = "0.1", verbose = FALSE),
               "single finite number")
  expect_error(filter_geno(d_maf, maf_above = 0.3, maf_below = 0.1,
                           verbose = FALSE), "larger than")
  # boundaries are inclusive and unchanged
  expect_equal(nrow(filter_geno(d_maf, maf_above = 0.25, maf_below = 0.25,
                                verbose = FALSE)), 2L)
})

test_that("LD specifications are range-checked", {
  expect_error(filter_geno(d_maf, indep_pairwise = c(2, 1, 1.5),
                           verbose = FALSE), "r2 threshold must be in")
  expect_error(filter_geno(d_maf, indep_pairphase = c(2, 1, 1.5),
                           verbose = FALSE), "r2 threshold must be in")
  expect_error(filter_geno(d_maf, indep = c(2, 1, 0.5), verbose = FALSE),
               "VIF threshold must be at least 1")
  expect_error(filter_geno(d_maf, indep_pairwise = c(2.5, 1.2, 0.5),
                           verbose = FALSE), "whole number")
  expect_error(filter_geno(d_maf, indep_pairwise = c(2, 1, NA),
                           verbose = FALSE), "three positive")
  # r2 = 1 is legal (a no-op), and a fractional kb window is legal
  expect_no_error(filter_geno(d_maf, indep_pairwise = c(2, 1, 1),
                              verbose = FALSE))
  expect_no_error(filter_geno(d_maf, indep_pairwise = c(2.5, 1, 0.5),
                              window_unit = "kb", verbose = FALSE))
})

test_that("logical controls and genotype column types are validated", {
  d <- mk_num(rbind(c(1, 1, 1, 1), c(1, -1, 0, -1)))
  expect_error(filter_geno(d, remove_monomorphic = 1, verbose = FALSE),
               "remove_monomorphic. must be TRUE or FALSE")
  expect_error(filter_geno(d, blocks = 1, verbose = FALSE),
               "blocks. must be TRUE or FALSE")
  expect_error(filter_geno(d, verbose = NA), "verbose. must be TRUE or FALSE")
  d$i1 <- as.character(d$i1)
  expect_error(filter_geno(d, verbose = FALSE), "must be numeric")
})

test_that("kb windows and blocks refuse a missing chr/pos; variant windows do not", {
  m <- rbind(c(-1, 0, 1, 1, 0, -1), c(-1, 0, 1, 1, 0, -1), c(1, 0, -1, 0, 1, 0))
  d <- mk_num(m, chr = NA_integer_, pos = NA_integer_)
  expect_error(filter_geno(d, indep_pairwise = c(2, 1, 0.2),
                           window_unit = "kb", remove_monomorphic = FALSE,
                           verbose = FALSE), "physical map")
  expect_error(filter_geno(d, blocks = TRUE, verbose = FALSE), "physical map")
  kept <- filter_geno(d, indep_pairwise = c(3, 1, 0.5),
                      remove_monomorphic = FALSE, verbose = FALSE)
  expect_lt(nrow(kept), nrow(d))
  # partial NA in pos is caught the same way
  d2 <- mk_num(m); d2$pos[2] <- NA
  expect_error(filter_geno(d2, blocks = TRUE, verbose = FALSE), "physical map")
})

test_that("MAF is computed on called genotypes only and boundaries are inclusive", {
  # 0 0 0 2 2 NA NA NA in 0/1/2 dosage: 4 alt alleles over 2 * 5 called = 0.4
  d <- mk_num(rbind(c(0, 0, 0, 2, 2, NA, NA, NA) - 1,
                    c(0, 0, 0, 0, 0, 0, 0, 2) - 1))
  expect_equal(nrow(filter_geno(d, maf_above = 0.4, verbose = FALSE)), 1L)
  expect_error(filter_geno(d, maf_above = 0.45, verbose = FALSE),
               "no markers passed")
})

test_that(".plink_hap_rsq() is symmetric and bounded on every short genotype pair", {
  g <- as.matrix(expand.grid(rep(list(0:2), 4L)))
  hr <- simplePHENOTYPES:::.plink_hap_rsq
  ok <- TRUE
  for (i in seq_len(nrow(g))) for (j in seq(i, nrow(g), by = 7L)) {
    a <- hr(as.integer(g[i, ]), as.integer(g[j, ]))
    b <- hr(as.integer(g[j, ]), as.integer(g[i, ]))
    if (!is.finite(a) || a < -1e-12 || a > 1 + 1e-9 || abs(a - b) > 1e-9)
      ok <- FALSE
  }
  expect_true(ok)
})

test_that(".plink_cubic_roots() recovers the roots of monic cubics", {
  cr <- simplePHENOTYPES:::.plink_cubic_roots
  # (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
  r <- cr(-6, 11, -6)
  expect_identical(r$n, 3L)
  expect_equal(r$sol, c(1, 2, 3), tolerance = 1e-8)
  # (x-1)^3 = x^3 - 3x^2 + 3x - 1: one (triple) root
  r3 <- cr(-3, 3, -1)
  expect_identical(r3$n, 1L)
  expect_equal(r3$sol[1L], 1, tolerance = 1e-8)
})

# ---------------------------------------------------------------------------
# write_phenotypes(): round trip (IO-F14) and per-QTN share robustness
# ---------------------------------------------------------------------------

test_that("write_phenotypes() round trips names, order, types and NA pattern", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  ph <- simulate_phenotype(SNP55K_maize282_maf04, n_traits = 2, h2 = 0.5,
                           seed = 11) |> additive(n_qtn = 4)
  f <- tempfile(fileext = ".txt")
  write_phenotypes(ph, file = f)
  back <- data.table::fread(f, data.table = FALSE)
  long <- phenotypes_long(ph)
  expect_identical(names(back), names(long))
  expect_identical(as.character(back$id), as.character(long$id))
  expect_identical(back$trait, as.character(long$trait))
  expect_identical(is.na(back$value), is.na(long$value))
  expect_lt(max(abs(back$value - long$value)), 1e-12)
})

test_that("qtn_table() shares stay finite for a constant locus and a zero-variance layer", {
  data("SNP55K_maize282_maf04", package = "simplePHENOTYPES")
  g <- SNP55K_maize282_maf04[1:60, ]
  g[3, -(1:5)] <- 1L                              # constant (monomorphic) locus
  ph <- simulate_phenotype(g, h2 = 0.4, seed = 5) |>
    additive(n_qtn = 5)
  qt <- qtn_table(ph)
  expect_false(anyNA(qt$var_explained))
  expect_true(all(is.finite(qt$var_explained)))
  ly <- ph$layers[[1L]]
  qe <- simplePHENOTYPES:::.layer_qtn_effect(ly, 1L, 1L)
  expect_equal(simplePHENOTYPES:::.qtn_var(ph, ly, 1L, qe$qtn, qe$effect, 0),
               rep(0, length(qe$qtn)))
})
