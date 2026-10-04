# test-counted-persist.R
#
# Round 9, gap A: the counted (+1) allele is persisted in an optional `counted`
# column (as_numeric(counted_column = TRUE)) that survives numeric text files and
# row subsetting, is accepted by every reader of numeric-format data, and feeds
# the cross-pool orientation guard; the default output is unchanged.
# Gap B: the default-name overwrite warning is raised when the file is written.

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# an 11-metadata-column HapMap table; `calls` is markers x individuals
cp_hmp <- function(calls, alleles = "A/G", pos = NULL, chrom = 1L) {
  calls <- as.matrix(calls)
  m <- nrow(calls)
  pos <- if (is.null(pos)) seq_len(m) * 100L else pos
  colnames(calls) <- paste0("L", seq_len(ncol(calls)))
  meta <- data.frame(
    `rs#` = paste0("s", seq_len(m)), alleles = rep_len(alleles, m),
    chrom = rep_len(chrom, m), pos = pos, strand = "+", `assembly#` = NA,
    center = NA, protLSID = NA, assayLSID = NA, panelLSID = NA, QCcode = NA,
    check.names = FALSE, stringsAsFactors = FALSE)
  cbind(meta, as.data.frame(calls, stringsAsFactors = FALSE))
}

cp_num <- function(calls, ..., counted_column = FALSE) {
  suppressMessages(as_numeric(cp_hmp(calls, ...), to_r = TRUE, verbose = FALSE,
                              counted_column = counted_column))
}

# panels that share the label "A/G" at every marker but count opposite alleles
cp_calls_a <- function() matrix(c("AA", "AA", "AG", "AA"), 4, 4)  # A counted
cp_calls_g <- function() matrix(c("GG", "GG", "AG", "GG"), 4, 4)  # G counted

# write a panel to a numeric text file, with or without the counted column
cp_file <- function(calls, name, counted_column = TRUE) {
  f <- file.path(tempdir(), name)
  suppressMessages(as_numeric(cp_hmp(calls, chrom = seq_len(nrow(calls))),
                              to_file = TRUE, to_r = FALSE, file_name = f, verbose = FALSE,
                              counted_column = counted_column))
  f
}

cp_read <- function(f) {
  suppressMessages(as_numeric(f, to_r = TRUE, verbose = FALSE))
}

cp_pop <- function(d, ids, pool) {
  k <- if (identical(names(d)[6L], "counted")) 6L else 5L
  names(d)[-seq_len(k)] <- ids
  d$cm <- synthetic_map(d$chr, d$pos)
  as_population(d, pool = pool)
}

# ---------------------------------------------------------------------------
# (1) the column is opt-in and the default output is unchanged
# ---------------------------------------------------------------------------

test_that("default as_numeric() output is unchanged (values, file bytes)", {
  hm <- cp_hmp(rbind(c("AA", "AG", "GG"), c("CC", "CT", "CC")),
               alleles = c("A/G", "C/T"))
  n <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE))
  expect_identical(names(n), c("snp", "allele", "chr", "pos", "cm",
                               "L1", "L2", "L3"))
  expect_identical(attr(n, "counted_allele"), c("A", "C"))
  f <- tempfile(fileext = ".txt")
  suppressMessages(as_numeric(hm, to_file = TRUE, to_r = FALSE, file_name = f,
                              verbose = FALSE))
  expect_identical(readLines(f), c(
    "snp\tallele\tchr\tpos\tcm\tL1\tL2\tL3",
    "s1\tA/G\t1\t100\tNA\t1\t0\t-1",
    "s2\tC/T\t1\t200\tNA\t1\t0\t1"))
  # asking for the column changes nothing else
  nc <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE,
                                    counted_column = TRUE))
  expect_identical(names(nc), c("snp", "allele", "chr", "pos", "cm",
                                "counted", "L1", "L2", "L3"))
  expect_identical(nc$counted, c("A", "C"))
  n_plain <- n
  attr(n_plain, "counted_allele") <- NULL
  expect_identical(nc[-6L], n_plain)
  expect_identical(attr(nc, "counted_allele"), attr(n, "counted_allele"))
  expect_identical(vapply(nc, typeof, ""), c(
    snp = "character", allele = "character", chr = "character",
    pos = "integer", cm = "double", counted = "character",
    L1 = "integer", L2 = "integer", L3 = "integer"))
  g <- tempfile(fileext = ".txt")
  suppressMessages(as_numeric(hm, to_file = TRUE, to_r = FALSE, file_name = g,
                              verbose = FALSE, counted_column = TRUE))
  expect_identical(readLines(g), c(
    "snp\tallele\tchr\tpos\tcm\tcounted\tL1\tL2\tL3",
    "s1\tA/G\t1\t100\tNA\tA\t1\t0\t-1",
    "s2\tC/T\t1\t200\tNA\tC\t1\t0\t1"))
  # to_r and to_file together: the returned table and the file agree
  both <- suppressMessages(as_numeric(hm, to_r = TRUE, to_file = TRUE,
                                      file_name = f, verbose = FALSE))
  expect_identical(both, n)
})

test_that("counted_column is validated and has a meaning only where an allele is counted", {
  hm <- cp_hmp(rbind(c("AA", "AG", "GG")))
  expect_error(as_numeric(hm, counted_column = NA, verbose = FALSE),
               "counted_column")
  expect_error(as_numeric(hm, counted_column = "yes", verbose = FALSE),
               "counted_column")
  expect_error(as_numeric(hm, model = "Dom", counted_column = TRUE,
                          verbose = FALSE), "counts no allele")
  # code_as = "012" counts the same allele
  n012 <- as_numeric(hm, code_as = "012", counted_column = TRUE, verbose = FALSE)
  expect_identical(n012$counted, "A")
  # numeric input: the attribute is promoted to the column; a table that has
  # neither the column nor the attribute is refused
  n <- as_numeric(hm, verbose = FALSE)
  again <- as_numeric(n, counted_column = TRUE, verbose = FALSE)
  expect_identical(again$counted, "A")
  stripped <- n
  attr(stripped, "counted_allele") <- NULL
  expect_error(as_numeric(stripped, counted_column = TRUE, verbose = FALSE),
               "needs a record of the counted allele")
  expect_identical(as_numeric(stripped, verbose = FALSE), stripped)
})

# ---------------------------------------------------------------------------
# (2) round trip through a numeric text file, and through row subsetting
# ---------------------------------------------------------------------------

test_that("the counted record round-trips through a text file into map$counted", {
  hm <- cp_hmp(rbind(c("AA", "AA", "AG", "GG"), c("AA", "AG", "AA", "GG"),
                     c("GG", "GG", "GG", "AA"), c("GG", "GG", "AG", "AA")))
  num <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE,
                                     counted_column = TRUE))
  expect_identical(num$counted, c("A", "A", "G", "G"))
  f <- tempfile(fileext = ".txt")
  suppressMessages(as_numeric(hm, to_file = TRUE, to_r = FALSE, file_name = f,
                              verbose = FALSE, counted_column = TRUE))
  # the file is detected as numeric format and read back with the column
  expect_identical(detect_format(f), "numeric")
  expect_identical(detect_format(num), "numeric")
  back <- cp_read(f)
  expect_identical(names(back)[6L], "counted")
  expect_identical(back$counted, num$counted)
  expect_identical(attr(back, "counted_allele"), num$counted)
  # converting the written file again reproduces it, byte for byte
  f2 <- tempfile(fileext = ".txt")
  suppressMessages(as_numeric(f, to_file = TRUE, to_r = FALSE, file_name = f2,
                              verbose = FALSE))
  expect_identical(readLines(f2), readLines(f))
  # as_population keeps it, identical to the in-memory (attribute) route
  num$cm <- synthetic_map(num$chr, num$pos)
  back$cm <- num$cm
  via_attr <- num[-6L]
  attr(via_attr, "counted_allele") <- num$counted
  pop_attr <- as_population(via_attr)
  pop_file <- as_population(back)
  pop_mem <- as_population(num)
  expect_identical(pop_file$map$counted, c("A", "A", "G", "G"))
  expect_identical(pop_file$map, pop_attr$map)
  expect_identical(pop_mem$map, pop_attr$map)
  expect_identical(colnames(dosages(pop_file)), paste0("L", 1:4))  # no extra id
  expect_identical(dosages(pop_file), dosages(pop_attr))
})

test_that("a plain fread() of the file works as well", {
  hm <- cp_hmp(rbind(c("AA", "AG"), c("GG", "AG")))
  f <- tempfile(fileext = ".txt")
  suppressMessages(as_numeric(hm, to_file = TRUE, to_r = FALSE, file_name = f,
                              verbose = FALSE, counted_column = TRUE))
  raw <- data.table::fread(f, data.table = FALSE)
  raw$cm <- c(0, 1)
  expect_identical(as_population(raw)$map$counted, c("A", "G"))
  # an all-unknown record reads back from text as a logical NA column: unknown
  raw$counted <- NA
  expect_null(as_population(raw)$map$counted)
})

test_that("the record survives row subsetting and reordering", {
  hm <- cp_hmp(rbind(c("AA", "AA", "AG"), c("GG", "GG", "AG"),
                     c("AA", "AG", "AA"), c("GG", "AG", "GG")), chrom = 1:4)
  num <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE,
                                     counted_column = TRUE))
  expect_identical(num$counted, c("A", "G", "A", "G"))
  num$cm <- synthetic_map(num$chr, num$pos)
  sub <- num[c(4, 2, 1), ]
  expect_identical(sub$counted, c("G", "G", "A"))
  p <- as_population(sub)
  expect_identical(p$map$snp, c("s4", "s2", "s1"))
  expect_identical(p$map$counted, c("G", "G", "A"))
  # column subsetting (individuals) keeps it too
  p2 <- as_population(num[, c(1:6, 8)])
  expect_identical(p2$map$counted, num$counted)
  expect_identical(colnames(dosages(p2)), "L2")
  # the attribute, in contrast, does not follow the rows: `[` leaves it in its
  # original order (documented; the column is the durable form)
  expect_identical(attr(num[c(4, 2, 1), ], "counted_allele"),
                   attr(num, "counted_allele"))
})

# ---------------------------------------------------------------------------
# (3) the cross-pool guard uses the column
# ---------------------------------------------------------------------------

test_that("panels from files whose counted columns disagree are refused, even after row subsetting", {
  da <- cp_read(cp_file(cp_calls_a(), "pa.txt"))
  dg <- cp_read(cp_file(cp_calls_g(), "pg.txt"))
  # the label alone cannot tell them apart (same order, same alleles)
  expect_identical(da$allele, dg$allele)
  ids_a <- paste0("a", 1:4)
  ids_g <- paste0("g", 1:4)
  pa <- cp_pop(da, ids_a, "P1")
  pg <- cp_pop(dg, ids_g, "P2")
  expect_error(cross(pa[1], pg[1], n = 1, seed = 1),
               "count different alleles as \\+1")
  expect_error(c(pa, pg), "count different alleles as \\+1")
  # after row subsetting and reordering (the same rows, a new order, on each
  # side so that the maps match) the record still follows its markers
  rows <- c(4, 1, 3)
  pa2 <- cp_pop(da[rows, ], ids_a, "P1")
  pg2 <- cp_pop(dg[rows, ], ids_g, "P2")
  expect_identical(pa2$map$counted, c("A", "A", "A"))
  expect_identical(pg2$map$counted, c("G", "G", "A"))   # s3 is all-het: a tie
  expect_error(c(pa2[1], pg2[1]), "count different alleles as \\+1")
  expect_error(cross(pa2[1], pg2[1], n = 1, seed = 1),
               "count different alleles as \\+1")
  # the same record on both sides: no complaint
  pa3 <- cp_pop(cp_read(cp_file(cp_calls_a(), "pa2.txt")), paste0("b", 1:4),
                "P3")
  expect_no_warning(cross(pa[1], pa3[1], n = 1, seed = 1))
  # without the column the same files fall back to the weaker label check, which
  # sees nothing (this is the gap the column closes)
  la <- cp_pop(cp_read(cp_file(cp_calls_a(), "la.txt", FALSE)), ids_a, "P1")
  lg <- cp_pop(cp_read(cp_file(cp_calls_g(), "lg.txt", FALSE)), ids_g, "P2")
  expect_null(la$map$counted)
  expect_no_error(suppressWarnings(cross(la[1], lg[1], n = 1, seed = 1)))
})

# ---------------------------------------------------------------------------
# (4) invalid records are rejected
# ---------------------------------------------------------------------------

test_that("an invalid counted column is rejected by as_population() and as_numeric()", {
  num <- cp_num(rbind(c("AA", "AG"), c("GG", "AG")), counted_column = TRUE)
  num$cm <- c(0, 1)
  values <- list(c("A/G", "A"),            # not one allele symbol
                 c("", "G"),               # empty
                 c("A B", "G"),            # whitespace
                 c("T", "G"),              # not in the label "A/G"
                 c(TRUE, FALSE),           # logical, not all NA
                 factor(c("A", "G")))      # factor
  for (i in seq_along(values)) {
    d <- num
    d[[6L]] <- values[[i]]
    expect_error(as_population(d), "counted", info = i)
    expect_error(suppressMessages(as_numeric(d, verbose = FALSE)), "counted",
                 info = i)
  }
  # a numeric sixth column is an individual, not a record (older behaviour)
  d <- num[-6L]
  d$counted <- c(1L, -1L)
  p <- as_population(d)
  expect_identical(colnames(dosages(p)), c("L1", "L2", "counted"))
  expect_null(p$map$counted)
  # a table with the record and no individual is refused
  only <- num[1:6]
  expect_error(as_population(only), "individual")
  expect_error(suppressMessages(as_numeric(only, verbose = FALSE)), "individual")
})

# ---------------------------------------------------------------------------
# (5) the other readers of numeric-format data accept the column
# ---------------------------------------------------------------------------

test_that("filter_geno() accepts the column and keeps it, and the attribute, aligned", {
  hm <- cp_hmp(rbind(c("AA", "AA", "AA", "AA"),     # monomorphic
                     c("GG", "GG", "AG", "AA"),
                     c("AA", "AA", "AG", "GG"),
                     c("GG", "GG", "GG", "AA")))
  num <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE,
                                     counted_column = TRUE))
  expect_identical(num$counted, c("A", "G", "A", "G"))
  fil <- suppressMessages(filter_geno(num, verbose = FALSE))
  expect_identical(fil$snp, c("s2", "s3", "s4"))
  expect_identical(names(fil)[6L], "counted")
  expect_identical(fil$counted, c("G", "A", "G"))
  expect_identical(attr(fil, "counted_allele"), c("G", "A", "G"))
  expect_identical(colnames(fil)[-(1:6)], paste0("L", 1:4))
  # without the column the attribute is subset with the rows (it used to keep
  # its full length and no longer match)
  plain <- suppressMessages(as_numeric(hm, to_r = TRUE, verbose = FALSE))
  fp <- suppressMessages(filter_geno(plain, verbose = FALSE))
  expect_identical(attr(fp, "counted_allele"), c("G", "A", "G"))
  fp$cm <- synthetic_map(fp$chr, fp$pos)
  expect_identical(as_population(fp)$map$counted, c("G", "A", "G"))
  # a non-numeric individual column is still reported as such
  bad <- num
  bad$L2 <- as.character(bad$L2)
  expect_error(suppressMessages(filter_geno(bad)), "non-numeric")
})

test_that("a numeric column named 'counted' is still an individual", {
  num <- cp_num(rbind(c("AA", "AG"), c("GG", "AG")))
  names(num)[6L] <- "counted"
  attr(num, "counted_allele") <- NULL
  num$cm <- c(0, 1)
  expect_null(as_population(num)$map$counted)
  expect_identical(colnames(dosages(as_population(num))), c("counted", "L2"))
  expect_identical(nrow(suppressMessages(filter_geno(num))), 2L)
})

# ---------------------------------------------------------------------------
# gap B: the default-name overwrite warning is raised when writing
# ---------------------------------------------------------------------------

test_that("a failing conversion does not warn about an existing default-named file", {
  dir <- tempfile("cpB")
  dir.create(dir)
  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)
  hm <- cp_hmp(rbind(c("AA", "AG"), c("GG", "AG")))
  writeLines("previous", "hm_numeric.txt")
  # method = "reference" without ref_allele fails after the default name exists
  expect_no_warning(expect_error(
    as_numeric(hm, to_file = TRUE, to_r = FALSE, method = "reference",
               verbose = FALSE),
    "ref_allele"))
  expect_identical(readLines("hm_numeric.txt"), "previous")   # untouched
  # a conversion that succeeds replaces the file and warns, once
  w <- NULL
  withCallingHandlers(
    suppressMessages(as_numeric(hm, to_file = TRUE, to_r = FALSE,
                                verbose = FALSE)),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  expect_length(w, 1L)
  expect_match(w, "default output file hm_numeric.txt already exists")
  expect_false(identical(readLines("hm_numeric.txt"), "previous"))
  # an explicit file name never warns, even over an existing file
  expect_no_warning(suppressMessages(
    as_numeric(hm, to_file = TRUE, to_r = FALSE, file_name = "hm_numeric.txt",
               verbose = FALSE)))
  # no default-named file yet: no warning
  file.remove("hm_numeric.txt")
  expect_no_warning(suppressMessages(
    as_numeric(hm, to_file = TRUE, to_r = FALSE, verbose = FALSE)))
  expect_true(file.exists("hm_numeric.txt"))
  # in-memory result only: no file is touched
  file.remove("hm_numeric.txt")
  expect_no_warning(suppressMessages(as_numeric(hm, verbose = FALSE)))
  expect_false(file.exists("hm_numeric.txt"))
})

test_that("a counted column that contradicts the attribute is an error, NA entries are filled (round 9)", {
  nc <- cp_num(cp_calls_a(), counted_column = TRUE)
  nc$cm <- synthetic_map(nc$chr, nc$pos)
  att <- attr(nc, "counted_allele")
  expect_false(is.null(att))
  # agreeing records convert as before
  expect_s3_class(as_population(nc), "Population")
  # a stale, contradicting attribute is no longer silently ignored
  bad <- nc
  attr(bad, "counted_allele") <- ifelse(att == "A", "G", "A")
  expect_error(as_population(bad), "disagree")
  # an all-unknown column takes the attribute instead of dropping it
  unk <- nc
  unk$counted <- NA
  expect_identical(as_population(unk)$map$counted, toupper(att))
  # an attribute of another length (stale after row subsetting) is ignored
  sub <- nc[1:2, ]
  attr(sub, "counted_allele") <- att
  expect_s3_class(as_population(sub), "Population")
})
