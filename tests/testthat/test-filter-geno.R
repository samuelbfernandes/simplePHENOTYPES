data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

test_that("filter_geno applies MAF, heterozygosity and monomorphic filters", {
  # this panel's MAF spans ~0.40-0.50, so filter within that range
  f_maf <- filter_geno(G, maf_below = 0.45, verbose = FALSE)
  expect_lt(nrow(f_maf), nrow(G))

  f_inc <- filter_geno(G, hets = "include", verbose = FALSE)
  expect_true(all(rowSums(as.matrix(f_inc[, -(1:5)]) == 0) > 0))

  f_rem <- filter_geno(G, hets = "remove", verbose = FALSE)
  expect_true(all(rowSums(as.matrix(f_rem[, -(1:5)]) == 0) == 0))

  # every returned frame keeps the numeric-format metadata columns
  expect_identical(tolower(names(f_maf)[1:5]),
                   c("snp", "allele", "chr", "pos", "cm"))
})

test_that("filtering to heterozygous markers lets dominance run", {
  g <- filter_geno(G, maf_above = 0.1, hets = "include", verbose = FALSE)
  ph <- simulate_phenotype(g, seed = 60) |>
    additive(prop = 0.4, n_qtn = 5) |>
    dominance(prop = 0.1)
  expect_s3_class(ph, "phenotype_sim")
})

test_that("PLINK-style LD pruning reduces markers (indep-pairwise and VIF)", {
  pw <- filter_geno(G, indep_pairwise = c(50, 5, 0.2), verbose = FALSE)
  expect_lt(nrow(pw), nrow(G))
  vf <- filter_geno(G, indep = c(50, 5, 2), verbose = FALSE)
  expect_lt(nrow(vf), nrow(G))
  # a stricter r2 keeps fewer markers than a looser one
  loose <- filter_geno(G, indep_pairwise = c(50, 5, 0.8), verbose = FALSE)
  expect_gte(nrow(loose), nrow(pw))
  # malformed spec is rejected
  expect_error(filter_geno(G, indep_pairwise = c(50, 5), verbose = FALSE),
               "window, step")
})

test_that("phased pruning and Gabriel blocks reduce markers", {
  sub <- G[G$chr == G$chr[1], ][1:200, ]           # keep it quick
  pp <- filter_geno(sub, indep_pairphase = c(50, 5, 0.2), verbose = FALSE)
  expect_lt(nrow(pp), nrow(sub))
  bl <- filter_geno(sub, blocks = TRUE, block_max_kb = 2000, verbose = FALSE)
  expect_lt(nrow(bl), nrow(sub))
  # two-locus EM r2 is 1 for a marker against itself, in [0, 1] generally
  g <- as.numeric(sub[1, -(1:5)]) + 1
  expect_equal(simplePHENOTYPES:::.hap_r2(g, g), 1)
})

test_that("kb windows and blocks need positions; maf/het still work on a matrix", {
  m <- matrix(sample(c(-1L, 0L, 1L), 600L, replace = TRUE), nrow = 20)
  expect_error(
    filter_geno(m, indep_pairwise = c(50, 5, 0.2), window_unit = "kb",
                verbose = FALSE),
    "pos")
  expect_error(filter_geno(m, blocks = TRUE, verbose = FALSE), "pos")
  expect_true(is.matrix(filter_geno(m, maf_above = 0, verbose = FALSE)))
})

test_that("filter_geno rejects a non-numeric-format data frame", {
  bad <- data.frame(a = 1:3, b = 4:6)
  expect_error(filter_geno(bad, verbose = FALSE), "numeric format")
})
