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
  # haplotypic r2 is 1 for a marker against itself
  g <- as.integer(as.numeric(sub[1, -(1:5)]) + 1)
  expect_equal(simplePHENOTYPES:::.plink_hap_rsq(g, g), 1)
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

test_that("filter_geno filters a Population and keeps map, ids and pedigree", {
  sub <- G[G$chr %in% unique(G$chr)[1:2], ]
  pop <- as_population(sub, individuals = 1:12, pool = "A")
  pf <- filter_geno(pop, maf_above = 0.42, verbose = FALSE)
  expect_s3_class(pf, "Population")
  expect_lt(nrow(pf$map), nrow(pop$map))
  expect_identical(pf$ids, pop$ids)
  expect_identical(pf$keys, pop$keys)
  expect_identical(pf$pedigree, pop$pedigree)
  expect_identical(nrow(pf$cis), nrow(pf$map))
  expect_identical(nrow(pf$trans), nrow(pf$map))
  expect_identical(rownames(dosages(pf)), pf$map$snp)
  expect_true(all(c("snp", "chr", "pos", "cm") %in% names(pf$map)))
  # same marker set as filtering the data frame it came from
  df_f <- filter_geno(sub[, c(1:5, 5 + 1:12)], maf_above = 0.42, verbose = FALSE)
  expect_identical(pf$map$snp, as.character(df_f$snp))
  # the retained markers all satisfy the threshold on the Population's dosages
  d <- dosages(pf)
  p <- rowMeans((d + 1) / 2)
  expect_true(all(pmin(p, 1 - p) >= 0.42))
  # hets filter on a Population
  ph <- filter_geno(pop, hets = "include", verbose = FALSE)
  expect_true(all(rowSums(dosages(ph) == 0) > 0))
  # a bare -1/0/1 Population cannot be declared 0/1/2
  expect_error(filter_geno(pop, code_as = "012", verbose = FALSE), "-1/0/1")
  # a factor `chr` drops the levels of chromosomes that lost every marker
  tiny <- data.frame(snp = c("m1", "m2"), allele = c("A/G", "A/G"),
                     chr = factor(c("1", "2")), pos = c(1, 1), cm = c(0, 0),
                     i1 = c(1, 1), i2 = c(-1, 1), i3 = c(0, 1), i4 = c(-1, 1),
                     stringsAsFactors = FALSE)
  pf1 <- filter_geno(as_population(tiny), remove_monomorphic = TRUE,
                     verbose = FALSE)
  expect_identical(pf1$map$snp, "m1")
  expect_identical(levels(pf1$map$chr), "1")
  expect_false(grepl("Inf", paste(capture.output(print(pf1)), collapse = "")))
})

test_that("filter_geno keeps a mate() result's plan attribute", {
  sub <- G[G$chr %in% unique(G$chr)[1:2], ]
  pop <- as_population(sub, individuals = 1:6)
  prog <- mate(mating_design(pop, design = "random", n_crosses = 3,
                            progeny_per_cross = 2, seed = 2), pop, seed = 7)
  skip_if(is.null(attr(prog, "plan")), "mate() result carries no plan here")
  pf <- filter_geno(prog, maf_above = 0.1, verbose = FALSE)
  expect_identical(attr(pf, "plan"), attr(prog, "plan"))
})

test_that("a filtered Population still crosses and simulates", {
  sub <- G[G$chr %in% unique(G$chr)[1:2], ]
  pop <- filter_geno(as_population(sub, individuals = 1:6),
                     maf_above = 0.42, indep_pairwise = c(50, 5, 0.5),
                     verbose = FALSE)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 3)
  expect_s3_class(f1, "Population")
  expect_identical(nrow(f1$map), nrow(pop$map))
  # progeny of a cross is itself filterable: an F2 segregates -1/0/1
  f2 <- selfcross(f1[1], n = 30, seed = 4)
  f2f <- filter_geno(f2, maf_above = 0.2, hets = "include", verbose = FALSE)
  expect_s3_class(f2f, "Population")
  expect_lte(nrow(f2f$map), nrow(f2$map))
  ph <- simulate_phenotype(f2f, seed = 5) |>
    additive(prop = 0.4, n_qtn = 3) |>
    dominance(prop = 0.1)
  expect_s3_class(ph, "phenotype_sim")
  expect_true(all(qtn_table(ph)$snp %in% f2f$map$snp))
})
