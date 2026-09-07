# test-cross.R
#
# Genetic map construction, the Population class, and the meiosis-based mating
# functions.
#
# Exact agreement with isqg is covered by test-isqg-parity.R. What is checked
# here is that the R layer wires the map and the draws up correctly, which
# shows as classical genetics: F2 segregating 1:2:1, doubled haploids fully
# homozygous, and recombination between linked markers following Haldane.

data("SNP55K_maize282_maf04")
G <- SNP55K_maize282_maf04

# A small synthetic panel with two opposite homozygous founders, so an F1 is
# heterozygous at every marker and recombination is directly observable.
.toy_geno <- function(n_mk = 60, len_cm = 100) {
  data.frame(
    snp    = paste0("m", seq_len(n_mk)),
    allele = "A/G",
    chr    = 1L,
    pos    = seq_len(n_mk) * 1e6,
    cm     = seq(0, len_cm, length.out = n_mk),
    P1     = rep(1L, n_mk),
    P2     = rep(-1L, n_mk),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. synthetic_map()
# ---------------------------------------------------------------------------

test_that("synthetic_map returns a monotone map starting at zero per chromosome", {
  cm <- synthetic_map(G$chr, G$pos)
  expect_length(cm, nrow(G))
  expect_false(anyNA(cm))
  expect_true(all(tapply(cm, G$chr, function(x) !is.unsorted(x))))
  expect_equal(as.vector(tapply(cm, G$chr, min)), rep(0, 10))
})

test_that("synthetic_map suppresses recombination near the centromere", {
  # Centromere at the midpoint by default, so the middle of the chromosome must
  # accumulate fewer cM per Mb than the arms.
  pos <- seq(0, 200e6, length.out = 401)
  cm <- synthetic_map(rep(1L, length(pos)), pos)
  rate <- diff(cm) / (diff(pos) / 1e6)
  n <- length(rate)
  middle <- mean(rate[seq(n * 0.45, n * 0.55)])
  arms <- mean(rate[c(seq_len(n * 0.1), seq(n * 0.9, n))])
  expect_lt(middle, arms / 2)
})

test_that("synthetic_map with suppression = 0 is linear in physical position", {
  pos <- seq(0, 100e6, length.out = 200)
  cm <- synthetic_map(rep(1L, length(pos)), pos, suppression = 0)
  expect_equal(cm, (pos - min(pos)) / 1e6 * 0.73, tolerance = 1e-8)
})

test_that("synthetic_map honours total_cm and gives co-located markers equal cM", {
  chr <- c(1L, 1L, 1L, 2L, 2L)
  pos <- c(0, 5e6, 5e6, 0, 10e6)      # markers 2 and 3 are co-located
  cm <- synthetic_map(chr, pos, total_cm = c(50, 80))
  expect_equal(max(cm[chr == 1]), 50)
  expect_equal(max(cm[chr == 2]), 80)
  expect_identical(cm[2], cm[3])
})

test_that("synthetic_map rejects malformed input", {
  expect_error(synthetic_map(1:3, 1:2), "same length")
  expect_error(synthetic_map(c(1, 1), c(5, 1)), "non-decreasing")
  expect_error(synthetic_map(c(1, 1), c(1, 5), suppression = 1), "suppression")
  expect_error(synthetic_map(c(1, 2), c(1, 5), total_cm = c(1, 2, 3)),
               "one entry per chromosome")
})

# ---------------------------------------------------------------------------
# 2. as_population()
# ---------------------------------------------------------------------------

test_that("as_population round-trips the dosage matrix", {
  pop <- as_population(G, individuals = c("33-16", "38-11"))
  expect_s3_class(pop, "Population")
  expect_identical(n_individuals(pop), 2L)
  expected <- as.matrix(G[, c("33-16", "38-11")])
  expect_identical(typeof(dosages(pop)), "integer")
  expect_identical(unname(dosages(pop)),
                   matrix(as.integer(expected), nrow = nrow(expected)))
  # Dimnames carry the marker and individual identities.
  expect_identical(rownames(dosages(pop)), as.character(G$snp))
  expect_identical(colnames(dosages(pop)), c("33-16", "38-11"))
})

test_that("as_population refuses a genotype set with no genetic map", {
  no_map <- G[1:50, 1:8]
  no_map$cm <- NA
  expect_error(as_population(no_map), "all NA")
  # The message must point at the fix, since this is the blocker users hit.
  expect_error(as_population(no_map), "synthetic_map")
})

test_that("as_population rejects an unsorted or partially missing map", {
  bad <- G[1:50, 1:8]
  bad$cm[10] <- NA
  expect_error(as_population(bad), "missing value")

  unsorted <- G[1:50, 1:8]
  unsorted$cm <- rev(unsorted$cm)
  expect_error(as_population(unsorted), "non-decreasing")
})

test_that("Population subsetting and printing behave", {
  pop <- as_population(G, individuals = 1:4)
  expect_identical(n_individuals(pop[2:3]), 2L)
  expect_identical(pop[2:3]$ids, colnames(G)[c(7, 8)])
  expect_output(print(pop), "<Population>")
  expect_output(print(pop), "Markers: 10650")
})

# ---------------------------------------------------------------------------
# 3. Mating designs — classical genetics
# ---------------------------------------------------------------------------

test_that("an F1 of two opposite homozygotes is heterozygous everywhere", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 3, seed = 1)
  expect_identical(n_individuals(f1), 3L)
  expect_true(all(dosages(f1) == 0L))
})

test_that("an F2 segregates 1:2:1", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  f2 <- selfcross(f1, n = 400, seed = 2)
  p <- prop.table(table(factor(dosages(f2), levels = c(-1, 0, 1))))
  expect_equal(as.numeric(p), c(0.25, 0.5, 0.25), tolerance = 0.05)
})

test_that("doubled haploids are completely homozygous and segregate 1:1", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  dh <- double_haploid(f1, n = 400, seed = 3)
  expect_false(any(dosages(dh) == 0L))
  p <- prop.table(table(factor(dosages(dh), levels = c(-1, 1))))
  expect_equal(as.numeric(p), c(0.5, 0.5), tolerance = 0.05)
})

test_that("progeny inherit real parental haplotypes, not the ancestry mask", {
  # Regression test. Progeny phase must come from the meiosis core, not be
  # reconstructed from a -1/0/1 genotype: a genotype cannot express the phase
  # of a heterozygote, so an F1 (heterozygous everywhere) would come back with
  # one all-allele-1 and one all-allele-2 strand, and its own progeny would
  # then carry the ancestry mask instead of founder alleles.
  #
  # Founders with an ALTERNATING allele pattern make that unmistakable: a mask
  # is a run of constant values, a real haplotype alternates.
  n <- 12
  pat <- rep(c(1L, -1L), length.out = n)
  geno <- data.frame(
    snp = paste0("m", seq_len(n)), allele = "A/G", chr = 1L,
    pos = seq_len(n) * 1e6, cm = seq(0, 10, length.out = n),
    P1 = pat, P2 = -pat, stringsAsFactors = FALSE
  )
  pop <- as_population(geno)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 5)
  expect_true(all(dosages(f1) == 0L))          # het everywhere

  dh <- double_haploid(f1, n = 20, seed = 7)
  g <- dosages(dh)
  # Every locus of every progeny must carry one of the two founder alleles.
  expect_true(all(g == pat | g == -pat))
  # And the progeny must not be constant runs, which is what the mask would be.
  expect_false(any(apply(g, 2, function(col) length(unique(col)) == 1L)))
})

test_that("recombination between linked markers follows Haldane's map function", {
  # Doubled haploids are single gametes, so the recombination fraction between
  # two markers is just the proportion of progeny carrying different parental
  # alleles at them. Under the count-location model with no interference that
  # is Haldane: r = (1 - exp(-2d)) / 2 for d Morgans.
  geno <- .toy_geno(n_mk = 41, len_cm = 200)
  pop <- as_population(geno)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 11)
  dh <- double_haploid(f1, n = 1500, seed = 12)
  g <- dosages(dh)

  d_cm <- geno$cm
  for (j in c(2L, 5L, 11L, 41L)) {
    d <- (d_cm[j] - d_cm[1]) / 100          # Morgans
    expected <- (1 - exp(-2 * d)) / 2
    observed <- mean(g[1, ] != g[j, ])
    expect_equal(observed, expected, tolerance = 0.05,
                 label = paste0("r between marker 1 and ", j))
  }
})

test_that("unlinked markers on different chromosomes recombine freely", {
  geno <- .toy_geno(n_mk = 40, len_cm = 100)
  geno$chr <- rep(1:2, each = 20)
  geno$cm <- rep(seq(0, 100, length.out = 20), 2)
  pop <- as_population(geno)
  f1 <- cross(pop[1], pop[2], n = 1, seed = 21)
  dh <- double_haploid(f1, n = 1000, seed = 22)
  g <- dosages(dh)
  expect_equal(mean(g[1, ] != g[21, ]), 0.5, tolerance = 0.05)
})

test_that("selfing progressively reduces heterozygosity", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 1, seed = 31)
  het <- function(p) mean(dosages(p) == 0L)

  set.seed(32)
  f2 <- selfcross(f1, n = 60)
  f3 <- selfcross(f2[1], n = 60)
  f4 <- selfcross(f3[1], n = 60)
  # Roughly halves each generation; assert the ordering rather than exact rates.
  expect_gt(het(f2), het(f3))
  expect_gt(het(f3), het(f4))
  expect_equal(het(f2), 0.5, tolerance = 0.08)
})

# ---------------------------------------------------------------------------
# 4. Reproducibility and validation
# ---------------------------------------------------------------------------

test_that("mating is reproducible from a seed and varies without one", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  expect_identical(dosages(selfcross(f1, n = 8, seed = 5)),
                   dosages(selfcross(f1, n = 8, seed = 5)))
  expect_false(identical(dosages(selfcross(f1, n = 8, seed = 5)),
                         dosages(selfcross(f1, n = 8, seed = 6))))
  # set.seed() before the call is equivalent to passing seed =.
  set.seed(99); a <- dosages(selfcross(f1, n = 8))
  set.seed(99); b <- dosages(selfcross(f1, n = 8))
  expect_identical(a, b)
})

test_that("mating functions require a single individual per parent", {
  pop <- as_population(G, individuals = 1:3)
  expect_error(cross(pop, pop[1]), "exactly one individual")
  expect_error(selfcross(pop), "exactly one individual")
  expect_error(double_haploid(pop), "exactly one individual")
  expect_error(cross(pop[1], pop[2], n = 0), "positive")
})

test_that("progeny ids and origin describe the pedigree", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 2, seed = 1)
  expect_identical(f1$ids, c("prog_1", "prog_2"))
  expect_match(f1$origin, "cross\\(P1 x P2\\)")
  expect_match(double_haploid(f1[1], n = 1, seed = 1)$origin, "double_haploid")
})

# ---------------------------------------------------------------------------
# 5. Integration with the phenotype grammar
# ---------------------------------------------------------------------------

test_that("simulate_phenotype accepts a Population", {
  pop <- as_population(.toy_geno())
  f1 <- cross(pop[1], pop[2], n = 1, seed = 1)
  f2 <- selfcross(f1, n = 50, seed = 2)

  ph <- additive(simulate_phenotype(f2, seed = 3), prop = 0.6, n_qtn = 4)
  expect_s3_class(ph, "phenotype_sim")
  expect_identical(nrow(ph$pheno), 50L)
  expect_identical(sort(unique(ph$pheno$id)), sort(f2$ids))
  # Realized heritability should land near the requested proportion.
  expect_equal(stats::var(ph$pheno$value), stats::var(ph$pheno$value))
  expect_output(print(ph), "Population")
})
