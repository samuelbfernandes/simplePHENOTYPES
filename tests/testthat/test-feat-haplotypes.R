# population_from_haplotypes() / haplotypes() (breedingDesigner SPEC-0020 item 5)

hap_map <- function(n_per = 20, nchr = 2, cm_step = 5) {
  n <- n_per * nchr
  data.frame(snp = paste0("m", seq_len(n)),
             chr = rep(seq_len(nchr), each = n_per),
             pos = rep(seq_len(n_per) * 1000L, nchr),
             cm = rep((seq_len(n_per) - 1) * cm_step, nchr),
             stringsAsFactors = FALSE)
}
hap_mats <- function(n_mark, n_ind, seed = 1) {
  set.seed(seed)
  list(cis = matrix(rbinom(n_mark * n_ind, 1, 0.5), n_mark, n_ind,
                    dimnames = list(NULL, paste0("i", seq_len(n_ind)))),
       trans = matrix(rbinom(n_mark * n_ind, 1, 0.5), n_mark, n_ind,
                      dimnames = list(NULL, paste0("i", seq_len(n_ind)))))
}

test_that("round trip: as_population -> haplotypes -> population_from_haplotypes", {
  data("SNP55K_maize282_maf04")
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:6)
  h <- haplotypes(pop)
  expect_identical(h$cis + h$trans - 1L, dosages(pop))
  back <- population_from_haplotypes(h$cis, h$trans, pop$map)
  expect_s3_class(back, "Population")
  expect_identical(dosages(back), dosages(pop))
  expect_identical(back$ids, pop$ids)
  expect_identical(back$map, pop$map)
  expect_identical(back$keys, pop$keys)          # same founders
  expect_identical(back$pedigree, pop$pedigree)
  expect_identical(back$cis, pop$cis)
  expect_identical(back$trans, pop$trans)
})

test_that("dosage coding: 1 = counted allele, dosage = cis + trans - 1", {
  m <- hap_map(3, 1)
  cis <- cbind(a = c(1, 1, 0), b = c(0, 0, 1))
  trans <- cbind(a = c(1, 0, 0), b = c(0, 1, 1))
  pop <- population_from_haplotypes(cis, trans, m)
  expect_equal(unname(dosages(pop)), cbind(c(1, 0, -1), c(-1, 0, 1)),
               ignore_attr = TRUE)
  expect_equal(rownames(dosages(pop)), m$snp)
  expect_equal(pop$origin, "founder")
})

test_that("known phase is kept (unlike as_population's arbitrary phasing)", {
  m <- hap_map(4, 1)
  # coupling vs repulsion double heterozygote: same dosages, different phase
  cp <- population_from_haplotypes(cbind(x = c(1, 1, 0, 0)),
                                   cbind(x = c(0, 0, 1, 1)), m)
  rp <- population_from_haplotypes(cbind(x = c(1, 0, 1, 0)),
                                   cbind(x = c(0, 1, 0, 1)), m)
  expect_identical(dosages(cp), dosages(rp))
  expect_false(identical(cp$cis, rp$cis))
  expect_false(identical(cp$keys, rp$keys))
})

test_that("individuals_in_rows and ids/dimnames handling", {
  m <- hap_map(5, 2)
  h <- hap_mats(10, 4)
  a <- population_from_haplotypes(h$cis, h$trans, m)
  b <- population_from_haplotypes(t(h$cis), t(h$trans), m,
                                  individuals_in_rows = TRUE)
  expect_identical(a$cis, b$cis)
  expect_identical(a$ids, c("i1", "i2", "i3", "i4"))
  # explicit ids override names; unnamed matrices need ids
  cis <- unname(h$cis); trans <- unname(h$trans)
  expect_error(population_from_haplotypes(cis, trans, m), "ids are required")
  p <- population_from_haplotypes(cis, trans, m, ids = c("a", "b", "c", "d"))
  expect_identical(p$ids, c("a", "b", "c", "d"))
  p2 <- population_from_haplotypes(h$cis, h$trans, m, ids = 1:4)
  expect_identical(p2$ids, c("1", "2", "3", "4"))
  # pool label is recorded and distinguishes founders
  q <- population_from_haplotypes(h$cis, h$trans, m, pool = "A")
  expect_true(all(q$pedigree$pool == "A"))
  expect_false(any(q$keys %in% a$keys))
  expect_error(population_from_haplotypes(h$cis, h$trans, m, pool = ""),
               "pool")
})

test_that("invalid inputs are rejected with clear errors", {
  m <- hap_map(5, 2)
  h <- hap_mats(10, 3)
  f <- function(...) population_from_haplotypes(...)
  expect_error(f(as.data.frame(h$cis), h$trans, m), "matrix")
  expect_error(f(h$cis, h$trans[-1, ], m), "same dimensions")
  expect_error(f(h$cis, h$trans, m[-1, ]), "1 marker|9 marker")
  # wrong orientation: individuals in rows but flag not set
  expect_error(f(t(h$cis), t(h$trans), m), "individuals_in_rows")
  # 2 / NA / non-numeric values
  bad <- h$cis; bad[2, 2] <- 2
  expect_error(f(bad, h$trans, m), "0 or 1")
  bad <- h$cis; bad[3, 1] <- NA
  expect_error(f(bad, h$trans, m), "missing")
  bad <- h$cis; mode(bad) <- "character"
  expect_error(f(bad, h$trans, m), "numeric")
  # logical haplotypes are accepted
  expect_s3_class(f(h$cis == 1, h$trans == 1, m), "Population")
  # marker names out of order / wrong
  c2 <- h$cis; t2 <- h$trans
  rownames(c2) <- rownames(t2) <- rev(m$snp)
  expect_error(f(c2, t2, m), "marker names")
  rownames(c2) <- rownames(t2) <- m$snp
  expect_s3_class(f(c2, t2, m), "Population")
  # strands of different individual order
  t3 <- h$trans; colnames(t3) <- rev(colnames(t3))
  expect_error(f(h$cis, t3, m), "different individual names")
  # ids
  expect_error(f(h$cis, h$trans, m, ids = c("a", "a", "b")), "unique")
  expect_error(f(h$cis, h$trans, m, ids = c("a", "b")), "one id per individual")
  expect_error(f(h$cis, h$trans, m, ids = c("a", NA, "b")), "complete")
  expect_error(f(h$cis, h$trans, m, ids = c("a", "", "b")), "non-empty")
  # map problems (the same checks as as_population)
  expect_error(f(h$cis, h$trans, m[, c("snp", "chr", "pos")]), "columns")
  m2 <- m; m2$cm[3] <- 100
  expect_error(f(h$cis, h$trans, m2), "non-decreasing")
  m3 <- m; m3$snp[2] <- m3$snp[1]
  expect_error(f(h$cis, h$trans, m3), "unique")
  m4 <- m; m4$cm <- NA_real_
  expect_error(f(h$cis, h$trans, m4), "all NA")
  m5 <- hap_map(20, 2, cm_step = 0.01)
  h5 <- hap_mats(40, 2)
  expect_warning(f(h5$cis, h5$trans, m5), "centiMorgans")
  expect_error(f(h$cis, h$trans, m, individuals_in_rows = NA), "TRUE or FALSE")
  expect_error(haplotypes(1), "Population")
})

test_that("counted/allele columns of the map are kept; absent otherwise", {
  m <- hap_map(5, 2)
  h <- hap_mats(10, 3)
  p0 <- population_from_haplotypes(h$cis, h$trans, m)
  expect_null(p0$map$counted)
  expect_null(p0$map$allele)
  m$allele <- rep(c("A/G", "G/A"), 5)
  m$counted <- rep(c("a", NA), 5)
  p1 <- population_from_haplotypes(h$cis, h$trans, m)
  expect_identical(p1$map$allele, m$allele)
  expect_identical(p1$map$counted, toupper(m$counted))
  # an all-NA counted column is not recorded
  m$counted <- NA_character_
  expect_null(population_from_haplotypes(h$cis, h$trans, m)$map$counted)
  # the orientation guard sees a constructed population's counted allele
  m$counted <- rep("A", 10)
  pa <- population_from_haplotypes(h$cis, h$trans, m)
  m$counted <- rep("G", 10)
  pb <- population_from_haplotypes(h$cis, h$trans, m)
  expect_error(cross(pa[1], pb[1], n = 1, seed = 1), "different alleles as \\+1")
})

test_that("a constructed population works in the whole engine", {
  m <- hap_map(30, 3)
  h <- hap_mats(90, 10, seed = 4)
  pop <- population_from_haplotypes(h$cis, h$trans, m)
  prog <- cross(pop[1], pop[3], n = 6, seed = 11)
  expect_s3_class(prog, "Population")
  expect_equal(n_individuals(prog), 6)
  expect_true(all(dosages(prog) %in% c(-1, 0, 1)))
  s <- selfcross(prog[1], n = 4, seed = 2)
  expect_equal(n_individuals(s), 4)
  dh <- double_haploid(pop[1], n = 3, seed = 3)
  expect_true(all(dosages(dh) %in% c(-1, 1)))
  G <- g_matrix(c(pop, prog))
  expect_equal(dim(G), c(16, 16))
  A <- a_matrix(c(pop, prog))
  expect_equal(dim(A), c(16, 16))
  expect_equal(unname(diag(A)[1:10]), rep(1, 10))
  ph <- additive(simulate_phenotype(prog, seed = 3), prop = 0.6, n_qtn = 5)
  expect_equal(nrow(phenotypes_long(ph)), 6)
  # pedigree of progeny traces to the constructed founders
  expect_true(all(parentage(prog)$mother %in% pop$ids))
  # print and subset work
  expect_output(print(pop), "Individuals: 10")
  expect_identical(haplotypes(pop[2:3])$cis[, 1], haplotypes(pop)$cis[, 2])
})

test_that("a seeded cross is bit-identical when built from the same haplotypes", {
  data("SNP55K_maize282_maf04")
  pop <- as_population(SNP55K_maize282_maf04, individuals = 1:4)
  h <- haplotypes(pop)
  pop2 <- population_from_haplotypes(h$cis, h$trans, pop$map)
  a <- cross(pop[1], pop[3], n = 5, seed = 42)
  b <- cross(pop2[1], pop2[3], n = 5, seed = 42)
  expect_identical(a$cis, b$cis)
  expect_identical(a$trans, b$trans)
  expect_identical(a$keys, b$keys)
  expect_identical(a$pedigree, b$pedigree)
  # and so does a doubled haploid / selfing run
  expect_identical(double_haploid(pop[1], n = 3, seed = 7)$cis,
                   double_haploid(pop2[1], n = 3, seed = 7)$cis)
  expect_identical(selfcross(pop[2], n = 3, seed = 8)$trans,
                   selfcross(pop2[2], n = 3, seed = 8)$trans)
})

test_that("a constructed population is consistent with as_population of the same dosages", {
  m <- hap_map(30, 2)
  h <- hap_mats(60, 6, seed = 9)
  pop <- population_from_haplotypes(h$cis, h$trans, m)
  geno <- data.frame(snp = m$snp, allele = "A/G", chr = m$chr, pos = m$pos,
                     cm = m$cm, dosages(pop), check.names = FALSE,
                     stringsAsFactors = FALSE)
  ap <- as_population(geno)
  expect_identical(dosages(ap), dosages(pop))
  # same genetic state in distribution terms: identical dosage, so identical
  # additive values and genomic relationship matrix
  expect_equal(g_matrix(ap), g_matrix(pop))
  # unknown phase differs from the known one only in the strands of heterozygotes
  expect_equal(unname(ap$cis + ap$trans), unname(pop$cis + pop$trans))
})

test_that("population_from_haplotypes of an existing population's map reuses its orientation record", {
  data("SNP55K_maize282_maf04")
  num <- SNP55K_maize282_maf04
  attr(num, "counted_allele") <- sub("/.*$", "", num$allele)
  pop <- as_population(num, individuals = 1:3)
  expect_false(is.null(pop$map$counted))
  h <- haplotypes(pop)
  back <- population_from_haplotypes(h$cis, h$trans, pop$map)
  expect_identical(back$map$counted, pop$map$counted)
  expect_no_error(c(back, pop[1:2]))   # pooling passes the guard
})
