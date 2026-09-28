# test-pedigree.R
#
# DECISION-024: a Population records its pedigree -- founders from
# as_population(), progeny parents from cross() / selfcross() /
# double_haploid(), ancestors kept through `[`, pedigrees pooled by c(). Links use
# internal keys, so colliding display ids (renamed by c()) do not break them, and
# the bookkeeping draws nothing, so genotypes are unchanged.

.ped_geno <- function(n = 4, m = 30, seed = 1) {
  set.seed(seed)
  g <- matrix(sample(c(-1L, 1L), n * m, replace = TRUE), m, n)
  colnames(g) <- paste0("P", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:3, each = m / 3), pos = rep(seq_len(m / 3), 3),
                   cm = rep(seq(0, 90, length.out = m / 3), 3),
                   stringsAsFactors = FALSE), as.data.frame(g))
}
FND <- as_population(.ped_geno(), pool = "A")

test_that("founders and crosses record parents, generation, pool and design", {
  p <- parentage(FND)
  expect_equal(p$id, paste0("P", 1:4))
  expect_true(all(is.na(p$mother) & is.na(p$father)))
  expect_equal(p$generation, rep(0L, 4)); expect_equal(p$pool, rep("A", 4))
  f1 <- cross(FND[1], FND[2], n = 3, seed = 1)
  q <- parentage(f1)
  expect_equal(q$mother, rep("P1", 3)); expect_equal(q$father, rep("P2", 3))
  expect_equal(q$generation, rep(1L, 3)); expect_equal(q$design, rep("cross", 3))
  s2 <- selfcross(f1[1], n = 2, seed = 2)
  r <- parentage(s2)
  expect_equal(r$mother, rep("prog_1", 2)); expect_equal(r$father, r$mother)
  expect_equal(r$generation, rep(2L, 2)); expect_equal(r$design, rep("self", 2))
  dh <- double_haploid(f1[2], n = 2, seed = 3)
  expect_equal(parentage(dh)$design, rep("dh", 2))
  # ancestors are kept: F1 parents -> founders
  a <- parentage(s2, ancestors = TRUE)
  expect_setequal(a$id, c("P1", "P2", "prog_1", "self_1", "self_2"))
})

test_that("pedigree bookkeeping does not change genotypes or the RNG stream", {
  a <- cross(FND[1], FND[2], n = 5, seed = 7)
  b <- cross(FND[1], FND[2], n = 5, seed = 7)
  expect_identical(a, b)                          # keys are deterministic
  set.seed(3); x <- cross(FND[1], FND[2], n = 2); after <- stats::runif(1)
  set.seed(3); y <- cross(FND[1], FND[2], n = 2); expect_identical(stats::runif(1), after)
  expect_identical(dosages(x), dosages(y))
})

test_that("pooling keeps links even when display ids collide", {
  x1 <- cross(FND[1], FND[2], n = 2, seed = 1)   # prog_1, prog_2
  x2 <- cross(FND[3], FND[4], n = 2, seed = 2)   # prog_1, prog_2 again
  pool <- c(x1, x2)
  p <- parentage(pool)
  expect_equal(anyDuplicated(p$id), 0L)          # renamed by c()
  expect_equal(p$mother, c("P1", "P1", "P3", "P3"))
  expect_equal(p$father, c("P2", "P2", "P4", "P4"))
  # a renamed individual's progeny point at its new id
  g2 <- selfcross(pool[3], n = 1, seed = 5)
  expect_equal(parentage(g2)$mother, pool$ids[3])
  expect_equal(parentage(g2, ancestors = TRUE)$id[1:2], c("P3", "P4"))
})

test_that("subsetting keeps only the ancestors of the kept individuals", {
  pool <- c(cross(FND[1], FND[2], n = 2, seed = 1),
            cross(FND[3], FND[4], n = 2, seed = 2))
  sub <- pool[3]
  a <- parentage(sub, ancestors = TRUE)
  expect_setequal(a$id, c("P3", "P4", pool$ids[3]))
  expect_error(pool["nope"], "out of bounds")
})

test_that("families() groups by shared parents", {
  fs <- c(cross(FND[1], FND[2], n = 3, seed = 1),
          cross(FND[2], FND[1], n = 2, seed = 2),       # same pair, reversed
          cross(FND[1], FND[3], n = 2, seed = 3))
  f <- families(fs, "full_sib")
  expect_equal(as.integer(table(f)), c(5L, 2L))
  expect_equal(levels(f), c("P1 x P2", "P1 x P3"))
  m <- families(fs, "maternal_half_sib")
  expect_equal(as.character(m), c(rep("P1", 3), rep("P2", 2), rep("P1", 2)))
  s <- c(selfcross(FND[1], n = 2, seed = 4), double_haploid(FND[2], n = 2, seed = 5),
         cross(FND[3], FND[4], n = 1, seed = 6))
  sf <- families(s, "selfed")
  expect_equal(as.character(sf), c("P1", "P1", "P2", "P2", NA))
  expect_true(all(is.na(families(FND))))
})

test_that("a Population without a recorded pedigree is treated as founders", {
  legacy <- FND; legacy$keys <- NULL; legacy$pedigree <- NULL
  expect_true(all(parentage(legacy)$design == "founder"))
  f1 <- cross(legacy[1], legacy[2], n = 1, seed = 1)
  expect_equal(parentage(f1)$mother, "P1")
  expect_identical(dosages(f1), dosages(cross(FND[1], FND[2], n = 1, seed = 1)))
})

test_that("scheme wrappers carry the pedigree through relabelling", {
  ssd <- single_seed_descent(FND[1:2], generations = 2, seed = 1)
  p <- parentage(ssd, ancestors = TRUE)
  expect_equal(max(p$generation), 2L)
  expect_true(all(parentage(ssd)$design == "self"))
})

test_that("keys stay distinct when draws coincide and across pools (review O1)", {
  # a 0 cM map draws no crossovers, so different seeds can give identical
  # meiosis draws; the RNG state still separates the matings
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = 0,
                  P1 = c(1L, 1L), P2 = c(-1L, -1L), stringsAsFactors = FALSE)
  f <- as_population(g)
  x <- cross(f[1], f[2], n = 1, seed = 2)
  y <- cross(f[1], f[2], n = 1, seed = 3)
  expect_false(identical(x$keys, y$keys))
  pool <- c(x, y)
  expect_equal(parentage(pool)$id, pool$ids)
  # identical founders imported into two pools are two founders
  A <- as_population(g[, 1:6], pool = "A"); B <- as_population(g[, 1:6], pool = "B")
  expect_false(identical(A$keys, B$keys))
  plan <- data.frame(mother = "P1", father = "P1", n = 1,
                     mother_pool = "A", father_pool = "B")
  ab <- mate(plan, A = A, B = B, seed = 1)
  anc <- parentage(ab, ancestors = TRUE)
  expect_setequal(anc$pool[anc$design == "founder"], c("A", "B"))
  expect_equal(parentage(ab)$design, "cross")
  # an exact re-run is the same individual under the same key
  expect_identical(cross(f[1], f[2], n = 1, seed = 2)$keys, x$keys)
})

test_that("parentage exposes stable keys; keys do not depend on rlang (review r2 O3/O4)", {
  A <- as_population(.ped_geno()[, 1:6], pool = "A")
  B <- as_population(.ped_geno()[, 1:6], pool = "B")      # same id P1, other pool
  plan <- data.frame(mother = c("P1", "P1"), father = c("P1", "P1"), n = 1,
                     mother_pool = c("A", "B"), father_pool = c("B", "A"))
  x <- mate(plan, A = A, B = B, seed = 1)
  p <- parentage(x)
  expect_equal(p$mother, c("P1", "P1"))                   # display ids collide ...
  expect_false(p$mother_key[1] == p$mother_key[2])        # ... keys do not
  expect_equal(p$mother_key[1], p$father_key[2])
  # keys are the package's own FNV-1a-128 hash, reproducible by construction
  expect_identical(stable_hash_core(""), "6c62272e07bb014262b821756295c58d")
  expect_equal(nchar(FND$keys), rep(33L, 4))
})

test_that("key encoding is unambiguous (review r3 O1)", {
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = 0,
                  X = c(1L, -1L), stringsAsFactors = FALSE)
  g2 <- g; names(g2)[6] <- "B|X"
  a <- as_population(g, pool = "A|B")
  b <- as_population(g2, pool = "A")
  expect_false(identical(a$keys, b$keys))
  anc <- parentage(c(a, b), ancestors = TRUE)
  expect_equal(nrow(anc), 2L)
  expect_setequal(anc$pool, c("A|B", "A"))
})

test_that("NA pool != '<NA>' pool; unique ids; keys ignore the RNG-kind header (review r4)", {
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = 0,
                  X = c(1L, -1L), stringsAsFactors = FALSE)
  expect_false(identical(as_population(g)$keys,
                         as_population(g, pool = "<NA>")$keys))
  dupg <- cbind(g, X = c(-1L, 1L))
  names(dupg)[7] <- "X"
  expect_error(as_population(dupg), "unique")
  # the RNG-kind header (RNGversion / sample.kind) does not enter the keys
  f <- as_population(.ped_geno())
  a <- cross(f[1], f[2], n = 3, seed = 91)
  suppressWarnings(RNGversion("3.5.0"))
  b <- cross(f[1], f[2], n = 3, seed = 91)
  RNGkind(sample.kind = "Rejection")
  expect_identical(dosages(a), dosages(b))
  expect_identical(a$keys, b$keys)
})

test_that("encoding-equivalent ids give the same key (review r8)", {
  u <- "\u00e9"
  l <- iconv(u, "UTF-8", "latin1")
  skip_if(is.na(l) || Encoding(l) != "latin1", "no Latin-1 conversion here")
  expect_identical(u, l)
  expect_identical(.stable_key("x", u), .stable_key("x", l))
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = 0,
                  X = c(1L, -1L), stringsAsFactors = FALSE)
  g1 <- g; names(g1)[6] <- u
  g2 <- g; names(g2)[6] <- l
  a <- as_population(g1, pool = "A"); b <- as_population(g2, pool = "A")
  expect_identical(a$keys, b$keys)
  expect_error(mating_design(a, b, design = "factorial"), "no matings")
})

test_that("keys ignore the numeric locale; an unseeded RNG still separates matings (review r9)", {
  f <- as_population(.ped_geno())
  a <- cross(f[1], f[2], n = 3, seed = 91)
  old <- Sys.getlocale("LC_NUMERIC")
  loc <- suppressWarnings(Sys.setlocale("LC_NUMERIC", "de_DE.UTF-8"))
  if (nzchar(loc)) {
    b <- tryCatch(cross(f[1], f[2], n = 3, seed = 91),
                  finally = suppressWarnings(Sys.setlocale("LC_NUMERIC", old)))
    expect_identical(b$keys, a$keys)
  }
  expect_identical(.stable_key(0.1), .stable_key(0.1))
  expect_false(identical(.stable_key(0.1), .stable_key(0.1 + 2^-55)))
  # 0 cM map: every draw coincides; with no RNG state before each mating the
  # state R seeds itself with still gives five distinct individuals
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2, cm = 0,
                  P1 = c(1L, 1L), P2 = c(-1L, -1L), stringsAsFactors = FALSE)
  z <- as_population(g)
  saved <- .Random.seed_safe()
  keys <- vapply(1:5, function(i) {
    if (exists(".Random.seed", envir = globalenv())) {
      rm(".Random.seed", envir = globalenv())
    }
    cross(z[1], z[2], n = 1)$keys
  }, character(1))
  if (!is.null(saved)) assign(".Random.seed", saved, envir = globalenv())
  expect_equal(anyDuplicated(keys), 0L)
})
