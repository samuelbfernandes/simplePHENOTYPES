# test-mate.R -- DECISION-025: mating_design() writes a plan, mate() runs it.

.mate_geno <- function(n = 6, m = 30, seed = 1) {
  set.seed(seed)
  g <- matrix(sample(c(-1L, 1L), n * m, replace = TRUE), m, n)
  colnames(g) <- paste0("P", seq_len(n))
  cbind(data.frame(snp = paste0("m", seq_len(m)), allele = "A/G",
                   chr = rep(1:3, each = m / 3), pos = rep(seq_len(m / 3), 3),
                   cm = rep(seq(0, 90, length.out = m / 3), 3),
                   stringsAsFactors = FALSE), as.data.frame(g))
}
G6 <- .mate_geno()
POP <- as_population(G6)

test_that("a one-row plan draws what cross() draws (same genotypes, RNG stream)", {
  plan <- data.frame(mother = "P1", father = "P2", n = 4)
  m <- mate(plan, POP, seed = 9)
  x <- cross(POP[1], POP[2], n = 4, seed = 9)
  expect_identical(unname(dosages(m)), unname(dosages(x)))
  expect_equal(m$ids, paste0("prog_", 1:4))
  expect_equal(attr(m, "plan")$progeny[[1]], m$ids)
  # a row naming one individual twice is a self; design = "dh" a doubled haploid
  s <- mate(data.frame(mother = "P3", father = "P3", n = 2), POP, seed = 1)
  expect_identical(unname(dosages(s)),
                   unname(dosages(selfcross(POP[3], n = 2, seed = 1))))
  d <- mate(data.frame(mother = "P3", father = "P3", n = 2, design = "dh"),
            POP, seed = 1)
  expect_equal(parentage(d)$design, rep("dh", 2))
})

test_that("rows run in plan order with one seed; pedigree matches the plan", {
  plan <- mating_design(POP, design = "half_diallel", progeny_per_cross = 2)
  expect_equal(nrow(plan), choose(6, 2))
  prog <- mate(plan, POP, seed = 3)
  expect_equal(n_individuals(prog), 2 * choose(6, 2))
  p <- parentage(prog)
  expect_equal(p$mother, rep(plan$mother, each = 2))
  expect_equal(p$father, rep(plan$father, each = 2))
  expect_equal(as.integer(table(families(prog, "full_sib"))), rep(2L, 15))
  expect_identical(prog, mate(plan, POP, seed = 3))
})

test_that("designs have the textbook shapes", {
  a <- POP[1:2]; b <- POP[3:6]
  f <- mating_design(a, b, design = "factorial")
  expect_equal(nrow(f), 8L)
  expect_setequal(paste(f$mother, f$father),
                  as.vector(outer(a$ids, b$ids, paste)))
  nd <- mating_design(b, a, design = "nested", mothers_per_father = 2)
  expect_equal(as.integer(table(nd$father)), c(2L, 2L))
  expect_equal(anyDuplicated(nd$mother), 0L)
  expect_error(mating_design(b, a, design = "nested", mothers_per_father = 3),
               "no valid assignment")
  dl <- mating_design(POP, design = "diallel")
  expect_equal(nrow(dl), 6L * 5L)
  expect_equal(nrow(mating_design(POP, design = "diallel", allow_self = TRUE)), 36L)
  expect_equal(nrow(mating_design(POP, design = "half_diallel", allow_self = TRUE)),
               21L)
  r <- mating_design(POP, design = "random", n_crosses = 50, seed = 1)
  expect_true(all(r$mother != r$father))
  expect_identical(r, mating_design(POP, design = "random", n_crosses = 50, seed = 1))
  ro <- mating_design(POP[1:3], POP[2:4], design = "random", n_crosses = 40, seed = 2)
  expect_true(all(ro$mother != ro$father))            # overlapping sets, no selfs
})

test_that("mating across pools records both pools and names progeny", {
  A <- as_population(G6[, c(1:5, 6:8)], pool = "A")
  B <- as_population(G6[, c(1:5, 9:11)], pool = "B")
  plan <- mating_design(A, B, design = "factorial")
  plan$mother_pool <- "A"; plan$father_pool <- "B"
  f1 <- mate(plan, A = A, B = B, seed = 1)
  expect_equal(f1$ids[1:2], c("AxB_1", "AxB_2"))
  anc <- parentage(f1, ancestors = TRUE)
  expect_setequal(anc$pool[anc$design == "founder"], c("A", "B"))
  expect_error(mate(plan[, 1:3], A = A, B = B), "mother_pool and father_pool")
  expect_error(mate(transform(plan, father_pool = "C"), A = A, B = B),
               "not given in")
  expect_error(mate(plan, A, B), "name each one")
})

test_that("plan validation errors clearly", {
  expect_error(mate(data.frame(mother = "P1"), POP), "columns mother, father and n")
  expect_error(mate(data.frame(mother = "P1", father = "P2", n = 0), POP),
               "positive whole")
  expect_error(mate(data.frame(mother = "P1", father = "Px", n = 1), POP),
               "not found")
  expect_error(mate(data.frame(mother = "P1", father = "P2", n = 1, design = "dh"),
                    POP), "different individuals")
  expect_error(mate(data.frame(mother = "P1", father = "P1", n = 1, design = "cross"),
                    POP), "with itself")
})

test_that("equal ids in different pools are crosses; nested keeps m per father (review O2/O3)", {
  A <- as_population(G6[, 1:7], pool = "A")          # P1, P2
  B <- as_population(G6[, c(1:5, 8:9)], pool = "B")
  colnames_B <- B$ids
  B <- .relabel(B, c("P1", "P2"))                    # same ids as A
  f <- mating_design(A, B, design = "factorial")
  expect_equal(nrow(f), 4L)
  r <- mating_design(A, B, design = "random", n_crosses = 200, seed = 1)
  expect_true(any(r$mother == r$father))             # A:P1 x B:P1 is allowed
  # nested within one pool skips the father itself and gives m to each father
  nd <- mating_design(c("A", "B"), design = "nested", mothers_per_father = 1)
  expect_equal(nd$mother, c("B", "A")); expect_equal(nd$father, c("A", "B"))
  nd2 <- mating_design(c("A", "B", "C", "D"), c("A", "B"), design = "nested",
                       mothers_per_father = 2)
  expect_equal(as.integer(table(nd2$father)), c(2L, 2L))
  expect_true(all(nd2$mother != nd2$father))
  expect_equal(anyDuplicated(nd2$mother), 0L)
})

test_that("self exclusion follows individuals, a named single pool runs (review r2 O1/O2)", {
  dup <- c(POP[1], POP[1], POP[2])                   # P1 twice (P1, P1_1), P2
  hd <- mating_design(dup, design = "half_diallel")
  expect_false(any(hd$mother == "P1" & hd$father == "P1_1"))
  expect_equal(nrow(hd), 1L)                         # P1 x P2 once (P1_1 is P1)
  dl <- mating_design(dup, design = "diallel")
  expect_equal(nrow(dl), 2L)
  r <- mating_design(dup, design = "random", n_crosses = 40, seed = 1)
  keys <- dup$keys[match(c(r$mother, r$father), dup$ids)]
  expect_true(all(keys[1:40] != keys[41:80]))
  plan <- mating_design(POP, design = "half_diallel")
  named <- mate(plan, A = POP, seed = 2)
  expect_equal(n_individuals(named), nrow(plan))
  expect_equal(named$ids[1], "A_1")
})

test_that("identity, not ids, decides selfs; nested matching; fathers from fathers (review r3)", {
  dup <- c(POP[1], POP[1], POP[2])                     # P1, P1_1 (same key), P2
  s <- mate(data.frame(mother = "P1", father = "P1_1", n = 2), dup, seed = 1)
  expect_equal(parentage(s)$design, rep("self", 2))
  sf <- as.character(families(s, "selfed"))         # one family, labelled by the
  expect_length(unique(sf), 1L)                     # parent's last display id
  expect_true(sf[1] %in% c("P1", "P1_1"))
  expect_equal(attr(s, "plan")$design, "self")
  # a feasible no-self nested assignment is found (a greedy one would fail)
  nd <- mating_design(c("A", "B"), c("C", "B"), design = "nested",
                      mothers_per_father = 1)
  expect_setequal(paste(nd$mother, nd$father), c("B C", "A B"))
  expect_error(mating_design(c("A", "B"), c("A", "B", "C"), design = "nested",
                             mothers_per_father = 1), "no valid assignment")
  # same individuals under other ids: fathers come from `fathers`
  q <- c(POP, POP)                                     # P1..P6, P1_1..P6_1
  r <- mating_design(q[1:2], q[7:8], design = "random", n_crosses = 20, seed = 1)
  expect_true(all(r$father %in% q$ids[7:8]))
  expect_true(all(q$keys[match(r$mother, q$ids)] != q$keys[match(r$father, q$ids)]))
})

test_that("cross(x, x) is a self; deterministic designs leave the RNG alone; NA design rows (review r5)", {
  s1 <- cross(POP[1], POP[1], n = 5, seed = 7)
  s2 <- selfcross(POP[1], n = 5, seed = 7)
  expect_identical(unname(dosages(s1)), unname(dosages(s2)))
  expect_equal(parentage(s1)$design, rep("self", 5))
  expect_equal(as.character(families(s1, "selfed")), rep("P1", 5))
  expect_identical(s1$keys, s2$keys)                 # the same individuals
  # factorial / nested / diallels draw nothing, so `seed` must not reseed
  for (d in c("factorial", "diallel", "half_diallel")) {
    set.seed(99); ref <- stats::runif(1)
    set.seed(99); mating_design(POP, design = d, seed = 1)
    expect_identical(stats::runif(1), ref)
  }
  set.seed(99); ref <- stats::runif(1)
  set.seed(99); mating_design(POP[1:4], POP[5:6], design = "nested",
                              mothers_per_father = 2, seed = 1)
  expect_identical(stats::runif(1), ref)
  # an NA design entry takes the row default (cross, or self for one individual)
  plan <- data.frame(mother = c("P1", "P3"), father = c("P2", "P3"), n = 1,
                     design = c(NA, NA_character_))
  m <- mate(plan, POP, seed = 1)
  expect_equal(attr(m, "plan")$design, c("cross", "self"))
  expect_error(mate(transform(plan, design = "x"), POP), "or NA")
})

test_that("identity comes from founder pool labels, not mate()'s pool names (review r6)", {
  g <- data.frame(snp = c("a", "b"), allele = "A/G", chr = 1, pos = 1:2,
                  cm = c(0, 50), P1 = c(-1L, 1L), stringsAsFactors = FALSE)
  # unlabelled re-imports of one genotype are one individual ...
  A0 <- as_population(g); B0 <- as_population(g)
  expect_identical(A0$keys, B0$keys)
  expect_error(mating_design(A0, B0, design = "factorial"), "pool =")
  # ... labelled pools are distinct founders
  A <- as_population(g, pool = "A"); B <- as_population(g, pool = "B")
  plan <- mating_design(A, B, design = "factorial", progeny_per_cross = 2)
  expect_equal(nrow(plan), 1L)
  x <- mate(transform(plan, mother_pool = "A", father_pool = "B"), A = A, B = B,
            seed = 1)
  expect_equal(parentage(x)$design, rep("cross", 2))
  expect_equal(nrow(parentage(x, ancestors = TRUE)), 4L)
  expect_true(all(is.na(families(x, "selfed"))))
})

test_that("mixed inputs compare ids; random finds the admissible pair; one row per identity pair (review r7)", {
  f <- as_population(G6[, 1:7])                        # P1, P2
  fm <- mating_design(f, "P1", design = "factorial")
  expect_equal(paste(fm$mother, fm$father), "P2 P1")    # P1 x P1 excluded
  for (s in 1:5) {
    r <- mating_design(c("A", "B"), "A", design = "random", n_crosses = 3, seed = s)
    expect_true(all(r$mother == "B" & r$father == "A"))
  }
  expect_error(mating_design("A", "A", design = "random", n_crosses = 1),
               "no admissible")
  dup <- c(POP[1], POP[1], POP[2])
  expect_equal(nrow(mating_design(dup, design = "half_diallel", allow_self = TRUE)),
               3L)                                      # P1xP1, P1xP2, P2xP2
  fd <- mating_design(dup, POP[3], design = "factorial")
  expect_equal(fd$mother, c("P1", "P2"))
  nd <- mating_design(dup, POP[3:4], design = "nested", mothers_per_father = 1)
  expect_setequal(nd$mother, c("P1", "P2"))             # P1 not shared via P1_1
})

test_that("the diallels ignore `fathers` (review r10)", {
  dup <- c(POP[1], POP[1], POP[2])
  for (d in c("diallel", "half_diallel")) {
    expect_identical(mating_design(dup, fathers = "ignored", design = d),
                     mating_design(dup, design = d))
  }
  expect_equal(nrow(mating_design(dup, fathers = "x", design = "half_diallel")), 1L)
})
