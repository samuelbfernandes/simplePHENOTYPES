# test-mabc.R
#
# Marker-assisted backcrossing: mabc_select() and recurrent_parent_recovery().
# Validation targets follow the breedingDesigner backcross theory review:
# unselected recovery R_t = 1 - (1/2)^(t+1), the complement of the expected
# donor proportion 1/2^(t+1) (Frisch & Melchinger 2005) (0.5 / 0.75 / 0.875 / 0.9375 at
# F1 / BC1 / BC2 / BC3), a selfing negative control that stays near 0.5, and hard
# foreground invariants for the selected individuals.

# Two opposite homozygous founders, two 1-Morgan chromosomes, 41 markers each:
# every marker informative (recurrent P1 = +1, donor P2 = -1).
.bc_geno <- function(n_per_chr = 41) {
  cm <- seq(0, 100, length.out = n_per_chr)
  data.frame(snp = paste0("m", seq_len(2 * n_per_chr)), allele = "A/G",
             chr = rep(1:2, each = n_per_chr),
             pos = rep(seq_len(n_per_chr), 2) * 1e6, cm = rep(cm, 2),
             P1 = 1L, P2 = -1L, stringsAsFactors = FALSE)
}
FOUNDERS <- as_population(.bc_geno())
REC <- FOUNDERS[1]
DON <- FOUNDERS[2]

# Backcross every individual of `pop` once to the recurrent parent.
.backcross_all <- function(pop, seed0) {
  do.call(c, lapply(seq_len(n_individuals(pop)), function(i)
    cross(pop[i], REC, n = 1, seed = seed0 + i)))
}
F1  <- cross(REC, DON, seed = 1)
BC1 <- cross(F1, REC, n = 300, seed = 2)
BC2 <- .backcross_all(BC1, 1000)
BC3 <- .backcross_all(BC2, 2000)

test_that("unselected recovery follows 1 - (1/2)^(t+1)", {
  expect_equal(as.numeric(recurrent_parent_recovery(F1, REC, DON)), 0.5)
  expect_equal(mean(recurrent_parent_recovery(BC1, REC, DON)), 0.75,
               tolerance = 0.02 / 0.75)
  expect_equal(mean(recurrent_parent_recovery(BC2, REC, DON)), 0.875,
               tolerance = 0.02 / 0.875)
  expect_equal(mean(recurrent_parent_recovery(BC3, REC, DON)), 0.9375,
               tolerance = 0.02 / 0.9375)
  expect_equal(attr(recurrent_parent_recovery(BC1, REC, DON), "n_markers"), 82L)
})

test_that("a selfing chain is not a backcross: recovery stays near 0.5", {
  F2 <- selfcross(F1, n = 300, seed = 3)
  F3 <- do.call(c, lapply(seq_len(n_individuals(F2)), function(i)
    selfcross(F2[i], n = 1, seed = 3000 + i)))
  expect_equal(mean(recurrent_parent_recovery(F2, REC, DON)), 0.5,
               tolerance = 0.03 / 0.5)
  expect_equal(mean(recurrent_parent_recovery(F3, REC, DON)), 0.5,
               tolerance = 0.03 / 0.5)
})

test_that("every selected BC individual carries the donor target allele", {
  sel <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 10, seed = 4)
  d <- attr(sel, "mabc")
  expect_equal(n_individuals(sel), 10L)
  expect_true(all(d$target[d$selected] == "H"))    # BC1: carrier = heterozygote
  expect_true(all(d$feasible[d$selected]))
  expect_equal(sort(colnames(sel$cis)), sort(d$id[d$selected]))
  # a BC1 has no donor homozygotes, so feasibility is exactly the heterozygotes
  expect_true(all(d$target %in% c("R", "H")))
  expect_identical(d$feasible, d$target == "H")
})

test_that("background selection lifts recovery above the cohort mean", {
  sel <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 5, seed = 5)
  d <- attr(sel, "mabc")
  expect_gt(mean(d$background_recovery[d$selected]),
            mean(d$background_recovery[d$feasible]) + 0.05)
  # lexicographic: no unselected feasible candidate outranks a selected one
  expect_true(min(d$background_recovery[d$selected]) >=
                max(d$background_recovery[d$feasible & !d$selected]))
})

test_that("donor_homozygote selects only donor homozygotes after a self", {
  bc3_carrier <- mabc_select(BC3, REC, DON, target_markers = "m10", n = 1,
                             seed = 6)
  bc3f2 <- selfcross(bc3_carrier, n = 60, seed = 7)
  sel <- mabc_select(bc3f2, REC, DON, target_markers = "m10",
                     target_requirement = "donor_homozygote", n = 3, seed = 8)
  d <- attr(sel, "mabc")
  expect_true(all(d$target[d$selected] == "D"))
  expect_identical(d$feasible, d$target == "D")
})

test_that("multiple targets: every target must meet the requirement", {
  sel <- mabc_select(BC1, REC, DON, target_markers = c("m10", "m60"), n = 3,
                     seed = 9)
  d <- attr(sel, "mabc")
  expect_true(all(d$target[d$selected] == "H;H"))
  expect_identical(d$feasible, d$target == "H;H")
})

test_that("flanking markers take precedence over background recovery", {
  sel <- mabc_select(BC1, REC, DON, target_markers = "m10",
                     flanking_markers = c("m8", "m12"), n = 5, seed = 10)
  d <- attr(sel, "mabc")
  expect_true(min(d$flank_recurrent[d$selected]) >=
                max(d$flank_recurrent[d$feasible & !d$selected]))
  info <- attr(sel, "mabc_info")
  expect_false(any(c("m8", "m10", "m12") %in% info$background_markers))
  expect_length(info$background_markers, info$n_background)
  expect_equal(info$n_background, 82L - 3L)
})

test_that("flanking markers must be linked to, and should bracket, the target", {
  # review O2: an unlinked "flank" would otherwise outrank genuine recombinants
  expect_error(mabc_select(BC1, REC, DON, target_markers = "m10",
                           flanking_markers = c("m8", "m50")),
               "not on a chromosome carrying a target")
  expect_warning(mabc_select(BC1, REC, DON, target_markers = "m10",
                             flanking_markers = "m8", n = 1, seed = 1),
                 "do not bracket")
  expect_silent(suppressMessages(
    mabc_select(BC1, REC, DON, target_markers = "m10",
                flanking_markers = c("m8", "m12"), n = 1, seed = 1)))
})

test_that("exclude_interval removes the target region from the background", {
  sel <- mabc_select(BC1, REC, DON, target_markers = "m10",
                     exclude_interval = list(chr = 1, from = 15, to = 30),
                     n = 1, seed = 11)
  info <- attr(sel, "mabc_info")
  in_iv <- FOUNDERS$map$snp[FOUNDERS$map$chr == 1 & FOUNDERS$map$cm >= 15 &
                              FOUNDERS$map$cm <= 30]
  expect_setequal(info$excluded_markers, in_iv)
  expect_equal(info$n_background, 82L - length(union(in_iv, "m10")))
})

test_that("tie-breaking is seeded and reproducible", {
  a <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 4, seed = 12)
  b <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 4, seed = 12)
  expect_identical(attr(a, "mabc"), attr(b, "mabc"))
  set.seed(99); before <- stats::runif(1)
  set.seed(99)
  invisible(mabc_select(BC1, REC, DON, target_markers = "m10", n = 4, seed = 12))
  expect_identical(stats::runif(1), before)          # ambient RNG restored
})

test_that("interval weights track genome length, not marker density", {
  map <- data.frame(snp = paste0("x", 1:6), chr = 1,
                    pos = 1:6, cm = c(0, 1, 2, 3, 50, 100))
  expect_equal(.mabc_weights(map, 1:6, "interval"),
               c(0.5, 1, 1, 24, 48.5, 25))
  expect_equal(.mabc_weights(map, 1:6, NULL), rep(1, 6))
})

test_that("non-informative targets and infeasible requests error clearly", {
  g <- .bc_geno(); g$P2[g$snp == "m10"] <- 1L          # m10 not informative
  fnd <- as_population(g)
  bc <- cross(cross(fnd[1], fnd[2], seed = 1), fnd[1], n = 20, seed = 2)
  expect_error(mabc_select(bc, fnd[1], fnd[2], target_markers = "m10"),
               "not homozygous for alternate alleles at m10")
  expect_error(mabc_select(BC1, REC, DON, target_markers = "m10",
                           target_requirement = "donor_homozygote", n = 1),
               "only 0 of 300 candidates")
  expect_error(mabc_select(BC1, REC, DON, target_markers = "m10",
                           flanking_markers = "m10"), "both a target")
  expect_error(mabc_select(BC1, REC, DON, target_markers = "nope"),
               "not in the map")
  expect_error(mabc_select(BC1, REC, FOUNDERS, target_markers = "m10"),
               "exactly one individual|single")
})

test_that("explicit non-informative markers are dropped with a warning", {
  g <- .bc_geno(); g$P2[g$snp %in% c("m1", "m2")] <- 1L
  fnd <- as_population(g)
  bc <- cross(cross(fnd[1], fnd[2], seed = 1), fnd[1], n = 5, seed = 2)
  expect_warning(r <- recurrent_parent_recovery(bc, fnd[1], fnd[2],
                                                markers = c("m1", "m2", "m3")),
                 "2 marker")
  expect_equal(attr(r, "n_markers"), 1L)
  expect_equal(attr(recurrent_parent_recovery(bc, fnd[1], fnd[2]), "n_markers"),
               80L)
})

test_that("a reverse-coded marker swaps the homozygote calls (documented)", {
  # review round 2, O1: the coding cannot be verified from dosages, so it is a
  # documented precondition, not a check. Negating one marker's dosages (what a
  # separately numericalized panel can do) swaps R and D there, keeps H.
  G <- dosages(BC1)
  G["m10", ] <- -G["m10", ]
  m <- FOUNDERS$map
  rev_pop <- as_population(cbind(data.frame(snp = m$snp, allele = "A/G",
                                            chr = m$chr, pos = m$pos, cm = m$cm,
                                            stringsAsFactors = FALSE),
                                 as.data.frame(G)))
  code <- function(p) {
    attr(mabc_select(p, REC, DON, target_markers = "m10", n = 1, seed = 1),
         "mabc")$target
  }
  a <- code(BC1); b <- code(rev_pop)
  expect_identical(b, unname(c(R = "D", H = "H", D = "R")[a]))
})

test_that("exclude_interval's length leaves interval-weighted recovery", {
  # review round 2, O3: markers at 0, 10, 20, 100 cM with [10, 20] excluded.
  # Cells of the background markers 0 and 100 are [0, 50] and [50, 100]; the
  # excluded 10 cM is removed from the first, not bridged: weights 40 and 50.
  map <- data.frame(snp = paste0("x", 1:4), chr = 1, pos = 1:4,
                    cm = c(0, 10, 20, 100), stringsAsFactors = FALSE)
  iv <- .mabc_parse_interval(list(chr = 1, from = 10, to = 20), map)
  expect_equal(.mabc_weights(map, c(1L, 4L), "interval"), c(50, 50))
  expect_equal(.mabc_weights(map, c(1L, 4L), "interval", iv), c(40, 50))
  # overlapping intervals are merged, never subtracted twice
  iv2 <- .mabc_parse_interval(data.frame(chr = 1, from = c(10, 5),
                                         to = c(20, 15)), map)
  expect_equal(nrow(iv2), 1L)
  expect_equal(.mabc_weights(map, c(1L, 4L), "interval", iv2), c(35, 50))
  # end to end: two candidates, x2 the target (inside the interval)
  g <- data.frame(snp = map$snp, allele = "A/G", chr = 1, pos = map$pos,
                  cm = map$cm, P1 = 1L, P2 = -1L, C1 = c(1L, 0L, 0L, -1L),
                  C2 = c(-1L, 0L, 0L, 0L), stringsAsFactors = FALSE)
  p <- as_population(g)
  sel <- mabc_select(p[3:4], p[1], p[2], target_markers = "x2",
                     exclude_interval = list(chr = 1, from = 10, to = 20),
                     marker_weights = "interval", n = 1, seed = 1)
  d <- attr(sel, "mabc")
  expect_equal(d$background_recovery, c(40 / 90, 25 / 90))
  expect_identical(attr(sel, "mabc_info")$background_markers, c("x1", "x4"))
})

test_that("interval cells run to the chromosome's mapped ends (round-3 O1)", {
  # chr1: 0, 50, 100 with the terminal marker (a target) removed; chr2: 0, 50, 100.
  # Each chromosome still counts its full 100 cM.
  map <- data.frame(snp = paste0("y", 1:6), chr = rep(1:2, each = 3), pos = 1:6,
                    cm = rep(c(0, 50, 100), 2), stringsAsFactors = FALSE)
  w <- .mabc_weights(map, c(1L, 2L, 4L, 5L, 6L), "interval")
  expect_equal(w, c(25, 75, 25, 50, 25))
  expect_equal(sum(w[1:2]), sum(w[3:5]))
  # one scored marker on a chromosome carries that whole chromosome
  expect_equal(.mabc_weights(map, c(2L, 4L, 5L, 6L), "interval"),
               c(100, 25, 50, 25))
  # zero-length chromosome: zero weight, warned; all zero -> equal weights
  flat <- rbind(map, data.frame(snp = "y7", chr = 3, pos = 7, cm = 0))
  expect_warning(w3 <- .mabc_weights(flat, c(1L, 7L), "interval"), "span 0 cM")
  expect_equal(w3, c(100, 0))
  expect_equal(.mabc_weights(flat, 7L, "interval"), 1)
})

test_that("weights are validated even when no background remains (round-3 O2)", {
  g <- data.frame(snp = "z1", allele = "A/G", chr = 1, pos = 1, cm = 0,
                  P1 = 1L, P2 = -1L, C = 0L, stringsAsFactors = FALSE)
  p <- as_population(g)
  expect_error(suppressWarnings(
    mabc_select(p[3], p[1], p[2], target_markers = "z1",
                marker_weights = "bogus")), "must be NULL")
  expect_error(suppressWarnings(
    mabc_select(p[3], p[1], p[2], target_markers = "z1",
                marker_weights = c(1, 2))), "one value per map marker")
})

test_that("exclude_interval chromosomes are validated (round-4 O2)", {
  # a missing or unknown chr would silently exclude nothing
  expect_error(mabc_select(BC1, REC, DON, target_markers = "m10",
                           exclude_interval = list(chr = NA, from = 0, to = 10)),
               "`chr` must not be missing")
  expect_error(mabc_select(BC1, REC, DON, target_markers = "m10",
                           exclude_interval = list(chr = 7, from = 0, to = 10)),
               "chromosome\\(s\\) 7 not in the map")
})

test_that("founders on a different map are rejected", {
  g <- .bc_geno(); g$cm <- g$cm * 2
  other <- as_population(g)
  expect_error(recurrent_parent_recovery(BC1, other[1], other[2]),
               "same marker map")
})

# ---- audit additions (reconciliation v2-ocs-usefulness-marker) ---------------

test_that("ranking is lexicographic: flank count outranks background recovery", {
  # A: both flanks recurrent-homozygous but a donor background (recovery 0);
  # B: no recurrent flank but a fully recurrent background (recovery 1).
  g <- data.frame(snp = paste0("m", 1:6), allele = "A/G",
                  chr = c(1, 1, 1, 2, 2, 2), pos = 1:6, cm = c(0, 10, 20, 0, 10, 20),
                  REC = 1L, DON = -1L, A = c(1L, 0L, 1L, -1L, -1L, -1L),
                  B = c(0L, 0L, 0L, 1L, 1L, 1L), stringsAsFactors = FALSE)
  pop <- as_population(g)
  sel <- mabc_select(pop[3:4], pop[1], pop[2], target_markers = "m2",
                     flanking_markers = c("m1", "m3"), n = 1, seed = 1)
  d <- attr(sel, "mabc")
  expect_identical(sel$ids, "A")
  expect_equal(d$flank_recurrent[d$id == "A"], 2L)
  expect_equal(d$background_recovery[d$id == "B"], 1)
  expect_equal(d$rank[d$id == "A"], 1L)
})

test_that("recurrent_parent_recovery() warns and stays NA with no informative marker", {
  g <- data.frame(snp = paste0("m", 1:3), allele = "A/G", chr = 1:3, pos = 1,
                  cm = 0, P1 = 1L, P2 = 1L, stringsAsFactors = FALSE)
  fnd <- as_population(g)
  expect_warning(r <- recurrent_parent_recovery(fnd, fnd[1], fnd[2]), "informative")
  expect_true(all(is.na(r)))
})

test_that("a zero-weight marker is identical to dropping the marker", {
  w <- rep(1, 82); w[5] <- 0
  r1 <- recurrent_parent_recovery(BC1, REC, DON, weights = w)
  r2 <- recurrent_parent_recovery(BC1, REC, DON, markers = setdiff(1:82, 5))
  expect_equal(as.numeric(r1), as.numeric(r2))   # ignore the n_markers attribute
})

test_that("flanking selection shortens the donor segment around the target", {
  # length (cM) of the contiguous non-recurrent run around target m10 on chr 1
  seg_len <- function(pop, tgt = 10L) {
    S <- (dosages(pop)[1:41, , drop = FALSE] * dosages(REC)[1:41, 1] + 1) / 2
    cm <- FOUNDERS$map$cm[1:41]
    apply(S, 2L, function(s) {
      lo <- tgt; while (lo > 1L && s[lo - 1L] < 1) lo <- lo - 1L
      hi <- tgt; while (hi < 41L && s[hi + 1L] < 1) hi <- hi + 1L
      cm[hi] - cm[lo]
    })
  }
  iv <- list(chr = 1, from = 0, to = 50)
  bg <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 10, seed = 1,
                    exclude_interval = iv)
  fl <- mabc_select(BC1, REC, DON, target_markers = "m10", n = 10, seed = 1,
                    flanking_markers = c("m6", "m14"), exclude_interval = iv)
  feas <- BC1[which(attr(bg, "mabc")$feasible)]
  expect_lt(mean(seg_len(fl)), mean(seg_len(bg)))
  expect_lt(mean(seg_len(bg)), mean(seg_len(feas)))
})

test_that("two rounds of foreground + background selection beat the unselected cohort", {
  rr <- function(p) mean(recurrent_parent_recovery(p, REC, DON))
  sel1 <- mabc_select(BC1, REC, DON, target_markers = "m20", n = 10, seed = 5)
  BC2s <- do.call(c, lapply(1:10, function(i)
    cross(sel1[i], REC, n = 30, seed = 7000 + i)))
  sel2 <- mabc_select(BC2s, REC, DON, target_markers = "m20", n = 10, seed = 6)
  BC3s <- do.call(c, lapply(1:10, function(i)
    cross(sel2[i], REC, n = 30, seed = 8000 + i)))
  expect_gt(rr(BC2s), rr(BC2) + 0.05)
  expect_gt(rr(BC3s), rr(BC3) + 0.03)
})
