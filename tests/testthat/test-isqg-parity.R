# test-isqg-parity.R
#
# isqg parity gate (DECISION-012).
#
# The Rust meiosis core must reproduce isqg 1.4's output EXACTLY, not merely
# with the same distribution. Each fixture in inst/extdata/isqg_v1_outputs/
# carries both isqg's output and the random draws that produced it, so these
# tests feed the recorded draws straight into Rust with no RNG in play: a
# failure is unambiguously a bug in the port.
#
# Regenerate fixtures only via capture_isqg_references.R, and only deliberately.

# Resolve against the installed package first so the gate also runs under
# R CMD check, where inst/extdata/ has been flattened to extdata/.
ISQG_DIR <- local({
  installed <- system.file("extdata", "isqg_v1_outputs",
                           package = "simplePHENOTYPES")
  if (nzchar(installed) && dir.exists(installed)) {
    installed
  } else {
    testthat::test_path("..", "..", "inst", "extdata", "isqg_v1_outputs")
  }
})

.isqg_ref <- function(name) readRDS(file.path(ISQG_DIR, paste0(name, ".rds")))

.have <- function(name) file.exists(file.path(ISQG_DIR, paste0(name, ".rds")))

# --- Fixture -> Rust argument marshalling -----------------------------------

# Chromosomes ascending, positions ascending within chromosome.
.layout <- function(ref) {
  m <- ref$map[order(ref$map$chr, ref$map$pos), ]
  list(
    loci_per_chr = as.integer(unname(table(m$chr))),
    positions    = as.numeric(m$pos),
    snp          = as.character(m$snp)
  )
}

# Flatten the recorded draws into the (event, chromosome) order the Rust core
# expects: counts[event * n_chr + chr], chiasmata concatenated in that order.
.flatten_draws <- function(draws) {
  list(
    counts    = as.integer(unlist(lapply(draws, `[[`, "counts"))),
    flips     = as.integer(unlist(lapply(draws, `[[`, "flips"))),
    chiasmata = as.numeric(unlist(lapply(draws, function(e) unlist(e$chiasmata))))
  )
}

.bits <- function(x) paste(as.integer(x), collapse = "")

# Run the Rust kernel for one captured mating scenario.
.rust_genotype <- function(ref, design, p1, p2) {
  lay <- .layout(ref)
  dr  <- .flatten_draws(ref$draws)
  v <- simplePHENOTYPES:::meiosis_core(
    loci_per_chr = lay$loci_per_chr,
    positions    = lay$positions,
    p1_cis       = .bits(p1$cis),
    p1_trans     = .bits(p1$trans),
    p2_cis       = .bits(p2$cis),
    p2_trans     = .bits(p2$trans),
    chiasmata    = dr$chiasmata,
    counts       = dr$counts,
    flips        = dr$flips,
    design       = design,
    n_prog       = as.integer(ref$n_prog)
  )
  # Rust returns loci-major, so read back with byrow = TRUE.
  matrix(as.integer(v), nrow = length(lay$positions), ncol = ref$n_prog,
         byrow = TRUE)
}

# Shared assertions for every mating scenario.
.expect_matches_isqg <- function(scenario, design, p1_key, p2_key = p1_key) {
  skip_if_not(.have(scenario), paste0("isqg fixture '", scenario, "' not captured"))
  ref <- .isqg_ref(scenario)
  obs <- .rust_genotype(ref, design, ref$founders[[p1_key]], ref$founders[[p2_key]])

  # The map order we send to Rust must be the order isqg reported its rows in;
  # otherwise the comparison below would be self-consistent but wrong.
  expect_identical(rownames(ref$genotype), .layout(ref)$snp)

  expect_identical(dim(obs), dim(ref$genotype))
  expect_identical(obs, unname(ref$genotype))
  invisible(ref)
}

# ---------------------------------------------------------------------------
# 1. Gamete masks — the recombination algorithm in isolation
# ---------------------------------------------------------------------------
# Highest-signal target: a raw ancestry mask, with no haplotype merging and no
# lossy -1/0/1 projection. This is the test that can detect a cis/trans swap.

test_that("Rust reproduces isqg gamete masks exactly", {
  skip_if_not(.have("gamete_masks"), "isqg fixture 'gamete_masks' not captured")
  ref <- .isqg_ref("gamete_masks")
  lay <- .layout(ref)
  dr  <- .flatten_draws(ref$draws)

  obs <- simplePHENOTYPES:::gamete_masks_core(
    loci_per_chr = lay$loci_per_chr,
    positions    = lay$positions,
    chiasmata    = dr$chiasmata,
    counts       = dr$counts,
    flips        = dr$flips
  )

  expect_identical(obs, ref$masks)
})

test_that("the captured gamete fixture exercises the edge cases it was built for", {
  skip_if_not(.have("gamete_masks"), "isqg fixture 'gamete_masks' not captured")
  ref <- .isqg_ref("gamete_masks")
  counts <- do.call(rbind, lapply(ref$draws, `[[`, "counts"))

  # Zero-crossover chromosomes still consume a Bernoulli; if the fixture never
  # contained one, the most likely desynchronisation bug would go untested.
  expect_true(any(counts == 0))
  # And chiasmata upstream of chromosome 1's first marker (0.15) give
  # breaks == 0, the whole-chromosome toggle.
  first_pos <- min(ref$map$pos[ref$map$chr == 1])
  upstream <- unlist(lapply(ref$draws, function(e) e$chiasmata[[1]]))
  expect_true(any(upstream < first_pos))
})

# ---------------------------------------------------------------------------
# 2. Mating designs
# ---------------------------------------------------------------------------

test_that("Rust reproduces isqg cross() genotypes exactly", {
  .expect_matches_isqg("cross", "cross", "P1", "P2")
})

test_that("Rust reproduces isqg cross() with parents swapped", {
  # cross(P2, P1) is not cross(P1, P2): parents are consumed in order and cis
  # always comes from the first, so this is what makes a phase swap fail.
  .expect_matches_isqg("cross_swapped", "cross", "P2", "P1")
})

test_that("cross and cross_swapped really do differ", {
  skip_if_not(.have("cross") && .have("cross_swapped"), "isqg fixtures not captured")
  expect_false(identical(.isqg_ref("cross")$genotype,
                         .isqg_ref("cross_swapped")$genotype))
})

test_that("Rust reproduces isqg selfcross() genotypes exactly", {
  .expect_matches_isqg("selfcross", "selfcross", "P1", "P1")
})

test_that("Rust reproduces isqg dh() genotypes exactly", {
  ref <- .expect_matches_isqg("dh", "dh", "P1", "P1")
  # Doubling a single gamete leaves no heterozygotes at all.
  expect_false(any(ref$genotype == 0L))
})

# ---------------------------------------------------------------------------
# 2b. The production path: R's draw order and the haplotype kernel
# ---------------------------------------------------------------------------
# The tests above feed isqg's recorded draws into meiosis_core(), so they pin
# the Rust consumption of the draws but not the other half of DECISION-012:
# that R draws them in isqg's order. The first test below pins the draw order;
# the second pins mate_haplotypes_core(), the string-strand kernel, on isqg's
# phased output. Neither is what cross()/selfcross()/double_haploid()/mate()
# run today: they go through .mate_many() -> mate_many_core(), which section
# 2c pins directly.

test_that(".draw_meiosis() reproduces isqg's draws under the fixture seed", {
  for (scenario in c("gamete_masks", "cross", "cross_swapped", "selfcross",
                     "dh")) {
    skip_if_not(.have(scenario), paste0("isqg fixture '", scenario,
                                        "' not captured"))
    ref <- .isqg_ref(scenario)
    lay <- .layout(ref)
    by_chr <- split(lay$positions, rep(seq_along(lay$loci_per_chr),
                                       lay$loci_per_chr))
    obs <- withr::with_seed(ref$seed, simplePHENOTYPES:::.draw_meiosis(
      unname(by_chr), length(ref$draws)))
    exp <- .flatten_draws(ref$draws)
    expect_identical(obs$counts, exp$counts, info = scenario)
    expect_identical(obs$flips, exp$flips, info = scenario)
    expect_identical(obs$chiasmata, exp$chiasmata, info = scenario)
  }
})

test_that("mate_haplotypes_core() reproduces isqg's phased progeny exactly", {
  designs <- list(cross = c("cross", "P1", "P2"),
                  cross_swapped = c("cross", "P2", "P1"),
                  selfcross = c("selfcross", "P1", "P1"),
                  dh = c("dh", "P1", "P1"))
  for (scenario in names(designs)) {
    skip_if_not(.have(scenario), paste0("isqg fixture '", scenario,
                                        "' not captured"))
    ref <- .isqg_ref(scenario)
    d <- designs[[scenario]]
    p1 <- ref$founders[[d[2]]]
    p2 <- ref$founders[[d[3]]]
    lay <- .layout(ref)
    dr <- .flatten_draws(ref$draws)
    strands <- simplePHENOTYPES:::mate_haplotypes_core(
      loci_per_chr = lay$loci_per_chr,
      positions    = lay$positions,
      p1_cis       = .bits(p1$cis),
      p1_trans     = .bits(p1$trans),
      p2_cis       = .bits(p2$cis),
      p2_trans     = .bits(p2$trans),
      chiasmata    = dr$chiasmata,
      counts       = dr$counts,
      flips        = dr$flips,
      design       = d[1],
      n_prog       = as.integer(ref$n_prog)
    )
    # Element 2i-1 is progeny i's cis strand, 2i its trans strand.
    unpack <- function(codes) {
      matrix(as.integer(unlist(strsplit(codes, "", fixed = TRUE))),
             nrow = length(lay$positions))
    }
    cis <- unpack(strands[seq(1, length(strands), by = 2)])
    trans <- unpack(strands[seq(2, length(strands), by = 2)])
    geno <- ifelse(cis & trans, 1L, ifelse(xor(cis, trans), 0L, -1L))
    # isqg's phased code is "<cis> <trans>", allele 1 for bit 1: unlike the
    # -1/0/1 genotype, it tells a cis/trans swap apart.
    phased <- ifelse(cis & trans, "1 1", ifelse(!cis & !trans, "2 2",
                     ifelse(cis == 1L, "1 2", "2 1")))
    expect_identical(geno, unname(ref$genotype), info = scenario)
    expect_identical(phased, unname(ref$genotype_phased), info = scenario)
  }
})

# ---------------------------------------------------------------------------
# 2c. The production path: mate_many_core() and the exported functions
# ---------------------------------------------------------------------------
# cross(), selfcross(), double_haploid() and mate() all run .mate_many() ->
# mate_many_core() (DECISION-040), an integer-strand re-implementation of the
# gamete selection that until here was only equivalence-tested against
# mate_haplotypes_core() (test-feat-crossing.R). These tests pin it to isqg
# itself: first the Rust core fed the recorded draws, then the exported
# functions under the fixture seed, i.e. the whole chain
# .draw_meiosis() -> .mate_many() -> mate_many_core() -> Population.

# The four captured mating scenarios: design, p1, p2.
.isqg_designs <- list(cross = c("cross", "P1", "P2"),
                      cross_swapped = c("cross", "P2", "P1"),
                      selfcross = c("selfcross", "P1", "P1"),
                      dh = c("dh", "P1", "P1"))

# isqg's phased code is "<cis> <trans>", allele 1 for bit 1: unlike the
# -1/0/1 genotype, it tells a cis/trans swap apart.
.isqg_phased <- function(cis, trans) {
  ifelse(cis & trans, "1 1", ifelse(!cis & !trans, "2 2",
         ifelse(cis == 1L, "1 2", "2 1")))
}

# The parental strand matrix and mating row .mate() builds for a design.
.isqg_strands <- function(design, p1, p2) {
  if (design == "cross") {
    list(strands = cbind(p1$cis, p1$trans, p2$cis, p2$trans),
         mating = 1:4)
  } else {
    list(strands = cbind(p1$cis, p1$trans), mating = c(1L, 2L, 1L, 2L))
  }
}

test_that("mate_many_core() reproduces isqg's phased progeny exactly", {
  for (scenario in names(.isqg_designs)) {
    skip_if_not(.have(scenario), paste0("isqg fixture '", scenario,
                                        "' not captured"))
    ref <- .isqg_ref(scenario)
    d <- .isqg_designs[[scenario]]
    lay <- .layout(ref)
    dr <- .flatten_draws(ref$draws)
    n_loci <- length(lay$positions)
    st <- .isqg_strands(d[1], ref$founders[[d[2]]], ref$founders[[d[3]]])
    strands <- st$strands
    storage.mode(strands) <- "integer"

    # Markers in isqg's own (chr, pos) order: `order` is the identity.
    res <- simplePHENOTYPES:::mate_many_core(
      loci_per_chr = lay$loci_per_chr,
      positions    = lay$positions,
      order        = seq_len(n_loci),
      strands      = strands,
      n_strands    = ncol(strands),
      mating       = as.integer(st$mating),
      design       = d[1],
      n_prog       = as.integer(ref$n_prog),
      chiasmata    = dr$chiasmata,
      counts       = dr$counts,
      flips        = dr$flips
    )
    cis <- matrix(as.integer(res$cis), nrow = n_loci)
    trans <- matrix(as.integer(res$trans), nrow = n_loci)
    expect_identical(dim(cis), dim(ref$genotype), info = scenario)
    expect_identical(cis + trans - 1L, unname(ref$genotype), info = scenario)
    expect_identical(.isqg_phased(cis, trans), unname(ref$genotype_phased),
                     info = scenario)

    # The same mating with the markers handed over in a scrambled order and
    # `order` the permutation .meiosis_layout() computes for that map, as
    # .mate_many() does: the progeny must come back in the caller's order and
    # still be isqg's. This pins the permutation path on the fixtures.
    perm <- order(seq_len(n_loci) %% 7L, -seq_len(n_loci))
    map <- data.frame(snp = lay$snp[perm], chr = ref$map$chr[perm],
                      pos = ref$map$pos[perm], cm = lay$positions[perm] * 100)
    play <- simplePHENOTYPES:::.meiosis_layout(map)
    expect_identical(play$positions, lay$positions, info = scenario)
    expect_identical(play$loci_per_chr, lay$loci_per_chr, info = scenario)
    expect_identical(map$snp[play$ord], lay$snp, info = scenario)
    res2 <- simplePHENOTYPES:::mate_many_core(
      loci_per_chr = play$loci_per_chr,
      positions    = play$positions,
      order        = play$ord,
      strands      = strands[perm, , drop = FALSE],
      n_strands    = ncol(strands),
      mating       = as.integer(st$mating),
      design       = d[1],
      n_prog       = as.integer(ref$n_prog),
      chiasmata    = dr$chiasmata,
      counts       = dr$counts,
      flips        = dr$flips
    )
    cis2 <- matrix(as.integer(res2$cis), nrow = n_loci)
    trans2 <- matrix(as.integer(res2$trans), nrow = n_loci)
    expect_identical(cis2[play$ord, ], cis, info = scenario)
    expect_identical(trans2[play$ord, ], trans, info = scenario)
  }
})

# A Population of the fixture founders. The fixture map is in Morgans; the
# package map wants centiMorgans (`cm`, divided by 100 by .meiosis_layout())
# and a non-negative physical `pos`, which plays no part in meiosis.
.isqg_population <- function(ref) {
  m <- ref$map[order(ref$map$chr, ref$map$pos), ]
  map <- data.frame(snp = as.character(m$snp), chr = m$chr,
                    pos = round(m$pos * 1e6), cm = m$pos * 100)
  cis <- cbind(P1 = ref$founders$P1$cis, P2 = ref$founders$P2$cis)
  trans <- cbind(P1 = ref$founders$P1$trans, P2 = ref$founders$P2$trans)
  suppressMessages(population_from_haplotypes(cis, trans, map))
}

test_that("cross()/selfcross()/double_haploid() reproduce isqg under the fixture seed", {
  for (scenario in names(.isqg_designs)) {
    skip_if_not(.have(scenario), paste0("isqg fixture '", scenario,
                                        "' not captured"))
    ref <- .isqg_ref(scenario)
    d <- .isqg_designs[[scenario]]
    pop <- .isqg_population(ref)
    # cross()/selfcross() run two meioses per progeny, double_haploid() one;
    # the fixture holds exactly that many events, so one seed covers both.
    events_per <- if (d[1] == "dh") 1L else 2L
    expect_identical(length(ref$draws), events_per * ref$n_prog,
                     info = scenario)

    run <- switch(d[1],
      cross     = function() cross(pop[d[2]], pop[d[3]], n = ref$n_prog),
      selfcross = function() selfcross(pop[d[2]], n = ref$n_prog),
      dh        = function() double_haploid(pop[d[2]], n = ref$n_prog))
    prog <- suppressMessages(withr::with_seed(ref$seed, run()))

    geno <- dosages(prog)
    expect_identical(rownames(geno), rownames(ref$genotype), info = scenario)
    expect_identical(unname(geno), unname(ref$genotype), info = scenario)
    h <- haplotypes(prog)
    expect_identical(unname(.isqg_phased(h$cis, h$trans)),
                     unname(ref$genotype_phased), info = scenario)

    # The exported function consumes exactly isqg's draws and nothing else
    # before or after them: the RNG state it leaves equals the state after
    # .draw_meiosis() alone (set.seed() here, not with_seed(), because the
    # latter restores the state we want to read).
    lay <- .layout(ref)
    by_chr <- split(lay$positions, rep(seq_along(lay$loci_per_chr),
                                       lay$loci_per_chr))
    set.seed(ref$seed)
    suppressMessages(run())
    after_run <- .Random.seed
    set.seed(ref$seed)
    simplePHENOTYPES:::.draw_meiosis(unname(by_chr), length(ref$draws))
    expect_identical(after_run, .Random.seed, info = scenario)
  }
})

# ---------------------------------------------------------------------------
# 3. Structural guards
# ---------------------------------------------------------------------------
# Commit 4167402 was a row/column-major mixup masked by uniform test values.
# These assertions are chosen so a transposition cannot survive them: the
# fixture is 37 x 5 (coprime, so a transpose cannot even be reshaped), and a
# single asymmetric cell is checked directly.

test_that("genotype output is loci-major and cannot be silently transposed", {
  skip_if_not(.have("cross"), "isqg fixture 'cross' not captured")
  ref <- .isqg_ref("cross")
  obs <- .rust_genotype(ref, "cross", ref$founders$P1, ref$founders$P2)

  expect_identical(dim(obs), c(37L, 5L))
  # A specific off-diagonal cell, hand-tied to isqg's own matrix. Row sums,
  # column sums and sorted values would all survive a transpose; this does not.
  expect_identical(obs[3, 5], unname(ref$genotype)[3, 5])
  expect_identical(obs[30, 2], unname(ref$genotype)[30, 2])
})

test_that("founder encoding round-trips through the recorded bit patterns", {
  skip_if_not(.have("founders"), "isqg fixture 'founders' not captured")
  ref <- .isqg_ref("founders")
  decode <- function(cis, trans) {
    ifelse(cis & trans, 1L, ifelse(xor(cis, trans), 0L, -1L))
  }
  expected <- cbind(
    P1 = decode(ref$founders$P1$cis, ref$founders$P1$trans),
    P2 = decode(ref$founders$P2$cis, ref$founders$P2$trans)
  )
  expect_identical(unname(ref$genotype), unname(expected))
  # The founders must not be trivially homozygous-opposite, or every progeny
  # locus would be heterozygous and most indexing errors would be invisible.
  expect_true(any(ref$genotype[, "P1"] == 0L))
})
