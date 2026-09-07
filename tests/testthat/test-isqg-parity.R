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
