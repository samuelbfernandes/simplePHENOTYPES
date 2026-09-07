# capture_isqg_references.R
#
# Run this script ONCE before implementing the Rust meiosis port.
# It installs isqg 1.4 from the in-repo source into a dedicated local library,
# runs each mating scenario under a fixed seed, and saves — for every scenario —
# BOTH isqg's output AND the random draws that produced it, as RDS files in this
# directory.
#
# These RDS files gate the isqg parity tests (test-isqg-parity.R). Under
# DECISION-012 the gate is EXACT bit-parity, so the Rust core must reproduce
# these outputs byte-for-byte from the recorded draws.
#
# Why the draws are recorded alongside the outputs: without them a failing test
# cannot distinguish an error in R's replication of isqg's draw ORDER from an
# error in the Rust kernel. With them, the Rust side is a pure-function test.
#
# Why isqg is installed from context/isqg rather than with remotes: isqg was
# archived from CRAN, so remotes::install_version() cannot resolve it. The
# source is vendored in-repo (DECISION-002 — we own this algorithm).
#
# Usage (from package root):
#   Rscript inst/extdata/isqg_v1_outputs/capture_isqg_references.R
#
# Requirements: Rcpp, R6, Rdpack, BH, and a C++14 compiler.

stopifnot(
  requireNamespace("Rcpp",  quietly = TRUE),
  requireNamespace("R6",    quietly = TRUE),
  requireNamespace("Rdpack", quietly = TRUE),
  requireNamespace("BH",    quietly = TRUE)
)

# ---------------------------------------------------------------------------
# 1. Install isqg 1.4 from the in-repo source into a dedicated library
# ---------------------------------------------------------------------------

lib_isqg <- path.expand("~/.R/isqg_v14_lib")
dir.create(lib_isqg, recursive = TRUE, showWarnings = FALSE)

installed_ver <- tryCatch(
  packageDescription("isqg", lib.loc = lib_isqg)$Version,
  error = function(e) NA_character_
)

if (is.na(installed_ver) || installed_ver != "1.4") {
  src <- normalizePath(file.path("context", "isqg"), mustWork = FALSE)
  if (!dir.exists(src)) {
    stop("isqg source not found at context/isqg. Run from the package root.",
         call. = FALSE)
  }
  message("Installing isqg 1.4 from ", src, " into ", lib_isqg, " ...")
  res <- system2("R", c("CMD", "INSTALL", paste0("--library=", shQuote(lib_isqg)),
                        "--no-docs", shQuote(src)),
                 stdout = FALSE, stderr = FALSE)
  if (res != 0L) stop("R CMD INSTALL of context/isqg failed.", call. = FALSE)
} else {
  message("isqg 1.4 already present at ", lib_isqg)
}

isqg_ns <- loadNamespace("isqg", lib.loc = lib_isqg)
attach(isqg_ns, name = "isqg_capture", warn.conflicts = FALSE)
on.exit(try(detach("isqg_capture"), silent = TRUE), add = TRUE)

# ---------------------------------------------------------------------------
# 2. Resolve output directory (where RDS files are saved)
# ---------------------------------------------------------------------------

this_file <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  f    <- sub("--file=", "", args[grepl("--file=", args)])
  if (length(f) && nzchar(f)) normalizePath(f) else stop()
}, error = function(e) {
  normalizePath("inst/extdata/isqg_v1_outputs/capture_isqg_references.R")
})
out_dir <- dirname(this_file)

# ---------------------------------------------------------------------------
# 3. The asymmetric test map
# ---------------------------------------------------------------------------
#
# Deliberately built so structural bugs in the port cannot hide:
#   * 3 chromosomes with UNEQUAL loci counts (11, 1, 25) = 37 loci total, and
#     n_prog = 5. 37 and 5 are coprime, so a transposed result cannot be
#     reshaped into the right dimensions.
#   * 37 is not a multiple of 64, exercising the trailing-word tail of a
#     bit-packed representation.
#   * Chromosome 2 has a SINGLE locus (n = 1 edge case).
#   * Chromosome 1 starts at a NONZERO position (0.15), so chiasmata drawn on
#     (0, L) can land upstream of the first marker; that gives breaks == 0,
#     which toggles the entire chromosome.
#   * Chromosome 3 is long (L = 2.4 Morgans), so several crossovers per meiosis
#     are common.
# Positions are strictly increasing within each chromosome: isqg mis-assigns bit
# indices for TIED positions (Functions.R:85 computes an order, not a reversed
# rank), so tied maps must never be used as a parity reference.

amap <- rbind(
  data.frame(chr = 1L, pos = round(seq(0.15, 1.35, length.out = 11), 4)),
  data.frame(chr = 2L, pos = 0.8),
  data.frame(chr = 3L, pos = round(seq(0.00, 2.40, length.out = 25), 4))
)
amap$snp <- sprintf("m%02d_c%d", seq_len(nrow(amap)), amap$chr)
amap <- amap[, c("snp", "chr", "pos")]

N_PROG <- 5L

spc  <- set_specie(amap)
chrs <- split(amap$pos, amap$chr)   # ascending positions, per chromosome

# ---------------------------------------------------------------------------
# 4. Helpers
# ---------------------------------------------------------------------------

# --- 4a. Founder haplotypes -------------------------------------------------
#
# Non-periodic, non-palindromic bit patterns rather than isqg's "AA"/"aa"
# founders. Two opposite homozygotes would make EVERY progeny locus
# heterozygous, which hides both cis/trans swaps and most indexing errors.
# Codes follow isqg: 1 = allele A, "cis trans", so "1 2" is cis = A, trans = a.

.bits_to_code <- function(cis, trans) {
  paste(ifelse(cis == 1L, "1", "2"), ifelse(trans == 1L, "1", "2"))
}

# Digits of an irrational-ish constant give an aperiodic pattern with no
# symmetry a transposition or reversal could preserve.
.pattern <- function(n, offset) {
  digits <- as.integer(strsplit(
    "14159265358979323846264338327950288419716939937510582097494459230781640628",
    "")[[1]])
  as.integer(digits[((seq_len(n) + offset - 1L) %% length(digits)) + 1L] %% 2L)
}

P1_cis   <- .pattern(nrow(amap), 0L)
P1_trans <- .pattern(nrow(amap), 17L)
P2_cis   <- .pattern(nrow(amap), 31L)
P2_trans <- .pattern(nrow(amap), 53L)

make_founder <- function(cis, trans) {
  code <- .bits_to_code(cis, trans)
  names(code) <- amap$snp
  import(spc, code)
}

P1 <- make_founder(P1_cis, P1_trans)
P2 <- make_founder(P2_cis, P2_trans)

# --- 4b. Replication of isqg's RNG stream -----------------------------------
#
# isqg draws, per whole-genome meiosis event, per chromosome in ascending order
# (context/isqg/src/Genetics.cpp:54-101):
#     n_x       ~ rpois(1, L)              L = LAST map position
#     chiasmata ~ sort(runif(n_x, 0, L))   NOT drawn when n_x == 0
#     flip      ~ rbinom(1, 1, 0.5)        ALWAYS drawn, even when n_x == 0
# The flip is unconditional; skipping it when n_x == 0 desynchronises every
# subsequent draw.

draw_event <- function() {
  n_chr  <- length(chrs)
  counts <- integer(n_chr)
  flips  <- integer(n_chr)
  chias  <- vector("list", n_chr)
  for (i in seq_len(n_chr)) {
    pos <- chrs[[i]]
    L   <- pos[length(pos)]
    k   <- rpois(1, L)
    chias[[i]] <- if (k > 0) sort(runif(k, 0, L)) else numeric(0)
    counts[i]  <- k
    flips[i]   <- rbinom(1, 1, 0.5)
  }
  list(counts = counts, chiasmata = chias, flips = flips)
}

draw_events <- function(n) lapply(seq_len(n), function(i) draw_event())

# Ancestry mask for one event: bit 1 => take the locus from `cis`.
event_mask <- function(ev) {
  unlist(lapply(seq_along(chrs), function(i) {
    pos <- chrs[[i]]
    n   <- length(pos)
    m   <- integer(n)
    for (x in ev$chiasmata[[i]]) {
      b <- sum(pos <= x)            # upper_bound: count of positions <= chiasma
      if (b < n) m[(b + 1L):n] <- 1L - m[(b + 1L):n]
    }
    if (ev$flips[i] == 1L) m <- 1L - m
    m
  }), use.names = FALSE)
}

gamete  <- function(cis, trans, mask) ifelse(mask == 1L, cis, trans)
decode  <- function(cis, trans) ifelse(cis & trans, 1L,
                                       ifelse(xor(cis, trans), 0L, -1L))
decode_phased <- function(cis, trans) {
  ifelse(cis & trans, "1 1",
         ifelse(!cis & !trans, "2 2", ifelse(cis == 1L, "1 2", "2 1")))
}

# Predict a whole mating design from pre-drawn events, exactly as the Rust core
# will: cross/selfcross consume two events per progeny (parent 1 then parent 2,
# cis from parent 1); dh consumes one and duplicates the gamete.
predict_design <- function(design, events, p1, p2 = p1) {
  per <- if (design == "dh") 1L else 2L
  n   <- length(events) / per
  cis_l <- trans_l <- vector("list", n)
  for (i in seq_len(n)) {
    if (design == "dh") {
      g <- gamete(p1$cis, p1$trans, event_mask(events[[i]]))
      cis_l[[i]] <- g; trans_l[[i]] <- g
    } else {
      a <- gamete(p1$cis, p1$trans, event_mask(events[[2L * i - 1L]]))
      b <- gamete(p2$cis, p2$trans, event_mask(events[[2L * i]]))
      cis_l[[i]] <- a; trans_l[[i]] <- b
    }
  }
  list(
    num    = mapply(decode,        cis_l, trans_l),
    phased = mapply(decode_phased, cis_l, trans_l),
    cis    = cis_l,
    trans  = trans_l
  )
}

# --- 4c. Scenario runner ----------------------------------------------------
#
# Runs isqg under `seed`, then replays the SAME seed through the R replication
# and asserts the two agree before saving. A capture that cannot be reproduced
# from its own recorded draws is worthless as a parity fixture, so this fails
# loudly at capture time rather than at test time.

run_scenario <- function(scenario, call_text, isqg_fn, n_events, replicate_fn) {
  message("\n--- Scenario: ", scenario, " ---")

  set.seed(SEED)
  pop <- isqg_fn()
  obs_num    <- genotype(pop)
  obs_phased <- genotype(pop, phase = TRUE)

  set.seed(SEED)
  events <- draw_events(n_events)
  pred   <- replicate_fn(events)

  exp_num <- pred$num
  dimnames(exp_num) <- dimnames(obs_num)
  if (!identical(unname(obs_num), unname(exp_num))) {
    stop("Replication of isqg's RNG stream FAILED for scenario '", scenario,
         "'. Refusing to save a fixture that cannot be reproduced.",
         call. = FALSE)
  }

  list(
    isqg_version = "1.4",
    captured     = Sys.time(),
    scenario     = scenario,
    call         = call_text,
    seed         = SEED,
    map          = amap,
    n_prog       = ncol(obs_num),
    founders     = list(
      P1 = list(cis = P1_cis, trans = P1_trans),
      P2 = list(cis = P2_cis, trans = P2_trans)
    ),
    draws        = events,
    genotype     = obs_num,
    genotype_phased = obs_phased
  )
}

SEED <- 20260907L

# ---------------------------------------------------------------------------
# 5. Reference scenarios
# ---------------------------------------------------------------------------

# -- 5a. Bare gamete masks ---------------------------------------------------
# Highest-signal target: spc$gamete() runs a full meiosis and returns the raw
# ancestry mask as a 0/1 string, with no genotype involved. It therefore
# isolates the recombination algorithm from haplotype merging, and — unlike the
# -1/0/1 genotype, which is a lossy 3-valued projection — it CAN detect a
# cis/trans swap. The string reads chromosome 1 ascending, then 2, then 3.
message("\n--- Scenario: gamete_masks ---")
set.seed(SEED)
obs_masks <- spc$gamete(n = 25L)

set.seed(SEED)
mask_events <- draw_events(25L)
exp_masks   <- vapply(mask_events,
                      function(ev) paste(event_mask(ev), collapse = ""),
                      character(1))

if (!identical(as.character(obs_masks), exp_masks)) {
  stop("Gamete-mask replication FAILED; refusing to save.", call. = FALSE)
}

saveRDS(
  list(
    isqg_version = "1.4",
    captured     = Sys.time(),
    scenario     = "gamete_masks",
    call         = "spc$gamete(n = 25)",
    seed         = SEED,
    map          = amap,
    draws        = mask_events,
    masks        = as.character(obs_masks)
  ),
  file.path(out_dir, "gamete_masks.rds")
)
message("Saved gamete_masks.rds")

# -- 5b. Founders ------------------------------------------------------------
# No meiosis: pins the founder encoding and the (cis, trans) -> -1/0/1 and
# phased decodings on their own.
message("\n--- Scenario: founders ---")
founder_pop <- list(P1 = P1, P2 = P2)
saveRDS(
  list(
    isqg_version = "1.4",
    captured     = Sys.time(),
    scenario     = "founders",
    call         = "import(spc, <hardcoded aperiodic patterns>)",
    seed         = NA_integer_,
    map          = amap,
    founders     = list(
      P1 = list(cis = P1_cis, trans = P1_trans),
      P2 = list(cis = P2_cis, trans = P2_trans)
    ),
    genotype        = genotype(founder_pop),
    genotype_phased = genotype(founder_pop, phase = TRUE)
  ),
  file.path(out_dir, "founders.rds")
)
message("Saved founders.rds")

# -- 5c. Cross (P1 x P2) -----------------------------------------------------
# 2 meioses per progeny, parent 1 then parent 2, cis from parent 1.
ref <- run_scenario(
  "cross",
  "cross(n = 5, P1, P2)",
  isqg_fn      = function() cross(n = N_PROG, P1, P2),
  n_events     = 2L * N_PROG,
  replicate_fn = function(ev) predict_design(
    "cross", ev, list(cis = P1_cis, trans = P1_trans),
    list(cis = P2_cis, trans = P2_trans))
)
saveRDS(ref, file.path(out_dir, "cross.rds"))
message("Saved cross.rds")

# -- 5d. Cross with parents SWAPPED ------------------------------------------
# cross(P2, P1) is not cross(P1, P2): the parents are consumed in order and cis
# always comes from the first. This scenario is what makes a cis/trans swap in
# the port fail loudly rather than silently.
ref <- run_scenario(
  "cross_swapped",
  "cross(n = 5, P2, P1)",
  isqg_fn      = function() cross(n = N_PROG, P2, P1),
  n_events     = 2L * N_PROG,
  replicate_fn = function(ev) predict_design(
    "cross", ev, list(cis = P2_cis, trans = P2_trans),
    list(cis = P1_cis, trans = P1_trans))
)
saveRDS(ref, file.path(out_dir, "cross_swapped.rds"))
message("Saved cross_swapped.rds")

# -- 5e. Selfcross -----------------------------------------------------------
# Structurally identical to cross, with a deep copy of the one parent as the
# second: still 2 independent meioses per progeny.
ref <- run_scenario(
  "selfcross",
  "selfcross(n = 5, P1)",
  isqg_fn      = function() selfcross(n = N_PROG, P1),
  n_events     = 2L * N_PROG,
  replicate_fn = function(ev) predict_design(
    "selfcross", ev, list(cis = P1_cis, trans = P1_trans))
)
saveRDS(ref, file.path(out_dir, "selfcross.rds"))
message("Saved selfcross.rds")

# -- 5f. Double haploid ------------------------------------------------------
# 1 meiosis per progeny, gamete duplicated: every progeny fully homozygous, so
# the genotype matrix must contain no zeros.
ref <- run_scenario(
  "dh",
  "dh(n = 5, P1)",
  isqg_fn      = function() dh(n = N_PROG, P1),
  n_events     = N_PROG,
  replicate_fn = function(ev) predict_design(
    "dh", ev, list(cis = P1_cis, trans = P1_trans))
)
stopifnot(!any(ref$genotype == 0L))
saveRDS(ref, file.path(out_dir, "dh.rds"))
message("Saved dh.rds")

# ---------------------------------------------------------------------------
# 6. Summary
# ---------------------------------------------------------------------------

message("\n", strrep("-", 70))
for (f in sort(list.files(out_dir, pattern = "\\.rds$"))) {
  x <- readRDS(file.path(out_dir, f))
  g <- if (!is.null(x$genotype)) paste(dim(x$genotype), collapse = " x ") else
       paste0(length(x$masks), " masks")
  message(sprintf("%-22s %-14s seed=%-10s %s", f, x$scenario,
                  as.character(x$seed), g))
}
message(strrep("-", 70))
message("\nDone. Commit the .rds files to inst/extdata/isqg_v1_outputs/")
message("Do NOT re-run this script after the Rust port is implemented —")
message("these fixtures are the frozen parity gate (DECISION-012).")
