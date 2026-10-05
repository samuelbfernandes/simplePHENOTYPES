# SPEC — native coalescent founders (breedingDesigner SPEC-0020 item 9)

Status: **accepted design, 2026-10-03.** The maintainer chose a full MaCS-style (sequential
Markov) coalescent over a forward burn-in or import-only, and approved C1 (DECISION-050).

## 1. Goal

`founders_coalescent()` returns a `Population` of founders with historical LD, mutation-drift
allele frequencies and a genetic map, so a breedingDesigner run with `source = "macs"` does not
need AlphaSimR. Target: AlphaSimR `runMacs()` presets and summary-statistic agreement, not bit
parity (MaCS output is not even seed-reproducible in AlphaSimR 2.1.0; BD SPEC-0020 R18).

## 2. API (proposed)

```r
founders_coalescent(n_ind, n_chr = 1, seg_sites = NULL, inbred = FALSE,
                    species = c("GENERIC", "MAIZE", "WHEAT", "CATTLE"),
                    ne = NULL, history = NULL, theta = NULL, rho = NULL,
                    chr_length_bp = NULL, chr_length_morgan = NULL,
                    split = NULL, seed = NULL)
```

- `species` presets reproduce the demographic commands `runMacs()` builds (read from AlphaSimR
  2.1.0 source): e.g. GENERIC `Ne = 100`, 1e8 bp, `-t 1E-5 -r 4E-6`, `-eN` history
  0.25/5, 2.5/15, 25/60, 250/120, 2500/1000; map length 1 M. These are parameter values, cited to
  AlphaSimR, not claims about the species.
- `history` = data frame `(time_4N, size_rel)`: piecewise-constant population size (ms `-eN`).
- `split` = generations ago that two subpopulations of equal size split (ms `-I 2` + `-ej`).
- `seg_sites`: sample this many segregating sites per chromosome (uniformly, as `runMacs`);
  `NULL` = all.
- `inbred = TRUE`: each individual = one haplotype twice (fully homozygous).
- Output via `population_from_haplotypes()`: map in cM from the position x chromosome length in
  Morgans; `pool = "coalescent"`; attributes record every parameter and the seed.

## 3. Model

- **Coalescent with recombination, SMC' approximation** (Marjoram & Wall 2006; McVean & Cardin
  2005 for SMC): walk along the chromosome, keep the marginal genealogy, at each recombination
  point detach a lineage and re-coalesce it into the current tree (SMC' allows it to coalesce
  back to its own branch). Population-size changes rescale coalescence rates piecewise.
- **MaCS extension** (Chen, Marjoram & Wall 2009): keep the trees of a retained window behind the
  current position (history parameter `h`) so a detached lineage can join branches that are
  absent from the current marginal tree. `h = 0` gives SMC'. Default: MaCS's own default,
  to be read from its source or paper before coding (not assumed here).
- **Mutation:** infinite sites, biallelic. With time in 4 N0 units a pair coalesces at rate
  2 / lambda and, for a tree of total length L, mutations and recombination points arrive along
  the [0, 1) sequence at rates theta L and rho L (theta = 4 N0 mu, rho = 4 N0 r per chromosome),
  so E[S] = theta sum_{i<n} 1/i. `seg_sites` keeps a uniform subset (reservoir sampling).
- Time in units of `4 N_e` generations, as ms/MaCS.

## 4. Decisions to take

- **C1 — RNG boundary (amends Critical Rule 4 / DECISION-006 / 012).** A coalescent draws an a-priori
  unknown number of random numbers (millions per chromosome at realistic `rho`). Drawing them in R
  and passing them to Rust is not practical, and a pure-R walk is orders of magnitude too slow for
  BD sizes (10,000 x 14,000 markers). Proposal: **Rust core with its own seeded PRNG**
  (a hand-written xoshiro256++ seeded through SplitMix64, so no new crate in the vendored
  dependency set). R draws one 32-bit seed per chromosome from its own stream (so `set.seed()` /
  `seed =` make the result reproducible; AlphaSimR's `runMacs()` is not). Rationale for the
  exception: the founder generator has no isqg reference, so it is not on the parity-critical
  path that Rule 4 protects. **Approved 2026-10-03** (DECISION-050; Rule 4 amended in `AGENTS.md`).
- **C2 — window `h`.** Expose `history_window` (default MaCS's) or fix it. Proposal: expose it.
- **C3 — map.** Linear physical-to-genetic map per chromosome (as `runMacs`), start at 0 cM.

## 5. Validation gates (analytic, no AlphaSimR in the test suite)

1. Constant size, no recombination: E[S] = theta x sum_{i=1}^{n-1} 1/i (Watterson 1975);
   site frequency spectrum E[xi_i] = theta / i; E[T_MRCA] = 2(1 - 1/n) (units 2N).
2. With recombination: E[r^2] between sites vs `rho` against Hudson's ms run on the same command
   (dev script, not shipped), and against the SMC' expectation; number of marginal-tree changes
   per unit `rho` vs the exact coalescent (SMC' is known to be close).
3. `-eN` history: E[T_MRCA] of a pair under a piecewise size history, closed form.
4. Split: F_ST between the two subpopulations vs its coalescent expectation.
5. Cross-engine evidence (dev script): summary statistics (S, SFS, LD decay, MAF) against
   AlphaSimR `runMacs()` with the same preset.

Every equation above gets a verified page before it is coded (`docs/THEORY_REVIEW.md`); the
references below are not yet page-verified.

## 6. Status

- **Phase (a) done 2026-10-04; Codex theory items all PASS (equation pages of Watterson 1975 / Fu 1995 not inspected: abstracts only):**
  `src/rust/src/coalescent.rs` (xoshiro256++ / SplitMix64, Kingman tree under `-eN` history,
  SMC' prune-and-regraft with a Fenwick tree over branch lengths, reservoir `seg_sites`) and the
  internal R wrapper `.coalescent_chromosome()` (`R/founders_coalescent.R`). Gates passed: E[T_MRCA]
  = 1 - 1/n; E[L] = sum 1/i; pair T_MRCA under a size change (with and without recombination);
  the tree at a fixed position stays Kingman under SMC' (after a fixed number of *events* it is
  length-biased, E[L^2]/E[L], which is expected); E[S] = theta a1 for rho = 0 and 20; SFS
  E[xi_i] = theta / i; seeded reproducibility. Rust: `cargo test` 64/64, clippy clean; R:
  `test-feat-coalescent.R` 34/34 (incl. Codex fixes: `history` bounded to times in (0, 1e12] and relative sizes in [1e-9, 1e9], plus runtime guards, so extreme sizes can neither overflow into an index panic nor collapse waiting times to 0; full 32-bit replayable `seed`; factor `history` refused; classed numbers such as `integer64` passed by value; non-finite `theta + rho` refused and a 2e9-event budget, checked upfront and in the loop, so huge rates error instead of hanging; references with PMIDs).
- Timing (debug build, GENERIC history, theta 1000, rho 400, 1400 kept sites): 200 / 2000 / 20000
  haplotypes = 0.1 / 0.9 / 9.9 s per chromosome. The O(n) scan for the lineages at the
  re-coalescence time is the hotspot to remove in phase (b).

- **Phase (b) 2026-10-04.** (1) *MaCS window dropped by maintainer decision:* AlphaSimR's MaCS
  source (`src/simulator.cpp:44` `dBasesToTrack = 1`; `src/datastructures.cpp:328` trailing gap
  `dBasesToTrack / dSeqLength`; `src/algorithm.cpp:1263-1271` pruning) keeps a 1-base window, and
  `runMacs()` never passes `-h`, so AlphaSimR runs effectively SMC'. The phase (a) core is that
  model; a longer window stays a possible later option. (2) *Speed:* the re-coalescence lineage
  is drawn by rejection from the root and the children of internal nodes older than tau
  (acceptance k / (2k - 1) >= 1/2; an O(log n) search plus O(1) expected draws, while keeping the
  sorted node list still shifts O(n) entries per regraft) instead of an O(n) scan per draw: 20000 haplotypes 9.9 -> 2.2 s per
  chromosome (debug build). (3) *Recombination gate:* the pair-TMRCA correlation at the two ends
  of a sequence of scaled length rho lies above the SMC value 1/(1+rho) and at or just below the
  ARG value (rho+18)/(rho^2+13rho+18) (Wilton, Carmi & Hobolth 2015, PMID 25786855): rho = 0.5 /
  1 / 3 / 10 gave 0.746 / 0.580 / 0.297 / 0.101 vs ARG 0.748 / 0.594 / 0.318 / 0.113 and SMC
  0.667 / 0.500 / 0.250 / 0.091. (4) *Cross-engine evidence* (`dev/parity-coalescent-alphasimr.R`,
  GENERIC preset, 200 haplotypes, 20 replicates each): segregating sites, mean MAF, the folded
  SFS and r^2 in seven distance bins all agree with AlphaSimR `runMacs()`, every |z| < 2 (largest
  1.96 of 14 statistics; rerun after Codex found the script's per-replicate `set.seed()` made the
  replicates dependent, which is fixed: one seed per run).

- **Phases (c) and (d) 2026-10-04.** `founders_coalescent()` exported (API in section 2 except
  that `ne` is the preset's and the `history` time / size bounds of phase (a) apply; `pool` and
  `bp` added). Split implemented in the Rust core (deme labels; per-deme sorted node times below the
  join; SMC' re-coalescence restricted to the lineage's deme below the join). Gates: split
  invariants through 5000 recombinations; between-deme pair E[T] = J + 1/2 at both ends of a
  recombining sequence; deme-0 subtree Kingman below a late join; between-minus-within pairwise
  differences = 2 theta J. BD scale (debug build): 10000 inbred WHEAT founders x 14000 markers on
  10 chromosomes in 42 s, about half of it `population_from_haplotypes()` (separate follow-up).

## 7. Size and plan

XL. Phases: (a) Rust SMC' core + PRNG + gates 1 and 3; (b) recombination gate 2 and the MaCS
window; (c) demography presets, split, `seg_sites`, `inbred`, Population output; (d) R API, docs,
DECISION-050, NEWS, BACKEND_CONTRACT, Codex theory review per phase.

## References (to verify)

- Chen GK, Marjoram P, Wall JD (2009) Fast and flexible simulation of DNA sequence data.
  *Genome Research* 19:136-142.
- Marjoram P, Wall JD (2006) Fast "coalescent" simulation. *BMC Genetics* 7:16.
- McVean GAT, Cardin NJ (2005) Approximating the coalescent with recombination. *Phil Trans R
  Soc B* 360:1387-1393.
- Hudson RR (2002) Generating samples under a Wright-Fisher neutral model of genetic variation.
  *Bioinformatics* 18:337-338.
- Watterson GA (1975) On the number of segregating sites in genetical models without
  recombination. *Theoretical Population Biology* 7:256-276.
- Gaynor RC, Gorjanc G, Hickey JM (2021) AlphaSimR. *G3* 11(2):jkaa017 (`runMacs()` presets).
