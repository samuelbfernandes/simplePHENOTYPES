# SPEC.md — simplePHENOTYPES v2 Public API

> Status: DRAFT v2 (open questions O1–O5 resolved; O6 deferred to ARCHITECTURE.md)
> Audience: implementers (Rust core + R bindings), CRAN maintainer, package users
> Companion docs: ARCHITECTURE.md, DECISIONS.md, BUGS.md
> Describes behavior, signatures, and data contracts. No implementation code.
> Reference implementations: `context/isqg/` (C++ meiosis/cross/DH), `context/PleioArch-main/` (exact rhoG pleiotropy).

---

## 1. Goals and Non-Goals

### Goals
- Replace the monolithic `create_phenotypes()` with a composable, tidyverse-style grammar.
- Single function call for simple cases; readable pipes for complex architectures.
- Faithful variance partitioning: heritability emerges from variance proportions.
- **Behavioral parity with v1.3.0** under matching genetic architectures (see §8).
- Single CRAN package; all functionality under `simplePHENOTYPES::`.
- Performance-critical paths ported to Rust (shared by R and, later, Python/Shiny).

### Non-Goals (v2)
- Expression-based simulation (deferred to v3).
- Partial pleiotropy as a dedicated architecture (recreated via `complex_phenotypes()`).
- `big_add_QTN_effect` in the new grammar (kept only in the frozen legacy
  `create_phenotypes()`; see §8 for how the reference runs remove it).

---

## 2. Core Concepts

1. **Foundation** — `simulate_phenotype()` chooses one genetic architecture and sets
   the residual variance. With no genetic layers the trait is pure noise: h² = 0.
2. **Variance-partition layers** — `additive()`, `dominance()`, and `epistasis()`
   add mean-effect genetic components; `vqtl()` adds genotype-dependent residual
   heterogeneity. Each `prop` is a requested marginal proportion of phenotypic
   variance. Layers never change the architecture.
3. **Combination (optional)** — `complex_phenotypes()` merges complete single-
   architecture models under a common heritability.

The pipe runs **eagerly** — there is no `simulate()` terminal. The object returned by
`simulate_phenotype()` and by each layer already carries the realized phenotypes.

### Variance identity (single model)
```
Σ mean-effect genetic proportions (additive + dominance + epistasis) = requested H²
Σ all layer proportions (mean effects + vqtl)                         ≤ 1
homoskedastic residual proportion                       = 1 − Σ all layer proportions
```
The vQTL share is residual rather than genetic variance and is excluded from H².
Supplying layer proportions that sum to > 1 is an error (§7).

**Modeling convention.** Each component is centered and scaled to its target
proportion. The coding is a simulation convention -- additive uses
the -1/0/1 dosage, dominance a heterozygote-deviation indicator, epistasis a
centered additive-by-additive (or a x d / d x d) product -- **not** Fisher's orthogonal
decomposition into average effects and orthogonal dominance deviations. When
dominance or epistasis share loci with the additive layer, the scaled
components are not exactly uncorrelated at non-0.5 allele frequencies, so the
reported H² is computed from the realized genetic and phenotypic values rather
than asserted from the requested budget. It usually tracks the sum of genetic
proportions closely but not to machine precision,
and the per-component "variances" are the simulation's, not the classical
orthogonal Va/Vd. For an additive-only model the additive proportion is the
narrow-sense h² under Hardy-Weinberg.

---

## 3. Default Behavior

- Architecture: `"independent"`; `n_traits = 1` (single trait is the default).
- Within-layer effect distribution: **geometric** (effect_i = effect^i), matching v1
  `sim_method = "geometric"`. No big-QTN effect in the new grammar.
- `dominance(same_as_add = TRUE)` by default — dominance reuses the additive QTNs;
  `same_as_add = FALSE` draws fresh dominance-only loci. **`vqtl()` follows the same
  pattern** (O3): `same_as_add = TRUE` reuses additive QTNs by default.
- **Output format defaults to long** (O2): one row per individual × trait × rep, with
  `id`, `trait`, `rep`, `value` columns. Wide / multi-file / gemma / plink outputs are
  separate `write_*()` exporters, not a core argument.
- Reproducibility: `seed` is set in `simulate_phenotype()` and threaded to all layers
  (§6).
- `vary_qtn = TRUE` draws an independent QTN set per replication (per-rep loci
  and effects), so replications are distinct genetic architectures rather than
  the same one with fresh residuals. Layers with user-supplied `qtn` keep their
  fixed loci across reps.
- `mean` adds a per-trait intercept at the phenotype level; genetic values
  (`genetic_values()`) stay centered.
- `individuals` restricts the simulated set (marker MAF is recomputed on the
  subset). It never copies the genotypes -- the accessor just returns the
  selected rows.

---

## 4. Function Reference

### 4.1 `simulate_phenotype()` — the foundation

```
simulate_phenotype(
  geno,                          # genotype input OR a Population (cross/selfcross/dh)
  architecture = "independent",  # "independent" | "pleiotropy" | "ld"
  n_traits     = 1,
  n_qtn        = 0,              # baseline QTN count (see O1 override rule below)
  n_reps       = 1,
  vary_qtn     = FALSE,          # redraw QTNs each replication
  seed         = NULL,
  h2           = NULL,           # requested genetic variance share (see below)
  mean         = NULL,           # per-trait intercept (scalar or length n_traits)
  individuals  = NULL,           # subset of individuals to simulate
  model        = "A",           # one-call model: "A" | "AD" | "AE"
  ...                            # architecture-specific args
)
```

Returns a `phenotype_sim` object (h² = 0 until a layer is added).

**One-call vs piped (folds in the former `sim_phenotypes()` shortcut).** If the
call already carries a self-sufficient genetic spec — `h2` supplied together with
`n_qtn > 0` — `simulate_phenotype()` builds the implied `model` (additive by
default; "AD"/"AE" split `h2` across components) and realizes a complete
  phenotype in one call. Otherwise it returns the h² = 0 foundation for the user to
complete via the layer pipe. Either way the result is a realized `phenotype_sim`,
so piping more layers onto a one-call result still works.

**`n_qtn` override rule (O1):** `n_qtn` may be set both here (baseline) and per layer.
A per-layer `n_qtn` **overrides** the baseline and emits a warning naming both values.

Architecture-specific arguments:
- `"pleiotropy"`: exact genetic-correlation control via the PleioArch algorithm (§13).
  Additional args (all with sensible defaults):
  - `cor = 0` — target genetic correlation, realized exactly in expectation for **any
    number of traits** by the PleioArch multivariate-normal engine (DECISION-013). Either
    a scalar applied to every trait pair, or a full `n_traits × n_traits` matrix
    (negative correlations allowed). Must be attainable: the implied genetic covariance
    matrix has to be positive semi-definite, which for two traits is exactly
    `cor² ≤ pi_1 × pi_2` and for more traits additionally rules out mutually
    inconsistent requests (three traits cannot all be strongly negatively correlated).
    Unattainable requests are an error, not a silent approximation. Supplying `cor`
    emits a one-per-session citation message for the correlation-control algorithm
    (published separately); it is an [rlang::inform()] message, silenceable with
    `suppressMessages()`, and is not emitted when `cor` is absent.
  - `pi` — proportion of each trait's genetic variance explained by the pleiotropic QTNs
    (default 1, pure pleiotropy). Scalar or length-`n_traits`. `pi_target` /
    `pi_secondary` remain as the two-trait spelling.
  - `n_pleio_major = 0` — number of "major" pleiotropic QTNs (large-effect; rest are
    minor). **Defaults to none**, so every pleiotropic QTN is drawn from the same
    multivariate normal.
  - `prop_var_major = 0` — fraction of the pleiotropic variance captured by the major
    QTN(s).

  These two default to zero deliberately. Concentrating a large share of the genetic
  variance in one locus makes the realized correlation depend chiefly on whether that
  single bivariate draw happens to agree in sign, so it stops converging on `cor` as
  `n_qtn` grows — defeating the reason for adopting this engine (DECISION-007), and
  reintroducing by the back door the big-QTN behavior that §1 lists as a non-goal for the
  grammar. Measured on `SNP55K_maize282_maf04` at `cor = 0.6` over 40 seeds, the realized
  genetic correlation has SD 0.13 / 0.09 / 0.08 at `n_qtn` = 20 / 100 / 400 with the
  defaults, against 0.25 / 0.32 / 0.23 — no convergence — under the former
  `n_pleio_major = 1, prop_var_major = 0.5`. Set both explicitly when a major-effect
  pleiotropic locus is the object of study.
  QTN counts: `n_qtn` sets total QTNs per trait; the pleiotropic subset is inferred from
  `pi`. Trait-specific QTNs fill the remainder, drawn separately for each trait.
- `"ld"`: **linked-but-not-pleiotropic** — a spurious genetic correlation
  between **exactly two traits** (`n_traits = 2`, else an error) whose causal
  loci are *distinct* but sit in linkage disequilibrium. For each of `n_qtn`
  loci a linked pair is drawn: one SNP is causal for trait 1, the other for
  trait 2, with squared correlation r2 in `[r2_min, r2_max]`. No SNP is causal
  for both traits — the correlation comes entirely from the linkage between the
  two traits' separate loci (contrast `"pleiotropy"`, where a shared locus and
  correlated effects drive the correlation). Arguments (O5):
  `ld_type = c("direct", "indirect")`, `r2_max = 0.8`, `r2_min = 0.2`.
  - `"direct"` (default): trait 1's causal SNP and trait 2's causal SNP are
    directly in LD (r2 in the window).
  - `"indirect"`: both causal SNPs flank a shared, **non-causal** "cause-of-LD"
    locus and are each in LD with it — the correlation is mediated by an
    unobserved variant. One flanking SNP is taken upstream, one downstream.

  `qtn_table()` reports the linked pair on each row: `QTN_t1` and `QTN_t2` name
  the causal SNP for trait 1 and trait 2, and `ld_r2` is the squared correlation
  between them. (For `"indirect"` the hidden cause-of-LD locus is not shown as a
  QTN — it is not causal — but is kept on the layer object.) (v1 names mapped in
  §12; this is the grammar equivalent of the legacy `qtn_linkage` engine.)
- `"independent"`: trait-specific QTNs, no enforced cross-trait sharing. Pass
  `distinct_chr = TRUE` to force every trait's QTNs onto a disjoint set of
  chromosomes (round-robin), giving genetically unlinked traits.

`geno` accepts v1 inputs (HapMap/numeric object, or file/path) and a `Population`.

### 4.2 Variance-partition layers

All take `prop` (proportion of V_P) and return an updated `phenotype_sim`.

```
additive(sim,  prop, n_qtn = NULL, qtn = NULL, effect = NULL,
         phase = c("coupling","repulsion"), dist = "geometric")
dominance(sim, prop, same_as_add = TRUE, n_qtn = NULL, qtn = NULL, dist = "geometric")
epistasis(sim, prop, n_pairs = NULL, interaction = 2, interaction_type = "a", qtn = NULL, effect = NULL, dist = "geometric")
vqtl(sim,      prop, same_as_add = TRUE, n_qtn = NULL, qtn = NULL, dist = "geometric")
```

- `prop`: scalar or length-`n_traits` vector.
- `dist`: within-layer effect distribution; default geometric. `effect` overrides with
  an explicit series (v1 `sim_method = "custom"`).
- `dominance()`: variance set by `prop`; `same_as_add` per §3. No separate
  degree-of-dominance argument (it would be washed out by variance scaling; the
  dominance/additive variance ratio is `prop_dom/prop_add`).
- `epistasis(interaction = 2)`: pairwise (2-way) by default.
- `epistasis(interaction_type=)`: `"a"` (additive; centered dosage) or `"d"`
  (dominance; centered het indicator) per interacting position -- `c("a","a")`
  is a x a (default), `c("a","d")` a x d, `c("d","d")` d x d. Each term is
  centered per locus, subtracting out the lower-order marginals. `"d"` terms
  need heterozygotes and are near-degenerate on an inbred panel.
- `vqtl(same_as_add=)`: reuse additive QTNs by default (O3). Its `prop` is an
  explicit marginal share for a heterogeneous residual component and is never
  drawn from `h2`. Conditional variance uses a log link, so it stays positive.
- `qtn`: fix this layer's causal loci by marker name or index (a vector for all
  traits, or a length-`n_traits` list); the other layers stay random. Epistasis
  takes an `n_pairs x interaction` matrix. `n_qtn`/`n_pairs` follow from it, and
  fixed loci are exempt from `vary_qtn`.
- `additive(phase=)`: `"coupling"` (default) or `"repulsion"`. Repulsion
  alternates effect signs so linked increasing/decreasing alleles oppose; the
  realized genetic variance is still `prop` (pinned by the partition), only the
  sign structure and cross-locus covariance change.

### 4.3 `complex_phenotypes()` — combine architectures

```
complex_phenotypes(..., h2)        # ... = two or more phenotype_sim objects
```

- Combines the inputs' genetic values, **weighted by their genetic variances**.
- Requires identical genotype data, marker order, individual subset, trait count,
  replication count, and trait means.
- Preserves every replication, scales the combined genetic value to requested
  variance `h2`, and adds a common residual; inputs' individual residuals are discarded.
- **Seed handling (O4):** if inputs were built with different seeds, the **first
  input's seed is used** and a warning is emitted.
- Recreates partial pleiotropy: combine a `"pleiotropy"` model with an `"independent"`
  model.

### 4.4 `create_phenotypes()` — frozen legacy function (DECISION-008)

- Full v1 signature preserved, including `big_add_QTN_effect`, `architecture =
  "partially"`, `sim_method`, `vary_QTN`, `cor`/`cor_res`, etc.
- **Frozen legacy: bugfix-only**, marked `lifecycle::badge("superseded")`. It retains the
  original v1 code paths, seed arithmetic, and `RNGversion('3.5.1')`. It does **not**
  delegate to the new grammar internals — the legacy engine and the new grammar are
  separate implementations that coexist.
- Its output is pinned by the `inst/extdata/v1_3_0_reference/*.rds` regression guard
  (§8.1): bug fixes must not change its output unintentionally.

### 4.5 One-call simulation (folded into `simulate_phenotype()`)

The beginner one-call shortcut is **not** a separate function (the name
`sim_phenotypes()` would confuse). It is folded into `simulate_phenotype()`:
supply `h2` together with `n_qtn > 0` and it targets that heritability and
realizes a model in a single call (see §4.1, "One-call vs piped"). The maintained
frozen `create_phenotypes()` (§4.4) also remains for v1-style one-call users.

### 4.6 Multi-generation crossing (DECISION-002, DECISION-004, DECISION-012)

Meiosis needs recombination distances, which a physical (`chr`, `pos`) map does not
supply. `synthetic_map()` builds a modelled cM map from physical positions when no real
one is available; `as_population()` phases founders from a `-1/0/1` dosage matrix; the
three mating functions consume and produce `Population`s.

```
synthetic_map(chr, pos, total_cm = NULL, cm_per_mb = 0.73, centromere = NULL,
             suppression = 0.85, width = 0.15)                    # -> numeric cM vector

as_population(geno, individuals = NULL)                           # -> Population
n_individuals(x)                                                  # -> integer
dosages(x)                                                        # -> -1/0/1 matrix
x[i]                                                              # subset individuals

cross(mother, father, n = 1, seed = NULL)                         # -> Population, n progeny
selfcross(parent, n = 1, seed = NULL)                             # -> Population, n progeny
double_haploid(parent, n = 1, seed = NULL)                        # -> Population, n progeny
```

**`synthetic_map()`** models pericentromeric suppression rather than a uniform rate:
each physical position gets a local rate `r(p) = 1 - s·exp(-0.5((p-c)/(wL))²)` (centromere
`c`, span `L`, suppression `s`, width `w`), integrated across marker intervals and
rescaled to `total_cm` (or `cm_per_mb × span` when `total_cm` is `NULL`). `suppression =
0` gives a uniform map. The result is a **model, not a measured linkage map** — it must
not be quoted as an estimate of real recombination distance. The bundled
[SNP55K_maize282_maf04] carries such a map in its `cm` column (~1,497 cM total, ~0.73
cM/Mb average); earlier package versions shipped `cm` as an all-`NA` placeholder.

**`as_population()`** requires a fully-resolved genetic map (`cm` non-missing and
non-decreasing within each chromosome) and dosage coded `-1/0/1`. Homozygotes phase
exactly; heterozygotes are phased arbitrarily as allele-1 on the first strand. This is
near-lossless for an inbred panel (rare heterozygotes) but means first-generation linkage
between heterozygous sites in an outbred sample is not realistic. The current public API
does not import external phase, so substantially heterozygous unphased founders are not
supported for realistic multi-generation recombination studies.

**`cross()` / `selfcross()` / `double_haploid()`** take single-individual `Population`s
(`mother`, `father`, or `parent`; subset a larger one with `x[i]`) and return `n` progeny
as a new `Population`. Recombination follows the Karlin & Liberman count-location
process — crossover count Poisson with mean equal to chromosome length in Morgans,
positions uniform along it, chromosomes independent — matching isqg exactly
(DECISION-012). `double_haploid()` progeny are fully homozygous by construction (a single
gamete, duplicated); `selfcross()` halves heterozygosity per generation; `cross()`
combines one gamete from each parent. All three draw every random quantity in R, in
isqg's exact order, before calling the Rust core, which is a pure function of those draws
and never calls an RNG (DECISION-012) — `seed` (or an ambient `set.seed()`) fully
determines the outcome.

A `Population` is accepted anywhere `simulate_phenotype()` takes `geno` (§4.1), so a
simulated pedigree can be phenotyped directly without converting back to a dosage matrix.

---

## 5. Data Structures

- `phenotype_sim`: S3 object holding the `geno` reference, `architecture`, `n_traits`,
  `seed`, an ordered list of layers (type, prop, QTN indices, effects), the realized
  long-format phenotype table, and a computed variance budget. `print()` shows a
  plain-language variance budget (§9).
- `Population`: S3 object (`as_population()`, `cross()`, `selfcross()`,
  `double_haploid()`) holding a marker `map` (`snp`, `chr`, `pos`, `cm`), phased `cis` /
  `trans` haplotype matrices (markers × individuals), individual `ids`, and an `origin`
  label. `n_individuals()`, `dosages()` (collapses to `-1/0/1`), `x[i]` (subset) and
  `print()` operate on it; accepted anywhere `geno` is (§4.1, §4.6).

---

## 6. Reproducibility / Seed Threading

- `seed` stored on the `phenotype_sim` at creation.
- Each layer derives a deterministic sub-seed from `(seed, layer_index, layer_type)`.
- Adding/removing/reordering a layer must not change other layers' realized values;
  tests assert this.
- This clean scheme is the grammar's own (DECISION-009 — no v1 bit-parity obligation).
  The frozen `create_phenotypes()` keeps its *own* legacy seed arithmetic
  (`(seed + z) * round(h2 * 10)`, `seed + i`, etc.) and `RNGversion('3.5.1')`; the two
  schemes are independent and need not agree.

---

## 7. Validation and Errors

- Σ all layer `prop` ≤ 1, else error reporting the running total and offending layer.
- `prop` must be finite, within [0,1], and have length 1 or `n_traits`.
- Counts are positive whole numbers (except baseline `n_qtn`, which may be zero).
- `"pleiotropy"` with `n_traits = 1` errors.
- `complex_phenotypes()` requires identical genotype data, individuals, traits,
  replications, and means.
- Architecture-specific arguments supplied to another architecture error.
- Correlation inputs must be finite, bounded, symmetric with unit diagonal, and
  imply a positive-semidefinite covariance model.
- Random QTNs are drawn only from polymorphic markers; a requested nonzero layer
  that has no usable design variance errors.
- Per-layer `n_qtn` overriding baseline `n_qtn` warns (O1).
- Differing seeds in `complex_phenotypes()` warn; first seed used (O4).

---

## 8. Testing Requirements

### 8.1 v1.3.0 regression guard + grammar validation (DECISION-009)
The new grammar owes v1 **no bit-for-bit parity**. The captured references serve two
distinct purposes:

**(a) Regression guard on the frozen `create_phenotypes()` (§4.4).**
1. For each README/vignette example, the v1.3.0 run (with `big_add_QTN_effect` removed
   where present) is captured: the **selected QTNs** and the **phenotype values**, stored
   as RDS under `inst/extdata/v1_3_0_reference/`.
2. `test-v130-parity.R` calls **`create_phenotypes()`** (the legacy fn, not the grammar)
   and asserts it still reproduces those references. Bug fixes that change the output must
   be deliberate and re-bless the references. These tests gate the legacy path.

**(b) Statistical / structural validation of the new grammar.**
The grammar is checked against properties, not bit-identity:
   - variance-partition identity (realized h² ≈ Σ genetic props within tolerance),
   - realized `cor` ≈ target for pleiotropy (PleioArch, 2 traits),
   - correct QTN-count structure per trait,
   - seed-threading invariance (§6).

`test-v130-parity.R` calls **`create_phenotypes()`** (not the bare grammar), matching
this section.

Coverage: single trait, pleiotropy (`seed = 10`), LD/spurious (`seed = 200`, both
`ld_type`), and the partial-pleiotropy example.

**Implementation constraint (resolved, ARCHITECTURE.md / DECISION-006):** stochastic steps
(QTN sampling, effect series, residual draws) stay on R's RNG; only deterministic numeric
work (genetic-value assembly, matrix ops) is ported to Rust. This keeps both the frozen
legacy path and the grammar's own seed reproducibility intact.

### 8.2 Internal consistency tests
- Variance-partition identity (realized h² ≈ Σ genetic props within tolerance).
- Seed-threading invariance (§6).
- `complex_phenotypes()` genetic-variance weighting and first-seed rule.

### 8.3 isqg parity and crossing tests (DECISION-012)
`test-isqg-parity.R` gates `genome.rs`/`meiosis.rs` against `inst/extdata/isqg_v1_outputs/`
(captured by `capture_isqg_references.R`) with **exact** bit-equality — not statistical
equivalence — because R replicates isqg's draw order exactly and the Rust core never
calls an RNG, so nothing but a real transcription bug can cause disagreement. Fixtures
store the drawn `counts`/`chiasmata`/`flips` alongside isqg's output so a failure can be
localized to the R draw order or the Rust core.

`test-cross.R` covers `as_population()`/`cross()`/`selfcross()`/`double_haploid()`:
correct phasing of homozygotes and heterozygotes, `dosages()` round-tripping, F2/DH
segregation ratios, and recombination fraction against Haldane's map function at known
genetic distances — the check that caught the phase-reconstruction bug described in
`docs/NEXT_STEPS.md` (segregation ratios alone did not).

### 8.4 Dataset
- `SNP55K_maize282_maf04` for all worked examples and tests, including its bundled
  synthetic `cm` map (§4.6) for crossing examples.

---

## 9. Example `phenotype_sim` print output (illustrative)

```
<phenotype_sim>  (realized · long format)
  Genotypes: SNP55K_maize282_maf04   Traits: 3   Architecture: pleiotropy   Seed: 10
  Variance partition (proportions of V_P):
    additive    0.40   (3 QTNs, geometric)
    dominance   0.10   (same QTNs as additive)
    epistasis   0.10   (2 pairs, 2-way)
    residual    0.40
  Implied h² (broad) = 0.60
```

---

## 10. Worked Examples (new grammar)

```r
# Single trait (default): additive only, h² = 0.5
ph <- simulate_phenotype(SNP55K_maize282_maf04, seed = 1) |>
  additive(prop = 0.5, n_qtn = 3)

# Pleiotropy: 3 traits, additive + dominance on the SAME QTNs
ph <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                         n_traits = 3, seed = 10) |>
  additive(prop = 0.4, n_qtn = 3) |>
  dominance(prop = 0.1, same_as_add = TRUE)

# Pleiotropy + LD combined under a common h²
pleio <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                            n_traits = 3, seed = 10) |> additive(0.4, n_qtn = 3)
ld    <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                            n_traits = 3, ld_type = "indirect", seed = 11) |>
         additive(0.3, n_qtn = 3)
both  <- complex_phenotypes(pleio, ld, h2 = 0.5)   # warns: differing seeds, uses 10

# Multi-generation cross, phenotyped directly (§4.6)
pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
f1  <- cross(pop[1], pop[2], n = 1, seed = 1)
f2  <- selfcross(f1, n = 200, seed = 2)
ph  <- simulate_phenotype(f2, seed = 3) |> additive(prop = 0.5, n_qtn = 3)
```

---

## 11. v1 → v2 Argument Mapping (user migration reference)

> These are equivalences for users porting v1 scripts to the v2 grammar. The frozen
> `create_phenotypes()` does NOT translate to the grammar at runtime (DECISION-008) — it
> keeps its own v1 code paths. `big_add_QTN_effect`, `cor`/`cor_res`, and
> `architecture = "partially"` live only in the legacy function.

| v1 `create_phenotypes()` | v2 grammar | Notes |
|---|---|---|
| `geno_obj` / `geno_file` / `geno_path` | `simulate_phenotype(geno=)` | unified |
| `ntraits` | `simulate_phenotype(n_traits=)` | |
| `architecture = "pleiotropic"` | `architecture = "pleiotropy"` | |
| `architecture = "partially"` | `complex_phenotypes()` | combine models |
| `architecture = "LD"` | `architecture = "ld"` | |
| `add_QTN_num` | `additive(n_qtn=)` | |
| `dom_QTN_num` | `dominance(n_qtn=)` | |
| `epi_QTN_num` | `epistasis(n_pairs=)` | |
| `add_effect` / `dom_effect` / `epi_effect` | `*(effect=)` | geometric default |
| `big_add_QTN_effect` | — (shim only) | removed from grammar; removed in parity runs |
| `h2` | Σ layer `prop` (single) / `complex_phenotypes(h2=)` | emerges |
| `same_add_dom_QTN` | `dominance(same_as_add=)` | |
| `degree_of_dom` | `dominance(degree=)` | |
| `epi_interaction` | `epistasis(interaction=)` | 2-way default |
| `sim_method = "geometric"` | `dist = "geometric"` | default |
| `sim_method = "custom"` | `effect = <series>` | |
| `type_of_ld` | `ld_type` | modernized (O5) |
| `ld_max` / `ld_min` | `r2_max` / `r2_min` | modernized (O5) |
| `cor` | `cor` (PleioArch, §13) | v1's buggy Cholesky `cor` replaced by the PleioArch engine for any number of traits (DECISION-013); the name `cor` is kept. Scalar or full matrix. v1 `cor` also survives in the frozen legacy fn |
| `cor_res` | residual-correlation arg | retained in legacy fn; grammar equivalent TBD |
| `rep` | `n_reps` | |
| `vary_QTN` | `vary_qtn` | |
| `seed` | `seed` | |
| `output_format` / `out_geno` / `output_dir` / `home_dir` / `to_r` | `write_*()` exporters | long is default; object always returned |

---

## 12. Remaining Open Question

- **O6 — RESOLVED** (DECISION-006 + DECISION-009): stochastic draws stay in R; only
  deterministic computation is ported to Rust. There is no need to replicate v1's RNG in
  Rust, because the new grammar owes v1 no bit-parity and the frozen `create_phenotypes()`
  keeps its own R-side seeding. No open questions remain blocking.

(O1–O5 resolved: n_qtn on both with per-layer override + warning; long default output;
vqtl uses same_as_add; complex_phenotypes uses first seed + warning; LD names modernized.)

---

## 13. PleioArch Algorithm — Exact Genetic Correlation for Pleiotropy

> Reference implementation: `context/PleioArch-main/Functions/simulateEffects.R`
> DECISION-007: adopted as the effect-generation engine for `architecture = "pleiotropy"`.

### Motivation
v1 "pleiotropy" merely shares QTN loci across traits; the resulting genetic correlation
is an emergent by-product of allele-frequency differences and cannot be set precisely.
The PleioArch algorithm draws allelic effects from correlated distributions so that the
realized genetic correlation equals the user's `cor` exactly in expectation (for
`n_traits = 2`).

### QTN classification
Each pleiotropy simulation partitions QTNs into three classes:

| Class | Description |
|---|---|
| Pleiotropic Major | One (or few) large-effect loci shared across traits |
| Pleiotropic Minor | Many small-effect loci shared across traits |
| Trait-Specific | Loci affecting only one trait (zero effect on others) |

### Variance budget
Given `n_traits = 2` (pairwise; extended per-pair for >2):
```
V_G_target    = h2_target      (genetic variance of target trait)
V_G_secondary = h2_secondary

V_pleio_target    = pi_target    × V_G_target     (pleiotropic share, target)
V_pleio_secondary = pi_secondary × V_G_secondary  (pleiotropic share, secondary)
Cov_pleio         = cor × sqrt(V_G_target × V_G_secondary)

V_spec_target    = (1 − pi_target)    × V_G_target
V_spec_secondary = (1 − pi_secondary) × V_G_secondary
```
Biological constraint (enforced, raises error if violated):
```
cor² ≤ pi_target × pi_secondary
```

### Effect draws
Pleiotropic effects are drawn from a **multivariate normal** (Cholesky decomposition of
the per-SNP covariance matrix scaled by QTN count; bivariate when `n_traits = 2`):
```
Sigma_major = Sigma_pleio × prop_var_major   / n_pleio_major
Sigma_minor = Sigma_pleio × (1−prop_var_major) / n_pleio_minor
```
Trait-specific effects are drawn from a **univariate normal** with variance
`V_spec / n_spec`.

Effects are then scaled by `1 / sqrt(2 × MAF × (1 − MAF))` to convert from allele
to per-genotype scale (matching the `scaleQTNEffects()` step in the reference code).

### Integration with the v2 grammar
- Effect draws (bivariate/univariate normals) **remain in R** (DECISION-006 — parity RNG
  constraint). Deterministic assembly (Cholesky, matrix multiply) may move to Rust.
- `additive(prop=, ...)` calls this engine when the parent `phenotype_sim` carries
  `architecture = "pleiotropy"`, passing the resolved QTN effects directly.
- v1 `create_phenotypes(architecture = "pleiotropic")` does **not** map to PleioArch — it
  is frozen legacy (§4.4, DECISION-008) and keeps v1's shared-loci behavior. PleioArch is
  the new grammar's pleiotropy engine, controlled by `cor` (DECISION-010, amended).
- **Any `n_traits`**: the exact engine generalizes directly (DECISION-013). The
  bivariate normal becomes an `n × n` multivariate normal with

  ```
  Sigma[i,i] = pi_i  × V_i                    (pleiotropic share of trait i)
  Sigma[i,j] = cor_ij × sqrt(V_i × V_j)       (whole covariance: specific loci are
                                               independent across traits)
  ```

  Trait-specific effects remain univariate normal with variance `(1 − pi_i) × V_i`, so
  each trait's total genetic variance is still `V_i` and every pair realizes `cor_ij`.
  For two traits this reduces algebraically to the bivariate reference implementation.
  The `cor² ≤ pi_i × pi_j` constraint generalizes to "`Sigma` must be positive
  semi-definite", checked by eigenvalue and raised as an error naming the smallest
  eigenvalue. The Cholesky-on-genetic-values fallback is **removed**: it never
  controlled individual QTN effects and is no longer needed.
