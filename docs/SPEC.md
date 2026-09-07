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
2. **Variance-partition layers** — `additive()`, `dominance()`, `epistasis()`, `vqtl()`
   each add a genetic variance component as a **proportion of total phenotypic
   variance**. Layers never change the architecture.
3. **Combination (optional)** — `complex_phenotypes()` merges complete single-
   architecture models under a common heritability.

The pipe runs **eagerly** — there is no `simulate()` terminal. The object returned by
`simulate_phenotype()` and by each layer already carries the realized phenotypes.

### Variance identity (single model)
```
Σ genetic proportions (additive + dominance + epistasis + vqtl) = h²   (≤ 1)
residual proportion (V_E / V_P)                                 = 1 − h²
```
Heritability is the sum of the genetic proportions; supplying proportions that sum to
> 1 is an error (§7).

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
  vary_qtn     = FALSE,
  seed         = NULL,
  h2           = NULL,           # one-call target heritability (see below)
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
  - `cor = 0` — target genetic correlation. For `n_traits = 2` the PleioArch
    bivariate-normal engine realizes it exactly in expectation; must satisfy
    `cor² ≤ pi_target × pi_secondary`. For `n_traits > 2`, `cor` is a scalar
    applied to every trait pair (or an `n_traits × n_traits` matrix) imposed via
    the v1.3 Cholesky decomposition, which **warns** that the phenotypes are
    correlated as requested but individual effect sizes are not guaranteed
    (DECISION-010, amended).
  - `pi_target` / `pi_secondary` — proportion of each trait's genetic variance explained
    by the pleiotropic QTNs (default: 1, i.e., pure pleiotropy).
  - `n_pleio_major = 1` — number of "major" pleiotropic QTNs (large-effect; rest are minor).
  - `prop_var_major = 0.5` — fraction of the pleiotropic variance captured by the major QTN(s).
  QTN counts: `n_qtn` sets total QTNs per trait; the pleiotropic subset is inferred from
  `pi_target`/`pi_secondary`. Trait-specific QTNs fill the remainder.
- `"ld"`: modern LD arguments (O5) — `ld_type = c("indirect", "direct")`,
  `r2_max = 0.8`, `r2_min = 0.2`, `r2_method = "composite"`. Both indirect and direct
  LD retained. (v1 names mapped in §12.)
- `"independent"`: trait-specific QTNs, no enforced cross-trait sharing.

`geno` accepts v1 inputs (HapMap/numeric object, or file/path) and a `Population`.

### 4.2 Variance-partition layers

All take `prop` (proportion of V_P) and return an updated `phenotype_sim`.

```
additive(sim,  prop, n_qtn = NULL, effect = NULL, dist = "geometric")
dominance(sim, prop, same_as_add = TRUE, n_qtn = NULL, degree = NULL, dist = "geometric")
epistasis(sim, prop, n_pairs = NULL, interaction = 2, effect = NULL, dist = "geometric")
vqtl(sim,      prop, same_as_add = TRUE, n_qtn = NULL, dist = "geometric")
```

- `prop`: scalar or length-`n_traits` vector.
- `dist`: within-layer effect distribution; default geometric. `effect` overrides with
  an explicit series (v1 `sim_method = "custom"`).
- `dominance(degree=)`: degree of dominance; `same_as_add` per §3.
- `epistasis(interaction = 2)`: pairwise (2-way) by default.
- `vqtl(same_as_add=)`: reuse additive QTNs by default (O3).

### 4.3 `complex_phenotypes()` — combine architectures

```
complex_phenotypes(..., h2)        # ... = two or more phenotype_sim objects
```

- Combines the inputs' genetic values, **weighted by their genetic variances**.
- Applies a common residual to reach `h2`; inputs' individual residuals are discarded.
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

---

## 5. Data Structures

- `phenotype_sim`: S3 object holding the `geno` reference, `architecture`, `n_traits`,
  `seed`, an ordered list of layers (type, prop, QTN indices, effects), the realized
  long-format phenotype table, and a computed variance budget. `print()` shows a
  plain-language variance budget (§9).
- `Population`: from the multi-generation functions; accepted anywhere `geno` is.

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

- Σ genetic `prop` ≤ 1, else error reporting the running total and offending layer.
- `prop` length must be 1 or `n_traits`.
- `"pleiotropy"` with `n_traits = 1` warns and proceeds as independent.
- `complex_phenotypes()` requires inputs on the same `geno` and identical `n_traits`.
- LD arguments warn and are ignored for non-`"ld"` architectures.
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

> Note: the current `test-v130-parity.R` calls the bare grammar
> (`simulate_phenotype() |> additive()`) and asserts bit-identity — that is the wrong
> target under DECISION-009 and must be repointed at `create_phenotypes()`.

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

### 8.3 Dataset
- `SNP55K_maize282_maf04` for all worked examples and tests.

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
  dominance(prop = 0.1, same_as_add = TRUE, degree = 0.5)

# Pleiotropy + LD combined under a common h²
pleio <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                            n_traits = 3, seed = 10) |> additive(0.4, n_qtn = 3)
ld    <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                            n_traits = 3, ld_type = "indirect", seed = 11) |>
         additive(0.3, n_qtn = 3)
both  <- complex_phenotypes(pleio, ld, h2 = 0.5)   # warns: differing seeds, uses 10
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
| `ld_method` | `r2_method` | modernized (O5) |
| `cor` | `cor` (PleioArch, §13) | v1's buggy Cholesky `cor` replaced by the PleioArch engine for `n_traits = 2`; the name `cor` is kept. For `n_traits > 2` the grammar falls back to the v1.3 Cholesky path with a warning (DECISION-010, amended). v1 `cor` also survives in the frozen legacy fn |
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
Pleiotropic effects are drawn from a **bivariate normal** (Cholesky decomposition of the
per-SNP covariance matrix scaled by QTN count):
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
- `n_traits = 2`: the exact bivariate-normal PleioArch engine above.
- `n_traits > 2`: **resolved** — the exact engine is not used; the grammar draws
  independent per-trait additive values and imposes the target `cor` via the v1.3
  Cholesky decomposition (see `base_line_multi_traits.R`), emitting a warning that
  the phenotypes are correlated as requested but individual effect sizes are not
  guaranteed. `cor` may be a scalar (every pair) or an `n_traits × n_traits`
  matrix (DECISION-010, amended).
