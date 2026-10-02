# SPEC.md — simplePHENOTYPES v2 Public API

> Status: DRAFT v2 (open questions O1–O5 resolved; O6 deferred to ARCHITECTURE.md)
> Audience: implementers (Rust core + R bindings), CRAN maintainer, package users
> Companion docs: ARCHITECTURE.md, DECISIONS.md, BUGS.md
> Describes behavior, signatures, and data contracts. No implementation code.
> Reference implementations: `context/isqg/` (C++ meiosis/cross/DH), `context/PleioArch-main/` (rhoG-controlled pleiotropy).

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
   variance. Layers never change the architecture. (Exception: `additive(orthogonal
   = TRUE, a =, d =)` switches that layer to the orthogonal genotypic model —
   per-locus `a`/`d`, with the additive/dominance split emerging from the effects
   and allele frequencies; see the Modeling convention below and DECISION-020.)
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
components are generally correlated, so the reported H² is computed from the
realized genetic and phenotypic values rather than asserted from the requested
budget, and it can differ from the sum of the genetic proportions. For
`additive()` + `dominance()` on the **same loci** the difference is structural,
not a finite-sample effect: for one additive and one dominance layer
`Var(g) = prop_A + prop_D + 2Cov(c_A, c_D)` (several layers of a type add their
mutual covariance to `Var(c_A)` / `Var(c_D)`; the general identity is
`Var(g) = Var(c_A) + Var(c_D) + 2Cov(c_A, c_D)`), and with
`x` the −1/0/1 dosage and `h` the heterozygote indicator, `Cov(x, h) =
−(2p − 1) 2pq` at each locus under Hardy-Weinberg (`p` = frequency of the +1
allele). With the default geometric series (one effect sign at every locus) each locus's
cross term has the sign of `1 - 2p`, so it does not average out when the
counted-allele frequencies are on one side of 0.5 (it partly cancels when they
straddle 0.5), and its sign flips with the allele coding
(audit 2026-09-29: an outbred panel with MAF 0.10–0.20, `n_qtn = 20`, A 0.4 + D
0.1, realized H² ≈ 0.63 with the minor allele coded +1 and ≈ 0.25 with the same
loci and seeds coded the other way, for a requested 0.5). The per-component
"variances" are the simulation's, not the classical orthogonal Va/Vd.
When additive and dominance layers share loci the `phenotype_sim` therefore
reports, per trait, the requested share, the realized share and the realized
`Var(A)`, `Var(D)` and `2Cov(A,D)` of the block, and `Var(c_A)`, `Var(c_D)` and
`2Cov(c_A,c_D)` (fractions of V_P; `$ad_report`,
and a note in `print()`; DECISION-034), and recommends
`additive(orthogonal = TRUE, a =, d =)` for a Fisher-orthogonal partition. For an
additive-only model the additive proportion is the narrow-sense h² under
Hardy-Weinberg.

**Orthogonal genotypic model (exception).** `additive(orthogonal = TRUE, a =, d =)`
opts a single additive layer out of the variance-partition convention above. It
builds each locus's genotypic value from a per-locus additive effect `a` and
dominance deviation `d` (values `-a`/`+d`/`+a` for gene content 0/1/2), scales the
**whole** genotypic value to the layer `prop`, and decomposes it into the
average-effect breeding value `A = Σ αⱼ(xⱼ − 2pⱼ)`, `αⱼ = aⱼ + dⱼ(1 − 2pⱼ)`
(Falconer 1985; the transmitting-ability form, consistent with the `on = "bv"`
breeding value, DECISION-019) and the realized dominance deviation `D = g − A`.
`A` and `D` are orthogonal (`Cov(A, D) = 0`) in expectation under random mating,
which puts each locus in Hardy-Weinberg proportions (e.g. an F2). This does *not*
require linkage equilibrium — between-locus LD is compatible with orthogonality in
this no-epistasis model — but per-locus HWE alone does not suffice under arbitrary
nonrandom multilocus genotype association (selection, structure). Whatever
covariance remains is reported, not normalized away, as an `add_dom_cov`
variance-budget row so `additive + dominance + add_dom_cov = prop` (≈0 for a large
random-mating / F2 sample). This mode is incompatible with
`vary_qtn`, a separate `dominance()` layer, and the `"pleiotropy"` (multi-trait) /
`"ld"` architectures. See DECISION-020.

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
  reps         = 1,             # records per entry; residual variance V_E / reps
  resid_cor    = NULL,          # residual correlation between traits (scalar or matrix); NULL = independent
  ...                            # architecture-specific args
)
```

Returns a `phenotype_sim` object (h² = 0 until a layer is added).

**Residual correlation (`resid_cor`, DECISION-042).** `cor` is the *genetic* correlation;
`resid_cor` is the target correlation between the traits' *residuals*: `NULL` (independent, as before), one
value for every pair (>= `-1/(n_traits - 1)`), or an `n_traits x n_traits` symmetric PSD matrix with unit
diagonal. Per-trait residual variance is still `1 - sum(prop)` exactly (each column is re-standardized), so
realized H2 and `var_budget` do not change; the realized sample correlation matches the target up to
sampling error of order `1/sqrt(n)`. A `vqtl()` heterogeneity component is drawn independently per trait and
dilutes it; `reps` scales each trait's residual and leaves the correlation as is. The phenotypic correlation
is a variance-weighted mix of the genetic and residual correlations.

**Entry-mean replication (`reps`, DECISION-038).** `reps` is a positive whole number (scalar or
one per trait): the phenotype is the mean of `reps` independent records of the same genotype
(AlphaSimR `setPheno(varE, reps)` semantics). The realized residual is the `reps = 1` realized residual
divided by `sqrt(reps)`, so its realized variance is the `reps = 1` residual variance / `reps` (for vqtl:
`[V0 + Vv + 2Cov(e0, ev)] / reps`, not the nominal `V_E / reps`, because the two standardized components
have non-zero sample covariance). Target (expected) heritabilities: single-record V_G/(V_G+V_E) (what `h2`
requests) and entry-mean V_G/(V_G+V_E/reps). `h2`, layer `prop` and `var_budget` stay on the
single-record scale. Realized heritabilities are Var(G)/Var(y) from the realized values: entry-mean
Var(y_bar) = V_G + V_E/reps + 2Cov(G,e)/sqrt(reps); record Var(y) = V_G + V_E + 2Cov(G,e); the allocation
formula is the realized value only when the sample Cov(G,e) = 0. The printed realized H2 is the
entry-mean value and `print()` also shows the single-record value when any `reps > 1` (and
`reps (per trait) = [...]` when `reps` varies). Records are iid given the genotype and, with a derived
transcriptome layer, the fixed transcriptome component: the environmental transcriptome part is a
persistent entry-level quantity that is not redrawn per record, so replication is conditional on the
fixed transcriptome covariate. `reps = 1` is bit-identical to the output before the argument existed.

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
- `"pleiotropy"`: genetic-correlation control via the PleioArch algorithm (§13), in
  every mean-effect layer (DECISION-023).
  Additional args (all with sensible defaults):
  - `cor = 0` — target genetic correlation for **any number of traits**, set by the
    PleioArch multivariate-normal engine (DECISION-013): the effect draw's cross-trait
    covariance and variances equal their targets in expectation, so `cor` is the ratio
    of expected moments. The realized correlation is a random ratio that converges to
    `cor` as the shared QTNs and the individuals grow (with a fixed sample it levels
    off at the sampling spread of a correlation over n) when the causal loci are in approximate linkage
    equilibrium and no major QTN keeps a fixed variance share (see `n_pleio_major`
    below); with few shared QTNs it is attenuated toward 0 on average, and strong
    LD can prevent convergence (DECISION-023, Scope). Either
    a scalar applied to every trait pair, or a full `n_traits × n_traits` matrix
    (negative correlations allowed). Must be attainable: the implied genetic covariance
    matrix has to be positive semi-definite, which for two traits is exactly
    `cor² ≤ pi_1 × pi_2` and for more traits additionally rules out mutually
    inconsistent requests (three traits cannot all be strongly negatively correlated).
    Unattainable requests are an error, not a silent approximation. Supplying `cor`
    emits a one-per-session citation message for the correlation-control algorithm
    (Prado et al., in preparation); it is an [rlang::inform()] message, silenceable
    with `suppressMessages()`, and is not emitted when `cor` is absent.
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
  `ld_type = c("direct", "indirect")`, `r2_max = 0.8`, `r2_min = 0.2`,
  `partner = c("strongest", "random")`. The window must satisfy `r2_max > 0` and
  `r2_min < 1`, and markers with r2 = 0 or r2 = 1 (identical genotype columns)
  are never partners. The default `partner = "strongest"` takes, for a randomly
  drawn anchor SNP, its **highest-r2** in-window partner (so the realized r2 skews
  toward `r2_max`; flanks are searched strongest-first for `"indirect"`);
  `"random"` takes a uniformly random in-window partner.
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
         phase = c("coupling","repulsion"), dist = "geometric",
         orthogonal = FALSE, a = NULL, d = NULL)
dominance(sim, prop, same_as_add = TRUE, n_qtn = NULL, qtn = NULL, dist = "geometric", effect = NULL)
epistasis(sim, prop, n_pairs = NULL, interaction = 2, interaction_type = "a", qtn = NULL, effect = NULL, dist = "geometric")
vqtl(sim,      prop, same_as_add = TRUE, n_qtn = NULL, qtn = NULL, dist = "geometric")
```

- `prop`: scalar or length-`n_traits` vector.
- `dist`: within-layer effect distribution; default geometric. `effect` overrides with
  a geometric base or an explicit series (v1 `sim_method = "custom"`), or a
  length-`n_traits` list of these (one per trait) in `additive()`, `dominance()`
  and `epistasis()`; rejected under multi-trait `"pleiotropy"` (DECISION-023).
- `additive(orthogonal = TRUE, a =, d =)`: the orthogonal genotypic model
  (Modeling convention, §2; DECISION-020). `a`/`d` are per-locus additive effects
  and dominance deviations (scalar or length-`n_qtn`); `effect` is rejected in this
  mode and `a`/`d` require `orthogonal = TRUE`. Unlike `dominance()` below, the
  degree of dominance `d/abs(a)` is meaningful here because the whole genotypic
  value is scaled together. Every `d != 0` locus must have a heterozygote.
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
  takes an `n_pairs x interaction` matrix (or a vector, filled by row) shared by
  every trait, or a length-`n_traits` **list of such matrices**, one per trait with
  the same number of sets; a list of sets is *not* accepted (it used to be read
  that way and silently replicated to every trait, giving a genetic correlation of
  1). `n_qtn`/`n_pairs` follow from it (a supplied `n_qtn`/`n_pairs` is ignored
  with a warning), fixed loci are exempt from `vary_qtn`, and `dominance()` /
  `vqtl()` record `same_as_add = FALSE` when `qtn` is given.
- `additive(phase=)`: `"coupling"` (default) or `"repulsion"`. Repulsion
  alternates effect signs so linked increasing/decreasing alleles oppose; the
  realized genetic variance is still `prop` (pinned by the partition), only the
  sign structure and cross-locus covariance change.

### 4.3 `complex_phenotypes()` — combine architectures

```
complex_phenotypes(..., h2, reps = 1, resid_cor = NULL)   # ... = two or more phenotype_sim objects
```

- Combines the inputs' genetic values, **weighted by their genetic variances**.
- Requires identical genotype data, marker order, individual subset, trait count,
  replication count, and trait means.
- Preserves every replication, scales the combined genetic value to requested
  variance `h2`, and adds a common residual; inputs' individual residuals are discarded.
  `reps` (default 1) divides that common residual by `sqrt(reps)` as in §4.1; the inputs' own
  `reps` are ignored with their residuals.
- **Seed handling (O4):** if inputs were built with different seeds, the **first
  input's seed is used** and a warning is emitted.
- Recreates partial pleiotropy: combine a `"pleiotropy"` model with an `"independent"`
  model.
- Inputs must be **complete** models (an input whose `h2` is not filled by its layers
  errors, exactly as `phenotypes_long()` does). The result is **terminal**: a layer
  added afterwards errors instead of being ignored, and it carries no per-input
  state (`mediation_split()` is `NULL`, no one-call hint).

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
between heterozygous sites in an outbred sample is not realistic. For known phase use
`population_from_haplotypes()` (below); substantially heterozygous *unphased* founders are
not supported for realistic multi-generation recombination studies.

**`population_from_haplotypes()` / `haplotypes()` (DECISION-039).**
`population_from_haplotypes(cis, trans, map, ids = NULL, pool = NA_character_,
individuals_in_rows = FALSE)` builds a `Population` from two known-phase 0/1 haplotype matrices
(markers x individuals by default, the `Population` layout; `individuals_in_rows = TRUE` takes the
transposed export layout). Entry 1 is the counted (+1) allele and dosage = `cis + trans - 1`, the
`as_population()` encoding; `cis`/`trans` carry no maternal/paternal meaning. The map is validated as in
`as_population()` (`snp`, `chr`, `pos`, `cm`; `cm` ordered, in centiMorgans); `map$counted` is recorded
only if the supplied map carries it (never inferred from `allele`). `haplotypes(pop)` returns
`list(cis, trans)`, integer markers x individuals with dimnames (`snp`, ids); `as_population()` ->
`haplotypes()` -> `population_from_haplotypes()` reproduces the `Population` exactly.

**`cross()` / `selfcross()` / `double_haploid()`** take single-individual `Population`s
(`mother`, `father`, or `parent`; subset a larger one with `x[i]`) and return `n` progeny
as a new `Population`. Recombination follows the Karlin & Liberman count-location
process — with `interference = NULL` (the default) the crossover count is Poisson with mean equal to
chromosome length in Morgans, positions uniform along it, chromosomes independent — matching isqg exactly
(DECISION-012). `double_haploid()` progeny are fully homozygous by construction (a single
gamete, duplicated); `selfcross()` halves heterozygosity per generation; `cross()`
combines one gamete from each parent. All three draw every random quantity in R, in
isqg's exact order, before calling the Rust core, which is a pure function of those draws
and never calls an RNG (DECISION-012) — `seed` (or an ambient `set.seed()`) fully
determines the outcome.

**Crossover interference (optional, DECISION-041).** `cross()`, `selfcross()`, `double_haploid()`,
`mate()` and `crossbreed()` take a trailing `interference = NULL`, as do (by propagation) every other
function that runs meiosis: `single_seed_descent()`, `bulk()`, `pedigree()`, `recurrent_selection()`,
`cross_usefulness()`, `combining_ability(method = "simulated")` (an error with `method = "expected"`) and
`progeny_test()`; `NULL` is the Poisson model and the
isqg stream above, bit-identical to versions without the argument. `interference = list(nu = , p = )`
(`1 <= nu <= 1e6`, `p` in [0, 1], `p` default 0) selects the two-pathway gamma model. Model: bivalent chiasmata
of intensity 2 per Morgan = a non-interfering Poisson pathway (share `p`) + a stationary renewal pathway
with Gamma(shape `nu`, rate `2 nu (1-p)`) gaps (share `1-p`); a gamete keeps each chiasma with
probability 1/2 (no chromatid interference), so the expected number of crossovers per Morgan stays 1 for
every `nu`, `p`. Then `r(d) = (1 - P0(d))/2`, `P0(d) = exp(-2 p d) [1 - F_e(2 nu (1-p) d)]`,
`F_e(y) = F_{nu+1}(y) + (y/nu)(1 - F_nu(y))` (`F_a`: Gamma(a, 1) cdf, `d` in Morgans);
`r(d)` is Haldane for `nu = 1` or `p = 1`, and `nu = 2.6, p = 0` is within 0.001 of Kosambi over
0-1 M. Tests: the recombination fraction and the pair correlation against these closed forms, and the
expected number of crossovers per Morgan unchanged. The draws are made in R (their own stream) and the
sorted chiasma positions go to the unchanged Rust core.

**Per-call cost and batching (DECISION-040).** Every crossing function runs through one batched
integer-strand call into the Rust core (`mate_many_core()`); `mate()` executes all rows of a plan in one
call, drawing row by row in plan order, so a plan equals the same sequence of `cross()` /
`selfcross()` / `double_haploid()` calls bit for bit. Many doubled-haploid families at once:
`mate(data.frame(mother = ids, father = ids, n = 100, design = "dh"), pop)`.

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
- Each layer derives a deterministic sub-seed from `(seed, layer type, occurrence of
  that type)`, where *occurrence* is the number of earlier layers of the same type
  (DECISION-033). Replication `r` of a `vary_qtn` layer and each trait's residual
  use the labels `"<type>_rep<r>"` and `"residual_t<t>"` under the same rule. The
  label is reduced with a position-sensitive rolling hash, so labels that differ
  only by a permutation of characters (`..._rep12` / `..._rep21`, `residual_t12` /
  `residual_t21`) get different sub-seeds; the earlier character-code sum made those
  replications and traits identical. The layer *index* is deliberately not used: it
  would break the invariance below.
- Adding/removing/reordering a layer must not change the realized values of layers
  **of other types**; tests assert this. Inserting another layer of the *same*
  type before an existing one shifts that layer's occurrence index, so its draws
  change (the second `additive()` is keyed by occurrence 1 whatever precedes it).
- The caller's RNG state is left untouched (sub-seeds are set and restored).
- `resid_cor = NULL` leaves the residual draw untouched (bit-identical). Otherwise each trait's unit
  residual is drawn under its own `residual_t<t>` sub-seed exactly as before, the columns are mixed
  through `chol(R)`, each column is re-standardized (mean 0, exact unit variance) and scaled by
  `sqrt(1 - sum(prop))`, so per-trait residual variance and the realized h2 / `var_budget` are
  unchanged and only correlation is induced (DECISION-042).
- `reps` does not alter the RNG stream: the residual of each trait is drawn under the same sub-seed
  `residual_t<t>` with the same number of draws, and is scaled by `1/sqrt(reps[t])` afterwards, so a
  seeded `reps = r` phenotype equals the `reps = 1` genetic value plus the `reps = 1` residual divided
  by `sqrt(r)`.
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
- Random QTNs are drawn only from markers with a non-constant dosage column
  (monomorphic and all-heterozygous markers are skipped); a requested nonzero layer
  that has no usable design variance errors.
- At least three individuals are required (with two, the exact-variance
  standardization makes the genetic value and residual collinear).
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
# ("ld" is a two-trait architecture, so both inputs use n_traits = 2)
pleio <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "pleiotropy",
                            n_traits = 2, seed = 10) |> additive(0.4, n_qtn = 3)
ld    <- simulate_phenotype(SNP55K_maize282_maf04, architecture = "ld",
                            n_traits = 2, ld_type = "indirect", seed = 11) |>
         additive(0.3, n_qtn = 3)
both  <- complex_phenotypes(pleio, ld, h2 = 0.5)   # warns: differing seeds, uses 10

# Multi-generation cross, phenotyped directly (§4.6)
pop <- as_population(SNP55K_maize282_maf04, individuals = c("33-16", "38-11"))
f1  <- cross(pop[1], pop[2], n = 1, seed = 1)
f2  <- selfcross(f1, n = 200, seed = 2)
ph  <- simulate_phenotype(f2, seed = 3) |> additive(prop = 0.5, n_qtn = 3)

# 5 records per entry: target residual variance V_E/5, single-record h2 = 0.4
ph <- simulate_phenotype(SNP55K_maize282_maf04, h2 = 0.4, n_qtn = 10, seed = 1, reps = 5)
ph   # prints the entry-mean realized H2 (target 0.4/(0.4+0.6/5) = 0.77) and the single-record H2 (target 0.4)
```

`h2` is the single-record target V_G/(V_G+V_E); with `reps` records per entry the entry-mean target is
V_G/(V_G+V_E/reps). The printed values are realized, Var(G)/Var(y) with
Var(y_bar) = V_G + V_E/reps + 2Cov(G,e)/sqrt(reps), so they match the targets up to sampling covariance.

---

## 11. v1 → v2 Argument Mapping (user migration reference)

> These are equivalences for users porting v1 scripts to the v2 grammar. The frozen
> `create_phenotypes()` does NOT translate to the grammar at runtime (DECISION-008) — it
> keeps its own v1 code paths. `big_add_QTN_effect`, `cor_res`, and
> `architecture = "partially"` live only in the legacy function; `cor` exists in
> both (in the grammar it is the PleioArch correlation control, §13).

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
| `degree_of_dom` | — (no grammar argument) | washed out by per-component variance scaling; the dominance/additive variance ratio is `prop_dom / prop_add`, and a per-locus degree `d/abs(a)` is available in `additive(orthogonal = TRUE, a =, d =)` (DECISION-020) |
| `epi_interaction` | `epistasis(interaction=)` | 2-way default |
| `sim_method = "geometric"` | `dist = "geometric"` | default |
| `sim_method = "custom"` | `effect = <series>` | |
| `type_of_ld` | `ld_type` | modernized (O5) |
| `ld_max` / `ld_min` | `r2_max` / `r2_min` | modernized (O5) |
| `cor` | `cor` (PleioArch, §13) | v1's buggy Cholesky `cor` replaced by the PleioArch engine for any number of traits (DECISION-013); the name `cor` is kept. Scalar or full matrix. v1 `cor` also survives in the frozen legacy fn |
| `cor_res` | `resid_cor` | residual (environmental) correlation between traits, `NULL` = independent (DECISION-042); scalar or full matrix, also on `complex_phenotypes()`. v1 `cor_res` scales the residual covariance as `sqrt(V) R sqrt(V)`; the grammar mixes standardized draws through `chol(R)` and re-standardizes each trait, so per-trait variance stays its `h2` target |
| `rep` | `n_reps` | |
| `vary_QTN` | `vary_qtn` | |
| `constraints = list(maf_above, maf_below, hets)` | `filter_geno(maf_above=, maf_below=, hets=)` | applied to the genotype once, upstream of every architecture and layer (additive/dominance/epistasis/vqtl); takes the data frame, matrix or a `Population`. v1 filtered only the randomly drawn QTNs, not the LD partner markers |
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

## 13. PleioArch Algorithm — Controlled Genetic Correlation for Pleiotropy

> Reference implementation: `context/PleioArch-main/Functions/simulateEffects.R`
> DECISION-007: adopted as the effect-generation engine for `architecture = "pleiotropy"`.

### Motivation
v1 "pleiotropy" merely shares QTN loci across traits; the resulting genetic correlation
is an emergent by-product of allele-frequency differences and cannot be set precisely.
The PleioArch algorithm draws allelic effects from correlated distributions whose
cross-trait covariance and variances equal their targets in expectation, so the user's
`cor` is the ratio of expected moments. The realized genetic correlation is a random
ratio: it converges to `cor` as the shared QTNs and the individuals grow (a sample
correlation over n individuals cannot beat its sampling spread) when the causal loci are in
approximate linkage equilibrium, is attenuated toward 0 on average with few shared QTNs,
and under strong LD need not converge (complete LD: every realized value is ±1). See
DECISION-023, Scope / limits, which also extends this to dominance and epistasis.

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
Pleiotropic effects are drawn from a **multivariate normal** (eigen (symmetric)
square root of the per-SNP covariance matrix scaled by QTN count, so exactly singular
but feasible covariances are still sampled; bivariate when `n_traits = 2`):
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
  constraint). Deterministic assembly (eigen square root, matrix multiply) may move to Rust.
- `additive(prop=, ...)` calls this engine when the parent `phenotype_sim` carries
  `architecture = "pleiotropy"`, passing the resolved QTN effects directly.
- v1 `create_phenotypes(architecture = "pleiotropic")` does **not** map to PleioArch — it
  is frozen legacy (§4.4, DECISION-008) and keeps v1's shared-loci behavior. PleioArch is
  the new grammar's pleiotropy engine, controlled by `cor` (DECISION-010, amended).
- **Any `n_traits`**: the engine generalizes directly (DECISION-013). The
  bivariate normal becomes an `n × n` multivariate normal with

  ```
  Sigma[i,i] = pi_i  × V_i                    (pleiotropic share of trait i)
  Sigma[i,j] = cor_ij × sqrt(V_i × V_j)       (whole covariance: specific loci are
                                               independent across traits)
  ```

  Trait-specific effects remain univariate normal with variance `(1 − pi_i) × V_i`, so
  each trait's expected genetic variance is still `V_i` and every pair targets `cor_ij`
  (in the sense above).
  For two traits this reduces algebraically to the bivariate reference implementation.
  The `cor² ≤ pi_i × pi_j` constraint generalizes to "`Sigma` must be positive
  semi-definite", checked by eigenvalue and raised as an error naming the smallest
  eigenvalue. The Cholesky-on-genetic-values fallback is **removed**: it never
  controlled individual QTN effects and is no longer needed.
