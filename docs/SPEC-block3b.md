# SPEC-block3b.md — engine dependencies in breedingDesigner's methods catalog

> **Status: DRAFT — scope for maintainer decision.** Date 2026-09-27. Package
> version 2.0.0.9000 (`master`). Scopes the eight items under `TODO.md` → "Block 3B —
> Engine work requested by breedingDesigner" → "Later — engine dependencies in BD's
> methods catalog (not yet specced; confirm scope)". Source of the requests:
> breedingDesigner (BD, `../breeding_designer`) `docs/BREEDING_METHODS_CATALOG.md`
> ("Engine:" notes and the "Engine (simplePHENOTYPES) dependencies, consolidated" list),
> with BD `specs/SPEC-0006` (computed criteria), `SPEC-0007` (mate allocation),
> `SPEC-0011` (reciprocal recurrent selection) and `R/rrs.R` (how BD composes GCA today).
> Companion docs here: `docs/BACKEND_CONTRACT.md` (the exported surface BD calls),
> `docs/THEORY_REVIEW.md` (the rubric every item must pass), `docs/DECISIONS.md`
> (latest DECISION-023), `docs/ROADMAP.md` §4 (polyploids) and §8b (Python).
> Nothing in this file is implemented; no code, tests or decisions are changed by it.
> Every equation below is either this package's own derivation from the stated model
> or is attributed to a source verified this session (§12); formulas are **not**
> attributed to papers whose text was not checked.

## 0. Summary

Six of the eight items are quantitative-genetics verbs that compose from primitives the
engine already exports — `cross()`, `genotypic_value()` (the fixed-scale A + D total
that closed BD's RRS-1b blocker), `additive_value()`, `phenotype_value()`,
`select_ind(on = <vector>)`, the family methods, `g_matrix()` and the average-effect
machinery behind `on = "bv"` (DECISION-019) — plus one piece of infrastructure the
engine lacks and three of the items need: a **pedigree carried by `Population`**
(mother/father per individual, propagated through mating, subsetting and pooling).
With that foundation, combining ability (item 1), progeny/family testing (item 4) and
crossbreeding bookkeeping (item 6) are mostly deterministic scoring around the existing
Rust-backed crossing core, and pedigree BLUP (item 3) becomes possible without new
dependencies because a simulation *knows* its variance components (no REML). Two items
are recommended for deferral: the polyploid model (item 7) is a v3-sized change and —
importantly for BD — is **not** needed for its allotetraploid `Wheat_div`, whose
diploid-per-subgenome coding is the correct model of disomic inheritance; and the Python
package (item 8) is not a genetics item, is already sequenced after CRAN (ROADMAP §8b),
and the only near-term action is to stop BD presenting generated Python as runnable.
Recommended order: pedigree foundation → `mate()`/`mating_design()` → combining ability →
families/progeny test → tandem/culling → marker selection → BLUP/EBV → crossbreeding
wrapper; defer 7 and 8. All additions are additive exports (a minor bump under
`BACKEND_CONTRACT.md`'s rules; maintainer decision 2026-09-28: no `2.1.0` — the
user-facing release stays `2.0`, and Block 3B ships in the development version
`2.0.0.9001`), RNG stays in R, nothing new moves to Rust.

## 1. Shared foundation

### F1. Pedigree in `Population` — new

**Gap.** `Population` holds `map`, `cis`/`trans`, `ids`, `origin` (SPEC §5;
`R/cross_population.R:103`). Parentage survives only as a batch-level `origin` string
(`"cross(A x B)"`) and as id prefixes the scheme wrappers embed (`"P1_ped_g1_3"`,
`"cyc1_x2_5"`, `R/select_schemes.R:312-344`), which `c.Population()` may rename
(`make.unique`). So there is no reliable way to recover which individuals are full sibs,
which sire a half-sib family shares, an individual's breed composition, or a numerator
relationship matrix. Items 3, 4 and 6 all need one of those; BD's `R/rrs.R` tracks
"which candidate produced which hybrid" by loop structure only.

**Proposal.** Add an optional `pedigree` element: a data frame with one row per
individual — `id`, `mother`, `father` (ids, or `NA` for founders), `generation`
(integer, founders 0) and `pool` (character, the founder pool a founder came from;
`NA` for progeny). Semantics:

- `as_population()` writes founders with `NA` parents (`pool` from a new optional
  `pool =` label, default `NA`).
- `.mate()` (`R/cross_mating.R:69`) appends the progeny rows; `selfcross()` writes the
  parent in both slots, `double_haploid()` likewise (a DH is a selfed gamete, doubled —
  the A-matrix rule for a DH is F = 1, see item 3).
- `[.Population` subsets rows but **keeps ancestor rows** (a pedigree is only useful
  with its ancestors), `c.Population()` unions rows (ids must be unique across the
  pooled set — today's `make.unique` renaming would silently break parent links, so
  `c()` should error on id collisions when pedigrees are present, or `mate()` should
  namespace ids; decision D2).
- Backward compatible: a `Population` without `pedigree` behaves exactly as now; every
  accessor that needs it errors with a clear message rather than guessing from ids.

**Accessors.** `pedigree(pop)` (the frame), `families(pop, by = c("full_sib",
"maternal_half_sib", "paternal_half_sib", "selfed"))` → a factor aligned with `pop$ids`
for `select_ind(family =)`, `a_matrix(pop, ids = NULL)` (item 3),
`breed_composition(pop)` (item 6).

**Theory.** None beyond bookkeeping; but the pedigree is what makes the A-matrix,
family structure and expected breed fractions *exact* rather than string-parsed.

**Validation.** Parents recorded by `cross()`/`selfcross()`/`double_haploid()` match
the `origin` string and the argument objects; `families()` on a `mating_design()` output
(item 6) reproduces the design table exactly; subsetting/pooling round-trips; a
`Population` built by `as_population()` from a dosage frame has no non-founder rows.
`test-isqg-parity.R` and every scheme test stay bit-identical (the pedigree is metadata;
no draw order changes).

**Effort** M. **Risk** low (additive field). **Decision D1** (adopt), **D2** (id
collision policy).

### F2. One mating executor — new (specified under item 6, used by items 1, 4, 6)

`.intermate()` (`R/select_schemes.R:328`), `cross_usefulness()`'s `.make_family()`, BD's
`.rrs_intermate()`, `.rrs_hybrid_pop()` and `.rrs_gca_select()` all hand-roll "loop over
pairs, `cross()`, relabel, `c()`". A single exported `mate(plan, ...)` that executes a
`{mother, father, n}` plan (BD SPEC-0007's output shape) over one or several
`Population`s, writes the pedigree, and prefixes ids deterministically, removes those
reimplementations and is the primitive the crossbreeding systems, testcross designs and
progeny tests are built from. Specified in §7.

---

## 2. Item 1 — combining-ability scorer (GCA / SCA, testcross merit)

### What BD needs
Catalog §B: "Half-sib recurrent selection (GCA, w/ tester) … **Engine:** combining-
ability score"; "Hybrid development (single / 3-way / double cross) … **Engine:**
SCA/heterosis (A+D value)"; consolidated dependency 2 "Combining-ability (GCA/SCA) &
testcross scoring — half-sib/full-sib recurrent selection, hybrid development". BD
SPEC-0011 RRS-1b: "score GCA on that total value so dominance/SCA drives selection and
`rrs_metrics()` reports the heterosis trend". Today BD composes GCA in
`R/rrs.R::.rrs_gca_select()`: sample `n_testers` from the opposite pool, `cross()` each
candidate to each tester (`n_progeny` each), GCA = mean `additive_value()` of the
hybrids; additive models only, dominance refused; no residual (an "environment-free"
estimate, as its header says). The engine side of the RRS-1b blocker —
`genotypic_value(x, qtn, a, d)` — shipped (TODO Block 3B "Next"); what is still missing
is the scorer that uses it.

### What already exists
- **Exists:** `cross()` between any two `Population`s sharing a map (cross-pool mating
  verified by BD); `genotypic_value()` — fixed-scale per-se `G = A + D`, whose docs
  already state that per-se value is *not* testcross merit and that hybrid programs must
  score the realized progeny; `additive_value()`; `phenotype_value()` (fixed residual,
  but **additive genetic part only**); `select_ind(on = <numeric>)` and
  `optimum_contribution(merit = <numeric>)` accept any external score; the average-
  effect helper `.avg_effect(a, d, p) = a + d(1 − 2p)` and `.layer_scaled_effects()`
  (reconstructs realized-scale `a_j`, `d_j` from a `phenotype_sim`; internal).
- **Partially exists:** the *simulated* GCA (BD's composition), additive only.
- **New:** an exported scorer; the *expected* (analytic) testcross value; the tester-
  frequency average effect; SCA; a residual on the *total* value; design bookkeeping
  (which candidate each hybrid came from → F1/F2).

### Proposed API
```r
combining_ability(candidates, testers, qtn, a, d,
                  design  = c("topcross", "factorial", "diallel"),
                  method  = c("expected", "simulated"),
                  n_progeny = NULL, h2 = NULL, var_e = NULL, ref = NULL,
                  seed = NULL)
```
- `candidates`, `testers`: `Population`s on one map (for `"diallel"`, `testers` is
  `NULL` and all candidate pairs are crossed). `qtn`, `a`, `d`: the frozen architecture
  on the realized scale, exactly as `genotypic_value()` takes them (`d = 0` for a purely
  additive locus). Epistatic architectures are refused (see Theory).
- `method = "expected"`: no simulation, no RNG — the expected progeny value of every
  candidate × tester (or pair) from the known effects and the parents' genotypes.
- `method = "simulated"`: realize each family with the crossing core (`mate()`),
  score with `genotypic_value()`, optionally add an independent residual
  `e ~ N(0, var_e)` per hybrid individual (exactly one of `h2`/`var_e`, `ref` as in
  `phenotype_value()`); this is a finite-sample *estimate* (tester + Mendelian +
  environmental sampling).
- Returns a `combining_ability` object: `gca` (named vector over candidates; for
  `"factorial"`/`"diallel"` also over testers), `sca` (candidates × testers matrix, or
  the symmetric pair matrix for a diallel), `testcross` (the per-cross means and, for
  `"simulated"`, sizes), `grand_mean`, attributes `method`, `design`, `var_e`, `seed`.
  `ca$gca` feeds `select_ind(sim, on = ca$gca)` / `optimum_contribution(merit =
  ca$gca)`; the hybrid `Population`s (with pedigree) are returned for `"simulated"` so
  `rrs_metrics()`-style heterosis series can score them.
- Lives in a new `R/select_combining.R`; contract addition under "Modern methods".
- Optional companion (**D5**): `template_effects(sim, trait = 1)` → `data.frame(snp, a,
  d)` on the realized scale (a thin export over `.layer_scaled_effects()`), so BD's
  `.freeze_from_template()` stops re-deriving effects from `qtn_table()`.

### Genetic theory
*Definitions.* Following the classical two-way model — general combining ability as the
average merit of a parent in cross combinations and specific combining ability as the
deviation of a particular cross from what its parents' GCAs predict (Sprague & Tatum
1942; Griffing 1956 for the diallel formulation) — the package states, as its own
convention, for the set of crosses the design produces:
```
E[Y_ik] = μ + g_i + g_k + s_ik,   Σ_i g_i = 0 over candidates (and over testers),
                                   Σ_k s_ik = 0 for every i
```
with `μ` the mean over the design's crosses, `g` the row/column means minus `μ`, and
`s` the residual. For a topcross to one pooled tester set, `sca` is per cross relative
to the candidate's GCA and the tester mean.

*Expected testcross value (package derivation).* At one locus with genotypic values
`−a / d / +a` for gene content 0/1/2 (the `genotypic_value()` convention), let a
candidate have gene content `x_i` and the tester *gametes* carry the counted allele with
probability `p_T` (a single tester's `x_T/2`; a tester population's allele frequency
under random mating). With `g = x_i/2`, the expected progeny value is
```
E[G] = g p_T a + [g(1 − p_T) + (1 − g) p_T] d − (1 − g)(1 − p_T) a
     = (x_i / 2) · α_T + c_T,     α_T = a + d(1 − 2 p_T),   c_T = p_T d − (1 − p_T) a.
```
So the candidate's expected testcross merit is linear in its gene content with slope
half the **tester-referenced average effect** `α_T` — DECISION-019's
`α = a + d(q − p)` with the tester's allele frequency in place of the population's — and
`GCA_i = ½ Σ_j α_{T,j} (x_ij − x̄_j)` over the candidate set. Consequences the review
should confirm: (i) a candidate's GCA is *not* intrinsic — it changes with the tester
(`p_T`); (ii) for an inbred tester `p_T ∈ {0, 1}`, `α_T = a ± d`; (iii) when the testers
are a random sample of the candidates' own random-mating population (`p_T = p`),
`GCA_i = ½ A_i`, half the DECISION-019 breeding value — the textbook half-sib result
and an exact known answer (§ Validation, 2); (iv) for `d = 0` the ranking equals the
`additive_value()` ranking (mid-parent argument), which is the "honest caveat" BD's
`rrs.R` already records; (v) per-locus expectations depend only on marginal gamete
frequencies, so **linkage does not enter the expected value** under the A + D model —
it enters the *variance* of a realized family, which `"simulated"` captures through the
meiosis core; (vi) an epistasis layer has no per-locus `a`/`d` and its expectation needs
the joint gamete distribution, so epistatic architectures are refused, as
`.breeding_value_matrix()` refuses them (DECISION-019).

*Pairwise expectation (SCA).* For two parents with `g_i = x_i/2`, `g_k = x_k/2`:
`E[G_ik] = g_i g_k a + [g_i(1 − g_k) + (1 − g_i) g_k] d − (1 − g_i)(1 − g_k) a`, summed
over loci; exact (not just in expectation) for inbred parents (`x ∈ {0, 2}`). With
`d = 0` everywhere `s_ik ≡ 0` — a purely additive model has no SCA by construction
(BD's `rrs.R` states this; it follows from the mid-parent identity).

*What is exact vs estimated.* `"expected"` is a true, idealized quantity (the
conditional expectation given the parents' genotypes and the frozen effects) —
legitimate in a simulation and labelled as such (rubric S5). `"simulated"` is an
estimate whose expectation is the `"expected"` value and whose error has tester,
Mendelian and (if `var_e > 0`) environmental components; Monte-Carlo agreement between
the two is itself an acceptance test. Neither is a per-se value; neither is the
transmissible breeding value unless `p_T = p`.

*Pitfalls the rubric would check.* Scale consistency of `a`, `d` (realized scale, as
`genotypic_value()`; raw `qtn_table()` effects are not on it — hence D5); the residual is
per hybrid *plot* (individual), `var_e` frozen across cycles as in DECISION-021;
monomorphic loci contribute constants only; `direction`; that the SCA sum-to-zero
convention is stated and tested; that RRS is documented as improving GCA — BD's
`rrs.R` reads Melchinger & Frisch (2023) that way, and this spec adopts that reading
without adding claims about SCA under GCA selection (not checked here).

### Validation plan
1. **Single-locus known answer.** Inbred candidates `AA`/`aa`, inbred tester `aa`,
   `a = 1`, `d = 0.5`: expected testcross values `d = 0.5` and `−a = −1`;
   `gca_AA − gca_aa = a + d = 1.5` exactly (`p_T = 0 → α_T = a + d`).
2. **Half the breeding value.** `method = "expected"` with `testers = candidates` on a
   random-mating (HWE) panel: `gca == 0.5 * bv` from `.breeding_value_matrix()` to
   floating-point tolerance (same centering, same `p`), for an A + D architecture.
3. **Additive ⇒ mid-parent.** With `d = 0`, `rank(gca) == rank(additive_value())` and
   `sca == 0` identically, for any tester.
4. **Simulated → expected.** Over ≥ 30 seeds, the mean `"simulated"` GCA per candidate
   is within `3·SE` of the `"expected"` GCA; the spread shrinks with `n_progeny`.
5. **Design identities.** `sum(gca) == 0`; `rowSums(sca) == 0`; diallel `sca` symmetric.
6. **RNG.** `"expected"` consumes no RNG (`.Random.seed` unchanged); `"simulated"`
   is reproducible from `seed` and identical to a hand-built `cross()` loop with the
   same draw order (exact).
7. **Independent review must confirm:** the `α_T` derivation and its sign convention;
   the SCA centering; the epistasis refusal; the S5 labelling ("expected" is idealized
   truth, "simulated" is an estimate); no per-se value is presented as testcross merit.

### Dependencies / risks / effort
F2 (`mate()`) for the simulated design; F1 for family provenance of the hybrids (nice,
not required). Risk: users pass raw `qtn_table()` effects (wrong scale) — mitigated by
D5. **Effort M** (expected: S; simulated + residual: S; docs/tests/review: M).

### Decisions needed
- **D3** Accept the `(qtn, a, d)` frozen-architecture interface (consistent with
  `genotypic_value()`) rather than a `phenotype_sim` template — yes/no.
- **D4** Ship `"expected"` and `"simulated"` together, or `"expected"` first.
- **D5** Export `template_effects()` so BD stops re-deriving realized-scale effects.
- **D6** Should `phenotype_value()` gain `d =` (fixed-scale A + D phenotype) so a single
  residual model serves items 1, 3 and 4, instead of a private residual here?

---

## 3. Item 2 — marker-based selection (MAS / gene pyramiding) and a marker index (MARS)

### What BD needs
Catalog §D: "Marker-assisted selection (MAS) / gene pyramiding — Major-gene tracking.
**Engine:** marker-based selection"; §B: "Marker-assisted recurrent selection (MARS) —
Marker-index recurrent selection. **Engine:** marker index"; consolidated dependency 6
"Marker-based selection (MAS/MABC/MARS) — foreground/background & marker-index".

### What already exists
- **Exists:** `mabc_select()` — foreground hard filter (`donor_carrier` /
  `donor_homozygote`), flanking/recombinant selection, background recovery ranking,
  seeded tie-break, diagnostic table — but tied to a `recurrent`/`donor` founder pair
  (informativeness is defined by the founders); `recurrent_parent_recovery()`;
  `select_ind(on = <numeric>)`; and, decisively, **`additive_value(x, markers,
  weights)` already *is* a marker index** `Σ_j w_j x_ij` when given marker names and
  weights — the MARS score needs no new arithmetic.
- **New:** a foreground filter whose favourable allele is given directly (no founders),
  per-marker requirements, staged pyramiding (`min_markers`), and a secondary ranking
  among feasible candidates. The *estimation* of marker effects for a MARS index is
  BD's SPEC-0006 (`criterion_gblup()`) or item 3's `predict_ebv()` (which back-solves
  marker effects); the engine's own effects are the truth, so an index on them is an
  *oracle* MARS.

### Proposed API
```r
marker_select(pop, markers, favorable = 1L,
              requirement = c("carrier", "homozygote"),   # scalar or per marker
              min_markers = length(markers),              # pyramiding stage
              n = NULL, prop = NULL,
              rank_on = NULL,                             # NULL | numeric | function(pop)
              direction = c("high", "low"), seed = NULL)
```
- `favorable`: `+1` or `−1` per marker — which homozygote (in `pop`'s −1/0/1 coding)
  carries the favourable allele; the same allele-coding caveat as `mabc_select()`
  (documented; not checkable from dosages).
- Feasible = satisfies `requirement` at ≥ `min_markers` of `markers`; among feasible,
  rank on `rank_on` (a numeric vector, e.g. `additive_value()` on an index or a GEBV;
  a phenotype; `recurrent_parent_recovery()`), then a seeded tie-break; error if fewer
  than `n` feasible. Returns a `Population` with a diagnostic attribute (per-marker
  state, `n_favorable`, feasible, rank, selected), mirroring `mabc_select()`.
- No new `marker_index()` function: document `additive_value()` as the marker index
  (D8 asks whether a documented alias is wanted for discoverability / BD's DAG).
- Implementation: refactor `.mabc_markers()`, `.mabc_recovery()` and the scoring into
  shared marker utilities; `mabc_select()` becomes the founder-informed special case.

### Genetic theory
With known major genes the process is Mendelian and mostly deterministic: after selfing
a heterozygote, `P(homozygous favourable) = 1/4` per unlinked locus in the F2 and
`1 − (1/2)^k` homozygosity after `k` selfing generations; `m` unlinked targets need
`n = ln(0.05) / ln(1 − (1/4)^m)` F2 plants for a 95 % chance of one full pyramid
(package derivation). For **linked** targets the joint probabilities depend on the
recombination fraction, which the meiosis core realizes under the count-location
(no-interference) model, so Haldane's map function `r = ½(1 − e^{−2d})` (`d` in
Morgans) gives the exact expected recombinant fractions to test against. The marker
index on estimated effects is the Lande & Thompson (1990) programme — indices combining
marker and phenotypic information, whose efficiency depends on the genetic parameters
and the scheme (their abstract; the specific formulas are not reproduced here) — and the
engine's oracle version (true effects on the causal loci) is its upper bound: ranking
on `additive_value()` with all causal effects equals ranking on `on = "bv"` for an
additive model. Pitfalls: label true-effect indices as oracle (S5); never call
`marker_select()` MABC (background recovery is `mabc_select()`'s job); `favorable`
orientation is unverifiable from dosages; `min_markers < length(markers)` must be
reported as a staged pyramid, not a full one.

### Validation plan
1. **Unlinked Mendelian ratios.** F2 (`selfcross`) of a double heterozygote, 2 000
   plants × 20 seeds: feasible fraction under `"homozygote"` ≈ 1/16 and `"carrier"`
   ≈ 9/16 within binomial `3·SE`.
2. **Linked targets.** Two targets 10 cM apart in coupling: `P(F2 double favourable
   homozygote) = ((1 − r)/2)^2` with `r` from Haldane's function; compare at 5, 10,
   50 cM.
3. **Determinism.** Same seed → same selection; `n` honoured; error when < `n`
   feasible; `min_markers` semantics on a 3-of-4 pyramid.
4. **Equivalence.** With `favorable` set to the donor allele, `marker_select()`'s
   feasible set equals `mabc_select()`'s at the same targets; with `rank_on =
   recurrent_parent_recovery()` the ranking equals `mabc_select()`'s background rank
   (given no flanking markers).
5. **Oracle MARS.** Index on all causal loci reproduces `select_ind(on = "bv")`'s
   selected set exactly (additive model); an index on a random half of the loci gives
   ≤ the oracle response over 30 seeds.
6. **Review must confirm:** the Mendelian expectations; no MABC wording; the oracle
   labelling.

### Dependencies / risks / effort
None hard; shares helpers with `select_mabc.R`. **Effort S–M.**

### Decisions needed
- **D7** Include staged pyramiding (`min_markers`) now or later.
- **D8** Add a documented alias `marker_index()` over `additive_value()` for
  discoverability, or documentation only.
- **D9** Should MARS's estimated weights come from item 3's `predict_ebv()` (engine) or
  BD's SPEC-0006 criteria (BD, rrBLUP in Suggests)? Ties to D12.

---

## 4. Item 3 — BLUP / EBV prediction verb (pedigree + own + relatives)

### What BD needs
Catalog §E: "Combined / index selection (BLUP-EBV) — Pedigree+own+relatives index.
**Engine:** BLUP/EBV prediction verb (GS accuracy is currently approximated)";
consolidated dependency 3 "Prediction verbs (BLUP-EBV, selection index, GEBV/ssGBLUP)
— animal index/progeny/family selection; today GS accuracy is a modeled
approximation". Concretely, BD's canvas projects response with a **user-entered**
accuracy — `inst/app/index.html:1971`: the OCS/GS node's `accuracy` parameter (default
0.6), `0.95` for selection on `gv`, else `sqrt(h2)` — rather than a realized
`cor(prediction, true breeding value)`. BD SPEC-0006 plans `criterion_gblup()` (rrBLUP)
and `criterion_rf()` in BD with per-cycle accuracy reporting, and asks the engine for a
`selection_methods()` manifest.

### What already exists
- **Exists:** `g_matrix()` (VanRaden 2008 method 1, moving or fixed base); the
  `select_ind(on = <vector/function>)` and `optimum_contribution(merit =)` hooks, which
  DECISION-015/016 designate as *the* GS/PS entry point ("none — plugs into the `on`
  hook"); the true breeding value `"bv"` to compute accuracy against (DECISION-019);
  `phenotype_value()` for a frozen-variance observable phenotype; the Lush
  `method = "combined"` (own record + family mean — a two-source selection index, not
  BLUP); the Smith–Hazel index on true BVs.
- **New:** any estimator; the numerator relationship matrix **A** (needs F1);
  reliability/accuracy reporting; multi-trait or single-step (H-matrix) variants.

### Scope question first (D12)
The engine "simulates; it should not take on the analysis side" (ROADMAP §5, on GWAS
wrappers), yet `BACKEND_CONTRACT.md` requires that every *selection computation* live
here, and BD's need is precisely a computation whose result feeds `select_ind()`. The
reconciling observation: in a simulation the variance components are **known** (`h2`,
`var_e` are inputs; `Var(A)` is computable from the frozen architecture), so BLUP with
known parameters is a deterministic linear solve of the mixed-model equations — no
REML, no dependency, no optimizer — and its output is an honest *estimated* criterion
for the `on` hook. That is a simulation verb, not an analysis package. What stays out:
variance-component estimation, Bayesian marker models, random forests (BD SPEC-0006
keeps those behind Suggests).

### Proposed API
```r
predict_ebv(x, pheno, method = c("gblup", "pedigree"),
            h2 = NULL, var_e = NULL, var_a = NULL, ref = NULL,
            K = NULL, candidates = NULL, fixed = NULL, base_freq = NULL,
            ridge = 0)
a_matrix(pop, ids = NULL)                                   # needs F1
prediction_accuracy(ebv, truth)                             # cor(), with a bias slope
```
- `x`: a `Population` (all individuals — phenotyped training set and unphenotyped
  candidates); `pheno`: a named numeric vector over the phenotyped subset (e.g. the
  output of `phenotype_value()`, or a `phenotype_sim`'s trait column). `K`: optional
  precomputed relationship matrix (`g_matrix()` or `a_matrix()`); built from `x`
  otherwise (`base_freq`/`ridge` forwarded to `g_matrix()`).
- Variance components: exactly one of `h2` (with `ref` — the base population — to set
  `var_a` as `phenotype_value()` sets `var_e`) or `var_a`+`var_e`. Frozen across
  cycles by passing the base values, as DECISION-021 recommends.
- Returns a named `ebv` vector over **all** individuals in `x` (candidates included),
  with attributes `reliability` (per individual), `lambda`, `method`, and — for
  `"gblup"` — `marker_effects` (back-solved). `select_ind(sim, on = ebv[sim$ids])`,
  `optimum_contribution(pop, merit = ebv)`, `marker_select(rank_on = ebv)`.
- `selection_methods()` — a small manifest (`id`, `label`, `params`, `citation`) of the
  engine's selection operators for BD's SPEC-0006 UI generation (**D13**).
- Lives in a new `R/select_blup.R`; `a_matrix()` in `R/cross_pedigree.R` (F1).

### Genetic theory
*Model.* `y = Xb + Zu + e`, `u ~ N(0, K σ²_A)`, `e ~ N(0, I σ²_e)`; Henderson's mixed-
model equations with known variances,
```
[ X'X      X'Z          ] [ b ]   [ X'y ]
[ Z'X   Z'Z + λ K⁻¹     ] [ u ] = [ Z'y ],      λ = σ²_e / σ²_A = (1 − h²)/h²,
```
`K = G` (VanRaden 2008, already cited by `g_matrix()`) gives GBLUP, `K = A` pedigree
BLUP. BLUP is the best linear unbiased predictor under the model, and its treatment
under selection is Henderson (1975) (abstract verified: methods for data arising from
selection). Unphenotyped candidates receive predictions through `K` (the `Z` rows are
zero for them; the solution is the conditional expectation given relatives' records).
Reliability `r²_i = 1 − PEV_i / (K_ii σ²_A)` with `PEV_i = C^{uu}_{ii} σ²_e` from the
inverse coefficient matrix (package statement of the standard quantity; fine for the
population sizes the engine handles — a direct solve, no sparse machinery). GBLUP marker
effects are back-solved as `û_m = M' G⁻¹ ĝ / (2 Σ_j p_j(1 − p_j))` — the package's
derivation from `g = M u_m`, `u_m ~ N(0, I σ²_m)`, `G = M M' / (2Σp(1−p))`.

*Numerator relationship matrix.* Coefficients of relationship and inbreeding are
Wright's (1922); `A` is built recursively from the pedigree by the tabular method
(Emik & Terrill 1949), `A_ii = 1 + F_i`, `A_ij = ½(A_{i,m(j)} + A_{i,f(j)})`, with
`A⁻¹` obtainable directly (Henderson 1976) — sizes here allow `solve()`. Selfed progeny:
`F = ½(1 + F_parent)`; a doubled haploid: `F = 1` (`A_ii = 2`); the engine's
`selfcross()`/`double_haploid()` pedigree rows must encode exactly this.

*What is exact / expected / approximate.* With known variance components BLUP is the
exact conditional mean under the (infinitesimal, multivariate-normal) model; the
simulation's truth is a finite set of QTN with possible dominance, so the predictor is
an **estimate** and its realized accuracy `cor(ebv, bv)` is an empirical quantity —
which is exactly what BD needs in place of a constant. Deterministic accuracy formulae
exist (Daetwyler et al. 2008 derive them as functions of the number of records per
effective locus and heritability, per their abstract); the validation checks
*monotone trends* against realized accuracy, not a formula. Assumptions: `K` captures the
causal loci (markers in LD with QTN for `G`; correct pedigree for `A`); a fixed base for
`σ²_A` and `G` across cycles (moving-base `G` and per-cycle `h2` inflate apparent
accuracy — the same trap DECISION-021 closed for phenotypes); records included for all
individuals selection was based on (Henderson 1975), else EBVs are biased.

*Pitfalls the rubric would check.* S5: the EBV is an estimated criterion and must be
computed from an **observable** phenotype (`pheno`/`phenotype_value()`), never from
`"gv"`/`"bv"` (that would launder truth into an "estimate"); the comparator for accuracy
is the DECISION-019 breeding value (under dominance), and with epistasis `"bv"` is
unavailable — say so rather than substitute; singular `G` for genotype-identical
individuals (DH copies, `sample_parents()` duplicates) → `ridge`; `λ` at `h2 → 0`;
alignment of `pheno` names to `x`; that `var_a` is the *base* additive variance (state
which).

### Validation plan
1. **Own-record known answer.** `K = I`, one record each: `ebv_i = h² (y_i − ȳ)` exactly,
   and accuracy of an own-record EBV equals `h` in expectation (`cor(ebv, y) = 1`,
   `cor(y, A) → h`). Executable to floating-point tolerance for the identity part.
2. **A-matrix known answers.** Full sibs 0.5, half sibs 0.25, parent–offspring 0.5,
   selfed progeny of a non-inbred parent `A_ii = 1.5`, DH `A_ii = 2`; and, on a
   pedigree simulated by `mate()`, `mean(G[class])` ≈ `A[class]` per relative class with
   a fixed founder base (tolerance; VanRaden's `G` is the realized counterpart of `A`).
3. **Textbook-size MME.** A five-individual pedigree/phenotype example solved by hand
   (recorded in the test) reproduces `ebv` and `b` to 1e-10.
4. **Unbiasedness under the model.** Regression of true `bv` on `ebv` across ≥ 30
   seeds has slope ≈ 1 (BLUP property `E[u | û] = û`), for GBLUP on an additive trait
   simulated with `phenotype_value()` and `G` from all markers.
5. **Trends.** Realized accuracy increases with training size and with `h2`; pedigree
   BLUP of an unphenotyped candidate with `n` half-sib progeny records approaches the
   progeny-test accuracy `sqrt(n h² / (4 + (n − 1) h²))` (item 4, derivation there).
6. **External cross-check (dev-only, not a test dependency).** `predict_ebv(method =
   "gblup")` equals `rrBLUP::mixed.solve()` with the variance components fixed, on the
   same `G`, to 1e-8 — recorded in `dev/`.
7. **Review must confirm:** MME assembly (which block is which, `λ` orientation); the
   `A` recursions incl. selfing/DH; that no truth leaks into the estimate; frozen-base
   wording; the reliability formula.

### Dependencies / risks / effort
F1 for `a_matrix()`; `"gblup"` needs nothing else. No new `Imports` (base `solve()`/
`chol()`). Risk: scope creep toward an analysis package — held by D12's "known variance
components only" rule. **Effort M** (GBLUP + accuracy + manifest) **+ M** (A-matrix +
pedigree BLUP, after F1).

### Decisions needed
- **D12** Is a known-variance-component BLUP in scope for the engine (recommended), or
  does prediction stay in BD (SPEC-0006, rrBLUP in Suggests) with the engine supplying
  only `g_matrix()`/`a_matrix()`/`prediction_accuracy()`?
- **D13** Ship `selection_methods()` (S) as SPEC-0006 requests.
- **D14** Multi-trait BLUP and single-step (`H`, Legarra et al. 2009) — defer (recommended).

---

## 5. Item 4 — progeny-mean scorer and family-structured phenotyping

### What BD needs
Catalog §E: "Family selection (between-family) — Select on family mean. **Engine:**
family-structured phenotyping + family-mean select"; "Within-family selection";
"Sib selection"; "Progeny testing — Judge parent by progeny mean (dairy sires).
**Engine:** progeny-mean scorer"; consolidated dependency 4. `TODO.md`: "`select_ind()`
already has `within_family` / `among_family` / `combined` methods; the gap is building
families in a design". BD's designer and canvas use no family method at all today.

### What already exists
- **Exists:** `select_ind(method = "within_family" | "among_family" | "combined",
  family = <vector>, h2 =, family_relationship =)` with the Lush index weights
  (`R/select_ind.R:472`); `.self_each()`/`.intermate()` produce family-structured progeny
  (but only as id strings); `cross()` for any mating; `genotypic_value()`/
  `phenotype_value()` for fixed-scale scoring.
- **New:** `families()` (F1); an explicit mating design that *creates* half-/full-sib
  structure (`mating_design()`, F2/§7); the progeny-mean scorer — which is item 1's
  `"simulated"` method with random mates instead of a tester panel.

### Proposed API
```r
mating_design(mothers, fathers, design = c("random", "nested", "factorial"),
              n_crosses = NULL, dams_per_sire = NULL, progeny_per_cross = 1L,
              seed = NULL)                     # -> plan {mother, father, n}
progeny_test(parents, mates, qtn, a, d, n_progeny, h2 = NULL, var_e = NULL,
             ref = NULL, seed = NULL)          # -> data.frame(id, progeny_mean, n) +
                                               #    the progeny Population (pedigree)
families(pop, by = c("full_sib", "maternal_half_sib", "paternal_half_sib", "selfed"))
```
`progeny_test()` is `combining_ability(method = "simulated", design = "topcross")` with
`mates` drawn at random from a population (or given); **D15** asks whether it is a
separate export or a documented call pattern. Family-structured phenotyping is then
`simulate_phenotype(progeny, ...)` or `phenotype_value(progeny, ...)` plus
`select_ind(..., family = families(progeny, "full_sib"))`; the `family_relationship`
for `"combined"` follows from `by` (0.5 full sibs, 0.25 half sibs — see the selfed
caveat below).

### Genetic theory
*Progeny mean as a predictor.* With mates a random sample of a random-mating population,
the expected half-sib progeny mean is `½ A_i` plus a constant (item 1, `p_T = p`);
dominance enters the parent-dependent part of a half-sib mean (the part that ranks
parents) only through the average effect
`α = a + d(1 − 2p)` (= `a + d(q − p)`, at the mates' frequencies; the breeding value `A_i`
is written with these `α`), while the parent's own dominance deviation does not, and
the full-sib family mean additionally carries the cross's dominance deviation.
Dominance also remains in the constant common to all parents, `a(p − q) + 2dpq` per
locus under HWE, so at `p = ½` the `d(1 − 2p)` term of `α` vanishes but the common
mean does not (a = 0, d = 1, p = ½: every parent's expected mean is 0.5). The accuracy of a progeny test with `n` half-sib progeny, from
`Cov(A_i, ȳ) = ½ σ²_A` and `Var(ȳ) = σ²_P [1 + (n − 1) t]/n` with `t = h²/4`, is
```
r_PT = sqrt( n h² / (4 + (n − 1) h²) )        (package derivation)
```
which rises to 1 as `n → ∞` — the reason progeny testing beats own performance for
low-`h²` traits. Family (among-family) selection response and its comparison with mass
and within-family selection follow Falconer & Mackay (1996) and Lynch & Walsh (1998)
(author-year only, C2); the engine already implements the Lush weights (rubric S3) and
this item adds only the *structure*.

*Selfed families — an existing-documentation check.* `select_ind()` documented
(until the 2026-09-27 fix below) `family_relationship = 0.5` "for full-sibs or selfed
families". For S1 sibs of a
non-inbred S0 plant the tabular rule gives `A_ij = 2Θ = 2·Θ_PP = 1` (with
`Θ_PP = ½(1 + F_P) = ½`), and the sibs are themselves inbred (`F = ½`), so the
correlation of their breeding values is `A_ij / sqrt(A_ii A_jj) = 1/1.5 = 2/3`, not 0.5.
*Resolved 2026-09-27:* `.combined_score()` uses `t = r·h²` and
`Cov(A, family mean) = V_A(1 + (n − 1)r)/n` with `h²` the candidates' own heritability,
so `r` is the correlation of breeding values `A_ij/A_ii` — 2/3 for these S1 sibs. The
`select_ind()` docs (and DECISION-015's note) now say so. The same review found the
index applied size-dependent weights to raw records (fixed: deviations from the mean).
The original flag read: which quantity the Lush derivation in `.combined_score()` needs (`r` as the additive
relationship vs. as the intraclass correlation of breeding values) should be re-derived
in this item's review and the docs/default corrected if needed — flagged here, not
asserted.

*Pitfalls.* Half-sib mean ≈ `½ A` only with random mates (a fixed tester makes it GCA
toward that tester — item 1); unequal family sizes (handled by `.sel_among_family()`);
family means from `simulate_phenotype()` are rescaled per population (use
`phenotype_value()` across cycles, DECISION-021); `families()` must never be inferred
from id strings.

### Validation plan
1. `families()` reproduces the `mating_design()` table exactly (full-sib groups =
   crosses; maternal half-sib groups = dams).
2. **Progeny-test accuracy.** Additive trait, HWE base, `n ∈ {2, 5, 20}` half-sib
   progeny with `phenotype_value()` records at `h2 ∈ {0.2, 0.5}`: realized
   `cor(progeny_mean, bv_sire)` within `3·SE` of `r_PT` over ≥ 30 seeds.
3. **Half the breeding value.** `progeny_test()` with `var_e = 0` and large `n`:
   regression of progeny mean on sire `bv` has slope → 0.5.
4. **Among-family vs mass.** At low `h2` (0.1) and family size 20, among-family
   selection on `phenotype_value()` records gives a higher realized `bv` response than
   mass selection at the same proportion over 30 seeds (Falconer & Mackay's qualitative
   ordering); at `h2 = 0.9` the ordering reverses.
5. **Review must confirm:** `r_PT`; the selfed-family `family_relationship` question
   above; that dominance enters a half-sib mean only in the parent-dependent part through
   `α = a + d(1 − 2p)` (a
   half-sib mean regressed on a dominance-free breeding value shows no separate
   parental-dominance term, and the `d(1 − 2p)` contribution vanishes at `p = ½`; the
   common mean constant `a(p − q) + 2dpq` does not), and
   that a full-sib mean also carries the cross's own dominance deviation, in the
   simulated data.

### Dependencies / risks / effort
F1, F2, item 1. **Effort S** (`families()`) **+ M** (`mating_design()` + `progeny_test()`).

### Decisions needed
- **D15** `progeny_test()` as its own export, or a documented `combining_ability()` call.
- **D16** Resolve the selfed-family `family_relationship` derivation (and fix docs/default
  if the review finds 0.5 wrong).

---

## 6. Item 5 — multi-trait sequential rules: tandem selection and independent culling

### What BD needs
Catalog §E: "Multi-trait rules: tandem · independent culling · selection index (Hazel
1943) — Economic-weight index vs. sequential culling. **Engine:** multi-trait selection".

### What already exists
- **Exists:** `select_ind(method = "index")` (Smith–Hazel `b = P⁻¹Ga` on true BVs),
  `"quadratic_index"`; single-trait selection via `trait =`; schemes take one `trait`.
- **New:** independent culling levels; a tandem *schedule* across generations; a per-
  trait external-prediction input (the index methods ignore `on`, DECISION-019 — a
  culling matrix is the natural place to accept per-trait GEBVs).

### Proposed API
```r
select_ind(sim, method = "culling", culling = c(0.3, 0.5),      # per-trait proportions
           on = "pheno" | "bv" | "gv" | <n x T numeric matrix>,
           direction = c("high", "low") | per-trait vector, ...)
pedigree(..., trait = c(1, 1, 2, 2))        # tandem: trait recycled over generations
recurrent_selection(..., trait = c(1, 2))   # likewise
```
- `culling`: one proportion per trait (or a per-trait threshold list); an individual is
  kept iff it is in the top `culling[t]` fraction of every trait (simultaneous
  independent culling). The number kept is then *emergent*: `n`/`prop`/`intensity`
  must be `NULL` (error otherwise), and the realized count and per-trait differentials
  are returned as attributes. A `sequential = TRUE` variant (cull trait 1, then trait 2
  among survivors) is the multi-stage form; **D17** picks which ships.
- Tandem needs no new operator: only the scheme wrappers' `trait` becomes a vector.

### Genetic theory
The classical comparison of tandem, independent culling and index selection is Hazel &
Lush (1942), extended by Young (1961) (titles verified; their derivations are not
reproduced here). This spec's own derivation under *their* idealized conditions —
`T` uncorrelated traits, equal `σ_A`, equal `h²`, equal economic weights, aggregate
`H = Σ_t A_t`, normality — gives, per generation with intensity `i(p)`:
```
index    R_H = i(p) · sqrt(T) · h σ_A         (b ∝ P⁻¹Ga = h² a, equal weights)
culling  R_H = T · i(p^{1/T}) · h σ_A          (each trait truncated at p^{1/T})
tandem   R_H = i(p) · h σ_A                    (one trait per generation)
```
so `tandem / index = 1/sqrt(T)` and `culling / index = sqrt(T) · i(p^{1/T}) / i(p)`;
for `T = 2`, `p = 0.1`: `i(0.1) = 1.755`, `i(0.316) = 1.125`, ratios 0.907 (culling)
and 0.707 (tandem) — index ≥ culling ≥ tandem, the ordering those papers are known
for, here as an executable expectation rather than a cited number. Pitfalls: the
normal-theory `i(p)` fails for few-QTN discrete distributions; culling on `"bv"` is
idealized (S5); with correlated traits the per-trait proportions interact and the
kept fraction is not the product; `direction` per trait; zero survivors → error, not an
empty `Population`.

### Validation plan
1. **Set semantics.** Kept set == intersection of per-trait top sets (deterministic);
   `T = 1` culling == mass selection at `prop`; error on `n`+`culling`.
2. **Ordering / ratios.** Two uncorrelated additive traits (`architecture =
   "independent"`, equal `h2 = 0.5`, 500 individuals, 30 seeds): realized aggregate-BV
   gains ordered index > culling > tandem, with culling/index and tandem/index within
   tolerance of 0.907 and 0.707 at `p = 0.1`.
3. **Tandem schedule.** A constant `trait` vector reproduces today's output bit-for-bit
   under a fixed seed (regression); alternating traits produce response in each.
4. **Review must confirm:** the three ratios' derivation and its assumptions; the
   emergent-count semantics; per-trait `direction`.

### Dependencies / risks / effort
None. **Effort S** (tandem) **+ S/M** (culling). Independent of the other items — a
good first PR if a quick win is wanted.

### Decisions needed
- **D17** Simultaneous culling (proportions, emergent count) vs sequential multi-stage
  (fixed `n`), or both.
- **D18** Accept an `n × T` matrix `on` for per-trait external predictions (closes the
  "no per-trait GEBV hook" gap of DECISION-019 for culling; the index methods stay on
  true BVs).

---

## 7. Item 6 — multi-population mating for crossbreeding

### What BD needs
Catalog §F: "Crossbreeding: two-way (F1) · three-way · rotational · terminal sire —
Exploit heterosis + breed complementarity. **Engine:** multi-population mating (ties
SPEC-0007)"; consolidated dependency 5 "Multi-population / multi-parent mating —
crossbreeding, rotational & terminal-sire, two-parent-cross-across-pools (RRS-3; ties
SPEC-0007)". BD SPEC-0007 (mate allocation) outputs "a mating list `{mother, father,
n}` executed via `simplePHENOTYPES::cross()`"; SPEC-0011 RRS-3 needs the canvas to draw
a cross across pools (BD SPEC-0012 two-source inputs is done, so the DAG can express it).

### What already exists
- **Exists:** `cross(mother, father)` mates individuals from *different* `Population`s
  sharing a map (BD's review verified this; `.mate()` checks map identity only);
  `c.Population()`; `.intermate()` (random pairs within one population);
  `sample_parents()`; `genotypic_value()` for realized heterosis on a frozen A + D
  architecture; `as_population()` (with its outbred phase caveat, SPEC §4.6).
- **New:** a plan executor `mate()` (F2), design generators, breed-composition
  bookkeeping (needs F1), heterosis reporting, and the crossbreeding systems.

### Proposed API
```r
mate(plan, ..., seed = NULL)
# plan: data.frame(mother, father, n); `...` named Population(s) whose ids the plan
# references (one pool: bare ids; several: "pool:id" or a `pool` column). Executes
# every row with cross() (parent == parent -> selfcross), writes the pedigree, ids
# "<pool>_<k>" deterministic. Replaces .intermate() and BD's .rrs_intermate()/
# .rrs_hybrid_pop() loops.

mating_design(mothers, fathers, design = c("random", "factorial", "nested",
              "diallel", "half_diallel"), n_crosses, progeny_per_cross, ...)  # -> plan

crossbreed(breeds = list(A = popA, B = popB, C = popC),
           system = c("two_way", "three_way", "backcross", "rotational", "terminal"),
           generations = 1L, n_progeny, sire_breed = NULL, seed = NULL)
# -> the final Population (pedigree; breed pools recorded) + attribute `history`
#    (per generation: expected breed composition, realized heterosis if `qtn,a,d` given)

breed_composition(pop)                        # expected breed fractions (exact, pedigree)
heterosis(f1, parents, qtn, a, d)             # realized mid-parent heterosis + expectation
```
**D19** decides whether `crossbreed()` (orchestration) lives here as a scheme wrapper
like `recurrent_selection()`, or in BD composed from `mate()` + `breed_composition()`
+ `heterosis()` (BD's rule: compose engine verbs). The primitives belong here either way.

### Genetic theory
*Breed composition* is exact from the pedigree: an individual's expected fraction from
pool `P` is the mean of its parents' fractions (founders 1/0). Realized per-locus
breed origin is **not** tracked (no founder-origin haplotypes — the same limit
`mabc_select()` documents); a marker-observed version exists only at breed-diagnostic
markers. Say "expected", never "realized", for pedigree fractions.

*Heterosis (package derivation, A + D model, no epistasis).* At one locus with
within-breed frequencies `p_A`, `p_B` of the counted allele and random mating within
breeds (per-locus HWE), the F1 mean minus the mid-parent is
```
H_F1 = d (p_A − p_B)²,        summed over loci;
F2 (F1 × F1, random):  ½ H_F1;   backcross F1 × A:  ½ H_F1;
three-way (A×B) × C:   ½ (H_AC + H_BC);
n-breed rotation, equilibrium:  (2ⁿ − 2)/(2ⁿ − 1) of the pairwise mean  (2/3 for n = 2,
                                6/7 for n = 3).
```
The fractions follow from the probability that an individual's two alleles at a locus
derive from *different* breeds under the breed-origin (dominance) model of heterosis,
conventionally attributed to Dickerson (1973; title verified, text not checked — treat
the attribution of the model as "source to be verified", the fractions as this spec's
derivation): e.g. in a two-breed rotation at equilibrium the sire is purebred and the
dam's composition alternates (2/3, 1/3), so `P(different breed) = 2/3`. Expectations are
linear, so within-breed HWE per locus suffices; linkage disequilibrium changes variances,
not these means. Additive-only architectures have `H = 0` exactly in expectation.
Epistasis is refused for the analytic expectation (no per-locus decomposition) but
realized heterosis via `genotypic_value()` is still reported.

*Pitfalls.* Sex is not modelled — a "sire line" is a designated pool, `mother`/`father`
are symmetric (documented at `cross()`); the outbred-phase guess of `as_population()`
affects multi-generation recombinant progeny (SPEC-0011's phase caveat applies to every
breed pool); id collisions across pools (D2); rotational equilibrium is an asymptote
(state the generation count the test uses); the terminal system needs a maintained dam
source — bookkeeping, not genetics.

### Validation plan
1. **Exactness of `mate()`.** A one-row plan equals `cross()` bit-for-bit under the same
   seed; `recurrent_selection()` refactored onto `mate()` stays bit-identical to the
   pre-change output under fixed seeds (regression); BD's `.rrs_intermate()` can then be
   retired (BD-side).
2. **Composition known answers.** F1 (½, ½); BC1 to A (¾, ¼); three-way (¼, ¼, ½);
   two-breed rotation converges to alternating (2/3, 1/3) within 1e-3 by generation 10.
3. **Heterosis known answers.** Two pools produced by random mating (`mate()` random
   design, so per-locus HWE within pools), frozen `a`, `d` with `d ≠ 0`, 30 seeds:
   realized F1 mid-parent heterosis (`genotypic_value()`) within `3·SE` of
   `Σ d_j (p_Aj − p_Bj)²`; F2 ≈ ½ H_F1; two-breed rotation → 2/3 H_F1 (generations
   6–10 averaged); additive-only → 0.
4. **Review must confirm:** the heterosis expectations and their HWE-within-breed
   condition; the rotation fractions; "expected" vs "realized" wording throughout;
   RNG in R via `cross()` only.

### Dependencies / risks / effort
F1 (composition/pedigree), item 1's scoring conventions. **Effort M** (`mate()` +
`mating_design()`) **+ M** (`crossbreed()` + composition + heterosis).

### Decisions needed
- **D19** Where `crossbreed()` lives (engine wrapper vs BD composition).
- **D20** Namespacing of ids across pools in `mate()` (`"pool_k"` prefix) — with D2.
- **D21** Whether `heterosis()` should also accept a `phenotype_sim` template (ties D5).

---

## 8. Item 7 — polyploid (tetrasomic) model

### What BD needs
Consolidated dependency 7: "`Wheat_div` (allotetraploid TGC 90K) is currently coded
diploid per subgenome-specific SNP (a standard but approximate treatment …). Faithful
polyploid simulation (homoeology, tetrasomic inheritance, dosage 0..2n) needs SP's
planned polyploid framework — `Wheat_div` migrates to it when it lands. Bread wheat
(hexaploid) and other polyploids would use the same." BD `TODO.md`: "Polyploid framework
(SP dependency) + migrate `Wheat_div`". This repo: ROADMAP §4 "Polyploids … Treat as v3."

### What already exists
- **Exists:** nothing polyploid. The Rust `Bits` strand representation generalizes to
  more than two strands (ROADMAP §4 notes this); everything else — `cis`/`trans`,
  −1/0/1 dosage, the variance formulas, `g_matrix()`, the average effect — assumes
  diploidy.

### The request should be re-scoped before any engine work
Durum/emmer (*T. turgidum*, AABB) is an **allo**tetraploid with **disomic** inheritance:
homologues pair within subgenome, homoeologues do not (in the presence of the *Ph1*
system). A subgenome-specific SNP on chromosome 1A therefore segregates as a diploid
locus linked only to other 1A markers — which is exactly what "coded diploid per
subgenome, one map per chromosome 1A…7B" simulates. For `Wheat_div` the diploid model is
**not an approximation of the inheritance**; the only approximation is assay-level
(markers that co-amplify homoeologues would carry composite 0..4 dosage), which is a
data-curation matter, not a meiosis model. Hexaploid bread wheat (AABBDD) is likewise
disomic. Recommendation: BD corrects the `Wheat_div` caveat; no engine change is needed
for allopolyploids with disomic inheritance.

What a polyploid framework *would* add is **autopolyploid, polysomic** inheritance
(potato, alfalfa, blueberry, some forages, sugarcane): 2k homologues, dosage 0..2k,
random pairing among homologues (bivalents, or multivalents with double reduction), and
a quantitative-genetics model in which a locus has more than one dominance term.

### Proposed design outline (v3; for the record, not for Block 3B)
- `Population` gains `ploidy` (even, default 2) and `strands` (a list of `ploidy`
  marker × individual matrices, replacing the `cis`/`trans` pair when `ploidy > 2`);
  `dosages()` returns 0..`ploidy` (or centred); `as_numeric()` learns a dosage input
  (0..2k calls) — deterministic, a Rust candidate under DECISION-006.
- Meiosis (Rust, pure; R draws): random-bivalent model first (draw the pairing of the
  2k homologues into k bivalents uniformly; each bivalent recombines under the same
  count-location process; the gamete takes one chromatid per bivalent) — no double
  reduction. Multivalent pairing with a double-reduction coefficient `α` later. There is
  no isqg reference for polyploids, so the gate is analytic segregation ratios, not
  bit-parity (DECISION-012's parity standard does not apply; a new decision is needed).
- Grammar: additive value linear in dosage; dominance with digenic (and optionally
  trigenic/quadrigenic) terms; average effects and the genotypic-value partition for a
  random-mating autotetraploid follow Kempthorne (1955, title verified); marker
  relationship matrices for additive and digenic-dominance effects from dosages as in
  Endelman et al. (2018, abstract verified: `G`, `D`, and `G#G` matrices from
  tetraploid dosage). Linkage under polysomic inheritance: Fisher (1947), Mather (1936),
  Haldane (1930) (titles verified) — the sources the review would work from.
- Known answers for the meiosis gate (package derivations under random chromosome
  segregation, no double reduction): duplex `AAaa` gametes `AA : Aa : aa = 1 : 4 : 1`,
  selfed → `1 : 8 : 18 : 8 : 1`; simplex `Aaaa` selfed → `AAaa : Aaaa : aaaa = 1 : 2 : 1`;
  heterozygosity decay under selfing slower than diploid (`1 − ⅔·…` per generation for
  tetrasomic vs ½), to be derived and tested.

### Effort / risk
**XL** (representation, Rust core, dosage I/O, every variance formula, `g_matrix()`,
breeding value, docs) — the ROADMAP's "deepest change on the list". Risk of a partial
model shipped as "polyploid support" (X1). **Recommendation: defer to v3 (D22), and
correct the BD-side claim about `Wheat_div` now.**

### Decisions needed
- **D22** Defer autopolyploid support to v3; ask BD to reword the `Wheat_div` caveat
  (disomic inheritance is modelled faithfully by the diploid-per-subgenome coding).
- **D23** Is there concrete autotetraploid demand (a potato/alfalfa dataset and scheme)
  that would justify pulling it forward?

---

## 9. Item 8 — Python package

### What BD needs
Catalog/TODO: "Python package — BD's canvas already generates Python for the planned
`simplephenotypes` API (already `docs/ROADMAP.md` §8b)". BD's codegen audit lists open
Python issues (nested parameters). DECISION-017 recorded that the Python transpilation
targets a mirror API whose header "notes the package is forthcoming".

### What already exists
- **Exists:** the deterministic Rust core (`numeric.rs`, `genome.rs`, `meiosis.rs`) that
  a PyO3 build could share; ROADMAP §8b's plan (PyPI + conda-forge/Bioconda; the
  stochastic grammar must be reimplemented; timing after CRAN); DECISION-011 (single
  rextendr package; the Cargo-workspace/maturin bootstrap was dropped).
- **New:** everything user-facing.

### Assessment
Not a genetics item and not a blocker for any BD *simulation* (BD runs the R engine via
webR/plumber). The live risk is presentational: BD shows generated Python as if runnable
(rubric X1 for BD's manuscript claims). Reimplementing the stochastic grammar in Python
doubles every grammar change and cannot share R's RNG stream, so parity would be
statistical, not bit-exact — a weaker standard than the package holds itself to
elsewhere (DECISION-012), and every one of the eight items above would then need a
second implementation and a second theory review.

### Options
1. **Defer (recommended).** Keep ROADMAP §8b's sequencing (after the CRAN release and
   a frozen grammar). Ask BD to label the Python tab "planned API — not runnable" or
   hide it until a package exists.
2. **Stopgap wrapper.** A thin `simplephenotypes` PyPI package that drives the installed
   R package through `rpy2` — runnable today, bit-identical results, no reimplementation;
   requires R. It does not deliver the shared-Rust-core vision and adds a support
   surface; only worth it if BD needs "Copy Python" to run soon.
3. **Full port** (PyO3/maturin over the Rust core; grammar in Python). **XL**; two
   implementations to keep in theory-review sync.

### Decisions needed
- **D24** Defer (option 1) for Block 3B; decide whether BD keeps emitting Python.
- **D25** If a runnable Python entry point is wanted before the port, accept the `rpy2`
  stopgap (option 2) as an explicitly labelled interim.

---

## 10. Recommended order

| # | Work | Effort | Why here |
|---|------|--------|----------|
| 1 | **F1 pedigree in `Population`** (+ `families()`, `pedigree()`) | M | Unblocks items 3, 4, 6 and tidies 1; purely additive metadata; no draw-order change. |
| 2 | **`mate()` + `mating_design()`** (item 6a, F2) | M | Retires four reimplementations (engine `.intermate()`, `.make_family()`; BD `.rrs_intermate()`, `.rrs_hybrid_pop()`); gives BD SPEC-0007 its executor and RRS-3 its cross-pool primitive; exact regression tests against `cross()`. |
| 3 | **`combining_ability()`** (item 1) | M | Highest BD value (RRS-1b heterosis RRS and hybrid schemes are waiting); reuses `genotypic_value()`, `.avg_effect()`, `cross()`; the expected method is deterministic and small; the simulated method is what BD already does, now exported. |
| 4 | **`progeny_test()` + family phenotyping** (item 4) | S–M | Falls out of 1–3; opens the animal-breeding rows of the catalog (family, sib, progeny testing). |
| 5 | **Tandem + culling** (item 5) | S + S/M | Self-contained; can be slotted anywhere (first, if a quick win is wanted). |
| 6 | **`marker_select()`** (item 2) | S–M | Refactor of `select_mabc.R` helpers; MARS index already covered by `additive_value()`; estimated weights wait for 7. |
| 7 | **`predict_ebv()` + `a_matrix()` + `selection_methods()`** (item 3) | M + M | Largest scope decision (D12); GBLUP can start any time, pedigree BLUP after 1; replaces BD's constant accuracy with a realized one. |
| 8 | **`crossbreed()` + `breed_composition()` + `heterosis()`** (item 6b) | M | After 1–3; mostly bookkeeping over `mate()` and `genotypic_value()`. |

Ordering principle: unblock the most BD rows per unit of new genetics, reuse the most
existing tested code, and keep the theory surface per PR small enough for the
independent review (each item is one `dev/dual.sh` loop, one DECISION entry, one
contract addition). Version: each landing is an additive export; per the maintainer
(2026-09-28) the release tag stays `2.0` and only the development version moves
(`2.0.0.9001`); add every new export to `tests/testthat/test-backend-contract.R`.

## 11. Deferrals and declines

- **Item 7 (polyploid) — defer to v3 (D22).** XL, touches the representation, the Rust
  core and every variance formula; and BD's stated need (`Wheat_div`) is not a need —
  disomic allotetraploids are modelled faithfully by the diploid-per-subgenome coding.
  Re-open when an autopolyploid dataset/scheme is on BD's catalog.
- **Item 8 (Python) — decline for Block 3B (D24).** Not genetics; sequenced after CRAN by
  ROADMAP §8b; a port would double the review burden of every item above. Near-term
  action is BD-side labelling; `rpy2` stopgap only on request.
- **Within items:** multi-trait BLUP and single-step `H` (D14), double reduction, and a
  realized (haplotype-origin) breed composition are out of scope; each is named in its
  item so it is a deliberate omission, not an overclaim.

## 12. Consolidated maintainer decisions

| ID | Item | Question | Recommendation |
|----|------|----------|----------------|
| D1 | F1 | Add an optional `pedigree` element to `Population`, written by every mating function and preserved by `[`/`c()`? | Yes |
| D2 | F1/6 | Id-collision policy when pooling/mating across populations with pedigrees: error, or namespace ids by pool in `mate()`? | Namespace in `mate()`; `c()` errors on collisions |
| D3 | 1 | `combining_ability()` takes the frozen `(qtn, a, d)` (as `genotypic_value()`), not a `phenotype_sim`? | Yes |
| D4 | 1 | Ship `"expected"` and `"simulated"` together? | Together (the Monte-Carlo agreement is the test) |
| D5 | 1 | Export `template_effects(sim, trait)` (realized-scale `a`, `d` per QTN) so BD stops re-deriving? | Yes (S) |
| D6 | 1/3/4 | Extend `phenotype_value()` with `d =` (A + D fixed-scale phenotype) as the one residual model? | Yes |
| D7 | 2 | Include staged pyramiding (`min_markers`) now? | Yes (trivial once the filter exists) |
| D8 | 2 | `marker_index()` alias over `additive_value()`? | Documentation only, unless BD's DAG needs a node name |
| D9 | 2/3 | MARS estimated weights from the engine (`predict_ebv()`) or BD (SPEC-0006)? | Engine, if D12 = yes |
| D12 | 3 | Known-variance-component BLUP (`predict_ebv()`) in the engine? | Yes — deterministic, dependency-free, an estimate for the `on` hook; no REML/Bayes |
| D13 | 3 | Ship `selection_methods()` manifest for BD SPEC-0006? | Yes (S) |
| D14 | 3 | Multi-trait / single-step BLUP? | Defer |
| D15 | 4 | `progeny_test()` as its own export? | Yes (animal-breeding vocabulary; thin over `combining_ability()`) |
| D16 | 4 | Re-derive `family_relationship` for selfed families (doc says 0.5; tabular rule gives `A = 1`, BV correlation 2/3)? | Re-derive in the item's review; fix docs/default if confirmed |
| D17 | 5 | Culling semantics: simultaneous proportions (emergent count) vs sequential with `n`? | Simultaneous first; sequential as an option |
| D18 | 5 | Accept an `n × T` matrix `on` for per-trait external predictions in culling? | Yes |
| D19 | 6 | `crossbreed()` in the engine (scheme wrapper) or composed in BD? | Primitives here; wrapper here too, for parity with `recurrent_selection()` |
| D20 | 6 | Id prefix scheme for `mate()` output | `"<pool>_<k>"` |
| D21 | 6 | `heterosis()` also from a `phenotype_sim` template? | Via D5's `template_effects()` |
| D22 | 7 | Defer autopolyploid support to v3; ask BD to reword the `Wheat_div` caveat? | Yes |
| D23 | 7 | Any concrete autotetraploid demand to pull it forward? | Maintainer to say |
| D24 | 8 | Defer Python; BD labels generated Python as planned/not runnable? | Yes |
| D25 | 8 | Accept an `rpy2` stopgap if a runnable Python entry point is needed sooner? | Only on request |

(D10–D11 intentionally unused.) Each accepted item gets its own DECISION entry
(DECISION-024 onward), a `docs/SPEC-<item>.md` in the format of
`SPEC-nonadditive-correlation.md` with executed acceptance criteria, and an independent
theory review before commit (`dev/dual.sh`), per `AGENTS.md`.

## 13. Sources verified this session

All entries below were checked on 2026-09-27 against the Crossref API (title, authors,
journal, volume, pages, DOI); where noted, the abstract was also read via Europe PMC.
They are cited above for what their titles/abstracts establish; no equation is
attributed to them unless stated as such.

- Sprague GF, Tatum LA (1942) General vs. specific combining ability in single crosses
  of corn. *Agronomy Journal* 34:923–932. doi:10.2134/agronj1942.00021962003400100008x
  — GCA/SCA concept (Crossref).
- Griffing B (1956) Concept of general and specific combining ability in relation to
  diallel crossing systems. *Australian Journal of Biological Sciences* 9:463–493.
  doi:10.1071/BI9560463 — diallel formulation (Crossref).
- Comstock RE, Robinson HF, Harvey PH (1949) A breeding procedure designed to make
  maximum use of both general and specific combining ability. *Agronomy Journal*
  41:360–367. doi:10.2134/agronj1949.00021962004100080006x — RRS origin (Crossref).
- Hallauer AR, Eberhart SA (1970) Reciprocal full-sib selection. *Crop Science*
  10:315–316. doi:10.2135/cropsci1970.0011183x001000030033x — catalog row (Crossref).
- Melchinger AE, Frisch M (2023) Genomic prediction in hybrid breeding: II. Reciprocal
  recurrent genomic selection with full-sib and half-sib families. *Theoretical and
  Applied Genetics* 136:203. doi:10.1007/s00122-023-04446-3 — as read by BD's `rrs.R`
  (Crossref; abstract via Europe PMC).
- Hazel LN, Lush JL (1942) The efficiency of three methods of selection. *Journal of
  Heredity* 33:393–399. doi:10.1093/oxfordjournals.jhered.a105102 — tandem / culling /
  index comparison (Crossref).
- Young SSY (1961) A further examination of the relative efficiency of three methods of
  selection for genetic gains under less-restricted conditions. *Genetical Research*
  2:106–121. doi:10.1017/S0016672300000598 (Crossref).
- Lande R, Thompson R (1990) Efficiency of marker-assisted selection in the improvement
  of quantitative traits. *Genetics* 124:743–756. doi:10.1093/genetics/124.3.743
  (Crossref; abstract via Europe PMC).
- Henderson CR (1975) Best linear unbiased estimation and prediction under a selection
  model. *Biometrics* 31:423–447. doi:10.2307/2529430 (Crossref; abstract via Europe PMC).
- Henderson CR (1976) A simple method for computing the inverse of a numerator
  relationship matrix used in prediction of breeding values. *Biometrics* 32:69–83.
  doi:10.2307/2529339 (Crossref).
- Wright S (1922) Coefficients of inbreeding and relationship. *The American
  Naturalist* 56:330–338. doi:10.1086/279872 (Crossref).
- Emik LO, Terrill CE (1949) Systematic procedures for calculating inbreeding
  coefficients. *Journal of Heredity* 40:51–55. doi:10.1093/oxfordjournals.jhered.a105986
  (Crossref).
- Meuwissen THE, Hayes BJ, Goddard ME (2001) Prediction of total genetic value using
  genome-wide dense marker maps. *Genetics* 157:1819–1829.
  doi:10.1093/genetics/157.4.1819 (Crossref; abstract via Europe PMC).
- Daetwyler HD, Villanueva B, Woolliams JA (2008) Accuracy of predicting the genetic
  risk of disease using a genome-wide approach. *PLoS ONE* 3:e3395.
  doi:10.1371/journal.pone.0003395 (Crossref; abstract via Europe PMC).
- Legarra A, Aguilar I, Misztal I (2009) A relationship matrix including full pedigree
  and genomic information. *Journal of Dairy Science* 92:4656–4663.
  doi:10.3168/jds.2009-2061 (Crossref; abstract via Europe PMC).
- Dickerson GE (1973) Inbreeding and heterosis in animals. *Journal of Animal Science*
  1973(Symposium):54–77. doi:10.1093/ansci/1973.symposium.54 — existence verified
  (Crossref); the breed-origin heterosis model's attribution to it is **to be verified
  against the text**.
- Kempthorne O (1955) The correlation between relatives in a simple autotetraploid
  population. *Genetics* 40:168–174. doi:10.1093/genetics/40.2.168 (Crossref).
- Endelman JB, Schmitz Carley CA, Bethke PC, et al. (2018) Genetic variance
  partitioning and genome-wide prediction with allele dosage information in
  autotetraploid potato. *Genetics* 209:77–87. doi:10.1534/genetics.118.300685
  (Crossref; abstract via Europe PMC).
- Fisher RA (1947) The theory of linkage in polysomic inheritance. *Philosophical
  Transactions of the Royal Society of London B* 233:55–87. doi:10.1098/rstb.1947.0006
  (Crossref).
- Mather K (1936) Segregation and linkage in autotetraploids. *Journal of Genetics*
  32:287–314. doi:10.1007/BF02982683 (Crossref).
- Haldane JBS (1930) Theoretical genetics of autopolyploids. *Journal of Genetics*
  22:359–372. doi:10.1007/BF02984197 (Crossref).
- Already cited by the package and re-used here without re-verification: VanRaden 2008
  (`g_matrix()`), Meuwissen 1997 (`optimum_contribution()`), Frisch & Melchinger 2001/2005
  (`mabc_select()`), Fisher 1918 / Falconer & Mackay 1996 / Lynch & Walsh 1998 (author-
  year only, rubric C2), Bernardo 2020 (author-year only), Toledo et al. 2019 (isqg).
