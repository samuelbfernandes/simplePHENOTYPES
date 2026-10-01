# SPEC-nonadditive-correlation.md — genetic-correlation control for non-additive layers

> **IMPLEMENTED (2026-09-25) — recorded as DECISION-023.** Companion to `SPEC.md`
> (§13 PleioArch) and `SPEC-transcriptome.md`. Extends `cor` / pleiotropy / LD
> architecture control from the **additive** layer to **dominance** and
> **epistasis** layers, per the 2026-09-17 audit (grammar P1/P4, effects-arch
> O2/X1) and the maintainer's decision to *extend* control (vs. document-only or
> reject). Code: `.pleio_nonadditive_draw()` / `.pleio_units()` /
> `.pleio_unit_effects()` (`R/effects_pleioarch.R`), `.epi_unit_column()`
> (`R/grammar_realize.R`); tests: `tests/testthat/test-nonadditive-cor.R`.
> **§5.3 outcome: the restriction fallback was taken** — under "ld", dominance must
> reuse the additive layer's linked loci (`same_as_add = TRUE`) and `epistasis()`
> errors (see §5.3). §6 acceptance criteria met after the independent review's
> findings were fixed — round 1: total-vs-component claim, constant units in the
> variance allocation, single-shared-unit warning text, component named in
> feasibility errors; round 2: "in expectation" wording for a random-ratio
> correlation, "ld" epistatic markers shared across traits, the attenuation
> multiplier; round 3: the total-correlation warning's "expected" label, "ld"
> epistasis outside the r² window / reusing other layers' loci, this file's stale
> sections; round 4: raw-draw vs rescaled covariance, the LE condition, one live
> shared unit after dropping hetless units; round 5: `cor = 0` in the single-unit
> guards, remaining "equals"/"exact" wording and the complete-LD case (±1 per
> replicate, not a limit), the sign of the constant-specific-unit warning; round
> 6: dated corrections to DECISION-007/010/013's "in expectation" claims and the
> major-QTN exception to convergence in the public `cor` docs; round 7: derived
> transcriptome signal outside `cor` (warned), a second additive layer under
> "ld" rejected, the dead-shared-unit warning's "near 0", stale README/vignette
> wording; round 8: single-draw direction wording in two warnings, the one-unit
> and `cor = 0` conditions, the direction of the `pi` advice. See DECISION-023 and
> the §3 correction. Round 9: the derived-transcriptome warning extended to "ld". Round 10: the
> multi-trait feasibility check made variance-free; no transcriptome warning at
> `prop = 0`. Round 11: one-trait transcriptome under "ld" warns; no inflation
> warning at `cor = 0`. Round 12: single-unit ±1 vs noisy classified per pair;
> relative (1%) threshold for the non-proportional-total warning. Round 13:
> zero-variance pairs excluded from the single-unit warning. Round 14: the v2
> vignette's convergence sentence, §6's stale "current code" line and its
> golden-set claim (the acceptance run is now `dev/accept-nonadditive-cor.R`).
> Round 15: dead-shared-unit warning at `cor = 0`; §2's "MAF scaling" for every
> component. Round 16: effective per-component targets in the total check,
> exact-±1 reclassification after constant specifics, no transcriptome warning at
> genetic h² = 0. Round 17: the transcriptome gene check covers every
> `vary_qtn` replication. Round 18: gene check limited to contributing traits;
> per-replication effective targets; the acceptance script gates on components.
> Round 19: hetless checks over every replication; transcriptome warning on the
> realized genetic signal. Round 20: convergence needs the individuals to grow
> too (a sample correlation over n); §6 item 6 no longer names the review tool.
> Round 21: TODO wording, covariance-vs-correlation in the complete reference,
> small-value formatting in the total warning.
>
> Audit 2026-09-29 corrections: the raw component's expected variance over effect
> draws is `E[Var(c_t)] = V_t` -- the shared units contribute `Σ_tt = π_t V_t` and
> the independent trait-specific units the remaining `(1 − π_t) V_t`; earlier
> wording (here, in DECISION-023 and in the `.pleio_nonadditive_draw()` comment)
> wrote `E[Var(c_t)] = Σ_tt`, which holds only for `π_t = 1`. Only the covariance
> `E[Cov(c_1, c_2)] = Σ_12` comes exclusively from the shared units, so the
> correlation target is unchanged. `epistasis(qtn = list(...))` is now one element
> per trait (a list of sets is refused); a trait with `π_t = 1` (or `0`) is
> assigned trait-specific (shared) loci of exactly zero effect, which are no
> longer reported as its QTNs.
>
> §1 and §4 below describe the code **as it was before** DECISION-023.

## 1. Problem (before DECISION-023)

`cor` (and the `pleiotropy` / `ld` architectures) controlled only the **additive**
genetic correlation. Dominance and epistasis layers ignored it: under
`architecture = "pleiotropy"` they reused the **same loci** across traits
(`arch_independent.R`: `rep(list(shared), nt)`) and were given the **same
deterministic effect series** (`effects_series.R`: `base ^ seq_len(n)`, identical
per trait), so those components realized a genetic correlation of ≈ **+1**
irrespective of the requested `cor`.

**Audit evidence (executed):**
- Pleiotropic epistasis, target `cor = 0`, 5 000 ind × 100 pairs → realized total
  genetic correlation **1.0**.
- One-call `model = "AD"`, target `cor = 0` → **0.45**.
- `cor = 0.2`, additive `prop = 0.25` + epistasis `prop = 0.25`, 100 QTN/pairs, 30
  seeds → additive correlation **0.196** (correct), **total** genetic correlation
  **0.60** (range 0.44–0.71), not 0.2.

The additive machinery was correct (audit P1/P2/P3 PASS); the gap was that the
non-additive components were not brought under the same control, while the public
docs advertised "controlled genetic correlation" for models (`AD`, `AE`) and
architectures that included them (grammar X1, effects-arch X1).

## 2. Goals and non-goals

### Goals
- Every mean-effect **component** (additive, dominance, epistasis) targets the
  requested `cor` — `cor` sets the effect draw's covariance; after rescaling to
  `prop` the realized correlation converges to `cor` as the units and the
  individuals grow, given
  approximate linkage equilibrium among the causal loci (a random ratio,
  attenuated toward 0 on average with few units, and possibly not converging at
  all under strong LD, as for the additive layer) —
  so the **total genetic correlation** targets `cor`
  whenever the layers' per-trait `prop` profiles are proportional — always so for
  scalar `prop` and the one-call models — and is otherwise reported (warned) at
  its attainable, attenuated value (§3).
- One consistent mechanism: the same PleioArch covariance construction
  (`Σ_ij = cor_ij · √(V_i V_j)`, PSD feasibility, per-unit design scaling) applied
  per genetic **component**, not just to the additive effects. The scaling is
  MAF-based `1/√(2·MAF·(1−MAF))` for additive and the realized design-column sd
  for dominance/epistasis (DECISION-023: non-additive design variance is not a
  clean function of MAF).
- Honest, self-consistent public documentation once implemented; retire the
  additive-only caveats added for the audit.
- Backward compatible for additive-only models (identical realized values under a
  fixed seed).

### Non-goals (this version)
- Cross-**component** covariance targeting (e.g. deliberately correlating trait 1's
  additive value with trait 2's dominance value). Components are drawn
  independently; only within-component cross-trait covariance is controlled.
- Correlation control for **vQTL** (a variance category, not a mean genetic
  component) or **transcriptome** layers (own SPEC).
- Extending the `ld` architecture's distinct-but-linked-loci construction to
  epistasis pairs beyond what §5.3 scopes; if it proves out-of-budget it drops to
  a documented restriction (two-trait additive/dominance LD only).

## 3. Why per-component control yields the target total

Let the total genetic value of trait `t` be `g_t = Σ_c x_{c,t}` over components
`c ∈ {A, D, E}`, each centered. If components are drawn **independently** of one
another, cross-component covariances vanish in expectation, so

```
Cov(g_1, g_2) = Σ_c Cov(x_{c,1}, x_{c,2})
Var(g_t)      = Σ_c Var(x_{c,t})
```

If **each component** has expected cross-trait covariance `cor · √(V_{c,1} V_{c,2})`
with `V_{c,t} = Var(x_{c,t})`, the total's target — expected covariance over the
root of the expected variances, which the realized correlation converges to as the
units and the individuals grow when the causal loci are in approximate linkage
equilibrium — is

```
cor · Σ_c √(V_{c,1} V_{c,2}) / √(Σ_c V_{c,1} · Σ_c V_{c,2})
```

By Cauchy–Schwarz the ratio is ≤ 1, with equality **iff the per-trait variance
profiles are proportional across components** (`V_{c,1} / V_{c,2}` the same for
every c). That holds whenever each layer's `prop` is a scalar (the case in every
example and the one-call models): then **the total genetic correlation's target
equals `cor`** (as it trivially does for `cor = 0`), and per-component control suffices — no joint optimization is
needed.

*Correction (independent review, 2026-09-25, O1):* an earlier draft of this
section asserted that component variances are always equal across traits "by the
`prop` budget". They are not when `prop` is a per-trait vector. E.g. additive
`prop = c(0.49, 0.01)` with dominance `prop = c(0.01, 0.49)` gives a total of
`cor · 2·√(0.49·0.01) / 0.5 = 0.28 · cor` (0.14 at `cor = 0.5`) — and since each component's correlation
is bounded by 1, a total of 0.5 is **not attainable** there at all with independent
components. Layers are also added one at a time (a pipe), so re-targeting earlier
layers to compensate is neither well defined nor order-independent. The design
therefore controls `cor` **per component**, and `.pleio_total_cor_check()` warns
with the total's large-sample target whenever the profiles are not proportional.

## 4. The additive mechanism that was generalized (state before DECISION-023)

- `effects_pleioarch.R::.pleio_draw()` partitions QTN into shared (pleiotropic)
  and trait-specific classes from `pi`, builds `Σ` from `cor` and per-trait
  additive `prop`, checks PSD feasibility (`cor² ≤ π_i π_j`, eigenvalue), draws
  correlated **additive** effects on the shared loci (eigen square root of `Σ`),
  and applies MAF scaling `1/√(2·MAF·(1−MAF))`. (Unchanged.)
- It was called **only** from `additive()`.
- `dominance()` / `epistasis()` drew loci via `arch_independent.R` (shared under
  `pleiotropy`) but effects via `effects_series.R` (identical series per trait) →
  correlation ≈ 1. They now use `.pleio_nonadditive_draw()` under "pleiotropy".

## 5. Design

### 5.1 Dominance
Draw **correlated dominance deviations** across traits on the shared het-bearing
loci, using the same covariance construction as `.pleio_draw` but on the dominance
component's variance budget (`prop` of the dominance layer). The dominance design
column is the centered heterozygote indicator; its per-locus variance depends on
heterozygosity, so:
- Shared loci must be het-bearing in **both** traits (they are the same
  individuals, so a locus is het-bearing or not regardless of trait — the existing
  hetless guards from the audit already apply).
- Factor heterozygosity scaling into the effect draw so realized (not nominal)
  dominance variance matches the budget, mirroring how additive uses MAF scaling.
  *(Implemented as division by the realized sd of each unit's design column, not a
  MAF formula — see DECISION-023.)*

### 5.2 Epistasis
Draw **correlated interaction effects** across traits on the shared interacting
sets. Each pair contributes one effect; build the pair-effect cross-trait
covariance from `cor` and the epistasis `prop`. The centered interaction design
already exists (`grammar_realize.R`); only the effect **draw** changes from
identical-series to correlated. `interaction_type = "d"` positions inherit the
dominance heterozygosity caveat (existing partial-hetless warning applies).

### 5.3 `architecture = "ld"`
Additive LD uses distinct-but-linked causal loci per trait (audit P4 PASS). For
non-additive LD, the minimal target is: dominance under LD analogous to additive
(distinct linked het-bearing loci). Epistasis-under-LD is the highest-risk piece;
if it cannot be made faithful within budget, **restrict** `ld` to additive +
dominance and error clearly for epistasis under `ld` (a documented restriction,
not silent misbehavior).

**Outcome (implemented):** the restriction was taken. Independent review showed
that ld epistasis — whether documented as trait-specific or drawn disjoint within
the layer — breaks DECISION-014's contract (causal markers outside the r² window,
and markers another layer made causal for the other trait). `epistasis()` under
"ld" errors; `dominance()` under "ld" must reuse the additive linked loci
(`same_as_add = TRUE`), because a fresh dominance draw can collide across layers.

### 5.4 Feasibility and signs
- Reuse the PSD feasibility check per component; a target `cor` infeasible for a
  component's `pi`/`prop` errors with the component named.
- Preserve the audit guards: zero-variance-trait `cor` (warn/undefined),
  single-shared-QTN `cor` (±1 warning) apply per component.

## 6. Acceptance criteria (executed-R)

Criteria 1–2 are checked by the ≥30-seed acceptance run
`dev/accept-nonadditive-cor.R` (18/18 cells passed; results in DECISION-023). The
committed `tests/testthat/test-nonadditive-cor.R` is a lighter 4-seed regression of
the same behaviour plus the guard and restriction tests. (The `evals/` golden set
seeds theory bugs to test the *reviewer*; no DECISION-023 entry was added there.)

1. **Total correlation hits target.** For `model ∈ {AD, AE, ADE}` and
   `architecture = "pleiotropy"`, over ≥30 seeds with adequate QTN counts, the mean
   realized **total** genetic correlation is within a stated tolerance of `cor`,
   for `cor ∈ {−0.5, 0, 0.5}`. (The code before DECISION-023 failed this:
   `cor = 0` → ~1.0.)
2. **Per-component correlation** each ≈ `cor` (diagnostic that §3 holds).
3. **Additive-only backward compatibility:** identical realized values under a
   fixed seed vs. pre-change (golden snapshot).
4. **Feasibility errors** name the offending component; **audit guards**
   (zero-variance, single-QTN) still fire per component.
5. **Docs:** grammar X1 / effects-arch X1 overclaim caveats removed; `cor` docs
   describe total-genetic-correlation control across A/D/E.
6. Full `devtools::test()` + `cargo test` green; the covariance construction is
   checked against the PleioArch reference implementation and passes the
   independent theory review (`docs/THEORY_REVIEW.md`) before commit.

## 7. Risks / open questions

- **PleioArch source.** The additive construction traces to Prado et al. (in
  preparation); no public page. Extending it to dominance/epistasis is our design,
  not a published method — label it as such (Rule #7), and keep the additive part's
  existing attribution.
- **Realized vs nominal variance.** Dominance/epistasis component variance is not a
  clean function of allele frequency the way additive is; the effect covariance may
  need calibration against the *realized* component sd (as the layers already scale
  to `prop`) rather than a closed form. Validate empirically, not just algebraically.
- **Cross-component covariance** is assumed ~0 (independent draws); confirm it does
  not creep in under LD or away from HWE, which would bias the total (this is the
  same non-orthogonality caveat as the variance budget, SPEC §2 / audit V2).
- **Epistasis-under-LD** may be out of budget → documented restriction fallback.
- Interaction with `complex_phenotypes()` (which combines realized models): total
  correlation there is a separate, already-documented emergent quantity — out of
  scope, but add a note so it is not mistaken for a regression.

## 8. Files likely touched

`R/effects_pleioarch.R` (generalize `.pleio_draw` or factor a
`.pleio_effects(cov, ...)` core), `R/effects_series.R` /
`R/arch_independent.R` (correlated non-additive effect draws under `pleiotropy`),
`R/grammar_layers.R` (`dominance()` / `epistasis()` call the correlated draw),
`R/arch_ld.R` (non-additive LD or the restriction), `R/grammar_realize.R` (no
change expected; the design columns already exist), docs (`SPEC.md` §13,
`grammar_simulate_phenotype.R` roxygen, `DECISIONS.md` DECISION-023), and
tests for §6 (the acceptance run landed in `dev/accept-nonadditive-cor.R`, not in
`evals/mutations.json`).

## 9. Decision to record

**DECISION-023 — recorded** in `docs/DECISIONS.md` (2026-09-25): genetic-correlation
control extends to dominance and epistasis via per-component PleioArch covariance
draws; each component targets `cor`, and the total targets `cor` for proportional
per-trait `prop` (otherwise a warning gives its large-sample target); under "ld"
the restriction of §5.3 applies.
