# SPEC-nonadditive-correlation.md — genetic-correlation control for non-additive layers

> **DRAFT / scoping — not yet implemented.** Companion to `SPEC.md` (§13 PleioArch)
> and `SPEC-transcriptome.md`. Records the design for extending `cor` / pleiotropy /
> LD architecture control from the **additive** layer to **dominance** and
> **epistasis** layers, per the 2026-09-17 audit (grammar P1/P4, effects-arch O2/X1)
> and the user's decision to *extend* control (vs. document-only or reject).
> Needs **DECISION-023** recorded once the approach here is accepted.

## 1. Problem

`cor` (and the `pleiotropy` / `ld` architectures) controls only the **additive**
genetic correlation. Dominance and epistasis layers ignore it: under
`architecture = "pleiotropy"` they reuse the **same loci** across traits
(`arch_independent.R`: `rep(list(shared), nt)`) and are given the **same
deterministic effect series** (`effects_series.R`: `base ^ seq_len(n)`, identical
per trait), so those components realize a genetic correlation of ≈ **+1**
irrespective of the requested `cor`.

**Audit evidence (executed):**
- Pleiotropic epistasis, target `cor = 0`, 5 000 ind × 100 pairs → realized total
  genetic correlation **1.0**.
- One-call `model = "AD"`, target `cor = 0` → **0.45**.
- `cor = 0.2`, additive `prop = 0.25` + epistasis `prop = 0.25`, 100 QTN/pairs, 30
  seeds → additive correlation **0.196** (correct), **total** genetic correlation
  **0.60** (range 0.44–0.71), not 0.2.

The additive machinery is correct (audit P1/P2/P3 PASS); the gap is that the
non-additive components are not brought under the same control, and the public
docs advertise "controlled genetic correlation" for models (`AD`, `AE`) and
architectures that include them (grammar X1, effects-arch X1).

## 2. Goals and non-goals

### Goals
- The **realized total genetic correlation** between traits equals the requested
  `cor` (in expectation, up to finite-QTN sampling) for models that combine
  additive with **dominance** and/or **epistasis** pleiotropic layers.
- One consistent mechanism: the same PleioArch covariance construction
  (`Σ_ij = cor_ij · √(V_i V_j)`, PSD feasibility, MAF scaling) applied per genetic
  **component**, not just to the additive effects.
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

If **each component** realizes the target correlation — `Cov(x_{c,1}, x_{c,2}) =
cor · √(Var(x_{c,1}) Var(x_{c,2}))` — and per-trait component variances are equal
across traits (they are, by the `prop` budget), then `Cov(g_1,g_2) = cor · Σ_c
Var(x_{c,t}) = cor · Var(g_t)`, i.e. **the total genetic correlation equals `cor`.**
So it suffices to make every mean-effect component realize `cor` on its own; no
joint optimization across components is required. This is the design's core claim
and the first thing the acceptance tests must confirm empirically.

## 4. Current mechanism (additive), to be generalized

- `effects_pleioarch.R::.pleio_draw()` partitions QTN into shared (pleiotropic)
  and trait-specific classes from `pi`, builds `Σ` from `cor` and per-trait
  additive `prop`, checks PSD feasibility (`cor² ≤ π_i π_j`, eigenvalue), draws
  correlated **additive** effects on the shared loci (eigen square root of `Σ`),
  and applies MAF scaling `1/√(2·MAF·(1−MAF))`.
- It is called **only** from `additive()` (`grammar_layers.R:206`).
- `dominance()` / `epistasis()` draw loci via `arch_independent.R` (shared under
  `pleiotropy`) but effects via `effects_series.R` (identical series per trait) →
  correlation ≈ 1.

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
- Factor MAF/heterozygosity scaling into the effect covariance so realized (not
  nominal) dominance variance matches the budget, mirroring how additive uses MAF
  scaling.

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

### 5.4 Feasibility and signs
- Reuse the PSD feasibility check per component; a target `cor` infeasible for a
  component's `pi`/`prop` errors with the component named.
- Preserve the audit guards: zero-variance-trait `cor` (warn/undefined),
  single-shared-QTN `cor` (±1 warning) apply per component.

## 6. Acceptance criteria (executed-R, added to the eval golden set)

1. **Total correlation hits target.** For `model ∈ {AD, AE, ADE}` and
   `architecture = "pleiotropy"`, over ≥30 seeds with adequate QTN counts, the mean
   realized **total** genetic correlation equals `cor` within a stated tolerance,
   for `cor ∈ {−0.5, 0, 0.5}`. (The current code fails this: `cor=0` → ~1.0.)
2. **Per-component correlation** each ≈ `cor` (diagnostic that §3 holds).
3. **Additive-only backward compatibility:** identical realized values under a
   fixed seed vs. pre-change (golden snapshot).
4. **Feasibility errors** name the offending component; **audit guards**
   (zero-variance, single-QTN) still fire per component.
5. **Docs:** grammar X1 / effects-arch X1 overclaim caveats removed; `cor` docs
   describe total-genetic-correlation control across A/D/E.
6. Full `devtools::test()` + `cargo test` green; independent review (Codex, and the
   PleioArch reference implementation) agrees the covariance construction is correct.

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
`evals/mutations.json` + tests for §6.

## 9. Decision to record

**DECISION-023** (pending): genetic-correlation control extends to dominance and
epistasis via per-component PleioArch covariance draws; the realized **total**
genetic correlation is the controlled quantity; epistasis-under-LD is faithful or
an explicit restriction. Record in `docs/DECISIONS.md` once the approach is
accepted and the acceptance tests (§6) pass.
