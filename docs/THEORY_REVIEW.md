# THEORY_REVIEW.md — genetic-theory review rubric

> The independent reviewer model applies this rubric to a change **before it is
> committed**. Purpose: catch errors in the **quantitative-genetics theory** (and the
> code that implements it), not just failing tests. The implementer's model must not run
> its own review — hand the diff to the *other* model (`dev/dual.sh review`).
>
> **Do not invent.** Every equation must trace to a primary source with a page you can
> verify. If you cannot verify a page or a claim, mark it `UNVERIFIABLE` — never guess.

## How to review

1. Read the change (a diff, a file list, or `--staged`). Read the functions it touches
   and the `@references` in their roxygen.
2. For each applicable checklist item below, decide **PASS / FAIL / UNVERIFIABLE** and
   give: the `file:line`, the equation *as implemented*, the primary source it should
   match, and — for FAIL — a **concrete falsifying case** (inputs → wrong output).
3. Try to break it. A confirming-only review is a failed review. Prefer numerical
   counterexamples over prose.
4. End with a one-line verdict: `THEORY: PASS` only if no FAIL remains; otherwise
   `THEORY: FAIL (n items)`.

## Output format (required)

```
### <function / file>
- [PASS|FAIL|UNVERIFIABLE] <rubric id> — <claim>
  evidence: <file:line> — <implemented form>
  source:   <Author Year, Journal vol:pages> (verified? yes/no)
  break:    <falsifying input → expected vs actual>   # FAIL only
...
THEORY: PASS | FAIL (n)
```

---

## Checklist

### C. Citations & equations (Rule #7 — the hard one)
- **C1** Every cited source in roxygen `@references` exists and the **page range is
  correct** for the equation used. Known-good, already verified in this package:
  Smith 1936 *Ann. Eugen.* 7:240–250; Hazel 1943 *Genetics* 28:476–490;
  Lush 1947 *Am. Nat.* 81:241–261, 362–379. Others the code cites (VanRaden 2008
  *J. Dairy Sci.* 91:4414–4423; Meuwissen 1997 *J. Anim. Sci.*; Zhong & Jannink 2007
  *Genetics* 176:2453–2461; Lehermeier 2017 *Genetics*) — **verify the page against the
  paper**, do not assume.
- **C2** No fabricated book chapter/page (e.g. "Bernardo ch. 7"): author-year only unless
  a page is verified.
- **C3** The equation in the code matches the equation in the cited source (symbols,
  scaling, sign), not merely the same name.

### V. Variance partitioning & heritability (SPEC §2)
- **V1** Σ mean-effect genetic proportions (add+dom+epi) = requested H²; Σ all layer
  props (incl. vqtl) ≤ 1; homoskedastic residual = 1 − Σ all props. `prop` sum > 1 errors.
- **V2** Reported h² is computed from **realized** genetic/phenotypic values, not asserted
  from the budget (components not orthogonal at non-0.5 freq). vQTL share is *residual*,
  excluded from H².
- **V3** Coding convention is the simulation one (−1/0/1 additive dosage; het-deviation
  indicator for dominance; centered a×a / a×d / d×d for epistasis), **not** Fisher's
  orthogonal average-effects decomposition — and the docs say so where it matters.
- **V4** `d`-type epistasis / dominance degenerate on hetless (inbred) loci is handled
  (errors or warns), not silently NaN.

### S. Selection engine (`select.R`, DECISION-015)
- **S1** Truncation response tracks **R = i·h²·σ_P**; realized intensity
  **i(p) = φ(Φ⁻¹(1−p))/p**; both `direction`s handled; ties/edge p→0,1 sane.
- **S2** Smith–Hazel index weights **b = P⁻¹ G a** (P phenotypic, G genetic covariance,
  a economic weights). Dimensions and which matrix is which are correct.
- **S3** Lush combined index **b = V⁻¹c** with Var(own)=vP, Var(fam mean)=Cov(own,fam
  mean)=vP(1+(n−1)t)/n, Cov(A,own)=vA, Cov(A,fam mean)=vA(1+(n−1)r)/n, t=r·h². Both
  weights ≥ 0; family term → 0 as h²→1.
- **S4** QGSI quadratic index **Î = w′y + y′Wy** implemented as documented; W symmetric;
  reduces to the linear index when quad weights = 0.
- **S5** Exactly one of `n`/`prop`/`intensity`; family methods require `family`; the
  `on` hook (`"pheno"`/`"gv"`/vector/function) is the GS/PS extension point and does not
  smuggle in true BV.

### O. Relationship & optimum contribution (`ocs.R`, `g_matrix`, DECISION-016)
- **O1** G = **ZZ′ / (2 Σ pⱼ(1−pⱼ))** (VanRaden 2008 method 1); monomorphic markers
  dropped; genomic inbreeding **Fᵢ = Gᵢᵢ − 1**; Z centered by 2pⱼ.
- **O2** OCS maximizes **c′g − (λ/2) c′Gc** on the simplex (c ≥ 0, 1′c = 1); group
  coancestry = **½ c′Gc**; `target/max_coancestry` met by bisection on λ; Frank–Wolfe
  stays feasible (convex combos of vertices) and terminates.
- **O3** `sample_parents()` turns contributions into an integer parent set without
  distorting the intended contribution proportions.

### U. Cross usefulness (`usefulness.R`, DECISION-016)
- **U1** **U = μ + i·σ** over the selected fraction; family **simulated** with the
  crossing engine (linkage enters σ), scored on the template's **fixed additive effects**
  via dosage×effect — **not** re-`simulate_phenotype()` (which rescales variance and
  flattens between-cross σ). This is the exact bug already fixed once — guard it.
- **U2** Additive/breeding-value basis only; DH/inbred families carry no dominance.

### P. Pleiotropy / correlation (SPEC §13, DECISION-013)
- **P1** Σ[i,i] = πᵢ·Vᵢ; Σ[i,j] = cor_ij·√(Vᵢ·Vⱼ); trait-specific var = (1−πᵢ)·Vᵢ, drawn
  independently per trait, so each trait's total genetic var is Vᵢ and every pair realizes
  cor_ij in expectation.
- **P2** Attainability = **Σ positive semi-definite** (checked by eigenvalue; error names
  the smallest). Two-trait reduces to **cor² ≤ π₁·π₂**. No silent approximation.
- **P3** Allele→genotype scaling **1/√(2·MAF·(1−MAF))** applied (scaleQTNEffects step).
- **P4** `"ld"`: two traits only; distinct causal loci in LD with r² ∈ [r2_min, r2_max];
  no SNP causal for both; `direct` vs `indirect` (shared non-causal cause-of-LD) correct.

### M. Meiosis / isqg parity (DECISION-012)
- **M1** Karlin & Liberman count-location: per chromosome, ascending — n_x ~ Poisson(L),
  L = **last map position in Morgans**; chiasmata ~ sort(Uniform(0,L)); **flip ~
  Bernoulli(0.5) is ALWAYS drawn**, even when n_x = 0.
- **M2** Deterministic XOR toggle downstream of each breakpoint; 0-based `breaks..n`
  (not the 1-based `rank > breaks`); never sort/dedupe/filter chiasmata in Rust; never
  divide by L.
- **M3** All draws in R in isqg's exact order (parent-1 then parent-2, progeny-major);
  Rust core pure; test asserts **exact** bit-equality, not distributional.

### R. Reproducibility / RNG boundary (DECISION-006/009)
- **R1** Every stochastic draw is on R's RNG; nothing random added to the Rust
  parity-critical path.
- **R2** Seed threading: `(seed, layer_index, layer_type)` for the grammar; adding/
  reordering a layer does not change other layers' realized values. Frozen
  `create_phenotypes()` keeps its own legacy seed math + `RNGversion('3.5.1')`.
- **R3** A given `seed` (or ambient `set.seed()`) fully reproduces crossing/selection.

### X. Overclaiming
- **X1** Capabilities described at true state; `synthetic_map()` output labeled a *model*,
  not a measured linkage map. No result presented that a seed cannot reproduce.

---

## Anchors
Functions: `R/grammar_*.R`, `R/arch_*.R`, `R/effects_*.R`, `R/select.R`, `R/ocs.R`,
`R/usefulness.R`, `R/schemes.R`, `src/rust/src/{numeric,genome,meiosis}.rs`.
Specs: `docs/SPEC.md` (§2 variance, §4 API, §13 PleioArch), `docs/DECISIONS.md`.
