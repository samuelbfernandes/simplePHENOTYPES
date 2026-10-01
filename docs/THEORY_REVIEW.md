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
  *J. Dairy Sci.* 91:4414–4423; Meuwissen 1997 *J. Anim. Sci.* 75:934–940;
  Zhong & Jannink 2007 *Genetics* 177(1):567–576; Lehermeier et al. 2017 *Genetics*
  207(4):1651–1661) — **verify the page against the paper**, do not assume.
- **C2** No fabricated book chapter/page (e.g. "Bernardo ch. 7"): author-year only unless
  a page is verified.
- **C3** The equation in the code matches the equation in the cited source (symbols,
  scaling, sign), not merely the same name.

### V. Variance partitioning & heritability (SPEC §2)
- **V1** Σ mean-effect genetic proportions (add+dom+epi) = requested H²; Σ all layer
  props (incl. vqtl) ≤ 1; homoskedastic residual = 1 − Σ all props. `prop` sum > 1 errors.
- **V2** Reported h² is computed from **realized** genetic/phenotypic values, not asserted
  from the budget (components not orthogonal at non-0.5 freq). vQTL share is *residual*,
  excluded from H². The transcriptome table follows V2 (DECISION-035): `h2_realized` =
  realized Var(G)/Var(P) (includes 2Cov(G, R), not bounded by 1), `h2_var_ratio` its alias,
  and the bounded allocation Var(G)/(Var(G)+Var(R)) is named `h2_allocated` (not a
  heritability). The marginal epistasis share is
  `eps/(eps + (1 − eps)/s_ct²)` with `s_ct² = 1 + 2 sqrt(ω(1 − ω)) cor(c, t)` (algebra of
  the blend in the `.tx` generator, verified numerically: 0.5714286 for `eps = 0.4`); it holds
  for the nondegenerate blend, whereas at the exact-cancellation fallback (perfect negative
  cis/trans correlation) the cis part is dropped and the share is `eps`. `h2_realized =
  h2/(1 + gr_cov)` holds only when Var(G) + Var(R) = 1; the general form is
  Var(G)/(Var(G) + Var(R) + gr_cov).
- **V3** Coding convention is the simulation one (−1/0/1 additive dosage; het-deviation
  indicator for dominance; centered a×a / a×d / d×d for epistasis), **not** Fisher's
  orthogonal average-effects decomposition — and the docs say so where it matters.
  *Exception:* `additive(orthogonal = TRUE, a =, d =)` is the orthogonal genotypic
  model (DECISION-020) — per-locus a/d, whole value scaled to `prop`, additive
  part = transmitting average effect `α = a + d(1−2p)` (as `on = "bv"`, DECISION-019),
  budget = realized `Var(A)/Var(g)`, `Var(D)/Var(g)` + an `add_dom_cov` row
  `2Cov(A,D)/Var(g)` (=0 in expectation under random mating → per-locus HWE; LD
  does **not** break it in this no-epistasis model, but arbitrary nonrandom
  multilocus genotype association does, and finite samples leave a residual).
  For this layer, check the a/d partition and its random-mating-conditional
  orthogonality, not the −1/0/1 variance-partition convention.
- **V4** `d`-type epistasis / dominance degenerate on hetless (inbred) loci is handled
  (errors or warns), not silently NaN.

### S. Selection engine (`select_ind.R`, DECISION-015)
- **S1** Truncation response tracks **R = i·Cov(A, P)/σ_P** (= **i·h²·σ_P** for
  additive models; in general exact only if E[A | P] is linear in P and equal to
  i·h²·σ_P only if Cov(A, P − A) = 0, incl. epistasis under LD); realized intensity
  **i(p) = φ(Φ⁻¹(1−p))/p**; both `direction`s handled; ties/edge p→0,1 sane.
- **S2** Smith–Hazel index weights **b = P⁻¹ G a** (P phenotypic, G genetic covariance,
  a economic weights). Dimensions and which matrix is which are correct.
- **S3** Lush combined index **b = V⁻¹c** with Var(own)=vP, Var(fam mean)=Cov(own,fam
  mean)=vP(1+(n−1)t)/n, Cov(A,own)=vA, Cov(A,fam mean)=vA(1+(n−1)r)/n, t=r·h². Both
  weights ≥ 0; family term → 0 as h²→1.
- **S4** QGSI quadratic index **Î = w′y + y′Wy** implemented as documented; W symmetric;
  reduces to the linear index when quad weights = 0.
- **S5** Exactly one of `n`/`prop`/`intensity`; family methods require `family`. The
  `on` criterion is `"pheno"`, `"gv"`, `"bv"`, a numeric vector, or a function. The
  named `"gv"` and `"bv"` are *idealized true* simulated criteria (total genetic value
  and the transmissible average-effect breeding value, DECISION-019) — legitimate in a
  simulation where the truth is known by design and explicitly labeled as true; `"bv"`
  is the one-generation random-mating transmitting ability α = a + d(q−p), refused
  under epistasis/complex architectures. The vector/function form is the GS/PS
  extension point for *estimated/predicted* values and must not be presented as a true
  breeding value. `method = "combined"` requires `on = "pheno"`; the multi-trait index
  methods score on true breeding values and ignore `on`.

### O. Relationship & optimum contribution (`select_ocs.R`, `g_matrix`, DECISION-016)
- **O1** G = **ZZ′ / (2 Σ pⱼ(1−pⱼ))** (VanRaden 2008 method 1); monomorphic markers
  dropped; genomic inbreeding **Fᵢ = Gᵢᵢ − 1**; Z centered by 2pⱼ.
- **O2** OCS maximizes **c′g − (λ/2) c′Gc** on the simplex (c ≥ 0, 1′c = 1); group
  coancestry = **½ c′Gc**; `target/max_coancestry` met by bisection on λ; Frank–Wolfe
  stays feasible (convex combos of vertices) and terminates.
- **O3** `sample_parents()` turns contributions into an integer parent set without
  distorting the intended contribution proportions.

### U. Cross usefulness (`select_usefulness.R`, DECISION-016)
- **U1** **U = μ + i·σ** over the selected fraction; family **simulated** with the
  crossing engine (linkage enters σ), scored on the template's **fixed additive effects**
  via dosage×effect — **not** re-`simulate_phenotype()` (which rescales variance and
  flattens between-cross σ). This is the exact bug already fixed once — guard it.
- **U2** Additive/breeding-value basis only; DH/inbred families carry no dominance.

### H. Crossbreeding / heterosis (`heterosis()`, DECISION-031)
- **H1** Retention of F1 heterosis by an F2 or backcross is 1/2 only under a per-locus
  condition (pure dominance): a backcross to A retains exactly 1/2 iff `h_A = 2 p_A (1 − p_A)`
  (recurrent breed only; any B); an F2 iff `[h_A − 2 p_A(1−p_A)] + [h_B − 2 p_B(1−p_B)] = 0`
  (deviations cancel). The backcross needs only the recurrent breed A in Hardy-Weinberg
  proportions; the F2 needs only that the two HWE deviations sum to zero (e.g. deviations −0.12 and
  +0.12 with neither breed in HWE), so HWE in each breed is sufficient but not necessary. A single
  fixed inbred line (`p ∈ {0,1}`, `h = 0`) qualifies; a mixture of several inbred lines that
  differ at the locus is not sufficient (`h = 0 < 2p(1−p)`). Numerically checked on a
  4×4×3×3 grid of `(p_A, p_B, h_A, h_B)` against the package's `.expected_cross_means()`;
  the two-thirds rotation fraction is the classical HWE-model statement and was not
  re-derived.

### P. Pleiotropy / correlation (SPEC §13, DECISION-013)
- **P1** Σ[i,i] = πᵢ·Vᵢ; Σ[i,j] = cor_ij·√(Vᵢ·Vⱼ); trait-specific var = (1−πᵢ)·Vᵢ, drawn
  independently per trait, so the *effect draw* gives each trait variance Vᵢ and cross-trait
  covariance cor_ij·√(VᵢVⱼ) in expectation. Each layer is then rescaled to its `prop`, so
  the realized quantity is a correlation r (covariance r·√(VᵢVⱼ)): a random ratio that
  converges to cor_ij as shared loci **and individuals** grow (with n fixed it levels off at
  the sampling spread of a correlation over n) **only under approximate linkage equilibrium among
  the causal loci and with no unit keeping a non-vanishing variance share** (major QTNs via
  n_pleio_major / prop_var_major prevent convergence); with few loci it is attenuated toward 0 on average (the multiplier
  depends on cor, π and the designs), and strong LD can prevent convergence (complete LD:
  every realized r is ±1, ensemble mean 2·asin(cor)/π — a mean, not a limit). Claims must
  say "targets/converges (given LE)", not "realizes", "equals" or "exact covariance".
  Single-shared-unit guards must cover every target strictly inside (−1, 1), `cor = 0`
  included; direction words ("inflated", "attenuated") must hold for negative `cor`. **Applies to every mean-effect layer** (DECISION-023): dominance and
  epistasis use the same Σ per component, each unit's effect divided by the *realized* sd of
  its design column (het indicator / centered product) — check that this normalizer matches
  the realization's design exactly, that the shared units are the ones common to every trait,
  and that the **total** genetic correlation (not just the additive one) tracks `cor`.
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
- **R2** Seed threading: `(seed, layer_type, occurrence of that type)` for the grammar
  (DECISION-033; a position-sensitive hash of the label, so permuted labels such as
  `rep12` / `rep21` get different sub-seeds and distinct labels are collision-resistant
  over ordinary ranges; the 31-bit sub-seed is NOT injective -- a known collision is
  `.layer_seed(123, "transcriptome_rep106", 0) == .layer_seed(123, "residual_t40160", 3)`.
  Reviewers must not accept an absolute "no two labels share a sub-seed" claim); adding or reordering a layer of
  a *different* type does not change other layers' realized values, while inserting a
  same-type layer shifts the later same-type layers' occurrence index. Frozen
  `create_phenotypes()` keeps its own legacy seed math + `RNGversion('3.5.1')`.
- **R3** A given `seed` (or ambient `set.seed()`) fully reproduces crossing/selection.

### X. Overclaiming
- **X1** Capabilities described at true state; `synthetic_map()` output labeled a *model*,
  not a measured linkage map. No result presented that a seed cannot reproduce.

---

## Anchors
Functions: `R/grammar_*.R`, `R/arch_*.R`, `R/effects_*.R`, `R/select_ind.R`, `R/select_ocs.R`,
`R/select_usefulness.R`, `R/select_schemes.R`, `src/rust/src/{numeric,genome,meiosis}.rs`.
Specs: `docs/SPEC.md` (§2 variance, §4 API, §13 PleioArch), `docs/DECISIONS.md`.
