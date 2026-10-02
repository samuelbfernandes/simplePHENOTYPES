# DECISIONS.md — simplePHENOTYPES v2 Architectural Decisions

> Canonical decision log. CLAUDE.md and ARCHITECTURE.md summarize these; this file is the
> source of truth. Each decision records the call, the rationale, and what it supersedes.
> Decisions 008–011 (2026-06-10) resolve contradictions found in a pre-implementation
> docs review and refine the original 001–007.

---

## DECISION-001: Package decomposition

**Decision:** Single package, internal modularization only (no `genoUtils`/`simBreed` split).
**Rationale:** Lower CRAN maintenance burden; one install for users; modules give clean
internal separation without multi-package interdependencies.
**Date:** 2026-06 (locked)

---

## DECISION-002: isqg integration

**Decision:** Absorb isqg's published C++ meiosis/cross/DH algorithms into our own Rust
port (`genome.rs`, `meiosis.rs`). Do not depend on the archived isqg CRAN package.
**Rationale:** Eliminates the archived-dependency risk; we own and can maintain the code;
the bitwise chromosome representation maps cleanly to Rust.
**Date:** 2026-06 (locked)

---

## DECISION-003: Backward compatibility for create_phenotypes()

**Decision:** `create_phenotypes()` is preserved with its full v1 signature.
**Rationale:** The published API is contractual; existing user scripts must keep working.
**Superseded in part by DECISION-008** — the original "compatibility shim that maps old
calls onto new grammar internals" reading is replaced by "frozen legacy function."
**Date:** 2026-06 (locked; refined by 008)

---

## DECISION-004: Multi-generation scope

**Decision:** Multi-generation cross simulation stays inside simplePHENOTYPES (not a
separate `simBreed` package).
**Date:** 2026-06 (locked)

---

## DECISION-005: Expression-based simulation timeline

**Decision:** Deferred to v3.
**Rationale:** Different input semantics (transcript vs marker); designing it alongside the
marker-based redesign would muddy both.
**Date:** 2026-06 (locked)
**Status:** **superseded by DECISION-022** (2026-10-02). Once the grammar, the crossing
engine and the fixed-scale accessors were stable, expression simulation was brought
forward into v2 as a process *downstream of the genome* (`simulate_transcriptome()`,
`transcriptome()`, `observe_counts()`; `SPEC-transcriptome.md`) rather than as a separate
transcript input format, which removes the "different input semantics" objection.

---

## DECISION-006: Rust scope — surgical, not a rewrite

**Decision:** Only deterministic, profiled-slow code and C++-origin (isqg) code moves to
Rust: `as_numeric()` numericalization, isqg meiosis/cross/DH, and genetic-value assembly
*only if profiling proves it matters*. The stochastic core (QTN sampling, effect-series
and residual draws, architecture mapping) stays in R. Rust never calls an RNG on a
parity-critical path.
**Rationale:** R and Rust RNGs differ; keeping all random draws in R preserves
reproducibility. CRAN acceptance is non-negotiable and a small surgical Rust surface is
the safest path to it.
**Date:** 2026-06 (locked; reaffirmed by 011)

---

## DECISION-007: PleioArch for the "pleiotropy" architecture

**Decision:** Adopt the PleioArch algorithm (`context/PleioArch-main/`) as the
effect-generation engine for `architecture = "pleiotropy"` — genetic-correlation
control via bivariate-normal effect draws.
**Rationale:** v1 pleiotropy merely shares loci and cannot control rho_g; PleioArch draws
effects whose cross-trait covariance equals the target in expectation, so the realized
correlation targets `cor`. *[Wording corrected 2026-09-25 (DECISION-023): the effect draw's cross-trait covariance and variances equal their targets in expectation, so `cor` is the ratio of expected moments; the realized correlation is a random ratio, not equal to `cor` in expectation — see DECISION-023, Scope / limits.]*
**Refined by DECISION-010** (the user-facing control is `rho_g`, replacing v1's buggy
`cor`).
**Date:** 2026-06 (locked; refined by 010)

---

## DECISION-008: create_phenotypes() is frozen legacy, not a delegation shim

**Question:** Should `create_phenotypes()` be reimplemented as a thin wrapper that
delegates to the new grammar internals, or kept as its own code path?

**Decision:** Keep `create_phenotypes()` as a **frozen legacy function**: bugfix-only,
marked `lifecycle::badge("superseded")`, retaining its original v1 code paths, seed
arithmetic, and `RNGversion('3.5.1')`. It is **not** a wrapper over the new grammar; the
legacy engine and the new grammar coexist as separate implementations.

**Rationale:** Making it delegate would force the new grammar to reproduce v1's messy,
h2-dependent seed arithmetic and exact RNG call order bit-for-bit — re-importing the very
fragility the redesign exists to remove. Freezing the legacy path keeps v1 reproducible
for existing users while leaving the new grammar free to be clean.

**Supersedes:** the "compatibility shim delegates to new grammar internals" interpretation
of DECISION-003 (and TODO.md's "Update create_phenotypes() shim to delegate to new grammar
internals").
**Date:** 2026-06-10

**Addendum (2026-09-13) — loader unified onto `as_numeric()`:** the *frozen* part of
`create_phenotypes()` is its simulation engine (QTN/effect/seed math, RNG order,
`RNGversion('3.5.1')`), not the genotype-file reading. The loader (`genotypes()`, formerly
via `file_loader()`) now performs its **numericalization** through the shared
`as_numeric()` / `format_conversion()` pipeline instead of a duplicate coder, so the whole
package has one genotype-coding implementation. `file_loader()`, `numericalization()`, and
`table_to_numeric()` (the old duplicate coders) were removed. This is a bug-fix, not a
grammar delegation: the buggy legacy coder mis-coded some HapMap/table cases (het coded as
a homozygote, major allele coded −1, imputation skipping the genetic-model transform), and
the shared kernel is correct. Numeric-input normalization/imputation still happens in the
legacy path (it is not numericalization). v1 output is preserved where it was correct —
`test-v130-parity.R` (numeric `SNP55K` input) stays green — and the frozen simulation
engine is unchanged. "Frozen" therefore means the *engine*, and the reproducibility
guarantee is against the v1.3.0 reference outputs, not against the old coder's bugs.

---

## DECISION-009: New grammar has no bit-for-bit parity obligation to v1

**Question:** Must the new grammar reproduce v1.3.0 phenotypes bit-for-bit?

**Decision:** No. The new grammar is an independent, clean implementation. SPEC §6's clean
`(seed, layer_index, layer_type)` seed-threading **stands** — it need not replicate v1's
seed math. The captured `inst/extdata/v1_3_0_reference/*.rds` references change role:
- a **frozen-regression guard on `create_phenotypes()`** — bug fixes must not change its
  output unintentionally; any change must be deliberate and re-blessed;
- **statistical / structural** validation for the new grammar (variance-partition
  identity, realized rho_g ≈ target, QTN-count structure) — not bit equality.

`test-v130-parity.R` must be repointed at `create_phenotypes()` (it currently calls the
bare grammar and asserts bit-identity, which is now the wrong target).

**Rationale:** Bit-parity between a clean redesign and the messy v1 engine is achievable
only by preserving v1's exact draw order, which defeats the redesign. DECISION-008 keeps
v1 reproducible on its own frozen path, so the grammar does not need to.

**Supersedes:** the "v1.3.0 bit-identical parity is NON-NEGOTIABLE and gates the grammar"
framing in CLAUDE.md (Rule #5) and SPEC §8.1.
**Date:** 2026-06-10

---

## DECISION-010: genetic-correlation control via `cor` (PleioArch + Cholesky fallback)

**Question:** How does v2 expose genetic-correlation control?

**Decision (amended 2026-06-11):** The grammar's pleiotropy correlation argument is
named **`cor`** (the v1 name is reused, not a new `rho_g`). Its engine depends on the
number of traits:
- `n_traits = 2`: the **PleioArch** bivariate-normal engine (DECISION-007) targets
  `cor` (effect covariance equal to the target in expectation); enforces
  `cor² ≤ pi_target × pi_secondary`. *[Wording corrected 2026-09-25 (DECISION-023): the effect draw's cross-trait covariance and variances equal their targets in expectation, so `cor` is the ratio of expected moments; the realized correlation is a random ratio, not equal to `cor` in expectation — see DECISION-023, Scope / limits.]*
- `n_traits > 2`: fall back to v1.3's **Cholesky** decomposition
  (`base_line_multi_traits.R`) to impose `cor` (scalar per pair, or an
  `n_traits × n_traits` matrix), emitting a **warning** that the phenotypes are
  correlated as requested but individual QTN effect sizes are not guaranteed.

**Rationale:** v1's Cholesky `cor` is unreliable for >1 trait, so PleioArch is the
correct engine where it applies (2 traits). Extending the exact bivariate algorithm to
many traits simultaneously is non-trivial; reusing the (imperfect but serviceable)
Cholesky path for >2 traits, behind an explicit warning, keeps the feature available now
without overclaiming precision. Keeping the name `cor` eases v1→v2 migration.

**Supersedes:** the original 2026-06-10 decision to drop the name `cor` in favor of
`rho_g` and to drop the Cholesky path entirely.

**Implication:** SPEC §4.1/§11/§13 use `cor`; the Cholesky path is retained (not dropped)
for the >2-trait fallback. v1's `cor`/`cor_res` also survive inside the frozen
`create_phenotypes()`.
**Date:** 2026-06-10 (amended 2026-06-11)

---

## DECISION-011: Single rextendr package; Cargo-workspace plan dropped

**Question:** Single rextendr package, or a multi-crate Cargo workspace
(`core` + `r-pkg` + `py-pkg`)?

**Decision:** Single rextendr package with `src/rust/` (per ARCHITECTURE.md §4). The
Cargo-workspace bootstrap in `TODO.md` (appendix, legacy item 7) (root `Cargo.toml`, members
`core`/`py-pkg`, `git mv` into `r-pkg/`, `maturin new py-pkg`) is **dropped**.

**Rationale:** CRAN acceptance is non-negotiable; a single rextendr package is the
proven-acceptable pattern, while a multi-crate workspace with a maturin/PyO3 member adds
CRAN risk for no near-term benefit. Python/Shiny sharing can be revisited later without
blocking v2.

**Reaffirms:** DECISION-001 and DECISION-006.
**Date:** 2026-06-10

---

## DECISION-012: meiosis randomness stays in R; exact isqg bit-parity is the gate

**Question:** isqg's meiosis is irreducibly stochastic. Does the Rust port draw its own
randomness (accepting statistical-only parity), call back into R's RNG, or consume
random draws made in R? ARCHITECTURE.md §6 and DECISION-006 gave contradictory answers.

**Decision:** **R draws every random quantity; the Rust meiosis core is a pure function
of pre-drawn randomness; exact bit-for-bit parity with isqg is the test gate.**

R performs the draws in isqg's exact order — per whole-genome meiosis event, per
chromosome in ascending order:

```
n_x       ~ rpois(1, L)          L = LAST map position, in Morgans
chiasmata ~ sort(runif(n_x, 0, L))    (not drawn when n_x == 0)
flip      ~ rbinom(1, 1, 0.5)         (ALWAYS drawn, even when n_x == 0)
```

and passes `(chiasmata, counts, flips)` into Rust, which performs only the deterministic
XOR chain and haplotype assembly. Rust never calls an RNG.

**Rationale:** verified empirically before adopting. A pure-R reimplementation of the
draw order reproduces isqg's `spc$gamete()` bitstrings byte-for-byte across five seeds
plus a continued multi-gamete stream, and reproduces `selfcross()` and `dh()` genotype
matrices exactly across three more seeds. So the DECISION-006 constraint ("Rust never
calls an RNG on a parity-critical path") costs nothing here — it is compatible with the
*strongest* available parity standard rather than forcing a weaker one. Exact parity
also catches transcription bugs (an off-by-one in the breakpoint index, a dropped
Bernoulli) that distributional checks would pass.

isqg draws from R's own RNG via Rcpp's `Rmath.h` wrappers inside an `Rcpp::RNGScope`
(`context/isqg/src/Genetics.cpp:74,84,93`), with no C++ `<random>` engine anywhere, which
is why the streams can be made to coincide at all.

**Consequences:**
- `test-isqg-parity.R` asserts exact equality, not distributional agreement.
- The port is pinned to the Karlin & Liberman count-location process and to isqg's draw
  order (parent-1 then parent-2, progeny-major). A different recombination model would be
  added alongside, not substituted.
- Reference fixtures store the drawn `counts`/`chiasmata`/`flips` next to isqg's output,
  so a failure distinguishes an error in the R draw order from an error in the Rust core.

**Supersedes:** ARCHITECTURE.md §6's "validated for statistical equivalence ... exact
bit-match is a goal where feasible but the firm requirement is correct recombination
behavior."
**Reaffirms:** DECISION-002 and DECISION-006.
**Date:** 2026-09-07

---

## DECISION-013: PleioArch generalized to any number of traits

**Question:** DECISION-010 restricted exact correlation control to two traits and sent
`n_traits > 2` to a Cholesky fallback that warned individual QTN effects were not
guaranteed. Can the exact engine cover any number of traits?

**Decision:** **Yes — the PleioArch engine now handles any `n_traits`, and the Cholesky
fallback is removed.** The bivariate normal generalizes directly to an `n x n`
multivariate normal:

```
Sigma[i,i] = pi_i   * V_i                 pleiotropic share of trait i
Sigma[i,j] = cor_ij * sqrt(V_i * V_j)     the whole genetic covariance
```

Trait-specific effects stay univariate with variance `(1 - pi_i) * V_i`. Because
trait-specific loci are independent across traits, each trait's total genetic variance is
still `V_i` while the entire covariance comes from the shared loci, so every pair targets
`cor_ij` (covariance and variances equal their targets in expectation). *[Wording corrected 2026-09-25 (DECISION-023): the effect draw's cross-trait covariance and variances equal their targets in expectation, so `cor` is the ratio of expected moments; the realized correlation is a random ratio, not equal to `cor` in expectation — see DECISION-023, Scope / limits.]* For two traits this reduces algebraically to the bivariate
reference implementation, so DECISION-007 is preserved rather than replaced.

`cor` accepts a scalar (applied to all pairs) or a full `n_traits x n_traits` matrix,
negative correlations included. `pi` accepts a scalar or one value per trait;
`pi_target`/`pi_secondary` remain as the two-trait spelling.

**Feasibility:** the `cor^2 <= pi_i * pi_j` constraint generalizes to "`Sigma` must be
positive semi-definite", checked by eigenvalue. This is strictly more informative: with
three or more traits it also catches mutually inconsistent requests that no pairwise check
would (three traits cannot all be strongly negatively correlated). Unattainable requests
raise an error naming the smallest eigenvalue, rather than being silently approximated.

**Rationale:** the Cholesky fallback correlated the *genetic values* after the fact, so
the QTN effects no longer corresponded to the realized correlation — the very defect
DECISION-007 adopted PleioArch to remove. Verified on `SNP55K_maize282_maf04` at
`n_qtn = 300` over 25 seeds: realized correlations match target at 2, 3, 5 and 8 traits
with no loss of precision as traits are added (mean 0.69 for a 0.7 target, SD ~0.05
throughout), and an asymmetric target matrix `(0.8, -0.4, -0.2)` realizes as
`(0.78, -0.39, -0.18)`.

**Supersedes:** DECISION-010's two-trait restriction and its Cholesky fallback. The rest
of DECISION-010 stands: the argument keeps the v1 name `cor`, and v1's `cor`/`cor_res`
survive inside the frozen `create_phenotypes()`.
**Reaffirms:** DECISION-007.
**Date:** 2026-09-07

---

## DECISION-014: the `"ld"` architecture simulates linked *distinct* causal loci

**Question:** What does `architecture = "ld"` mean in the v2 grammar? The first
implementation drew each trait's QTNs independently (as for `"independent"`) and then, per
causal QTN, annotated a *companion* tag marker in an r2 window for **reporting only** —
under `ld_type = "indirect"` the tag was reported in place of the (hidden) causal marker.
That is a single-trait, GWAS-style "you observe a tag SNP" annotation: the companion never
enters any genetic value, and the two traits' causal loci are not linked to each other, so
the traits are effectively independent. It does not reproduce the classic simplePHENOTYPES
linkage architecture (`legacy_QTN_linkage.R`), whose point is a *spurious* cross-trait
correlation produced by linkage between distinct causal loci.

**Decision:** `"ld"` is a **two-trait** architecture (`n_traits = 2`, else an error) in
which the two traits have **distinct** causal loci that sit in linkage disequilibrium, so
the traits covary through linkage rather than through a shared locus. For each of `n_qtn`
loci a linked pair is drawn — one SNP causal for trait 1, the other for trait 2, with
squared correlation r2 in `[r2_min, r2_max]`. No SNP is causal for both traits; the genetic
correlation is a consequence of the linkage, not of shared effects. Two flavors:

- `ld_type = "direct"` (**default**): trait 1's and trait 2's causal SNPs are directly in
  LD.
- `ld_type = "indirect"`: both causal SNPs flank a shared, **non-causal** cause-of-LD locus
  and are each in LD with it (one flanking marker upstream, one downstream).

`qtn_table()` reports the linked pair on each row: `QTN_t1` / `QTN_t2` name the two traits'
causal SNPs and `ld_r2` is the squared correlation between them (corrected 2026-09-29: this
paragraph used to describe a `companion` column and a flank-to-cause r2 that the code no longer
produces; for `"indirect"` the hidden cause-of-LD locus is kept on the layer object but is not a
QTN and is not shown). Partner rule (documented 2026-09-29): for a randomly drawn anchor the
`"direct"` partner is its **strongest** in-window partner (highest r2), so the realized r2 skews
toward `r2_max`; the opt-in `partner = "random"` takes a uniformly random in-window partner;
markers with r2 = 0 or (numerically) 1 are never partners, and the window needs `r2_max > 0`,
`r2_min < 1`. The old reporting-only `.annotate_ld()` / `.ld_reported_qtn()` helpers are removed; the
linkage is now established during QTN sampling (`.draw_qtn_ld()`), RNG staying in R
(DECISION-006).

**Rationale:** the companion-annotation model answered a different question (tag-SNP
observation for one trait) than the one the architecture name implies and than the v1
engine implements (linked-but-not-pleiotropic cross-trait correlation). Verified on
`SNP55K_maize282_maf04` (`seed = 200`, `n_qtn = 3`): the realized genetic correlation
between traits is strong under `"ld"` (|r| ≈ 0.67 direct / 0.76 indirect) and ≈ 0 for the
same setup under `"independent"`, with the two traits' causal-locus sets disjoint. `direct`
is the default because it is the more common and more directly interpretable request (the
two observed causal SNPs are themselves the linked pair); `indirect` remains for the hidden
cause-of-LD scenario. This is the grammar counterpart of the frozen legacy `qtn_linkage`
engine, which keeps its own behavior (DECISION-008).

**Scope:** two traits only, by design (matching v1). Generalizing linked causal loci to
`n_traits > 2` (round-robin or star topologies) was considered and deferred — no v1
precedent and no current demand.

**Reaffirms:** DECISION-006 (RNG stays in R), DECISION-009 (grammar owes v1 no bit-parity —
this is a clean reimplementation of the *concept*, not the v1 code path).
**Date:** 2026-09-09

---

## DECISION-015: the selection engine and named breeding-scheme wrappers

**Question:** The breeding-program designer needs a "selection bucket" (SSD, bulk,
etc.), but the package had no selection or generation-advance code — only crossing.
What is the selection API, which methods, and how do schemes compose?

**Decision:** Add a headless, tested selection engine and scheme wrappers, built
ahead of any UI and aligned with the standard texts (Bernardo; Falconer & Mackay;
Lynch & Walsh).

- **`select_ind()`** (not `select()` — avoids masking `dplyr::select`, and reads as
  "select individuals"). Truncation selection returning the chosen individuals as a
  crossable `Population` subset (or their ids if the sim was not Population-backed),
  carrying attributes `selected`, `differential` (S), `intensity` (realized i),
  `criterion`, `method`.
  - **Criterion `on`:** `"pheno"` (observed phenotype — realistic mass selection,
    response tracks R = i·Cov(A, P)/σ_P, which is i·h²·σ_P for additive models or when
    Cov(A, P − A) = 0 and E[A | P] is linear; see the 2026-09-29 amendment), `"gv"` (true genetic value — idealized upper
    bound), a **numeric vector** (named by id or in population order), or a
    **function** `on(sim)`. The vector/function is the deliberate extension point
    for **genomic selection, phenomic selection, and other predictive criteria**
    computed externally. True breeding value (`on = "bv"`, the additive average
    effect α) is deferred to the orthogonal genotypic-model / average-effects
    feature (ROADMAP §3).
  - **Methods:** `mass` (individual truncation), `within_family`, `among_family`,
    `combined` (Lush combined index), `index` (Smith–Hazel multi-trait economic
    index, `b = P⁻¹ G a`), `random` (drift control). Family methods require a
    `family` grouping.
  - **`combined`** ranks on the selection-index prediction of breeding value from
    the individual's own record and its family mean, `b = V⁻¹ c` built from `h2`
    and `family_relationship`, the within-family correlation of breeding values
    A_ij/√(A_ii A_jj) (for families of non-inbred, unrelated parents: 0.25
    half-sibs default, 0.5 full-sibs, 0.5 doubled haploids, 2/3 S1 sibs).
    *[Corrected 2026-09-27: this was first recorded as the "additive relationship"
    with "0.5 full-sibs/selfed". Since `r` enters as t = r·h² with h² the
    candidates' own heritability, it is the correlation of breeding values, which
    equals A_ij only for non-inbred members; S1 sibs of a non-inbred plant have
    A_ij = 1, A_ii = 1.5, so r = 2/3. Inbred or related parents change all these
    values (the package's derivation: 2(1+F)/(3+F) for S1 sibs, (1+F)/2 for doubled
    haploids of a parent with inbreeding F). The code was right; the advice was not.
    Derivation (tabular A, parent P with inbreeding F): coancestry of two gametes of
    P is Θ_PP = (1+F)/2, so two S1 sibs have A_ij = 2Θ_PP = 1+F and each has
    inbreeding Θ_PP, A_ii = 1 + (1+F)/2, giving r = 2(1+F)/(3+F); two doubled
    haploids also have A_ij = 1+F but A_ii = 2, giving r = (1+F)/2. Same review:
    the index weights depend on family size but were applied to raw records, so with
    unequal families the ranking depended on the phenotype origin (adding a constant
    changed the selection); the score now uses deviations from the candidate mean
    (Hazel 1943), a family of one is scored h²·deviation, and the docs state that one
    h² and one r are shared by all families (mixing family types is out of scope).
    The singular-case test also used an absolute floor on the denominator, so tiny-unit
    records fell back to own-record scoring and changed the selection; it now tests
    the dimensionless 1 − r·h², and the weights are computed in their vP-free form
    b₁ = h²(1−r)/(1−rh²), b₂ = h²nr(1−h²)/{[1+(n−1)rh²](1−rh²)} (equal to V⁻¹c to
    ~1e-15; the vP² products under/overflowed at extreme scales), with
    1 − rh² = (1−r) + r(1−h²) so the own-record fallback is needed only at exactly
    r = h² = 1 (a tolerance there reversed valid near-singular selections). The docs
    no longer call the index "optimal" outright: it is optimal under the additive
    model t = r·h²; dominance / epistasis / family environment add to full-sib and
    selfed covariances. The Lush (1947) reference gained Part II's DOI.]* Derived from first principles
    (Var(own)=vP; Var(fam mean)=Cov(own,fam mean)=vP(1+(n−1)t)/n; Cov(A,own)=vA;
    Cov(A,fam mean)=vA(1+(n−1)r)/n; t=r·h²). Both weights are non-negative and the
    family term vanishes as h²→1 (own record becomes sufficient). Requires `h2`.
  - **Intensity:** exactly one of `n`, `prop`, or standardized `intensity` (mapped
    to a count via i(p)=φ(Φ⁻¹(1−p))/p). Both `direction`s (`high`/`low`).

- **Scheme wrappers** compose `select_ind()` with the crossing primitives, seeding
  the RNG once and threading the ambient stream through `seed = NULL` primitive
  calls (so one `seed` reproduces the whole scheme):
  - `single_seed_descent()` — one selfed seed per line per generation, no
    selection; line count preserved.
  - `bulk()` — mass selfing, progeny pooled and subsampled to a bulk size; no line
    identity.
  - `pedigree()` — self + select each generation; each selected line is selfed into
    an equal-sized family and the families are pooled to a fixed `pop_size`, so the
    population does not drift down as it inbreeds.
  - `recurrent_selection()` — select parents, intercross them (`n_crosses` random
    pairs × `progeny_per_cross`) to form the next cycle.
  Both selecting schemes take a **`phenotype` callback** (`Population → realized
  phenotype_sim`), applied each generation, so the caller controls the genetic
  model (and can fix causal loci with `additive(qtn = ...)`). Each returns a
  `Population` carrying a per-generation `history` (n selected, S, i).
  - **`c.Population()`** — pools populations sharing a marker map (column-binds
    haplotypes, uniquifies ids), the primitive the wrappers use to reassemble a
    generation from per-parent progeny.

**Efficiency (DECISION-006 preserved):** selection ranking and the stochastic
orchestration stay in R; the meiosis core stays in Rust; Populations are
reference-based, so selected individuals cross forward without copying genotypes.

**Rationale:** matches the designer's node buckets while keeping the primitives
usable and testable on their own; the method set is the standard textbook coverage;
the custom-criterion hook gives the maintainer the requested space for GS/PS and
other predictive approaches without committing the engine to any one estimator.
AlphaSimR (Gaynor et al. 2021) is the quality reference, not a competitor —
simplePHENOTYPES keeps its ease-of-use edge.

**Documentation constraint:** every selection/designer decision is mirrored in the
designer manuscript (`breeding_designer/manuscript.tex`, "Design decisions log" +
Methods), per the standing instruction, so the code and paper cannot drift.

**Amendment 2026-09-29 (selection audit; wording and contracts, no equation changed):**
- *Smith-Hazel `G`.* `method = "index"` keeps `G = Cov(A)` (the covariance of the
  true breeding values), the standard case. Hazel's (1943) g_ij is
  Cov(phenotype_i, breeding value_j); `Cov(P, A) = Cov(A)` needs Cov(A, D) = 0
  (random mating / Hardy-Weinberg). In selfed, inbred or previously selected
  populations with dominance the two differ, so the index is Smith-Hazel *under
  that assumption* (documented in `select_ind()`). The `on = "pheno"` response
  sentence likewise reads R = i·Cov(A, P)/σ_P, not "governed by V_A". This is the
  linear-regression prediction of the response of the mean breeding value: it is
  exact only if E[A | P] is linear in P (e.g. A and P jointly normal), not an
  identity for arbitrary joint distributions; it reduces to i·h²·σ_P (h² = V_A/V_P)
  only when Cov(A, P − A) = 0, i.e. every non-additive component of the phenotype
  (dominance, epistasis, residual) is uncorrelated with A -- including epistasis
  under linkage disequilibrium, where the centred epistatic products need not be
  orthogonal to A (mirrors the `select_ind()` roxygen). Phenotype rows are joined to individuals
  by id, not by position.
- *Contracts.* `trait` must be one in-range index (a vector is an error, never
  recycled; tandem schedules are checked against the callback's `n_traits`);
  `family` may not contain `NA`; `direction` is one value; `pedigree()` refuses
  `prop` and `n_select` together; `recurrent_selection()` refuses
  `n_parents >= N`; scheme callbacks must be backed by the population they are
  given; non-finite/overflowing scores or SDs are errors.
- *Statistics/semantics.* `method = "random"` reports S and i realized on `on`
  (expected 0) and validates `on`; `among_family` keeps whole families (`n` is a
  floor); `bulk()` draws its `n` seeds uniformly from the pooled progeny
  (multinomial parental contributions, random even for `n = k·N`); the scheme
  wrappers validate `seed` and restore the caller's RNG state.
- *Scale.* `pedigree()`/`recurrent_selection()` histories are on the scale of the
  callback's phenotype; a `simulate_phenotype()` callback re-standardises the
  genetic layer per generation (see DECISION-021), so the fixed-scale response is
  obtained by ranking on `phenotype_value()` through `on`.

**Addendum 2026-09-30 (SPEC-0020 items 7 and 8; no equation changed):** `within_family`
accepts an explicit per-family count (`n_per_family`, the last argument of `select_ind()`): one
whole number applied to every family, or a vector named by family label (or unnamed in sorted
family-label order). It is exclusive with `n`/`prop`/`intensity` and with every method other
than `within_family`, and it is never capped silently: a family smaller than the request is an
error. Within each family the top scores (after `direction`) are kept; S and i keep their
definitions on the pooled candidates (S = mean(selected) - mean(all), i = S / sd(score)), so with
unequal per-family counts S is the realized differential of the pooled selected set. The result
carries the realized counts as attribute `n_per_family` only when the argument is used.
`pedigree()`/`recurrent_selection()` stay mass-selection only (they have no `family`) and do not
forward it. `sample_parents()` returns a `source` attribute (data frame `slot`, `id`, `index`,
`name`, one row per slot in draw order): a pure record attached after the draw, no RNG consumed,
the draws unchanged for both `method = "allocate"` and `"multinomial"`.

**Reaffirms:** DECISION-006 (RNG in R, meiosis in Rust), DECISION-004
(multi-generation stays in-package).
**Date:** 2026-09-10

---

## DECISION-016: modern selection methods — OCS, genomic relationship, cross usefulness

**Question:** Beyond the textbook truncation/index methods (DECISION-015), which
*modern* methods are relevant enough to build into the selection engine, and how,
given the CRAN dependency constraints?

**Decision:** Survey of modern methods and what each needs:

| Method | Engine work | Call |
|--------|-------------|------|
| Genomic selection (GBLUP/rrBLUP, Bayes) | none — plugs into the `on = <vector/function>` hook | `select_ind(on = gebv)` |
| Phenomic selection (Rincent 2018) | none — same hook (NIRS-predicted values) | `select_ind(on = pred)` |
| **Optimum contribution selection** (Meuwissen 1997) | **new** — needs a relationship matrix + constrained optimizer | `optimum_contribution()` |
| **Genomic relationship matrix** (VanRaden 2008) | **new** — infrastructure for OCS + inbreeding reporting | `g_matrix()` |
| **Cross usefulness** (Zhong & Jannink 2007; Lehermeier 2017) | **new** — simulate a family, score μ + iσ | `cross_usefulness()` |
| Weighted/optimal GS, look-ahead mating (Moeinizade 2019) | large/advanced | deferred (roadmap) |

**Built this pass:**

- **`g_matrix()`** — VanRaden (2008) method-1 genomic relationship matrix
  \eqn{G = ZZ'/(2\sum p_j(1-p_j))} from the marker dosages already held; drops
  monomorphic markers; diagonal gives genomic inbreeding \eqn{F_i = G_{ii}-1}.
  Genomic **G**, not pedigree **A** (chosen: we have genotypes, no pedigree
  tracking needed; A deferred).
- **`optimum_contribution()`** — maximizes \eqn{c'g - \tfrac{\lambda}{2}c'Gc} over
  contributions on the simplex (`c ≥ 0`, `1'c = 1`); group coancestry is
  \eqn{\tfrac12 c'Gc}. Give `lambda` to trace the gain–diversity frontier, or a
  `target_coancestry`/`max_coancestry` met by bisection on `lambda`. The optimizer
  is a **dependency-free Frank–Wolfe active set** on the simplex (exact line search
  per step; no QP-solver dependency — chosen to respect the CRAN/`Imports` budget).
  Merit uses the same `on` hook (gv/pheno/custom), so GS/PS EBVs drive OCS.
  `sample_parents()` turns contributions into a parent set for mating and, by
  owner decision (2026-09-29: keep the OCS optimization), **allocates** slots
  (`method = "allocate"`, default): \eqn{\lfloor n c_i \rfloor} each, then the
  remaining slots to the largest fractional parts (ties: larger contribution, then
  random with `seed`), so every count is within one slot of \eqn{n c_i}. This
  controls the marginal count error; it is an approximation, not a guarantee on the
  quadratic group coancestry \eqn{\tfrac12 c'Gc}, and the realized-vs-optimized gap
  depends on `n` and `G` (n = 1, G = I, c = (0.5, 0.5): optimum 0.25, realized 0.5).
  Measured on an F2 example (target 0.0521): allocation 0.0518 at `n = 30`,
  independent multinomial draws 0.0672 on average (sd 0.015). Only the slot ORDER is
  shuffled.
  The former independent weighted draw with replacement remains as
  `method = "multinomial"` (bit-identical for the same seed) for users who want
  sampling variation; its realized coancestry is on average above the optimum, more so for
  small `n`. Also fixed (audit 2026-09): the length-1 `sample()` gotcha and the
  ambient-RNG handling. The merit
  `"bv"` is the *transmitting average effect* \eqn{\alpha = a + d(q-p)} of
  DECISION-019 (not a least-squares/Fisher projection), refused under an epistasis
  layer or `architecture = "complex"`. A supplied `G` must carry identical, unique
  individual ids as dimnames (audit: an unnamed `G` no longer yields a corrupt
  `ocs` object), and a target coancestry above the unconstrained optimum's is
  reported by a warning (`lambda = 0`) instead of being met silently; the penalty
  search is rescaled to unit merit spread for extreme merit scales. OCS is invariant
  to a constant added to all merits (merit centred internally; reported on the
  original scale); target-above-unconstrained detection uses a machine-epsilon-relative
  band (16*eps*max(1,|target|,|c0|)), no absolute 1e-6 band (review round 2,
  2026-09-30).
  **Note 2026-09-30 (audit-fix implementation record):** the audit-fix record's OCS-F2
  entry ("DOC-ONLY per D5 ... sampling scheme unchanged") was superseded by the owner's
  D5 decision: `sample_parents()` defaults to `method = "allocate"`; the independent
  weighted draw is kept as `method = "multinomial"`.
  **Note 2026-09-30 (audit-fix implementation record, round 3):** round 3 closes the
  Codex round-2 findings (R3-1 to R3-15): numerical fixes in OCS centring, supplied-`K`
  symmetrization, the transcriptome `h2_*` fields, allele-label and HapMap detection,
  integer genotype columns and the v1 seed-overflow message, plus wording corrections
  (allocation coancestry comparison, scheme accuracy, heterosis criteria, the one-A/one-D
  variance identity, `h2_realized` scaling, the epistasis-share fallback). See NEWS.md
  (Review round 3).
  **Note 2026-09-30 (audit-fix implementation record, round 4):** round 4 (R4-1 to R4-5)
  makes the OCS above-optimum band purely relative (no `max(1, .)` absolute floor; it
  supersedes the `16*eps*max(1,|target|,|c0|)` band in the round-2 note above), makes
  the v1 seed-overflow message state the accepted interval `[lo, hi]` (or "no seed
  accepted") when no interval centred at 0 is accepted (`abs(seed) <= N` only in the
  ordinary case), makes the transcriptome mimic rescale always hit the requested
  per-gene variance when the realized unit variance is finite and positive (warning
  when it is < 1e-12), and mirrors the `select_ind()` response assumption
  (`E[A | P]` linear; `Cov(A, P - A) = 0` for the `i*h2*sigma_P` reduction) in
  DECISIONS/THEORY_REVIEW/ROADMAP. See NEWS.md (Review round 4).
- **`cross_usefulness()`** — ranks candidate biparental crosses by
  \eqn{U = \mu + i\,\sigma}, the expected value of the best `select_top` fraction of
  progeny. The family is **simulated** with the crossing engine (so linkage enters
  the variance, not an approximation) and scored on the template's **fixed additive
  effects** via a plain dosage×effect dot product — deliberately **not** through a
  re-`simulate_phenotype()` call, which rescales each family's genetic values to a
  fixed variance and would flatten the between-cross σ the criterion depends on.
  Additive (breeding-value) basis only; DH/inbred families carry no dominance.
  **Heterozygous parents (audit USE-F1):** the `"dh"`/`"selfcross"` schemes derive
  each family from a *single* F1 individual. For inbred parents that F1 is unique
  and the family moments are exact; for a heterozygous parent the reported mean and
  sd are conditional on the one F1 drawn (they vary with the seed; sd is on average
  below the F1-averaged variance). This is documented in the help and a warning is
  emitted once per call; the sampling scheme is unchanged (`scheme = "cross"` draws
  every progeny from an independent meiosis). The intensity \eqn{i} is the
  infinite-population normal value (1.755 at `select_top = 0.1`), so \eqn{U} is a
  ranking score, and a family with `n_progeny * select_top < 1` is only a
  normal-theory extrapolation.

**Deferred to roadmap (user request):** a PopVar-style function that selects the
best crosses from **estimated marker effects on real training data** (Mohammadi,
Tiede & Smith 2015), i.e. the real-data counterpart of `cross_usefulness()`'s
known-effect simulation. Also deferred: weighted/optimal GS and multi-generation
look-ahead mating.

**Efficiency (DECISION-006 preserved):** the relationship matrix, the optimizer,
and the usefulness scoring are deterministic linear algebra in R; the family
simulation inside `cross_usefulness()` uses the Rust crossing core; Populations stay
reference-based.

**Reaffirms:** DECISION-006 (deterministic R, meiosis in Rust), DECISION-015 (the
`on` hook is the GS/PS extension point).
**Date:** 2026-09-11

---

## DECISION-017: the designer — DAG-JSON contract, headless executor, and codegen

**Question:** How is the breeding-program designer built so it is (a) promptly
available from the package, (b) able to Run a design and Copy an equivalent script
in R *and* Python, and (c) capable of "Apple-product" quality without risking the
CRAN tarball or the correctness of the genetics?

**Decision:** Build the designer as a **headless engine in the R package first**, and
ship the **canvas UI as a companion web app** later. The engine is the contract; the
UI only emits and displays it.

- **The contract is DAG-JSON**, defined before any UI (manuscript, Approach 3.3). A
  design is a list/JSON object `{version, seed, nodes:[{id, type, params, inputs}]}`.
  Node types each map to exactly one already-tested package function: `founders`
  (`as_population`), `cross`/`self`/`dh` (crossing core), `ssd`/`bulk` (schemes),
  `phenotype` (grammar `model` sub-spec), `select` (`select_ind`), `ocs`
  (`optimum_contribution` + `sample_parents`), `usefulness` (`cross_usefulness`),
  `pedigree`/`recurrent` (selecting schemes, additive loci frozen once so every
  generation scores on the same causal loci). So the designer adds **orchestration,
  not new genetics** — correctness rides on the tested primitives.
- **`run_design()`** (the Run button's backend) validates, topologically sorts, and
  executes a design with no server; **`validate_design()`** checks structure and
  acyclicity; **`design_breeding_program()`** is the exported launcher. One `seed`
  reproduces a whole run. (`R/designer.R`, `test-designer.R`.)
- **`design_script()`** (the Copy-script button's backend) transpiles a design into a
  **runnable R script** (verified: the generated R evaluates and reproduces
  `run_design()`'s output) and an equivalent **Python script** against the planned
  Python mirror API (Python is 0-based; header notes the package is forthcoming).
  Same topological walk as the executor, so script and Run agree. (`R/designer_codegen.R`.)
- **Delivery (recommended, logged for the UI pass):** the canvas is a **React Flow
  single-page app** hosted statically (e.g. GitHub Pages), **not** in the CRAN tarball
  (§5, 5 MB budget, bundled-JS scrutiny). Its **Run** button POSTs the DAG to a thin
  local `plumber` endpoint that `design_breeding_program()` launches — so genotypes
  never leave the user's machine — which calls `run_design()`. **Copy script** calls
  `design_script()`. Everything the UI produces must round-trip through
  `validate_design()`. This separation is what lets the UI pursue "Apple-quality"
  polish without touching correctness, which the theoretical-scrutiny gate (ROADMAP
  §8a) locks on the engine.
- **Dependency:** `jsonlite` added to **Suggests** (JSON parsing only; the native
  R-list interface needs nothing). No hard dependency added.

**Reaffirms:** DECISION-015/016 (the engine is the tested primitives), DECISION-006
(deterministic orchestration in R, meiosis in Rust), DECISION-011 (single package;
frontend weight stays out of the tarball).
**Date:** 2026-09-11

---

## DECISION-018: the designer becomes a separate product (`breedingDesigner`)

**Question:** DECISION-017 shipped the designer (DAG executor + code generator)
inside `simplePHENOTYPES`. The maintainer has since decided the designer should be
its own product — its own release cadence, UI, hosting, and G3 application-note
paper — relying on `simplePHENOTYPES` for the engine but open to other engines
(e.g. AlphaSimR for founder haplotypes with historical recombination). Where does
the designer code live?

**Decision:** Split the designer into a **separate package, `breedingDesigner`**
(repo `github.com/fernandes-lab/breeding_designer`, working copy at
`collaboration/Software/breeding_designer`), which `Imports: simplePHENOTYPES`.

- **Stays in simplePHENOTYPES (the engine):** phenotype grammar, PleioArch, the
  isqg-derived meiosis core, and the **selection engine** (`select_ind`, schemes,
  `g_matrix`, `optimum_contribution`, `sample_parents`, `cross_usefulness`). These
  are genetics and remain here, exported.
- **Moves to `breedingDesigner` (the designer layer):** the DAG-JSON contract, the
  headless executor (`run_design`, `validate_design`, `design_breeding_program`),
  the code generator (`design_script`), the engine adapter, the web canvas, and the
  execution back ends (webR + plumber). Files `R/designer.R`,
  `R/designer_codegen.R`, `inst/designer/plumber_api.R`, and
  `tests/testthat/test-designer.R` are copied there and re-namespaced to
  `simplePHENOTYPES::`.
- **Removal from simplePHENOTYPES happens only after `breedingDesigner` is green**
  (its `R CMD check` passes), so the engine is never broken mid-flight. Until then
  the files coexist; this is the one migration in flight.
- **Execution model (breedingDesigner, ADR-0002):** the canvas runs the R engine in
  the browser via **webR** (serverless; requires a pure-R engine path — the isqg
  parity reference and the pure-R stochastic core make this feasible), with a local
  **plumber** back end for large/Rust-native runs. This is the "click Run → results"
  goal, not just code generation.
- **Engine-agnostic (ADR-0003):** the executor calls engines through an adapter;
  default simplePHENOTYPES, pluggable AlphaSimR (historical-LD founders) and others.

**Consequences:** two citable units (engine + designer) with a citation flywheel;
the designer can pursue heavy UI/hosting without touching the CRAN engine; the
engine must keep exporting the functions the designer needs, and (for webR) a pure-R
execution path. Spec-driven development for the new product lives in
`breeding_designer/specs/` and `docs/adr/`.

**Supersedes:** DECISION-017's "designer ships inside simplePHENOTYPES." The DAG-JSON
contract, executor semantics, and code-generation behavior from DECISION-017 are
unchanged — only their home moves.
**Reaffirms:** DECISION-006 (RNG in R; enables the webR pure-R path), DECISION-015/016
(the selection engine stays in simplePHENOTYPES).
**Date:** 2026-09-11

---

## DECISION-019: breeding value = analytic transmissible average effect

**Question:** `select_ind(on = "bv")`, the OCS default merit, and the Smith–Hazel /
quadratic (QGSI) index merit all need a *breeding value*. What quantity, and how
computed? The first implementation used the sample least-squares projection of the
total genetic value onto the causal loci's gene content (the Fisher/NOIA
*statistical* additive value). The independent theory review showed this equals
the transmissible breeding value only under Hardy–Weinberg: on a deliberately
non-HWE locus it was ~64% off the Mendelian expected-offspring value, so OCS and
the indices would not optimize transmitted progeny merit after inbreeding or
selection.

**Decision:** use the **classical average-effect breeding value**
\(A_i = \sum_j \alpha_j (x_{ij} - 2p_j)\) with per-locus average effect of
substitution \(\alpha_j = a_j + d_j(q_j - p_j)\), and compute \(a_j\), \(d_j\)
**analytically from the simulation's own QTN effects** (rescaled by the same
`sqrt(prop)/sd` factor the realization applies), not by regression on the sample.
Because it is a simulation the per-locus effects are known exactly, so the breeding
value is exact and robust to **both** linkage disequilibrium (an F2/biparental
family is in strong LD yet returns the exact additive value — a per-locus *sample*
regression does not) **and** departures from HWE (the failure mode of the sample
projection). This is the merit `on = "bv"` (the OCS default) and the index methods
use.

**Rejected alternatives:** (a) the sample least-squares projection — LD-robust but
wrong off-HWE (the original defect); (b) a Hardy–Weinberg-weighted per-locus sample
regression — fixes HWE but is contaminated by LD, so it is only ~0.91-correlated
with the true additive value in the common F2 case; (c) a hybrid that also
reconstructs epistasis-induced marginal average effects — most complete but most
code/most review surface, deferred.

**Consequences / limitation:** the average effects come from **additive and
dominance** layers only. An **epistasis** layer has no single per-locus \(a\)/\(d\),
so its induced additive average effects cannot be reconstructed; rather than return
a breeding value that silently omits them (which misranks transmissible merit —
e.g. it is identically zero for a purely epistatic model whose progeny merit
varies), `.breeding_value_matrix()` **errors when any epistasis layer is present**,
exactly as it does for `architecture = "complex"`. Pass a custom `on` criterion
(externally supplied breeding values / merit) for epistatic models. (The fuller
hybrid — reconstructing epistasis-induced marginal average effects — was
deferred.) The index methods (`index`, `quadratic_index`) score on all traits'
true breeding values and therefore ignore `on` (they warn if one is supplied);
there is no external per-trait GEBV hook.

**Related hardening (same review round):** `method = "combined"` (Lush index) is
restricted to `on = "pheno"` (its covariance derivation assumes phenotypic
records); `g_matrix()` gains an optional fixed `base_freq` (the default remains the
current-sample "moving base", now documented as such); `optimum_contribution()`
validates that `G` is finite, symmetric and positive semi-definite before the
Frank–Wolfe optimizer (which assumes non-negative coancestry curvature).

**Rubric:** `docs/THEORY_REVIEW.md` S5 is reconciled with this decision — `"bv"` is
blessed as an idealized *true* simulated criterion parallel to the already-permitted
`"gv"` (both explicitly labeled true; the vector/function hook remains the
estimated/predicted extension point and must not present a value as true BV).

**Terminology:** the quantity is the *one-generation random-mating transmitting
ability* α = a + d(q−p), not the current-population Fisher/NOIA *statistical*
additive value (the sample-genotype-frequency regression slope), which it equals
only under HWE. Documentation is worded accordingly.

**Amendment 2026-09-29:** the breeding value above is the merit of the QGSI and the
Smith-Hazel index. Smith-Hazel uses `G = Cov(A)` of these breeding values; that
equals Hazel's Cov(phenotype, breeding value) only when Cov(A, D) = 0 (HWE /
random mating) -- an assumption documented in `select_ind()` and DECISION-015, not
a change of estimator.

**Reaffirms:** DECISION-015/016 (selection engine + modern methods stay here);
DECISION-006 (deterministic, in R).
**Date:** 2026-09-13

---

## DECISION-020: orthogonal genotypic model via `additive(orthogonal = TRUE)`

**Question:** the variance-partition grammar codes additive as -1/0/1 dosage and
dominance as a heterozygote indicator, each scaled independently to a chosen
`prop`. That is not Fisher's orthogonal decomposition: at non-0.5 allele
frequencies the additive and dominance components are correlated, realized H2 only
approximates the sum of props, and a "degree of dominance" scalar is washed out by
per-component scaling. How to offer the theoretically correct additive/dominance
partition?

**Decision:** add an orthogonal genotypic model as a mode of the additive layer,
`additive(orthogonal = TRUE, a =, d =)` (API chosen over a separate
`genotypic_model()` builder). Each locus's genotypic value is built from a
per-locus additive effect `a` and dominance deviation `d` (`-a`/`+d`/`+a` for gene
content 0/1/2), and the **whole** genotypic value is scaled to `prop` of the
phenotypic variance. The additive (breeding-value) component uses the average
effect \(\alpha_j = a_j + d_j(1 - 2 p_j)\) and the dominance deviation is the
realized residual \(D = g - A\). The additive and dominance components are
orthogonal (\(Cov(A, D) = 0\)) in expectation under **random mating**, which puts
each locus in Hardy–Weinberg proportions. This does *not* require linkage
equilibrium — between-locus LD is compatible with orthogonality in this
no-epistasis model — but per-locus HWE alone is not sufficient under arbitrary
nonrandom multilocus genotype association (selection, structure), which can
correlate one locus's gene content with another's heterozygosity. \(\alpha\) is
the one-generation transmitting-ability average effect (Falconer 1985). The name
"orthogonal" reflects that random-mating regime; it is not a claim of
orthogonality in every population. A finite or structured sample carries a
(usually small) non-zero \(Cov(A, D)\) that breeders customarily ignore — this
implementation does not, and normalizing it away would misreport the shares. The additive/dominance variances **emerge** from (a, d, p) and are
reported as the *realized* fractions of \(Var(g)\): an `additive` row
\(Var(A)/Var(g)\), a `dominance` row \(Var(D)/Var(g)\), and an `add_dom_cov` row
\(2\,Cov(A, D)/Var(g)\) so the three close to the layer `prop` (the covariance row
is \(\approx 0\) for a large random-mating / F2 sample and non-zero otherwise) —
rather than being set by separate layer props.
Because the whole value is scaled together (not per component), the **degree of
dominance** `d/abs(a)` is meaningful (magnitude 1 = complete, > 1 = overdominance;
`abs` because `a` may be negative; it drives the emergent Va:Vd ratio). The breeding
value (DECISION-019) picks up the dominance-induced average effect automatically.

**Scope / limits:** total genetic fraction (H2) is still the layer's `prop` (the
grammar's budget contract is kept); what emerges is the additive-vs-dominance
split within it. The mode fixes per-locus effects, so it is incompatible with
`vary_qtn`, with a separate `dominance()` layer (it already carries dominance), and
— like `qtn =` — with the correlation-controlling `"pleiotropy"` (multi-trait) and
`"ld"` architectures. A dominance deviation needs heterozygotes (errors on hetless
loci, as `dominance()` does).

**Reaffirms:** DECISION-019 (transmissible average-effect breeding value);
DECISION-006 (deterministic, in R).
**Date:** 2026-09-13

---

## DECISION-021: fixed-scale phenotype accessor `phenotype_value()`

**Context:** cross-generation `on = "pheno"` selection needs a phenotype whose
genetic scale and residual variance are *frozen*, so the parametric heritability
`Var(g)/(Var(g)+var_e)` declines as selection exhausts genetic variance. `simulate_phenotype()` / `genetic_values()`
re-scale the genetic layer to `prop` on every population, holding the genetic share
constant — correct for `on = "gv"` (rescaling is monotone, ranking unchanged) but
optimistic for `on = "pheno"` (additive-only accuracy stays near sqrt(h2) under re-standardization, 0.695-0.714 vs 0.707 in the review run, and declines only when dominance / non-orthogonal components are present, 0.564 -> 0.538 -> 0.477). `additive_value()`
(DECISION-020 companion) already froze the genetic value; the residual was missing.

**Decision:** add `phenotype_value(x, qtn, effect, h2 = NULL, var_e = NULL,
ref = NULL, seed = NULL)`. Phenotype = the fixed additive value
(`additive_value()`, −1/0/1 dosage × effect, no per-population rescale) + a normal
residual on a **fixed** variance. Exactly one of `h2` or `var_e` sets that variance:
`var_e` directly (the robust cross-generation choice — compute once, reuse), or `h2`
converted via `var_e = Var(g_ref)(1 − h2)/h2` from a reference population's genetic
variance (`ref`, default `x`). The residual is an **independent** draw
`e ~ N(0, var_e)` (a genuine normal, *not* sample-standardized to exactly `var_e` —
so it stays normal, works at `n = 1`, and has `Cov(g, e) ≈ 0` in expectation),
seed-restoring via `.Random.seed_safe()`/`.restore_seed()`; `seed` is checked with
`.validate_seed()` and an overflowing `var_e` (extreme `h2`) errors. Returns a named
phenotype vector with `var_e` and `genetic_value` attributes. Bumps the dev version
(1.4.0-9002) so downstream (breedingDesigner) can pin the minimum and drop its
`.bv_index` fallback.

**Scope / limits:** RNG stays in R (DECISION-006). `h2` with the default `ref = x`
re-derives `var_e` per population and so is *not* frozen across generations — the
docs flag this; cross-generation callers pass `var_e` or `ref = <base>`.

**Amendment 2026-09-29:** the multi-generation schemes `pedigree()` /
`recurrent_selection()` take their phenotype from a `simulate_phenotype()`
callback, which re-standardises the genetic layer to `prop` in every population;
their `history` S and i are therefore on a per-generation scale; additive-only selection
accuracy stays near sqrt(h2) under re-standardization (0.695-0.714 vs 0.707 in the
review run) and declines when dominance / non-orthogonal components are present
(0.564 -> 0.538 -> 0.477). The fixed-scale response of this decision is reached by
ranking on `phenotype_value()` through `on = function(s) phenotype_value(s$geno,
qtn, effect, var_e = ve)` (documented, with an example, in `pedigree()`).

**Reaffirms:** DECISION-020 (`additive_value()` fixed scale); DECISION-006/019.
**Date:** 2026-09-14

---

## DECISION-022: transcriptome simulation and the genome → transcriptome → phenotype model

**Question:** DECISION-005 deferred expression-based simulation to v3 ("different
input semantics; would muddy the marker redesign"). The v2 grammar, the crossing
engine, and the fixed-scale accessors are now stable. How do we add gene-expression
simulation and expression-mediated phenotypes — and how much of it now?

**Decision:** Add expression as a **downstream genetic process** and expose three
phenotype **bases** through one interface, reusing the existing machinery.

- **Expression generator — hybrid latent-factor eQTL, not a mechanistic GRN.**
  Normalized (~Gaussian) expression `E_g = mu_g + G_g + R_g` per gene:
  - `G_g` = a per-gene **cis** score (markers in a physical-bp window; the dosage ×
    effect model, localized) plus a **trans** score mediated by `Q << T` latent
    regulatory factors (each a genetic hub with a few QTL and sparse gene
    loadings). cis and trans covary (LD/structure), so the **combined** genetic
    score is scaled to the per-gene target `h2_g`, and `Cov(cis,trans)` is reported
    as its own budget row — the **same pattern as the orthogonal model's
    `add_dom_cov` row** (DECISION-020).
  - `R_g` = shared **non-genetic** module factors (same loadings) + gene noise, so
    co-expression is **decoupled from heritability** (genes co-express even at
    `h2_g = 0`). `omega_g` = cis fraction; `kappa_g` = residual module fraction;
    `h2_g` = total genetic variance — three non-competing controls.
  - **Fixed reference calibration:** all means/frequencies/sd's are frozen on the
    founder population and never re-estimated on descendants — the fixed-scale
    principle of `additive_value()`/`phenotype_value()` (DECISION-020/021), which
    is what makes a simulated transcriptome usable across a crossing/selection
    pipeline.

- **No annotation, no reference data required.** Defaults come from a named
  calibration **profile** (`generic_bulk`), presented as a benchmarking compromise
  **not** biological constants. Gene coordinates are synthesized on physical bp
  when no annotation is given (cis is never silently a cM window). A user may
  optionally supply an annotation, or supply expression to **mimic** (calibrate the
  generator's moments + low-rank co-expression, and h² if genotypes are paired,
  then regenerate) — mimic calibrates the **expression model only**, never the
  phenotype slopes. When paired genotypes are supplied, per-gene `h2` is calibrated
  by a **GREML-style** estimator (REML on a GRM). An **example dataset** (a
  synthetic gene annotation and an example expression matrix for the bundled SNP55K
  panel) ships so every basis runs out-of-the-box.

- **Phenotype bases inferred from inputs, one `transcriptome()` layer.** The basis
  follows from which inputs `simulate_phenotype()` receives (no explicit `basis=`);
  markers-only is the default. `expression =` is **real/observed** expression;
  `transcriptome =` (`TRUE` or a `transcriptome_sim`) is **genome-derived**
  expression. The continuous-predictor layer is `transcriptome()`, held to the
  **same variance-budget / realized-h² / effect-table / reproducibility guarantees
  as the genome path**. Four modes:
  - `geno` (default) → markers only, `y = Zδ + η` (unchanged);
  - `expression = E` (no `geno`) → transcriptome alone, `y = Ẽs + η`;
  - `geno` + `transcriptome =` → **derived** G→E→Y, `y = Z(δ + Bs) + Rs + η`,
    separating mediated-genetic (`Bs`), direct-genetic (`δ`), and env-mediated
    (`Rs`); the total genetic score is scaled jointly to the target h² and the
    mediated/direct split reported with its covariance;
  - `geno` + `expression = E` (**both real**) → `y = Zδ + Ẽs + η`, a direct marker
    effect plus an observed-expression effect with **no simulated mediation**; their
    realized covariance (expression is biologically downstream of the genome) is
    reported, not asserted.

- **Scope.** The normalized generator, the transcriptome and
  genome→transcriptome phenotype bases, `mimic` calibration, the RNA-seq **count**
  observation layer (`observe_counts()`, NB), the genotype-free generator and
  additive-by-additive epistatic expression (`simulate_transcriptome(epistasis =)`)
  are all implemented. Still explicit non-scope: directed regulatory networks and
  tissue specificity.

- **Reuse & boundary.** Reuse genotype ingestion, QTN sampling, the
  variance-realization ("scale-to-target") step, the map, and the RNG-in-R
  discipline. All stochastic draws stay in R (DECISION-006); the model is sparse /
  low-rank and scales to `T ~ 10^4` genes with sparse products; **no Rust** until
  profiling proves a deterministic bottleneck.

**Rationale:** Gene expression is itself a set of genetically controlled
quantitative traits, so the genome→transcriptome arrow is the existing additive
model applied per gene; only the transcriptome→phenotype arrow is new (continuous
predictors → generated slopes). A latent-factor model gives cis + trans + modules
from one low-rank construct that scales to whole transcriptomes, where a mechanistic
GRN would need wiring no annotation can justify and a dense trait-by-trait
correlation (PleioArch) cannot scale. Two of the design's load-bearing pieces — the
covariance budget row and fixed reference calibration — are patterns the package
already uses, which keeps the addition coherent with the codebase. Known ground
truth (which variant is a cis/trans-eQTL, which gene mediates) makes the generator a
benchmark for eQTL/TWAS/mediation methods.

**Supersedes:** DECISION-005 ("expression-based simulation deferred to v3") — brought
forward, with expression modeled as downstream of the genome rather than as a
separate transcript input format.
**Reaffirms:** DECISION-006 (RNG in R; surgical Rust), DECISION-020/021 (fixed-scale
values; realized variance reporting with a covariance row), DECISION-009 (realized,
not asserted, variance partition).
**Date:** drafted 2026-09-16; ratified 2026-10-02, when the full feature set of
`SPEC-transcriptome.md` was built and reviewed (the draft lived in
`DECISION-022-transcriptome-DRAFT.md` until then).

---

## DECISION-023: genetic-correlation control extends to dominance and epistasis

**Context:** the 2026-09-17 audit (grammar P1/P4, effects-arch O2/X1) found that
`cor` controlled only the **additive** layer. Under `architecture = "pleiotropy"`
`dominance()` / `epistasis()` reused shared loci and gave every trait one identical
effect series, so those components realized a correlation of ~ +1 whatever `cor`
was — target `cor = 0` realized **1.0** for pleiotropic epistasis and **0.45** for
the one-call `model = "AD"`; `cor = 0.2` with an epistatic layer realized ~0.60. A
fixed `qtn =` on a non-additive layer was also accepted under "pleiotropy"/"ld"
(replicated to every trait → correlation 1.0, and under "ld" the same set causal
for both traits). The maintainer chose to *extend* control (over document-only or
reject). Scoped in `docs/SPEC-nonadditive-correlation.md`.

**Decision:** the PleioArch covariance construction (DECISION-007/013) is applied
per mean-effect component. `.pleio_nonadditive_draw()` partitions a dominance or
epistasis layer's units (a locus / an interacting set) by the same `pi` into shared
units (common to every trait; effects jointly MVN(0, Σ/n), Σᵢᵢ = πᵢVᵢ,
Σᵢⱼ = corᵢⱼ√(VᵢVⱼ)) and trait-specific units (independent), with the same `cor`,
feasibility check and zero-variance / single-shared-unit guards. Each unit's effect
is divided by the **realized** standard deviation of its design column (heterozygote
indicator; centered product via `.epi_unit_column()`, now the single source of truth
for realization too), so every unit contributes equal design variance on the
simulated sample; then, over effect draws, E[Cov(c₁,c₂)] = Σ₁₂ and E[Var(c_t)] = V_t
(the shared units give Σ_tt = π_tV_t and the independent trait-specific units the
remaining (1−π_t)V_t; only the covariance comes from the shared units; corrected
2026-09-29, this line used to read E[Var(c_t)] = Σ_tt, which holds only for π_t = 1;
any LD: effects of different units are independent), so the component targets `cor`
in the same sense as the additive layer (see Scope for what that means for the
realized correlation). Constant design columns (e.g. hetless loci) get effect 0 and are
left out of the allocation, so Σ is split over the informative shared units and
(1−π)V over the informative specific ones (counting dead units would shift the
shared:specific ratio and bias the correlation); warnings fire when degeneracy
leaves no informative shared units (correlation cannot be carried) or no
informative specific units for a trait (targeted correlation inflated in magnitude,
i.e. further from 0 than a nonzero `cor`; a zero target stays 0 and does not warn). Independently
drawn components add no cross-component covariance in expectation, so the
**total** genetic correlation targets `cor · Σ_c√(V_c1V_c2) / √(ΣV_c1·ΣV_c2)` (a
ratio of expected moments): exactly `cor` when the layers' per-trait `prop` profiles are proportional or `cor = 0` (proportional always for
scalar `prop` and the one-call models), otherwise attenuated toward 0 (Cauchy–
Schwarz) — possibly below what any per-layer correlation could restore — and
`.pleio_total_cor_check()` warns with that large-sample target (a ratio of
expected moments, not the finite-unit mean correlation). The total is not
re-targeted (layers arrive one at a time). With `same_as_add = TRUE` dominance
reuses the additive layer's shared/trait-specific loci (shared = present for every
trait). The additive draw is unchanged. Fixed
`qtn =` on `dominance()`/`epistasis()` is rejected under "pleiotropy" and "ld", and
`effect =` / non-default `dist` under "pleiotropy" (as `additive()` already did).

**Why the realized sd, not MAF:** the additive `1/√(2·MAF·(1−MAF))` is the HWE
expectation of the same normalizer. Non-additive design variance is not a clean
function of MAF — a near-inbred panel has far fewer heterozygotes than 2pq — so the
realized sd is what equalizes unit contributions on the sample actually simulated.

**"ld" (SPEC §5.3 restriction taken):** the ld contract (DECISION-014) is that each
trait's causal loci are distinct markers in LD within the r² window, with no SNP
causal for both traits in the mean-effect (genetic-value) layers. It holds within
one additive draw only, so a **second `additive()` layer under "ld" now errors**
(round 7: two layers made six SNPs causal for both traits). A genome-derived
`transcriptome()` layer adds a shared genetic cause (shared genes / eQTL) outside
this design, so it warns under "ld" as under "pleiotropy" (round 9) — even when only
one trait receives it, since its eQTL markers can be the other trait's causal loci
(round 11). Dominance satisfies it only by reusing the additive layer's
linked loci (`same_as_add = TRUE`, the default); `dominance(same_as_add = FALSE)`
under "ld" now errors, because a fresh draw knows nothing of the additive loci and
can make a marker causal for both traits across layers. `epistasis()` under "ld"
now errors: epistatic sets have no linked-distinct construction — drawn anywhere
they are unlinked (outside the r² window) and can reuse another layer's causal
marker for the other trait. (Earlier drafts documented ld epistasis as
trait-specific and, in round 2, made its sets disjoint within the layer; round 3
showed both still break the contract, so the restriction fallback the SPEC had
anticipated was taken.) No test or vignette used either form.

**Verification:** 18/18 acceptance cells pass (outbred HWE panel, 30 seeds,
cor ∈ {−0.5, 0, 0.5} × {AD reused loci, AD fresh loci, AE a×a, AE a×d, A+D+E with
π = 0.7, one-call AD}): mean realized total correlation within max(0.06, 3·se) of
target, every component likewise (`dev/accept-nonadditive-cor.R`). A 150-seed independent re-check of the reused-loci
dominance case shows no bias (z = 0.62 at cor = 0, 0.16 at 0.5). Additive
pleiotropy (incl. major QTN, vary_qtn, one-call A), independent A/D/E and ld models
are bit-identical to the previous release. `tests/testthat/test-nonadditive-cor.R`
fails on the pre-change code and passes now. **Independent theory review
(THEORY_REVIEW rubric) BLOCKed the first version** on four points, all confirmed and
fixed before commit: (O1) the claim that per-component control always gives a total
of `cor` — false for non-proportional per-trait `prop` (0.14 at target 0.5 in the
reviewer's counterexample); (O2) constant units still counted in the variance
allocation, biasing the correlation (0.306 instead of 0.5); (O3) the single-shared-
unit warning claimed ±1 even with trait-specific units; (O4) infeasibility errors
did not name the component, although the SPEC's acceptance criteria required it.
Each has a regression test. **Round 2 BLOCKed** on three further points, also
confirmed and fixed: (O1) "realizes `cor` in expectation" is imprecise — the
raw draw's covariance is exact in expectation but the realized correlation is a random ratio,
attenuated toward 0 with few units (2 shared units: mean 0.41 at target 0.5); the
wording is now precise everywhere, with independently verified numbers; (O2) under
"ld", independent per-trait epistatic draws could make a marker (or a whole set)
causal for both traits — first fixed by drawing the sets jointly, then superseded in
round 3 by the "ld" restriction; (O3) the
attenuation example's multiplier was misstated as 0.14 × `cor` (it is 0.28 × `cor`,
i.e. 0.14 at `cor = 0.5`; the warning's computed value was always right).
**Round 3 BLOCKed** on four more, confirmed and fixed: (O1) the non-proportional
warning called the moment ratio the "expected" total correlation — it is the
large-sample target (reviewer: target 0.40, realized mean 0.29 with two shared
loci), now worded so; (O2, O3) ld epistasis sets, even disjoint within the layer,
were outside the r² window and could reuse another layer's causal marker for the
other trait — resolved by the "ld" restriction above; (O4) the SPEC still described
the pre-change mechanism as current and the test header overstated the total-
correlation claim. The attenuation multipliers are also now stated at `cor = 0.5`
(they grow toward 1 as |cor| → 1).
**Round 4 BLOCKed** on three more, confirmed and fixed: (O1) "the expected
covariance is exact" held only for the raw effect draw, not the component after
rescaling to `prop` (reviewer: realized 0.156 vs 0.2 with two shared loci) — the
wording now distinguishes them; (O2) convergence was claimed without its condition —
the docs now state that it needs approximate linkage equilibrium among causal loci
and that the attenuation depends on `cor`, π and the designs; (O3) with reused additive loci, hetless units could leave exactly one
informative shared unit without the single-unit warning — `.pleio_unit_effects()`
now warns in that case.
**Round 5 BLOCKed** on three more, confirmed and fixed: (O1) the single-shared-unit
guards skipped `cor = 0`, which one shared unit cannot realize either (reviewer:
+1.000 with no warning) — both guards now fire for any target strictly inside
(−1, 1); (O2) several passages (and canonical SPEC §4/§13) still said "realizes",
"equals", "exact in expectation" or "expected value", and the complete-LD case was
misdescribed as a limit — under complete LD every realized r is ±1 however many
units, and 2·asin(cor)/π is only the ensemble mean (reviewer: 300 seeds, all ±1,
mean 0.360); all now say "targets" / "converges (given approximate LE)"; (O3) the
constant-specific-unit warning said "inflated above `cor`", wrong in direction for
negative `cor` (reviewer: −0.83 at target −0.4) — it now says inflated in magnitude.
**Round 6 BLOCKed** on two documentation points, confirmed and fixed: (O1) DECISION-007,
-010 and -013 still said the realized correlation equals `cor` in expectation (reviewer:
mean 0.407 at target 0.5 with two shared units) — each now carries a dated correction;
(O2) the public `cor` docs stated convergence without the major-QTN exception — now
stated there, in the PleioArch header and in THEORY_REVIEW P1.
**Round 7 BLOCKed** on four more, confirmed and fixed: (O1) a `transcriptome()` layer
on a genome-derived source adds its genome-mediated signal to `genetic_values()`
outside `cor` (reviewer: `cor = 0`, total 0.64, no warning) — it now warns under
"pleiotropy" and the docs say the total is then not targeted; (O2) several
`additive()` layers under "ld" could make a SNP causal for both traits — a second
layer is now rejected; (O3) the all-constant-shared-unit warning promised a
correlation "near 0", false under LD among the specific units — it now says the
correlation is no longer controlled by `cor`; (O4) README, the complete reference,
the v2 vignette and ARCHITECTURE.md still said "in expectation" / "exact" — corrected
(the round-5 note's "all now say" was premature).
**Round 8 BLOCKed** on four wording/advice points, confirmed and fixed: (O1) two
warnings stated the direction of a single finite draw ("is inflated", "is smaller")
where only the target or the average moves (reviewer: 0.25 at target 0.4; 0.82 at
target 0.48) — they now describe the target and note the scatter; (O2) the partition
helper's doc said one shared unit gives ±1 without the no-trait-specific-units
condition; (O3) "exactly `cor` when proportional" omitted `cor = 0`, and a test comment
said "realizes"; (O4) the marker-shortage error advised lowering `pi`, which raises the
number of distinct markers needed (nt·n − (nt−1)·pleio_n) — it now says raise `pi`; the
additive draw had no such guard at all and failed inside `sample()`, so it now gets the
same message.
**Round 9 BLOCKed** on one point, confirmed and fixed: (O1) the derived-transcriptome
warning covered "pleiotropy" only, although under "ld" shared genes likewise add a
genetic cause outside the distinct-linked-loci design (reviewer: additive correlation
0.55 → total 0.94, no warning) — it now warns under "ld" too.
**Round 10 BLOCKed** on two, confirmed and fixed: (O1) the multi-trait feasibility check
tested the variance-scaled Σ with a tolerance set by its largest eigenvalue, so a block
of tiny-variance traits could pass while infeasible (reviewer: `prop = c(1e-16, 1e-16,
0.5)`, `pi = c(0.1, 0.1, 1)`, cor₁₂ = 0.5 > √(0.1·0.1) accepted) — it now tests the
variance-free M (Mᵢᵢ = πᵢ, Mᵢⱼ = corᵢⱼ; Σ = D^½ M D^½ is PSD exactly when M is);
(O2) `transcriptome(prop = 0)` emitted the correlation warning although it adds
nothing — the warning now requires a contributing layer.
**Round 11 BLOCKed** on two, confirmed and fixed: (O1) under "ld" that warning required
two contributing traits, but a one-trait layer can reuse the other trait's causal
marker as an eQTL (reviewer: total −0.48 → −0.32, no warning) — any contribution now
warns; (O2) the lost-trait-specific-units warning claimed inflation at `cor = 0`, where
the target stays 0 (reviewer: mean 0.002) — it now fires only for a nonzero `cor`.
**Round 12 BLOCKed** on two, confirmed and fixed: (O1) the single-shared-unit warning
chose "exactly ±1" vs "one noisy draw" globally; it is per pair — two traits with no
trait-specific variance (none drawn, or π = 1) are exactly ±1 even when other traits
have specific units (reviewer: π = (1, 1, 0.1) gave −1.000 but the "noisy" text) — the
warning now lists pairs under each case; (O2) the non-proportional-total warning used an
absolute 0.01 threshold, silent at `cor = 0.01` despite a 72% attenuation — it is now
relative (more than 1% of `cor`).
**Round 13 BLOCKed** on one, confirmed and fixed: (O1) with `prop = 0` for a trait the
single-unit warning still claimed "exactly ±1" for its pairs, whose correlation is
undefined (NA; `.pleio_check_zero_var()` already says so) — zero-variance pairs are
now skipped there.
**Round 14 BLOCKed** on three documentation/evidence points, confirmed and fixed: (O1)
the v2 vignette still said the mean "sits on the target at every size" with SD
falling as 1/√n_qtn (reviewer, bundled maize: 0.593 / 0.597 / 0.575, SD 0.146 / 0.077 /
0.080 at 20 / 100 / 400 QTNs) — reworded; (O2) SPEC §6 still said "the current code
fails this"; (O3) SPEC §6 claimed the ≥30-seed criteria were added to the `evals/`
golden set — they were run as a one-off script, now committed as
`dev/accept-nonadditive-cor.R`, and the SPEC says so.
**Round 15 BLOCKed** on two, confirmed and fixed: (O1) the all-constant-shared-unit
warning said the correlation is no longer controlled even at `cor = 0`, where the
independent specific effects still target 0 (reviewer: mean 0.0007) — at `cor = 0` it
now reports only the lost shared variance share; (O2) SPEC §2 said "MAF scaling" is
applied per component, whereas dominance/epistasis divide by the realized design-column
sd — corrected.
**Round 16 BLOCKed** on four, confirmed and fixed: (O1) this record named the review
tool, which AGENTS.md forbids in files — removed (also from DECISION-019's older
text); (O2) `.pleio_total_cor_check()` assumed every layer targets `cor`, false once
a component loses its shared or specific units (reviewer: dead shared dominance units,
reported total 0.320, correct 0.160) — each pleiotropic non-additive layer now records
its effective target (`layer$target_cor`) and the check sums those; (O3) a fresh draw
with one shared unit kept the partition's nominal "noisy draw" text when every
specific unit proved constant (realized −1.000) — the informative-unit recheck now runs
then and reports the exact ±1; (O4) the transcriptome warning fired for a derived
source at h² = 0, whose genetic expression is constant — it now requires genetic
variance in the layer's causal genes.
**Round 17 BLOCKed** on one, confirmed and fixed: (O1) that gene check looked only at
the canonical draw, so under `vary_qtn` a later replication's genome-mediated genes
went unwarned — it now checks every replication's genes. (The complete-LD
2·asin(cor)/π value is now given with its derivation.)
**Round 18 BLOCKed** on three, confirmed and fixed: (O1) the gene check pooled the genes
of every trait, including `prop = 0` traits that receive no signal — it now uses only
contributing traits; (O2) under `vary_qtn` only the canonical replication's effective
target was kept (reviewer: replication 5 targeted 0.231, reported 0.366) — every
replication's target is now stored (`target_cor_reps`) and checked, and a replication
that differs is reported by number; (O3) `dev/accept-nonadditive-cor.R` printed the
per-component means but gated only on the total — it now gates on both (still 18/18).
**Round 19 BLOCKed** on two, confirmed and fixed: (O1) the dominance/epistasis hetless
checks looked only at the canonical draw, so a `vary_qtn` replication could silently
lose units (reviewer: live counts 4, 3, 4, 4, 2, no warning) — they now check every
replication; (O2) the transcriptome warning tested gene-level genetic variance, not the
layer's realized genetic signal, so zero slopes still warned — it now tests the
realized genome-mediated component (`.tx_raw(..., "genetic")`) per contributing trait
and replication.
**Round 20 BLOCKed** on two, confirmed and fixed: (O1) convergence was stated as the
shared units grow, but r is a sample correlation over n individuals, so with n fixed
its spread levels off (SD 0.18 at n = 20 even with 5 000 QTNs) — every convergence
statement now requires the individuals to grow too; (O2) SPEC §6 item 6 (written
before implementation) named the review tool and pre-announced its agreement —
reworded.
**Round 21 BLOCKed** on three, confirmed and fixed: (O1) the TODO entry still said
"converges as units grow" without the individuals; (O2) the complete reference said
the effect covariance "matches `cor`" — it is `cor·√(V₁V₂)`, and it is the correlation
that targets `cor`; (O3) the total warning printed `%.3f`, so at `cor = 1e-4` a 72%
attenuation read "0.000 (cor 0.000)" — tiny values now print significant digits.

**Scope / limits:** what holds exactly in expectation (over effect draws, whatever
the LD) is the **raw effect draw's** cross-trait covariance cor·√(VᵢVⱼ) and variances
Vᵢ, so `cor` is the ratio of expected moments — the target. Each layer is then
rescaled to its `prop`, so what is realized is a correlation r (covariance
r·√(VᵢVⱼ)), a random ratio, and E[r] ≠ `cor` in general. r is a sample correlation
over the n individuals, so it converges to `cor` only as the shared units **and n**
grow — with n fixed, more units bring r to the finite-n sampling distribution of a
correlation, whose spread does not vanish (reviewer: n = 20, 5 000 shared QTNs, SD
0.18 at `cor = 0.5`) — and **only if the units' design columns are close to
uncorrelated** — linkage equilibrium among causal loci, disjoint loci — which the
draw does not enforce (it samples distinct markers) — and only if no unit keeps a
non-vanishing share of the variance: additive major QTNs (`n_pleio_major` /
`prop_var_major`) keep theirs however many QTNs there are, so r then does not
converge (reviewer: SD 0.31 / 0.28 / 0.28 at 20 / 100 / 400 shared QTNs with one
major unit holding half the shared variance; cf. DECISION-013's 0.25 / 0.32 / 0.23).
Under strong LD it need not converge at all either: with complete LD (identical design columns) every realized r is
±1 however many units — the sign of the product of the two summed effects, which are
bivariate normal with correlation `cor` — so its ensemble mean is
P(XY>0) − P(XY<0) = 2·asin(cor)/π, from the bivariate-normal orthant probability
P(XY>0) = 1/2 + asin(cor)/π (reviewer: 100 linked
units, 300 seeds, mean 0.360 at `cor = 0.5`). Under LE, with few units r is
attenuated toward 0 on average by an amount
that depends on `cor`, π and the designs: 0.68 / 0.82 / 0.92 / 0.97 / 0.98 / 0.99 ×
`cor` at 1 / 2 / 5 / 10 / 20 / 60 shared units for `cor = 0.5`, π = 1 and
orthogonal unit-variance designs (40 000-draw check); with trait-specific units
(π = 0.5) the reviewer measured 0.947 × `cor` at two shared units. This holds
equally for the additive PleioArch layer: DECISION-007/013's "realizes `cor` in
expectation" should be read this way (SPEC §4/§13 now say so). A single shared unit warns: with no
trait-specific units the realized correlation is exactly ±1; with them, the whole
covariance rests on one effect pair, so it is one noisy draw (any target strictly
inside (−1, 1), `cor = 0` included). `n_pleio_major` / `prop_var_major` shape
the additive layer only. The derivation assumes units uncorrelated with one another
(linkage equilibrium, disjoint loci); at a shared locus the additive and dominance
design columns correlate when p ≠ 0.5, but their effects are drawn independently,
so the cross-component covariance is zero only in expectation. Hetless units
contribute nothing (existing warnings apply). The PleioArch additive construction is
Prado et al. (in preparation); this non-additive extension is this package's design,
not a published method (Rule #7).

**Refines:** DECISION-007/013 (PleioArch, now all mean-effect layers); DECISION-014
("ld" non-additive behavior). **Date:** 2026-09-25

---

## DECISION-024: a Population records its pedigree

**Context:** `docs/SPEC-block3b.md` F1. Parentage survived only as a batch `origin`
string and as id prefixes that `c()` may rename, so full-sib / half-sib families,
breed composition and a relationship matrix could not be recovered. Block 3B items
(combining ability, progeny testing, BLUP, crossbreeding) all need them.

**Decision:** every `Population` carries `keys` (one per individual) and a `pedigree`
frame (`key`, `id`, `mother`, `father`, `generation`, `pool`, `design`), maintained by
`as_population()` (founders; new optional `pool =` label), `.mate()` (every
`cross()` / `selfcross()` / `double_haploid()` progeny), `[` (keeps the ancestors of
the individuals kept) and `c()` (union by key). **Links are keys, not display ids**
(maintainer decision on SPEC D2, replacing the SPEC's "c() errors on collisions",
which would have broken pooling crosses whose progeny all start at `prog_1`): `c()`
still renames colliding display ids and the links are untouched. Keys are
deterministic — a founder's is a hash of its pool, id and haplotypes, a progeny's a
hash of the mating (design, parent keys, the RNG generator's state words before and
after the draws — not its kind header, which `RNGversion()` changes; the state after
separates matings when no state existed before — and the drawn meiosis events, doubles
by their exact bits, so keys do not depend on the numeric locale)
plus its index — so seeded runs reproduce them, an exact re-run (same parents, same
seed) reproduces the same individuals under the same keys, and matings whose draws
merely coincide (e.g. on a 0 cM map) stay distinct. The hash is FNV-1a-128 in the
Rust core (`stable_hash_core()`, a deterministic kernel) over a canonical text
encoding, so for the same parents and random stream the keys do not change with R,
package or dependency versions (the first draft used `rlang::hash()`, whose values
rlang 1.3.0 changed). Founder ids must be unique; a missing pool (any `NA`/`NaN` label, normalized to one) is
encoded apart from any pool string; the bookkeeping draws nothing, so every genotype and the RNG
stream are unchanged (isqg parity and golden snapshots bit-identical). A `Population`
without a pedigree (an object from an older version) is treated as founders.
Accessors: `parentage(x, ancestors = FALSE)` (the SPEC's `pedigree()` name is taken
by the pedigree-selection scheme; it returns `key` / `mother_key` / `father_key`
beside the display ids, which need not be unique) and `families(x, by = c("full_sib",
"maternal_half_sib", "paternal_half_sib", "selfed"))`, a factor for
`select_ind(family =)`. The canonical encoding is length-prefixed (no value can
imitate a separator) and in UTF-8 (ids R treats as identical hash identically, whatever
their declared encoding). `mating_design()`, `mate()` and `cross()` treat the same individual (same
key) in two positions as a self (`cross(x, x)` draws what `selfcross(x)` draws and is
recorded and keyed as one), whatever its display ids; `"nested"` assigns mothers
by bipartite matching, so a feasible no-self assignment is always found.

**Date:** 2026-09-27

---

## DECISION-025: mating plans — `mating_design()` and `mate()`

**Context:** `docs/SPEC-block3b.md` F2 / item 6. The engine (`.intermate()`,
`.make_family()`) and breedingDesigner (`.rrs_intermate()`, `.rrs_hybrid_pop()`) each
hand-rolled "loop over pairs, `cross()`, relabel, `c()`"; BD SPEC-0007 outputs a
`{mother, father, n}` plan with no engine executor.

**Decision:** `mating_design(mothers, fathers, design = c("random", "factorial",
"nested", "diallel", "half_diallel"), ...)` writes a plan; `mate(plan, ..., seed,
prefix)` runs it — rows in plan order after one `set.seed()`, `cross()` per row, a
self when mother and father are the same individual, a doubled haploid for
`design = "dh"` — across one or several named pools sharing a map (plan columns
`mother_pool` / `father_pool`), and names progeny `<prefix>_<k>` (default the pool
names joined by `x`, SPEC D20). A one-row plan draws the same random numbers and
gives the same progeny genotypes as the equivalent `cross()`, `selfcross()` or
`double_haploid()` call (the object differs in ids, `origin` and the `plan`
attribute). Identity is the pedigree key, not the pool names given to
`mate()` (those only locate each parent): equal ids in founder pools imported with
different `as_population(pool =)` labels are different individuals, never a self, while
the same genotypes imported twice under one label (or none) are the same individuals,
as are those of `pop[1:3]` and `pop[2:4]`; `"nested"` gives each father exactly its mothers, skipping itself within
one pool. `mating_design(seed =)` seeds only `"random"`: the other designs draw nothing
and leave the RNG stream untouched. An `NA` in a plan's `design` column takes that
row's default. Each individual enters a parent set once (a second id of the same
individual is dropped), so no identity pair is listed twice; with character ids on
either side individuals are compared by id; `"random"` draws a mother among those with
an admissible father.
`recurrent_selection()` keeps its own `.intermate()`: it draws each random pair
between meioses, and moving it onto a pre-drawn plan would change the RNG stream,
so the SPEC's "refactor, bit-identical" was not attainable; it is not refactored.

**Date:** 2026-09-27

---

## DECISION-026: combining ability on a frozen A + D architecture

**Context:** `docs/SPEC-block3b.md` item 1 (BD catalog: half-sib RS with a tester,
hybrid development, RRS-1b). BD composed an additive-only GCA from `cross()` +
`additive_value()`; no engine scorer existed.

**Decision:** `combining_ability(candidates, testers, qtn, a, d, design =
c("topcross", "factorial", "diallel"), method = c("expected", "simulated"), ...)`
(SPEC D3: frozen `(qtn, a, d)` as `genotypic_value()`; D4: both methods together).
`"expected"` is the conditional expectation of every cross given the parents'
genotypes, `E[G_ik] = -2 d g_i g_k + (g_i + g_k)(a + d) - a` per locus with
`g = x/2` (package derivation), no RNG; `"simulated"` realizes `n_progeny` per cross
with `mate()`, scores `genotypic_value()`, optional residual via `phenotype_value()`
(D6: `phenotype_value()` gains `d =`, the fixed-scale `A + D` phenotype with a
broad-sense `h2`). Centering: factorial / topcross GCA = row (column) mean − μ; diallel
(Griffing's method-4 layout) `g_i = (m_i − μ)(p − 1)/(p − 2)`; `SCA = Y − μ − g_i − g_k`,
GCAs and each candidate's SCAs sum to zero. Candidates and testers must each list an
individual once (by pedigree key), so both methods score the same crosses. With testers = the candidates' own
population the *expected* GCA equals half the DECISION-019 breeding value exactly, and
with `d = 0` every *expected* SCA is zero; a *simulated* SCA carries Mendelian sampling
and is not zero for `d = 0` from heterozygous parents (zero only in expectation, or exactly from fully inbred parents when no residual is requested -- an independent
residual per progeny is environmental sampling error that inbred parents do not remove). `ref` is refused without `h2`, and `seed` restores the
caller's RNG state. Epistasis is outside the model. `template_effects(sim, trait, rep)` (D5)
exports the realized-scale `a`, `d` a simulation uses (the `.layer_scaled_effects()`
reconstruction behind `on = "bv"`), refusing epistasis / `"complex"` and a derived
`transcriptome()` layer with `prop > 0` (its heritable part is expression-mediated, not
per-locus, so a template would silently omit it). The diallel is Griffing's method 4
(no reciprocals); identity is the pedigree key, which includes the founder `pool` label.

**Date:** 2026-09-28

---

## DECISION-027: progeny testing

**Decision:** `progeny_test(parents, mates, qtn, a, d, n_progeny, h2 | var_e, ref,
seed)` (SPEC D15: its own export) mates each parent to `n_progeny` distinct random
mates, one progeny each — a genuine half-sib family: never the parent itself, each
individual of `mates` once (by pedigree key), and an error when fewer distinct mates
exist — scores on the frozen architecture with an optional residual and returns each
parent's progeny mean with the progeny `Population` (pedigree recorded). With mates drawn
from one population separate from the parents, the expected progeny mean is half the
breeding value in the mates' population plus a common constant, the breeding value using the average effects `α = a + d(q − p)` at
the mates' frequencies (so dominance enters the parent-dependent part of the mean through them, not as a parental dominance deviation;
dominance also stays in the constant common to all parents, `a(p - q) + 2dpq` per locus under HWE, so at `a = 0`, `p = 1/2` all parents
have the same expected mean `d/2`, not zero); each parent is listed once (by pedigree key); the accuracy on `n` half-sib
records of an additive trait, under the classical half-sib assumptions (large
random-mating, non-inbred mate population; a different, independent mate per
progeny), is `sqrt(n h² / (4 + (n − 1) h²))` (package derivation),
validated by simulation. `h2` converts to a residual variance from the *broad-sense*
variance of the genotypic value, so with `d != 0` it is not the `h2` in the accuracy formula
(stated for an additive trait); `ref` is refused without `h2`; `seed` restores the
caller's RNG state (mates, meioses and residuals are drawn from one stream in that order).

**Date:** 2026-09-28

---

## DECISION-028: independent culling and tandem selection

**Decision:** `select_ind(method = "culling", culling = <per-trait proportions>,
sequential = FALSE)` keeps individuals in the top `culling[t]` of every trait
(simultaneous; `sequential = TRUE` culls traits in order among survivors, SPEC D17);
the count kept is emergent (`n` / `prop` / `intensity` must be `NULL`), `direction` may
be per trait, and `on` may be an individuals × traits matrix of external predictions
(D18). Tandem selection is a `trait` vector on `pedigree()` / `recurrent_selection()`,
recycled over generations (a scalar `trait` is unchanged, bit-identical). Under Hazel &
Lush's (1942) idealized conditions the package derivation gives index : culling :
tandem = `sqrt(T) i(p)` : `T i(p^{1/T})` : `i(p)`, i.e. 1 : 0.907 : 0.707 at `T = 2`,
**Amendment 2026-09-29:** a tandem `trait` schedule is validated against the
phenotype callback's `n_traits` (clear error), `select_ind()` accepts only one
in-range `trait` outside `method = "culling"` (a vector was recycled silently), and
`direction` must be a single value in the schemes.

`p = 0.1`, reproduced by simulation.

**Date:** 2026-09-28

---

## DECISION-029: marker-assisted selection and gene pyramiding

**Decision:** `marker_select(pop, markers, favorable, requirement = c("carrier",
"homozygote"), min_markers, n | prop, rank_on, direction, seed)`: a founder-free
foreground filter (the favourable allele given per marker), staged pyramiding
(`min_markers`, D7), ranking of feasible candidates on any score, seeded tie-break.
`additive_value()` is documented as the MARS marker index (D8: no alias); an index on
the simulation's own causal loci weighted by their average effects (`a` for an additive
architecture, `a + d(q − p)` with dominance) is labelled an oracle and keeps the
individuals `on = "bv"` keeps, up to ties at the cut-off (random vs input-order
tie-breaks). `mabc_select()`'s `exclude_interval` bounds are inclusive; with equal
weights only the markers are removed (the interval's genome share stays in the
denominator unless `marker_weights = "interval"`), and the `1 - 2^-(t+1)` recovery
benchmark is Mendelian transmission (the F&M 2005 attribution carries no section
pointer). `recurrent_parent_recovery()` warns when no marker is informative.
Validated against Mendelian
F2 ratios (1/16, 9/16), Haldane's map function for linked targets, `mabc_select()`'s
foreground, and `on = "bv"` for the oracle index.

**Date:** 2026-09-28

---

## DECISION-030: known-variance BLUP, the A matrix and a methods manifest

**Decision (SPEC D12, D13, D14):** `predict_ebv(x, pheno, method = c("gblup",
"pedigree"), h2 | var_a + var_e, ref, K, base_freq, ridge)` — BLUP with the variance
components known (a simulation knows them; no REML, no new dependency), in the GLS form
`û = K Z'(Z K Z' + λI)⁻¹(y − 1μ̂)` (no `K⁻¹`, so singular `G` is fine), equal to
Henderson's MME to 1e-10; reliability `1 − PEV/(K_ii σ²_A)`; GBLUP marker effects
back-solved on `g_matrix()`'s scale when `G` is built internally (`ridge = 0`); a
supplied `K` must be finite, symmetric (used in its symmetric form) and positive
semidefinite up to floating-point rounding only (`n ε` times the largest eigenvalue on
the correlation scale), with `K + λI` over the phenotyped individuals positive
definite (its Cholesky factor) and every reliability in [0, 1]; it must name its
individuals identically on both axes (reordered to `x`), and carries no marker
scale; an unrepresentable `var_e / var_a` or a non-finite solution is an error, never
`NaN` EBVs. `a_matrix(pop, ids, founder_f = 0)` by the tabular method from the
recorded pedigree (selfs `F = (1 + F_P)/2`, doubled haploids `A_ii = 2`); founders are
unrelated with `A_ii = 1 + founder_f` (default `0`, non-inbred; inbred lines such as the
bundled maize panel need `founder_f = 1`, giving `A[P, F1] = 1` and `A[P, P] = 2`).
`predict_ebv()` is one trait, one record per individual, intercept only; with `h2` the
phenotypic variance is that of `ref` (default the *phenotyped subset*, truncated after
selection; because `var_a = h2 * V_P` and `var_e = (1 - h2) * V_P`, `V_P` (hence `ref`) rescales var_a and var_e together and leaves lambda, EBV and reliability unchanged -- only the reported absolute variance attributes change); one record gives EBV 0 and reliability 0 (NA where the prior variance K_ii sigma_A^2 is zero); `a_matrix(founder_f=)` names are pedigree keys or (unambiguous) founder ids, and an id shared by several founders is an error; `ref` without `h2`, duplicated names in `prediction_accuracy()`, or one
individual phenotyped under two ids are errors, and a supplied `K` whose correlation
scale overflows is refused as not positive semidefinite.
`prediction_accuracy(ebv, truth)` reports `cor` and the regression slope.
`selection_methods()` lists the engine's selection operators (tandem included) for BD's
SPEC-0006. The
EBV is an estimate from observable phenotypes; multi-trait and single-step BLUP are
deferred (D14, kept in TODO).

**Date:** 2026-09-28

---

## DECISION-031: crossbreeding

**Decision (SPEC D19–D21):** `breed_composition(pop)` — expected breed fractions from
the pedigree (founder pools; mean of parents); `heterosis(pop, breeds, qtn, a, d)` —
realized heterosis (mean minus the composition-weighted breed means) and the expected
pairwise F1 heterosis computed from the breeds' own genotypes, which equals
`Σ d (p_A − p_B)²` under within-breed HWE and `Σ d (p_A + p_B − 2 p_A p_B)` for
inbred-line breeds (package derivation; the SPEC's HWE-only form generalized);
`crossbreed(breeds, system = c("two_way", "backcross", "three_way", "terminal",
"rotational"), n_progeny, generations, sire_breed, seed)` as a scheme wrapper over
`mate()`. Each name in `breeds` must be the founder pool (`as_population(pool =)`) its
population traces to, checked, so compositions and breed means cannot be attributed
to the wrong breed. `expected_f1` covers distinct pairs (`NA` diagonal). The empty string and
`"<unassigned>"` (the composition column for unlabelled founders) are reserved pool
labels in `as_population()`, so an unlabelled founder cannot be read as a breed. The HWE case is Falconer & Mackay's `H_F1 = Σ d y²`; the general
form `d [h_AB − (h_A + h_B)/2]` is derived in the help. Validated: compositions ½:½, ¾:¼, ¼:¼:½, rotation → 2/3 : 1/3; F2 retains ½ only when the two Hardy-Weinberg deviations sum to zero (see THEORY_REVIEW H1) and
a two-breed rotation ≈ 2/3 of the F1 heterosis; additive architectures give 0.

**Date:** 2026-09-28

---

## DECISION-032: multi-trait known-covariance BLUP

**Decision (lifts the multi-trait half of SPEC D14; single-step `H` stays deferred):**
`predict_ebv()` takes an individuals x traits `pheno` matrix (row names ids, column
names traits, `NA` for a missing record) with `var_a` = `G0` and `var_e` = `R0`, the
known traits x traits base additive genetic and residual covariance matrices, and
returns an individuals x traits matrix of predictions. One verb, not a new export: a
vector `pheno` is single-trait BLUP exactly as before. The model is Henderson & Quaas
(1976): `y_t = 1 μ_t + Z_t u_t + e_t`, `Var(u) = G0 ⊗ K`, residuals correlated only
within an individual (`R0` between one individual's records, 0 between individuals), a
mean per trait. It is solved in the GLS form of DECISION-030 over the stacked records,
`û = Cov(u, y) V⁻¹ (y − X μ̂)`, so a singular `K` or `G0` (genetic correlation 1) is
fine; reliability is `1 − PEV / (G0_tt K_ii)` with PEV from
`P = V⁻¹ − V⁻¹X (X'V⁻¹X)⁻¹ X'V⁻¹`. Missing records are allowed (every trait needs
one). `h2` / `ref` cannot set between-trait covariances, so they are rejected for a
matrix `pheno`; the covariance matrices are validated as symmetric and PSD (the
DECISION-030 correlation-scale check) with positive variances, and named axes must be
the traits (reordered to `pheno`'s). GBLUP marker effects per trait come from the same
back-solve. The BLUP of an aggregate genotype `a'u` is `a'û`, ranked through
`select_ind(on =)`. Validated: equals Henderson's multi-trait MME (with
`(G0 ⊗ A)⁻¹`) to 1e-9 for `μ`, `û` and PEV with partly missing records; a one-column
matrix equals single-trait BLUP, and diagonal `G0`, `R0` separate into per-trait BLUP;
marker effects reproduce the GEBVs; under its own model a correlated recorded trait
lifts the accuracy of an unrecorded one (0.67 vs 0.45 single-trait at `r_g = 0.8`),
the slope of truth on prediction is ≈ 1, and mean reliability ≈ squared accuracy.
Review r1: one-marker marker effects kept as a markers x traits matrix; a
correlation-scale entry that overflows (e.g. `[[1e-320, 1e200], [1e200, 1e-320]]`) now
fails the PSD check directly instead of making its tolerance infinite (this also closes
the same hole for a supplied `K` in single-trait BLUP).

## DECISION-033: position-sensitive layer sub-seed (grammar seed rule)

**Context:** the grammar derived a layer's sub-seed from the *sum* of the character codes
of its draw label (`additive_rep<r>`, `residual_t<t>`, ...), which is permutation-invariant:
`additive_rep12` and `additive_rep21` (and `residual_t12` / `residual_t21`) got the same
sub-seed, so `vary_qtn` replications 12/21 and multi-trait residuals 12/21 were
byte-identical (84 of 150 realistic labels collided; audit 2026-09-29, GRAM-O1 / EFF-F1).

**Decision (owner D2):** the label is reduced with a position-sensitive polynomial rolling
hash (`h <- (257 h + code) mod 2147483629`, in double precision) and mixed with
`(seed, occurrence)` exactly as before (`seed*1009 + hash*7919 + occurrence*104729`, mod
`.Machine$integer.max`; `NULL` and the largest valid seed stay safe). The rule is
**(seed, layer type, occurrence of that type)**, *not* the layer index (SPEC §6 used to say
`layer_index`, which would break the invariance): reordering or adding layers of *other*
types never changes a layer's draws, but inserting a same-type layer before an existing one
shifts that layer's occurrence index and so changes its draws (documented, tested).
The property claimed is **collision resistance over ordinary ranges**, not injectivity: the
sub-seed is a 31-bit integer, so no such map can be injective. Absence of duplicates over
`<type>_rep1..500`, `residual_t1..100` (and the vQTL / complex families) and all layer types x
occurrence 0:20 is asserted in `tests/testthat/test-audit-grammar.R`. A real cross-family
collision exists far outside realistic use (independent review 2026-09-30): with `seed = 123`,
`.layer_seed(123, "transcriptome_rep106", 0)` and `.layer_seed(123, "residual_t40160", 3)` both
give 2116039371 (replication 106 of the first transcriptome layer vs the residual of trait 40160
in replication 4). Pinned in `tests/testthat/test-fix2-grammar.R`.

**Consequence:** every seeded grammar value changed; the grammar owes v1 no bit-parity
(DECISION-009), and the frozen `create_phenotypes()` keeps its own seed arithmetic.
Hard-coded seeded expectations in the grammar tests were re-derived.

**Date:** 2026-09-29

---

## DECISION-034: report the realized additive/dominance partition on shared loci

**Context:** in the variance-partition coding the additive component (dosage x effect) and the
dominance component (heterozygote indicator x effect) are scaled separately, so on shared loci
`Var(g) = prop_A + prop_D + 2Cov(c_A, c_D)` for one additive and one dominance layer (general
identity `Var(g) = Var(c_A) + Var(c_D) + 2Cov(c_A, c_D)`; e.g. two additive layers and one
dominance layer: 0.888 vs 0.588), and under Hardy-Weinberg
`Cov(dosage, het) = -(2p - 1) 2pq` at each locus. With the fixed-sign geometric series each locus's
cross term has the sign of `1 - 2p` (p = frequency of the counted +1 allele), so it is one-signed across
loci only when the counted-allele frequencies are all on one side of 0.5; with frequencies straddling 0.5
the per-locus terms partly cancel (exact HWE columns at p = 0.2 and 0.8 give Cov(x,h) = +0.194 and
-0.194). Cross-locus (LD) terms are additional. The effect is a **structural**, allele-coding-dependent bias
(outbred MAF 0.10-0.20, A 0.4 + D 0.1: realized H2 0.63 with the minor allele coded +1, 0.25
with the same loci/seeds coded the other way, requested 0.5), not the "finite-sample"
fluctuation the help and SPEC §2 claimed (audit 2026-09-29, GRAM-F1).

**Decision (owner D3):** document it and report it. When a non-orthogonal `additive()` layer
and a `dominance()` layer overlap in loci, the object carries `$ad_report` (per trait: `requested`
= summed `prop`s on the unit-variance scale, `realized` = `Var(g)/V_P`, the average-effect partition
`var_A`, `var_D`, `cov2_AD` = `2Cov(A, D)/V_P`, summing to `realized`, and the component partition
`var_cA` = `Var(c_A)/V_P`, `var_cD` = `Var(c_D)/V_P`, `cov2_comp` = `2Cov(c_A, c_D)/V_P`, also summing
exactly to `realized` for any number of layers; all as fractions of V_P). The exact link to the request
is `realized - requested/V_P = (var_cA + var_cD - requested/V_P) + cov2_comp` (NOT
`realized - requested`; the bracket is the covariance among same-type layers, zero for one additive +
one dominance layer),
`print()` shows it with a note recommending `additive(orthogonal = TRUE, a =, d =)` (DECISION-020),
and the "finite sample / usually tracks closely" wording is removed from the help and SPEC §2.
The requested budget (`var_budget`) is unchanged (it stays the request); the coding convention
itself is not changed. Epistasis products and LD-correlated causal loci can correlate components
too; no report is produced for those.

**Date:** 2026-09-29

---

## DECISION-035: transcriptome realized heritability naming (`h2_realized`, `h2_allocated`, `h2_var_ratio`)

**Context:** the audit-fix round 1 made `genes$h2_realized` the bounded allocation
`Var(G)/(Var(G)+Var(R))`, which equals the target on the reference panel by construction and so
is not a realized heritability (independent review 2026-09-30, round-2 A1: at `n = 3`,
`h2 = 0.5`, seed 1, the old `h2_realized` was 0.5 for every gene while `Var(G)/Var(E)` ranged
0.27-617.7).

**Decision (owner, round 2 A1):** `h2_realized` = the realized `Var(G)/Var(P)` computed from the
realized genetic values and the realized expression (it includes `2Cov(G, R)`, equals
`h2/(1 + gr_cov)` on the reference panel only when `Var(G) + Var(R) = 1` (general form
`Var(G)/(Var(G) + Var(R) + gr_cov)`), and is not bounded by 1; a warning is issued below 30
individuals). `h2_var_ratio` is kept as an identical alias for compatibility. The bounded
allocation `Var(G)/(Var(G)+Var(R))` is now `h2_allocated`: a variance allocation, not a
heritability, equal to the target on the reference population by construction and informative
under `predict()`. The mimic GREML identifiability guard tests the spectrum of `K` after
projecting out the intercept (`K = a(I - 11'/n) + b 11'` is flagged; h2 is set to 0 with a
warning). The marginal-epistasis-share text gives `epsilon/(epsilon + (1 - epsilon)/s_ct^2)`,
`s_ct^2 = 1 + 2 sqrt(omega(1 - omega)) cor(c, t)`; this holds for the nondegenerate blend, while at the
exact-cancellation fallback (perfect negative cis/trans correlation) the cis part is dropped and the
share is `epsilon`. Round 3 (R3-4): the `h2_*`, `cis_fraction_realized` and
`epistasis_realized` ratios are scale-free, reported whenever the denominator is finite and positive
(no absolute 1e-12 cutoff; an exactly zero `Var(P)` still reports 0).

**Date:** 2026-09-30

---

## DECISION-036: `as_numeric()` records the counted allele as an R attribute (`counted_allele`)

**Context:** the numeric `allele` label alone cannot reveal which allele was coded `+1`: an all-`AA`
and an all-`GG` panel converted separately both encode `+1`, so crossing them mis-reads
heterozygotes (independent review 2026-09-30, round-2 A2).

**Decision (owner, round 2 A2):** the allele counted as `+1` (`2` under `code_as = "012"`) is
recorded per marker in `attr(, "counted_allele")` of the `as_numeric()` result (absent under
`model = "Dom"`; `NA` where the label cannot name it), computed once in `.apply_coding()` so every
reader gets it. It is not a new column and does not change the `allele` label or any dosage.
`as_population()` keeps it as `map$counted`; the cross-pool orientation guard errors when both
populations carry it and disagree at a marker, and falls back to the label-only check (warning on
opposite order, error on disjoint alleles) when either side lacks it. Legacy numeric input without
the record behaves as before. Persistence through numeric text files is out of scope (a text file
cannot carry the attribute; row-subsetting a data frame drops it). Documented in `?as_numeric`,
`?as_population` and `docs/BACKEND_CONTRACT.md`.

**Addendum (round 6, Codex G5):** the default output file name of `as_numeric()` for an inline object
(DECISION-038 `.geno_label()`) is sanitized: characters outside `[A-Za-z0-9._-]` become `_` and the ends are
trimmed (a symbol `hmp` gives `hmp_numeric.txt`; an inline data frame
`inline_data.frame_1077_x_670_numeric.txt`), so the label's `<`/`>` no longer reach the file system. The
written content is unchanged. The `counted_allele` validation of the haplotype constructor is in DECISION-039.

**Date:** 2026-09-30

---

## DECISION-037: v1 duplicate `chr_pos` accepted in pleiotropic/partial, rejected in LD

**Decision (owner, round 2 A3):** the frozen `create_phenotypes()` accepts duplicated `chr_pos`
markers in the fully pleiotropic and partially pleiotropic architectures (no error; documented in
`?create_phenotypes`, section "Duplicated marker positions"); the `"LD"` architecture keeps its own
stricter check and rejects them with an "LD contract" error. The v2 grammar
(`simulate_phenotype()`, one-call and layered) also accepts duplicated `chr`/`pos`.

**Date:** 2026-09-30

---

## DECISION-038: entry-mean replication (`reps`) in the grammar

**Status:** accepted (round 5, breedingDesigner SPEC-0020 item 6).

**Context:** AlphaSimR `setPheno(varE, reps)` gives each individual a phenotype that is the mean of
`reps` independent records of the same genotype, so the error variance is `varE / reps`. The grammar
could only produce single-record phenotypes.

**Decision:** `simulate_phenotype(reps = 1)` and `complex_phenotypes(reps = 1)` take a positive whole
number (scalar or one per trait). The residual of trait t, including the vqtl heterogeneity component,
is drawn exactly as for `reps = 1` (same sub-seed, same number of draws) and then multiplied by
`1/sqrt(reps[t])`; `reps = 1` skips the multiplication (bit-identical output). The realized residual is
therefore the `reps = 1` realized residual divided by `sqrt(reps[t])`, so its realized variance is the
`reps = 1` residual variance divided by `reps[t]`: for a vqtl trait `[V0 + Vv + 2Cov(e0,ev)]/reps[t]`, not
the nominal `V_E/reps[t]`, because the two standardized components have a non-zero sample covariance. The
genetic value and the transcriptome component are not rescaled (the expression-mediated environmental part is a
persistent property of the entry, not record noise).

**Scales:** `h2`, every layer `prop`, `var_budget` and the printed residual share are single-record
(`V_G + V_E = 1`). Target (expected) heritabilities: single-record `h2 = V_G/(V_G+V_E)` (what `h2`
requests) and entry-mean `V_G/(V_G + V_E/reps)`; the allocation formula is a target. Realized
heritabilities are `Var(G)/Var(y)` from the realized values: entry-mean
`Var(y_bar) = V_G + V_E/reps + 2Cov(G,e)/sqrt(reps)`, record `Var(y) = V_G + V_E + 2Cov(G,e)`; the
allocation formula is the realized value only when the sample `Cov(G,e) = 0`. The stored phenotype is the
entry mean, so `.realized_h2()` and the shares "of V_P" (`$ad_report`, `mediation_split()`) are entry-mean;
the single-record realized value is reconstructed as `Var(G)/Var(y + (sqrt(reps)-1)e)` and printed
alongside when `reps > 1`. `print()` shows `reps (per trait) = [1, 4, 1]` when `reps` varies.
`complex_phenotypes()` ignores its inputs' own `reps` (their residuals are discarded); its own `reps`
scales the common residual.

**Records and the transcriptome:** without a transcriptome layer, records are independent given the
genotype (no shared permanent environment is modelled). With a derived transcriptome layer, the
environmental transcriptome component `Tx_env` is a persistent entry-level quantity that is NOT redrawn
per record: replication is conditional on the fixed transcriptome covariate, and only the phenotype
residual is rescaled (its value by `1/sqrt(reps)`, its variance by `1/reps`). The realized denominator is `Var` of the full stored phenotype
(`G_marker + Tx_total + e/sqrt(reps) + mu`). Alternative considered and rejected for now: redraw `Tx_env`
per record (would change numbers; revisit if a repeated-measures model is wanted).

**Rust boundary:** unaffected (the stochastic draw stays in R; the rescale is deterministic arithmetic in R).

**Naming cost (item 1, same round):** `simulate_phenotype()` (and `as_numeric()`) named the genotype input
with `deparse(substitute(geno))`, which renders a whole inline `Population` when called through
`do.call()` (9 s for 2000 x 14000). A constant-time `.geno_label()` replaces it: a symbol keeps its name,
a small call keeps its old first deparse line, anything carrying a large embedded object gets a type label
such as `"<inline Population>"`. No change to any simulated value.

**Date:** 2026-09-30

---

## DECISION-039: phased-haplotype constructor `population_from_haplotypes()` (SPEC-0020 item 5)

**Decision:** `population_from_haplotypes(cis, trans, map, ids, pool, individuals_in_rows)` is the public
constructor for known-phase founders and `haplotypes()` its inverse. Layout is markers x individuals by
default (the `Population` convention; `individuals_in_rows = TRUE` accepts the individuals x markers
layout, and a transposed input without the flag errors). Entries are 0/1 with 1 = the counted (+1)
allele, dosage = `cis + trans - 1` (the `as_population()` encoding; heterozygote phase is the caller's, not
guessed; `cis`/`trans` carry no maternal/paternal meaning). The map is validated by the helpers
`as_population()` uses (`.check_map`, the cM-units check, now also the internal `.make_map()` and
`.check_pool()`). `map$counted` is recorded only when the supplied map has a `counted` column (for
example `pop$map`); it is never inferred from the `allele` label, which does not say which allele is
counted (DECISION-036), so without it the cross-pool guard stays label-only. `as_population()` ->
`haplotypes()` -> `population_from_haplotypes()` is `identical()` to the original and seeded
`cross()`/`selfcross()`/`double_haploid()` are bit-identical. No RNG, no Rust. `as_population(haplotypes =)`
was not added (a separate function leaves the `as_population()` signature untouched).

**Addendum (round 6, Codex H1):** a supplied `map$counted` is validated by `.check_counted()`: one entry
per marker; character (a numeric column is rejected); each non-`NA` entry one non-empty allele symbol (`""`
is rejected) and, when the `allele` label names two alleles, one of them. Only `NA` means unknown; the
internal `.make_map()` and the orientation guard `.check_orientation()` also treat `""` as unknown, so a
malformed record can no longer silently disable the guard (DECISION-036).

**Date:** 2026-09-30

---

## DECISION-040: batched, integer-I/O crossing core (SPEC-0020 item 2)

**Decision:** the meiosis core gains `mate_many_core()`: parental strands as 0/1 integer vectors in the
caller's marker order, a marker-rank permutation, a table of matings (strand indices, design, n) and one
shared event stream (the draws of the matings concatenated in plan order); it returns the progeny strands
as integer vectors in the caller's marker order. Every crossing function goes through `.mate_many()`;
`mate()` executes all plan rows in one call. R still draws every random quantity, row by row in plan order,
exactly as the same sequence of `cross()` / `selfcross()` / `double_haploid()` calls would (the batch equals
the sequential run draw for draw, bit for bit, including pedigree keys); Rust draws nothing
(DECISION-006/012). The string kernel entry points (`mate_haplotypes_core`, `meiosis_core`,
`gamete_masks_core`) stay untouched for the parity tests.

**Measured** (14,000 markers, 20 chromosomes, one heterozygous parent, 100 calls of `double_haploid(parent,
n = 100)`): 14.1 s to 2.3 s through the single-call API (6.1x), 14.3 s to 2.1 s through one 100-row
`mate()` plan (6.9x); what remains per call is the R draws (must follow R's stream), the pedigree-key
hashing and a few ms of Rust.

**Rejected:** moving the draws to Rust (breaks R-seed reproducibility); vectorising the Poisson draw
(changes the isqg stream).

**Date:** 2026-09-30

---

## DECISION-041: crossover interference as an option of the meiosis draws (SPEC-0020 item 3)

**Decision:** `interference = NULL` (default: Poisson, isqg stream, bit-identical) or `list(nu, p)`
appended to `cross()`, `selfcross()`, `double_haploid()`, `mate()` and `crossbreed()`: the two-pathway
gamma model on the bivalent. Chiasmata live on the four-strand bivalent with intensity 2 per Morgan; a
share `p` is a non-interfering Poisson pathway, the rest a stationary renewal process with
Gamma(shape `nu`, rate `2 nu (1-p)`) gaps and equilibrium start `U * Gamma(nu+1, rate)`; a gamete keeps each
chiasma independently with probability 1/2 (no chromatid interference). The map keeps its meaning:
expected crossovers per Morgan = 1 for every `nu`, `p`. The draws are made in R and handed to the unchanged
Rust core as sorted chiasma positions in [0, L]. `nu = 1` or `p = 1` is Poisson in distribution but uses its
own stream, not isqg's. `1 <= nu <= 1e6` (negative interference is not modelled; larger `nu` made the renewal draw hang, so the
rate is computed as `nu * (2 * (1 - p))` and the renewal loop is guarded against zero progress), `p` in [0, 1]. The recombination
fraction is `r(d) = (1 - P0(d))/2` (closed forms in SPEC section 4.6), `nu = 2.6, p = 0` lies within 0.001 of
Kosambi over 0-1 M. Applies to all rows of a `mate()` plan. Citations in `?cross` (McPeek & Speed 1995, Genetics; Housworth & Stahl 2003, Am. J. Hum.
Genet.) are author/year/journal only and unverified; it is not claimed that AlphaSimR's code uses exactly
this construction.

**Propagation (round 6, maintainer note 2026-09-30):** the `interference` option is a property of the
meiosis model, not of one function; every function that draws meiosis forwards it: `single_seed_descent()`,
`bulk()`, `pedigree()`, `recurrent_selection()`, `cross_usefulness()`, `combining_ability(method =
"simulated")` and `progeny_test()`, in addition to `cross()`, `selfcross()`, `double_haploid()`, `mate()` and
`crossbreed()`. It is the last formal, default `NULL` = isqg stream (bit-identical), validated up front by
`.check_interference()` with the caller's name, and captured once so every generation or cycle of a scheme
uses it (`cross_usefulness()` also passes it to its nested `single_seed_descent()`).
`combining_ability(method = "expected")` errors on a non-`NULL` value rather than ignoring it. Functions that
do not run meiosis (`select_ind()`, `select_ocs()`, `optimum_contribution()`, `sample_parents()`,
`predict_ebv()`, `marker_select()`, `mabc_select()` scoring, `g_matrix()`, `a_matrix()`, the grammar) do not
take it. Rationale: a scheme with a single `interference` value must not silently mix models. The Rust
boundary is unchanged (R draws every event). The legacy `create_phenotypes()` is untouched (DECISION-008).

**Date:** 2026-09-30

---

## DECISION-042: residual correlation between traits (`resid_cor`)

**Decision:** `simulate_phenotype(resid_cor = NULL)` and `complex_phenotypes(resid_cor = NULL)`: the target
correlation of the traits' **residuals** (`cor` stays the genetic correlation, DECISION-010). `NULL` leaves
the draw path untouched (same sub-seeds, same draws, bit-identical output). A scalar (common pairwise value,
>= `-1/(n_traits - 1)`) or an `n_traits x n_traits` matrix, validated symmetric, unit diagonal, entries in
[-1, 1], positive semi-definite (the same style as `cor`); needs `n_traits >= 2`.

**Construction:** per trait the unit-variance residual is drawn under the unchanged `residual_t<t>` sub-seed
(`complex_resid_t<t>` for `complex_phenotypes()`), the n x n_traits matrix is multiplied by the upper Cholesky
factor of the target (eigen square root if singular), and each column is re-standardized (mean 0, exact unit
sample variance) and scaled by `sqrt(1 - sum(prop))`. Each trait therefore keeps exactly its `h2`-implied
residual variance; realized h2 and `var_budget` are unchanged and only the cross-trait correlation is
induced. RNG stays in R (DECISION-006).

**Scope / limits:** the realized sample correlation equals the target up to sampling error of order
`1/sqrt(n)` (it is not forced exact, unlike the variance). A `vqtl()` heterogeneity component is drawn
independently per trait after the mixing and dilutes the correlation; `reps` rescales each trait by its own
`1/sqrt(reps)` afterwards, which keeps the correlation. The phenotypic correlation mixes the genetic and
residual correlations by their variance shares. This is the grammar counterpart of v1 `cor_res` (frozen
legacy function unchanged, DECISION-008); no v1 bit-parity is owed (DECISION-009).

**Date:** 2026-10-02

---

## Decision Log Summary

| ID | Decision | Status |
|----|----------|--------|
| 001 | Single package, internal modularization | locked |
| 002 | Port isqg algorithms to Rust (own the code) | locked |
| 003 | `create_phenotypes()` v1 signature preserved | locked (refined by 008) |
| 004 | Multi-generation stays in simplePHENOTYPES | locked |
| 005 | Expression-based simulation → v3 | superseded by 022 (2026-10-02) |
| 006 | Rust surgical/bottleneck-only; stochastic core in R | locked (reaffirmed by 011) |
| 007 | PleioArch for `"pleiotropy"` | locked (refined by 010) |
| 008 | `create_phenotypes()` = frozen legacy, not a delegation shim | locked (2026-06-10) |
| 009 | New grammar has no bit-for-bit v1 parity; references guard the legacy fn | locked (2026-06-10) |
| 010 | Correlation via `cor`: PleioArch (2 traits) + Cholesky fallback (>2, warns) | superseded in part by 013 |
| 011 | Single rextendr package; Cargo-workspace plan dropped | locked (2026-06-10) |
| 012 | Meiosis randomness drawn in R; Rust core pure; exact isqg bit-parity is the gate | locked (2026-09-07) |
| 013 | PleioArch generalized to any n_traits; Cholesky fallback removed | locked (2026-09-07) |
| 014 | `"ld"` = two-trait linked *distinct* causal loci (`direct` default); companion-annotation model removed | locked (2026-09-09) |
| 015 | Selection engine `select_ind()` (mass/within/among/combined/Smith–Hazel/random; pheno/gv/custom criterion) + scheme wrappers (SSD/bulk/pedigree/recurrent) + `c.Population()` | locked (2026-09-10) |
| 016 | Modern methods: `g_matrix()` (VanRaden), `optimum_contribution()` (Meuwissen, dependency-free Frank–Wolfe) + `sample_parents()`, `cross_usefulness()` (Zhong & Jannink/Lehermeier); GS/PS via the `on` hook; PopVar real-data cross selection deferred | locked (2026-09-11) |
| 017 | Designer = DAG-JSON contract + headless `run_design()`/`validate_design()`/`design_breeding_program()` + `design_script()` (R runs today, Python mirror); React Flow canvas ships as a companion app, not in the tarball; `jsonlite` to Suggests | superseded by 018 (home moves) |
| 018 | Designer split into a separate product `breedingDesigner` (Imports simplePHENOTYPES); selection engine stays here; webR in-browser execution + plumber back end; engine adapter (AlphaSimR etc.); move after new pkg green | locked (2026-09-11) |
| 019 | Breeding value (`on = "bv"`, OCS default, index merit) = classical average-effect A = Σαⱼ(xⱼ−2pⱼ), αⱼ = aⱼ+dⱼ(qⱼ−pⱼ) reconstructed analytically from known QTN effects (LD- and HWE-robust); epistasis-induced marginals omitted; index methods ignore `on` | locked (2026-09-13) |
| 020 | orthogonal genotypic model as `additive(orthogonal = TRUE, a =, d =)` (orthogonal in expectation under random mating → per-locus HWE, e.g. F2; LD is fine, but nonrandom multilocus association breaks it; A is the transmitting average effect, Falconer 1985): per-locus a/d, whole value scaled to `prop`; budget reports realized Var(A)/Var(g), Var(D)/Var(g) + an `add_dom_cov` row 2Cov(A,D)/Var(g) (=0 in expectation under random mating; nonzero for structured / finite samples) closing to `prop`; degree of dominance d/abs(a) meaningful; d!=0 requires a het per locus (checked per-locus); `qtn_table()` gains a `d` column; new args appended to the signature (positional compat kept); incompatible with vary_qtn / dominance() / pleiotropy(multi) / ld | locked (2026-09-13) |
| 021 | `phenotype_value(x, qtn, effect, h2/var_e, ref, seed)` = fixed additive value (`additive_value()`) + independent residual on a **frozen** variance (no per-population rescale), so the parametric h²=Var(g)/(Var(g)+var_e) declines as variance is exhausted (faithful cross-gen `on="pheno"`); exactly one of h2/var_e (h2 → var_e=Var(g_ref)(1−h2)/h2); RNG in R; version bump 1.4.0-9002 so downstream can pin | locked (2026-09-14) |
| 022 | transcriptome simulation: hybrid latent-factor eQTL (cis in a physical-bp window + trans via Q≪T latent regulatory factors + non-genetic co-expression modules + gene noise), normalized Gaussian scale; per-gene `h2_g`, cis fraction `omega_g`, residual module fraction `kappa_g` as non-competing knobs; **joint** cis/trans genetic scaling with a reported `cis_trans_cov` budget row (cf. `add_dom_cov`, DECISION-020); **fixed reference calibration** (DECISION-020/021); no annotation/reference data required (named `generic_bulk` profile; synthetic physical-bp coords; optional `mimic=` calibrates the expression generator only, GREML per-gene h²; ships an example annotation + expression dataset for SNP55K); phenotype bases inferred from inputs via a `transcriptome()` layer — markers-only (default) / `expression=` real transcriptome alone / `geno`+`transcriptome=` derived G→E→Y (`y = Z(δ+Bs)+Rs+η`, mediated+direct) / `geno`+`expression=` both real (`y = Zδ+Ẽs+η`, no simulated mediation, realized G–E cov reported) — all with genome-path rigor; counts (`observe_counts()`) and epistatic expression implemented, GRN/tissue specificity remain non-scope; RNG in R (DECISION-006), no Rust yet | locked (2026-10-02; drafted 2026-09-16) |
| 023 | `cor` control extends to dominance + epistasis under "pleiotropy": per-component PleioArch covariance (shared units MVN-correlated, trait-specific independent, split by `pi`), each unit scaled by its realized design-column sd, constant units left out of the allocation; every component targets `cor` (realized correlation converges as units and individuals grow under approximate linkage equilibrium; attenuated on average with few units; strong LD can prevent convergence), and the total targets `cor` when layers' per-trait `prop` profiles are proportional (scalar `prop`), else attenuated with a warning giving its large-sample target; fixed `qtn=` rejected under pleiotropy/ld and `effect=`/`dist` under pleiotropy; under "ld" dominance must reuse the additive linked loci (`same_as_add = TRUE`), epistasis and a second additive layer are rejected (SPEC §5.3 restriction); a derived `transcriptome()` layer's genome-mediated signal is outside `cor` (warned under pleiotropy); additive draw unchanged (bit-identical) | locked (2026-09-25) |
| 024 | A `Population` records its pedigree (`keys` + `pedigree` frame; founders from `as_population(pool =)`, progeny from every mating, ancestors kept by `[`, pooled by key in `c()`); links are deterministic keys, not display ids, so colliding ids stay safe; bookkeeping draws nothing (genotypes / RNG bit-identical); accessors `parentage()`, `families()` | locked (2026-09-27) |
| 025 | `mating_design()` writes random / factorial / nested / diallel / half-diallel plans; `mate()` runs `{mother, father, n}` plans across named pools (one seed, plan order; self / DH rows; ids `<prefix>_<k>`); a one-row plan equals the equivalent `cross()` / `selfcross()` / `double_haploid()`; `recurrent_selection()` keeps `.intermate()` (RNG order) | locked (2026-09-27) |
| 026 | `combining_ability()` (topcross / factorial / diallel; `"expected"` = exact conditional cross means on a frozen `(qtn, a, d)`, `"simulated"` via `mate()`), GCA/SCA centred to sum zero, Griffing method-4 diallel GCA; `template_effects()` exports a simulation's realized `a`, `d`; `phenotype_value(d =)` scores `A + D` with a broad-sense `h2` | locked (2026-09-28) |
| 027 | `progeny_test()`: half-sib progeny of each parent on random mates, scored on the frozen architecture; accuracy √(n h² / (4 + (n − 1) h²)) | locked (2026-09-28) |
| 028 | `select_ind(method = "culling")` (simultaneous or sequential per-trait proportions; emergent count; external prediction matrix via `on`); tandem = `trait` vector on `pedigree()` / `recurrent_selection()` (scalar unchanged) | locked (2026-09-28) |
| 029 | `marker_select()`: foreground carrier / homozygote filter, staged pyramiding (`min_markers`), ranking on any score, seeded tie-break; `additive_value()` documented as the MARS index | locked (2026-09-28) |
| 030 | `predict_ebv()` known-variance BLUP (GBLUP / pedigree, GLS form = MME), `a_matrix()` tabular method, `prediction_accuracy()`, `selection_methods()` manifest; multi-trait / single-step deferred (D14, TODO) | locked (2026-09-28) |
| 031 | `breed_composition()`, `heterosis()` (realized; exact expected F1 from the breeds' genotypes), `crossbreed()` two-way / backcross / three-way / terminal / rotational over `mate()` | locked (2026-09-28) |
| 032 | Multi-trait known-covariance BLUP in `predict_ebv()` (individuals x traits `pheno` with `NA` for missing records; `var_a` = G0, `var_e` = R0; Henderson & Quaas 1976, GLS form = multi-trait MME); single-step `H` still deferred | locked (2026-09-29) |
| 033 | Layer sub-seed = position-sensitive rolling hash of the draw label, mixed with `(seed, occurrence of the layer type)`; replaces the permutation-invariant character-code sum (replications 12/21 and traits 12/21 were byte-identical); collision-resistant over ordinary ranges (31-bit, not injective; see the seed-123 collision in the body); adding/removing/reordering a layer never changes other-type layers, inserting a same-type layer shifts later same-type layers by design; every seeded grammar value changed (no v1 parity owed, DECISION-009) | locked (2026-09-29) |
| 034 | Additive + dominance on shared loci (variance-partition coding): the realized genetic variance is `prop_A + prop_D + 2Cov(c_A,c_D)` for one additive and one dominance layer (general: `Var(c_A) + Var(c_D) + 2Cov(c_A,c_D)`), a structural, allele-coding-dependent bias; report it (`$ad_report`: requested, realized, Var(A), Var(D), 2Cov(A,D), component `Var(c_A)`, `Var(c_D)`, `2Cov(c_A,c_D)` (exact closure for any number of layers), per trait, fractions of V_P) and note `additive(orthogonal = TRUE, ...)` in `print()`, help and SPEC §2; the "finite-sample / usually tracks closely" wording is removed | locked (2026-09-29) |
| 035 | Transcriptome `genes$h2_realized` = realized `Var(G)/Var(P)` (includes `2Cov(G,R)`, not bounded by 1); `h2_var_ratio` = identical alias; bounded allocation `Var(G)/(Var(G)+Var(R))` renamed `h2_allocated` (not a heritability); mimic GREML guard on the intercept-projected spectrum of K | locked (2026-09-30) |
| 036 | `as_numeric()` records the `+1`-counted allele per marker as the `counted_allele` attribute (no dosage, column or label change); `as_population()` keeps `map$counted`; cross-pool guard errors on disagreement, else label-only fallback; not persisted through text files | locked (2026-09-30) |
| 037 | v1 `create_phenotypes()`: duplicated `chr_pos` accepted in pleiotropic / partially pleiotropic, rejected ("LD contract") in `"LD"`; the grammar also accepts duplicated `chr`/`pos` | locked (2026-09-30) |
| 038 | Entry-mean replication in the grammar: `reps` divides the residual (incl. vqtl) by `sqrt(reps)` after the unchanged draw (target residual variance `V_E/reps`; realized: residual value = `reps = 1` realized residual / `sqrt(reps)`, so its variance is the `reps = 1` realized variance / `reps`; vqtl `[V0+Vv+2Cov]/reps`); `h2`/props stay single-record targets; realized H2 = `Var(G)/Var(y)` (includes `2Cov(G,e)/sqrt(reps)`) on the entry-mean scale with the single-record value alongside; replication is conditional on a fixed transcriptome covariate; `reps = 1` bit-identical; constant-time `.geno_label()` replaces `deparse(substitute(geno))` | locked (2026-09-30) |
| 039 | `population_from_haplotypes()` / `haplotypes()`: known-phase 0/1 constructor (1 = counted allele, dosage = cis + trans - 1, markers x individuals), `map$counted` only if supplied (validated: character, non-empty symbol, `NA` = unknown), same map validation as `as_population()`, no RNG, no Rust | locked (2026-09-30) |
| 040 | Batched integer-I/O crossing core `mate_many_core()`; `mate()` runs every plan row in one call; R draws in plan order so the batch equals the sequential run bit for bit; ~6x cheaper per call at 14,000 markers | locked (2026-09-30) |
| 041 | Crossover interference `interference = NULL \| list(nu, p)` (`1 <= nu <= 1e6`) on `cross`/`selfcross`/`double_haploid`/`mate`/`crossbreed` and, by propagation, `single_seed_descent`/`bulk`/`pedigree`/`recurrent_selection`/`cross_usefulness`/`combining_ability` (simulated only)/`progeny_test`: two-pathway gamma model on the bivalent, 1/2 thinning to the gamete, expected crossovers per Morgan unchanged, drawn in R, default NULL = Poisson/isqg stream bit-identical | locked (2026-09-30) |
| 042 | Residual correlation between traits: `resid_cor = NULL \| scalar \| matrix` on `simulate_phenotype()` / `complex_phenotypes()` (genetic stays `cor`); unit draws under the unchanged sub-seeds mixed through `chol(R)` and re-standardized per trait, so each trait's residual variance, realized h2 and `var_budget` are unchanged; `NULL` bit-identical; sample correlation = target up to `1/sqrt(n)`; vqtl dilutes it | locked (2026-10-02) |
