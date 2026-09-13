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
effect-generation engine for `architecture = "pleiotropy"` — exact genetic-correlation
control via bivariate-normal effect draws.
**Rationale:** v1 pleiotropy merely shares loci and cannot control rho_g; PleioArch makes
the realized correlation equal the target in expectation.
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
- `n_traits = 2`: the **PleioArch** bivariate-normal engine (DECISION-007) realizes the
  target `cor` exactly in expectation; enforces `cor² ≤ pi_target × pi_secondary`.
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
Cargo-workspace bootstrap in `TODO_newfeatures.md` item 7 (root `Cargo.toml`, members
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

**Decision:** **Yes — the exact engine now handles any `n_traits`, and the Cholesky
fallback is removed.** The bivariate normal generalizes directly to an `n x n`
multivariate normal:

```
Sigma[i,i] = pi_i   * V_i                 pleiotropic share of trait i
Sigma[i,j] = cor_ij * sqrt(V_i * V_j)     the whole genetic covariance
```

Trait-specific effects stay univariate with variance `(1 - pi_i) * V_i`. Because
trait-specific loci are independent across traits, each trait's total genetic variance is
still `V_i` while the entire covariance comes from the shared loci, so every pair realizes
`cor_ij` in expectation. For two traits this reduces algebraically to the bivariate
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

`qtn_table()` reports the linked marker (`companion` — the other trait's causal SNP for
`"direct"`, the shared cause-of-LD locus for `"indirect"`) and the r2 with it (`ld_r2`).
The old reporting-only `.annotate_ld()` / `.ld_reported_qtn()` helpers are removed; the
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
    response tracks R = i·h²·σ_P), `"gv"` (true genetic value — idealized upper
    bound), a **numeric vector** (named by id or in population order), or a
    **function** `on(sim)`. The vector/function is the deliberate extension point
    for **genomic selection, phenomic selection, and other predictive criteria**
    computed externally. True breeding value (`on = "bv"`, the additive average
    effect α) is deferred to the Fisher-orthogonal average-effects feature
    (ROADMAP §3).
  - **Methods:** `mass` (individual truncation), `within_family`, `among_family`,
    `combined` (Lush combined index), `index` (Smith–Hazel multi-trait economic
    index, `b = P⁻¹ G a`), `random` (drift control). Family methods require a
    `family` grouping.
  - **`combined`** ranks on the selection-index prediction of breeding value from
    the individual's own record and its family mean, `b = V⁻¹ c` built from `h2`
    and the within-family additive relationship `family_relationship` (0.25
    half-sibs default, 0.5 full-sibs/selfed). Derived from first principles
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
  `sample_parents()` turns contributions into a drawn parent set for mating.
- **`cross_usefulness()`** — ranks candidate biparental crosses by
  \eqn{U = \mu + i\,\sigma}, the expected value of the best `select_top` fraction of
  progeny. The family is **simulated** with the crossing engine (so linkage enters
  the variance, not an approximation) and scored on the template's **fixed additive
  effects** via a plain dosage×effect dot product — deliberately **not** through a
  re-`simulate_phenotype()` call, which rescales each family's genetic values to a
  fixed variance and would flatten the between-cross σ the criterion depends on.
  Additive (breeding-value) basis only; DH/inbred families carry no dominance.

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
*statistical* additive value). The independent review (Codex) showed this equals
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

**Reaffirms:** DECISION-015/016 (selection engine + modern methods stay here);
DECISION-006 (deterministic, in R).
**Date:** 2026-09-13

---

## Decision Log Summary

| ID | Decision | Status |
|----|----------|--------|
| 001 | Single package, internal modularization | locked |
| 002 | Port isqg algorithms to Rust (own the code) | locked |
| 003 | `create_phenotypes()` v1 signature preserved | locked (refined by 008) |
| 004 | Multi-generation stays in simplePHENOTYPES | locked |
| 005 | Expression-based simulation → v3 | locked |
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
