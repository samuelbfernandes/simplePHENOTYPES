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
