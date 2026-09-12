# AGENTS.md — simplePHENOTYPES v2

> **Canonical, cross-tool project instructions.** Read this at the start of every
> session — whether you are Claude Code, Codex, or any other agent. This file is
> **committed**; it is the single source of truth. `CLAUDE.md` (gitignored) only points
> here. Keep under ~220 lines.

---

## Critical Rules (read first, every session — human and AI alike)

1. **NEVER add `Co-Authored-By`, "Generated with …", 🤖, or any AI/assistant
   attribution** (Claude, Codex, GPT, Copilot, Gemini, …) to commits, PRs, or files.
   Commit messages read as a human wrote them. This is enforced by
   `.githooks/commit-msg`, which **rejects** such trailers — do not try to bypass it.
   (Rationale: an AI co-author trailer once registered the assistant as a GitHub
   *collaborator*; that must never recur.)
2. **NEVER commit `CLAUDE.md`, `.claude/`, or anything in `.gitignore`.**
3. **NEVER break the v1 `create_phenotypes()` signature** — backward-compat is
   contractual. It is **frozen legacy** (bugfix-only, `lifecycle::badge("superseded")`),
   NOT a wrapper over the new grammar (DECISION-008).
4. **Rust is surgical, not a rewrite.** Only deterministic, profiled bottlenecks and
   C++-origin (isqg) code move to Rust. **Stochastic steps stay in R** (see Rust Boundary).
5. **NEVER push code that fails `devtools::test()`.** Grammar owes v1 no bit-parity
   (DECISION-009); the isqg port owes **exact** parity (DECISION-012).
6. **Push to a remote ONLY with the maintainer's explicit OK**; never force-push
   without explicit authorization.
7. **Do not invent citations, page numbers, equations, or capabilities.** Every
   equation traces to a primary source with a *verified* page. Unsure → say so, don't
   fabricate. (See `docs/THEORY_REVIEW.md`.)

---

## Dual-model development (Claude + Codex)

Two independent models are used deliberately so that a second set of eyes catches
**implementation bugs and — most importantly — errors in the genetic theory** before
they land. The protocol:

- **Roles are split per task: the implementer is never its own reviewer.** One model
  writes the change; the *other* reviews it independently.
- **Every change touching genetics/quantitative-genetics gets a theory review** against
  `docs/THEORY_REVIEW.md` before commit. Genetics = the grammar, PleioArch, effects,
  selection engine (`select.R`, `ocs.R`, `usefulness.R`, `schemes.R`), meiosis/isqg,
  variance partitioning, any cited equation.
- **The reviewer is read-only.** It reports PASS / FAIL / UNVERIFIABLE per rubric item
  with file:line evidence and the primary source; it does not edit code.
- **Objective gate:** `devtools::test()` (and `rcmdcheck` for releases) must be green
  regardless of what either model says.
- Driver: `dev/dual.sh` (see `dev/README.md`). `dev/dual.sh review <paths|--staged>`
  runs a standalone theory review; `dev/dual.sh loop "<task>"` runs implement→review→fix.

The reviewer's job is **not** to agree. It should try to falsify: find the input that
breaks the equation, the frequency where the identity fails, the citation whose page is
wrong. A review that only confirms is a failed review.

---

## Project Overview

simplePHENOTYPES is a CRAN R package simulating pleiotropic, linked, and epistatic
phenotypes from real genomic data (Fernandes & Lipka, BMC Bioinformatics 2020). v2 adds
a composable tidyverse-style grammar (`docs/SPEC.md`), multi-generation cross simulation
(isqg algorithms ported to Rust), a **selection engine** (mass/index/OCS/usefulness),
and keeps one repository under `simplePHENOTYPES::`. The **breeding-program designer** is
now a separate product (`breedingDesigner`, DECISION-018); the *selection engine stays
here*.

### Context folder (`context/`)
- `context/isqg/` — C++ isqg source; meiosis/cross/DH being ported to Rust.
- `context/PleioArch-main/` — R prototype for exact-rhoG pleiotropy (`simulateEffects()`);
  the reference for the `"pleiotropy"` architecture (DECISION-007/013).

---

## The Rust Boundary (DECISION-006) — apply on every porting task

**Moves to Rust** (deterministic + profiled-slow only): `as_numeric()` numericalization;
isqg meiosis/recombination/cross/DH; genetic-value assembly *only if profiling proves it*.

**Stays in R** (parity-critical — never port): QTN sampling, effect-series & residual
draws, architecture QTN-to-trait mapping.

**Why:** R and Rust RNGs differ. R draws every random quantity and passes it into Rust;
Rust never calls an RNG on the parity-critical path (DECISION-012).

---

## Locked Decisions (source of truth: `docs/DECISIONS.md`)

| ID | Decision |
|----|----------|
| 001 | Single package, internal modularization |
| 002 | Port isqg algorithms to Rust (own the code) |
| 003/008 | `create_phenotypes()` = frozen legacy, not a delegation shim |
| 004 | Multi-generation stays in-package |
| 005 | Expression-based simulation → v3 |
| 006 | Rust surgical; stochastic core stays in R |
| 007/013 | PleioArch for `"pleiotropy"`, generalized to any `n_traits`; Cholesky fallback removed |
| 009 | Grammar owes v1 NO bit-parity; RDS references guard the legacy fn |
| 010 | Correlation control keeps the v1 name `cor` |
| 011 | Single rextendr package |
| 012 | Meiosis RNG drawn in R; Rust core pure; **exact** isqg bit-parity is the gate |
| 014 | `"ld"` = two-trait linked *distinct* causal loci (`direct` default) |
| 015 | Selection engine `select_ind()` + scheme wrappers + `c.Population()` |
| 016 | Modern methods: `g_matrix()`, `optimum_contribution()`+`sample_parents()`, `cross_usefulness()` |
| 017/018 | Designer split to separate product `breedingDesigner` |

---

## Public API (NAMESPACE exports)
`create_phenotypes()` (frozen) · `simulate_phenotype()` + `additive()`/`dominance()`/
`epistasis()`/`vqtl()` · `complex_phenotypes()` · crossing (`as_population`, `cross`,
`selfcross`, `double_haploid`, `synthetic_map`) · selection (`select_ind`,
`single_seed_descent`, `bulk`, `pedigree`, `recurrent_selection`, `c.Population`,
`g_matrix`, `optimum_contribution`, `sample_parents`, `cross_usefulness`) · `as_numeric()` ·
`format_conversion()` · `phenotypes_long/wide()` · `write_phenotypes()`.

---

## Build & Test

```bash
Rscript -e "devtools::load_all()"
Rscript -e "devtools::document()"
Rscript -e "devtools::test()"                                          # gate every commit
Rscript -e "testthat::test_file('tests/testthat/test-grammar.R')"      # single file
Rscript -e "testthat::test_file('tests/testthat/test-v130-parity.R')"  # gates releases
Rscript -e "rcmdcheck::rcmdcheck()"
cd src/rust && cargo test && cargo clippy -- -D warnings && cargo fmt  # Rust core
```
Deps: SNPRelate + gdsfmt are Bioconductor (`BiocManager::install()`).

---

## Conventions
- Git: no AI attribution; branches `<type>/<desc>`; commits
  `fix:`/`feat:`/`refactor:`/`port:`/`docs:`/`test:`.
- R: snake_case; internal funcs `@noRd`; v1 public params never renamed; every export tested.
- Rust: `thiserror`/`anyhow`; no `unwrap()` in lib code; no `unsafe` without `// SAFETY:`;
  no nightly; **no RNG on the parity-critical path**.
- Every selection/designer decision is mirrored in `docs/DECISIONS.md` and the designer
  manuscript so code and paper cannot drift.

---

## Avoid
- Don't port stochastic/RNG code to Rust (breaks parity — Rule #4).
- Don't add AI attribution; don't commit CLAUDE.md / .claude/.
- Don't break v1 `create_phenotypes()` signature.
- Don't invent citations/pages/equations; don't overclaim capabilities.
- Don't put logic in binding layers; no `unwrap()` in Rust lib code.
- Don't self-review a genetics change — hand it to the other model.
