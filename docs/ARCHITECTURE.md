# ARCHITECTURE.md — simplePHENOTYPES v2 Redesign

> Status: APPROVED — decisions finalized; Rust scope narrowed to bottlenecks only
> Working directory: ~/Library/CloudStorage/OneDrive-UniversityofArkansas/UARK/collaboration/Software/simplePHENOTYPES
> Strategic goals: maximize citations, maximize adoption, preserve v1.3.0 parity

---

## 1. Project Background

simplePHENOTYPES simulates pleiotropic, linked, and epistatic phenotypes from real
genomic data (Fernandes & Lipka, BMC Bioinformatics 2020). v2 introduces a composable,
tidyverse-style grammar (see SPEC.md), adds multi-generation cross simulation
(isqg-style algorithms), and keeps a single CRAN deliverable — while reproducing
v1.3.0 behavior exactly under matching architectures.

---

## 2. Locked Decisions

| ID | Decision |
|----|----------|
| 001 | Single package, internal modularization |
| 002 | Port isqg's algorithms to Rust (matching the published version) |
| 003 | `create_phenotypes()` preserved (frozen legacy — refined by 008) |
| 004 | Multi-generation stays in simplePHENOTYPES |
| 005 | Expression-based simulation deferred to v3 |
| 006 | **Rust is surgical, not a rewrite** — only C++-origin code and profiled bottlenecks move to Rust; the stochastic simulation core stays in R |
| 007 | PleioArch adopted for `"pleiotropy"` — exact rhoG control via bivariate-normal draws |
| 008 | `create_phenotypes()` = frozen legacy (bugfix-only, superseded), NOT a delegation shim |
| 009 | New grammar has NO bit-for-bit v1 parity; RDS references guard the legacy fn + statistical checks for the grammar |
| 010 | v1 `cor` dropped; genetic correlation reimplemented as `rho_g` via PleioArch |
| 011 | Single rextendr package; Cargo-workspace plan (TODO_newfeatures item 7) dropped |
| 012 | Meiosis randomness drawn in R; Rust core pure; **exact** isqg bit-parity is the gate |

> Full rationale for every decision lives in `DECISIONS.md` (canonical log).

---

## 3. The Rust Boundary (DECISION-006)

This is the most important architectural principle for v2 and it overrides the earlier,
broader "port function by function" framing.

### What moves to Rust
- **isqg-origin C++ code** — meiosis, recombination, crossing, double-haploid
  generation. This is already C++; porting to Rust removes the archived-package
  dependency and unifies the codebase.
- **Numericalization** — `as_numeric()` (the genotype → −1/0/1 conversion). This is a
  measured bottleneck that scales with dataset dimensions and is purely deterministic.
- **Profiled hot loops only** — any other function proven slow by benchmarking on a
  large dataset, and only if it is deterministic.

### What stays in R (parity-critical — must NOT move)
- **QTN sampling** — which markers become causal loci.
- **Effect-size generation** — the geometric series and any custom effect series.
- **Residual / environmental draws** — the noise added to reach the target variance.
- **Architecture logic** — pleiotropy / LD / independent QTN-to-trait mapping.

### Why this split
Exact v1.3.0 parity (SPEC.md §8) pins the random-number stream. R and Rust use
different RNGs, so any stochastic step performed in Rust would break bit-for-bit
reproducibility of QTN selection and phenotype values. Keeping all RNG-dependent steps
on R's generator guarantees parity; Rust handles only deterministic numeric work, where
results are identical regardless of language.

### The contract between layers
```
R (stochastic + orchestration)                Rust (deterministic, hot)
  ├─ select QTNs (set.seed)                      ├─ as_numeric(): geno → -1/0/1
  ├─ draw effect series (geometric/custom)       ├─ genetic-value assembly (matrix ops)
  ├─ draw residuals                              └─ meiosis / recombination / crossing
  └─ assemble phenotype_sim, call Rust for ─────►   (isqg port)
     the heavy deterministic pieces
```
R passes already-drawn random quantities (QTN indices, effects, residuals) into Rust;
Rust never calls a random generator on the parity-critical path.

---

## 4. Internal Package Structure (single package)

> **Module convention (implementation reality):** a standard R package build sources
> only top-level `R/*.R` files — it does **not** recurse into `R/grammar/`,
> `R/architectures/`, etc. The "modules" below are therefore expressed as **filename
> prefixes in a flat `R/`**, not subdirectories. Implemented mapping:
> `grammar_*.R` (simulate_phenotype, layers, complex, realize),
> `arch_*.R` (independent / ld / pleiotropy QTN logic),
> `effects_*.R` (series, residual, pleioarch),
> `io_*.R` (exporters). The legacy engine keeps its historical filenames
> (`create_phenotypes.R`, `base_line_*`, `QTN_*`, …). The tree below is the logical
> grouping; read each `←` as a filename prefix.

```
simplePHENOTYPES/
├── R/
│   ├── grammar_*            ← simulate_phenotype(), additive(), dominance(),
│   │                         epistasis(), vqtl(), complex_phenotypes() (one-call
│   │                         folded into simulate_phenotype(); no sim_phenotypes())
│   ├── arch_*              ← pleiotropy / ld / independent QTN-to-trait logic (R, parity-critical)
│   ├── effects_*           ← effect-series + residual + PleioArch draws (R, parity-critical)
│   ├── (legacy filenames)  ← create_phenotypes() frozen v1 engine (bugfix-only, superseded)
│   ├── cross_*             ← R wrappers around Rust meiosis/cross/dh (later)
│   ├── io_*                ← write_*() exporters (long default; wide/gemma/plink/multi-file)
│   └── (utils)             ← as_numeric() R wrapper, check_in, helpers
├── src/                    ← Rust via rextendr (CRAN-required location)
│   └── rust/src/
│       ├── lib.rs
│       ├── numeric.rs      ← as_numeric() core (deterministic)
│       ├── genome.rs       ← bitwise chromosome representation (isqg port)
│       ├── meiosis.rs      ← recombination / crossing / DH (isqg port)
│       └── gvalue.rs       ← genetic-value assembly (deterministic matrix ops)
├── inst/extdata/
│   ├── v1_3_0_reference/   ← frozen QTNs + phenotypes from v1.3.0 (big-QTN removed)
│   └── isqg_v1_outputs/    ← isqg parity references
├── tests/testthat/
│   ├── test-v130-parity.R  ← identical QTNs + identical phenotypes vs v1.3.0
│   ├── test-isqg-parity.R
│   ├── test-as-numeric.R   ← Rust as_numeric() == R reference
│   └── test-grammar.R      ← variance identity, seed threading, complex weighting
├── vendor/                 ← vendored Rust deps for CRAN
├── python/                 ← PyO3 bindings (later; in .Rbuildignore)
└── docs/
    ├── ARCHITECTURE.md      ├── DECISIONS.md      ├── SPEC.md
    ├── BUGS.md              └── NEXT_STEPS.md
```

---

## 5. Rust Port Candidates (revised — bottlenecks only)

| Priority | Target | Type | RNG? | Move to Rust? |
|----------|--------|------|------|---------------|
| 1 | `as_numeric()` numericalization | deterministic | no | **Yes** — measured bottleneck, already started |
| 2 | isqg meiosis / recombination / cross / DH | deterministic C++ | no | **Yes** — already C++, removes archived dep |
| 3 | genetic-value assembly (matrix ops) | deterministic | no | **Yes, if profiling shows it matters** |
| — | QTN sampling | stochastic | yes | **No** — parity-critical, stays in R |
| — | effect-series / residual draws | stochastic | yes | **No** — parity-critical, stays in R |
| — | architecture QTN-to-trait mapping | logic | partial | **No** unless deterministic & profiled |

Rule: a function moves to Rust only if it is (a) deterministic and (b) demonstrated slow
by profiling on a large dataset. Stochastic steps never move.

---

## 6. isqg Algorithm Port

Port isqg's published C++ meiosis/cross/DH algorithms to Rust (`genome.rs`,
`meiosis.rs`), validated against captured isqg reference outputs in
`inst/extdata/isqg_v1_outputs/`.

**Exact bit-for-bit parity with isqg is the test gate (DECISION-012)** — not statistical
equivalence. This is achievable because isqg draws from R's own RNG, so R can make every
draw in isqg's exact order and pass the results into a Rust core that never calls an RNG.
That satisfies DECISION-006 and the strongest parity standard simultaneously; see
DECISION-012 for the verification behind it.

isqg's C++ source is kept in-repo at `context/isqg/` for traceability (it is
`.Rbuildignore`d, so it does not affect the CRAN bundle).

---

## 7. Parity Strategy Summary

- v1.3.0 parity (SPEC.md §8) is guaranteed because all RNG-dependent steps stay in R.
- `as_numeric()` in Rust is validated against the R numericalization on multiple
  datasets (`test-as-numeric.R`); since it is deterministic, exact equality is required.
- The Rust meiosis port is validated against isqg references.
- Capture the v1.3.0 references (QTNs + phenotypes, big-QTN removed) BEFORE writing new
  code, so the parity harness exists from day one.

---

## 8. CRAN Constraints
- Rust deps vendored (`rextendr::vendor_pkgs()`); `SystemRequirements: Cargo, rustc`.
- Bundle under 5MB; builds offline; no nightly Rust.
- `python/`, `docs/`, `.claude/`, `CLAUDE.md` excluded via `.Rbuildignore`.

---

## 9. Open Questions

- None blocking. O6 (RNG strategy) is resolved by DECISION-006: stochastic steps stay
  in R. Revisit only if profiling reveals a stochastic step is itself a bottleneck — in
  which case the fix is a faster R approach or pre-drawing in bulk, not moving RNG to Rust.
