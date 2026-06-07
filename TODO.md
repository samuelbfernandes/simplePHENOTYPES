# TODO — simplePHENOTYPES

> Items 1–6 are complete. `rcmdcheck::rcmdcheck()` passes with 0 errors, 0 warnings, 0 notes.
> Phase 0 workspace bootstrap (item 7) and new features (item 8) are next.

---

## Status legend

- `[ ]` — not started
- `[~]` — in progress
- `[x]` — done
- `[?]` — needs investigation
- `[-]` — rejected / won't fix

---

## 7 · Phase 0 workspace bootstrap

> Start here next session.

- `[ ]` Update `.gitignore` (add Rust entries)
- `[ ]` Create `.claudeignore`
- `[ ]` Create root `Cargo.toml` (workspace members: `core`, `py-pkg`)
- `[ ]` Create `core/` Rust library scaffold
- `[ ]` `git mv` R package files into `r-pkg/`; run `rextendr::use_extendr()`
- `[ ]` Create `py-pkg/` via `maturin new --bindings pyo3 py-pkg`
- `[ ]` Create `docs/SPEC.md`, `docs/BUGS.md`, `docs/DECISIONS.md` stubs
- `[ ]` Add PostToolUse hooks to `~/.claude/settings.json`

---

## 8 · New features (Phase 4 — implement in Rust after port)

These are queued for after the R → Rust port is complete. Do not implement in R.

### Core simulation

- `[ ]` Simulate null trait directly (no QTN effect)
- `[ ]` Provide QTN list for only one effect type (add/dom/epi/var) and randomize the other
- `[ ]` Proportion of variance explained by a subset of QTNs
- `[ ]` Simulate alleles in repulsion (repulsion phase LD)
- `[ ]` Different means per sub-population (PCA-based)
- `[ ]` Simulate independent traits with QTNs from different chromosomes
- `[ ]` Fix h2 and vary number of QTNs (discuss with Kaio)
- `[ ]` Implement GxE via envirotyping
- `[ ]` Implement GxB — add extra random variable to total genetic value (without changing `pheno()`)
- `[ ]` Simulate HTP: allelic effects must be correlated across environments
- `[ ]` Categorical phenotypes via threshold model: y=1 if l<γ1, y=2 if γ1<l<γ2, …, y=C if l>γ(C-1)
- `[ ]` Create `create_complex_phenotype()` — combine multiple architecture types
- `[ ]` Co-heritability in the package (rGh1h2)
- `[ ]` Simulate eQTL
- `[ ]` Haplotype-based simulation
- `[ ]` Genetic map integration (Kaio)
- `[ ]` Polyploid support — encode each homologous chromosome separately (0/1 per chromosome)

### Input / output

- `[ ]` Read `.gz` compressed genotype files
- `[ ]` GDS to VCF conversion (package `SeqArray`)
- `[ ]` Filter non-polymorphic SNPs
- `[ ]` Subset individuals
- `[ ]` Relative path support for `home_dir`

### Diagnostics & tooling

- `[ ]` Diagnostic plots (QQ, Manhattan, trait distributions)
- `[ ]` Unit tests (testthat)
- `[ ]` Power calculation
- `[ ]` Wrapper for GWAS tools: GAPIT, GEMMA

### Interfaces

- `[ ]` Shiny app
