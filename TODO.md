# TODO — simplePHENOTYPES

> **Rule**: Fix all bugs and resolve branch decisions BEFORE running the Phase 0 workspace
> bootstrap. This file is the gate that keeps broken code from being ported to Rust.
> **Keep this file updated** — mark each item `[x]` as soon as it is done.

---

## Status legend

- `[ ]` — not started
- `[~]` — in progress
- `[x]` — done
- `[?]` — needs investigation
- `[-]` — rejected / won't fix

---

## 1 · Branch strategy

- **`restructure`** is the active development branch. All work happens here.
- **`master`** stays frozen as the stable CRAN release until `restructure` is ready.
- When `restructure` passes all checks (section 5), open a PR to merge into `master`.

### `vQTL` branch

Only `vignettes/.gitignore` changes and compiled vignette files — nothing is broken.

- `[x]` Merge `vQTL` into `restructure` (or apply the `.gitignore` change manually)

---

## 2 · Restructure tasks

Split `create_phenotypes.R` into focused functions and expose internal utilities properly.
The public API (`create_phenotypes()` signature) must stay backward-compatible throughout.

### Extraction complete

- `[x]` `R/as_numeric.R` — numeric coercion extracted (internal, not exported)
- `[x]` `R/check_in.R` — pre-flight parameter validation extracted (internal)
- `[x]` `R/format_conversion.R` — format conversion utilities extracted (internal)
- `[x]` `R/handling_input_formats.R` — genotype input dispatch extracted (internal)
- `[x]` `R/table_to_numeric.R` — character table → numeric matrix extracted (internal)
- `[x]` `R/vQTL.R` — variance QTL simulation extracted (exported)

### Still to do

- `[ ]` Verify `create_phenotypes.R` no longer contains logic that belongs in the above files
- `[ ]` Confirm all extracted functions have roxygen docs and are listed correctly in NAMESPACE (exported or not)
- `[ ]` `devtools::check()` passes with no ERRORs or WARNINGs after all extractions

---

## 3 · Pending changes to review

Changes already partially implemented (in `restructure` or locally) that need review before merge.

- `[ ]` Option to simulate a null trait directly
- `[ ]` Put all effects for each trait in the same table as the QTNs
- `[ ]` Added `mean = c(1, 1)` parameter to `create_phenotypes()` call — verify behavior
- `[ ]` Removed mean-shift lines from output:
  - `pheno_ade_simple$Trait_1_H2_1 <- pheno_ade_simple$Trait_1_H2_1 + gv$gv_mu[1]`
  - `pheno_ade_simple$Trait_2_H2_1 <- pheno_ade_simple$Trait_2_H2_1 + gv$gv_mu[2]`

---

## 4 · Bugs to fix in R source

> Use `/bug-create` → `/bug-analyze` → `/bug-fix` → `/bug-verify` for each.
> Do NOT port any function that still has an open bug.

| # | File / function | Description | Status |
| - | --------------- | ----------- | ------ |
| 1 | `R/QTN_linkage.R` | `ldmax`/`ldmin` bounds not enforced correctly | `[x]` fixed in `restructure` |
| 2 | `R/file_loader.R` | GDS file reading error (`snp.rs.id` vs `snp.id`) | `[x]` fixed in `restructure` |
| 3 | `R/Genotypes.R` | Genotype names incorrect in output | `[x]` fixed in `restructure` |
| 4 | `R/file_loader.R` | Numericalizing HapMap with numeric columns: skips numericalization when column 12 is already numeric | `[x]` |
| 5 | `R/Phenotypes.R` | Correlation not printed when multiple h2 simulated with ntraits > 1 | `[x]` |
| 6 | `R/file_loader.R` | Multi-file log always said "HapMap files" — BED/PED/GDS had no message | `[x]` |
| 7 | `R/file_loader.R` | `out_name` included full path when using `geno_path`, breaking output file names | `[x]` |
| 8 | `R/Phenotypes.R` | `phenotypes/` folder silently overwrites files on re-run | `[x]` |

---

## 5 · Housekeeping

- `[ ]` Add citation to GitHub repo (CITATION file / README badge)
- `[ ]` Include timestamp (date/time) in simulation log output

---

## 6 · Pre-port checklist (gate)

All items below must be checked before starting Phase 0 bootstrap.

- `[ ]` Branch decisions (section 1) are resolved
- `[ ]` Restructure tasks (section 2) are complete
- `[ ]` Pending changes (section 3) are reviewed and merged or rejected
- `[ ]` All bugs in section 4 are marked `[x]` or `[-]`
- `[ ]` `devtools::test()` passes on `restructure`
- `[ ]` `rcmdcheck::rcmdcheck()` reports no ERRORs or WARNINGs

---

## 7 · Phase 0 workspace bootstrap

> Only start after section 6 is fully checked.

- `[ ]` Update `.gitignore` (add Rust entries)
- `[ ]` Create `.claudeignore`
- `[ ]` Create root `Cargo.toml` (workspace members: `core`, `py-pkg`)
- `[ ]` Create `core/` Rust library scaffold
- `[ ]` `git mv` R package files into `r-pkg/`; run `rextendr::use_extendr()`
- `[ ]` Create `py-pkg/` via `maturin new --bindings pyo3 py-pkg`
- `[ ]` Create `docs/SPEC.md`, `docs/BUGS.md`, `docs/DECISIONS.md` stubs
- `[ ]` Add PostToolUse hooks to `~/.claude/settings.json`

---

## 7 · New features (Phase 4 — implement in Rust after port)

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
