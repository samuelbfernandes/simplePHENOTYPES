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
| 007 | PleioArch adopted for `"pleiotropy"` — rhoG control (target `cor`) via bivariate-normal draws |
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
> `io_*.R` (genotype readers, format detection/conversion, `as_numeric()`, exporters),
> `qc_*.R` (`filter_geno()` and its LD methods),
> `cross_*.R` (Population, map, meiosis wrappers, pedigree, mating plans, crossbreeding),
> `select_*.R` (selection engine: `select_ind()`, schemes, marker selection, MABC, OCS,
> usefulness, combining ability, progeny testing, BLUP),
> `transcriptome_*.R` (transcriptome generator, layer, counts, mimic),
> `legacy_*.R` (the frozen `create_phenotypes()` engine). Files that R tooling names by
> convention keep their names: `extendr-wrappers.R` (generated by rextendr), `data.R`,
> `simplePHENOTYPES-package.R`, `zzz.R`. The tree below is the logical grouping; read
> each `←` as a filename prefix.

```
simplePHENOTYPES/
├── R/
│   ├── grammar_*            ← simulate_phenotype(), additive(), dominance(),
│   │                         epistasis(), vqtl(), complex_phenotypes() (one-call
│   │                         folded into simulate_phenotype(); no sim_phenotypes())
│   ├── arch_*              ← pleiotropy / ld / independent QTN-to-trait logic (R, parity-critical)
│   ├── effects_*           ← effect-series + residual + PleioArch draws (R, parity-critical)
│   ├── legacy_*            ← create_phenotypes() frozen v1 engine (bugfix-only, superseded)
│   ├── cross_*             ← Population, genetic map, R wrappers around Rust meiosis/cross/dh,
│   │                         pedigree, mating plans, crossbreeding
│   ├── select_*            ← selection engine (select_ind, schemes, marker_select, MABC,
│   │                         OCS/g_matrix, usefulness, combining ability, progeny test, BLUP)
│   ├── transcriptome_*     ← genome → transcriptome → phenotype layer
│   ├── io_*                ← readers + format detection/conversion + as_numeric() wrapper;
│   │                         write_phenotypes() exporter (long default, or wide)
│   ├── qc_*                ← filter_geno() and its PLINK-parity LD methods
│   └── (conventional)      ← extendr-wrappers.R, data.R, simplePHENOTYPES-package.R, zzz.R
├── src/                    ← Rust via rextendr (CRAN-required location)
│   └── rust/src/
│       ├── lib.rs
│       ├── numeric.rs      ← as_numeric() core (deterministic)
│       ├── genome.rs       ← bitwise chromosome representation (isqg port)
│       ├── meiosis.rs      ← recombination / crossing / DH (isqg port)
│       └── hash.rs         ← FNV-1a-128 content hash for pedigree keys (DECISION-024)
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

### 5.1 Function inventory

Audit of every top-level function in `R/`, 2026-09-28. The table was generated by
parsing the sources, not written by hand: 274 functions (54 exported, 9 registered S3
methods, 211 internal). Call edges come from a static call graph
(`codetools::findGlobals`).

Column meanings:

- **Exported**: `yes` means `export()` in `NAMESPACE`, `S3` means a registered S3 method,
  and `no` means internal.
- **Tests**: `direct` means `tests/testthat/` calls the function by name (including through
  `:::`). `indirect` means no test names it, but a directly tested function reaches it; an S3
  method counts as reached when a test calls its generic. `no` means nothing reaches it.
- **Module**: the §4 filename prefix. `rust (generated)` is `extendr-wrappers.R`, and
  `package` covers the conventionally named package files.
- **Rust candidate**, per the §5 rule:
  - `is Rust` is the generated extendr wrapper.
  - `done` means the function already calls the named Rust core for its deterministic
    step. Any random draws stay in R.
  - `no — RNG` means it draws random numbers itself or through a callee (DECISION-006).
  - `no — frozen legacy` is the `create_phenotypes()` engine (DECISION-008).
  - `profile first` is a deterministic numeric kernel (matrix algebra, LD/EM, pedigree
    recursion, dosage transforms). It moves only if profiling on a large dataset shows a
    bottleneck.
  - `no — glue` covers validation, argument resolution, printing and bookkeeping.
- **Parity-critical** means an exact-parity test pins the output:
  - `v1.3.0`: on the `create_phenotypes()` call path (`test-v130-parity.R`).
  - `PLINK`: on the `filter_geno()` LD path (`test-filter-geno-plink-parity.R`).
  - `isqg`: pinned to isqg 1.4 by `test-isqg-parity.R`. This covers R's draw order
    (`.draw_meiosis`) and the Rust cores fed the recorded draws: `meiosis_core`,
    `gamete_masks_core`, and `mate_haplotypes_core`, which `cross()`, `selfcross()` and
    `double_haploid()` run. `mate_haplotypes_core` is checked on phased haplotypes, so a
    cis/trans swap cannot pass.
  - `—`: no bit-parity obligation. The v2 grammar owes v1 none (DECISION-009).

Summary:

| | Count |
|---|---|
| Tests: direct / indirect / none | 96 / 177 / 1 (`.onAttach`, a package hook) |
| Rust: is Rust / done / RNG / frozen legacy / profile first / glue | 5 / 3 / 52 / 14 / 41 / 159 |
| Parity: v1.3.0 / PLINK / isqg / none | 36 / 12 / 4 / 222 |

Notes:

- The audit found one dead function, `.causal_loci()`, and removed it. The only functions
  that no export or S3 method reaches are now `meiosis_core` and `gamete_masks_core`.
  They are kept as the isqg parity gate's entry points.
- Only three R functions call a Rust core: `.mate`, `.apply_coding` and `.stable_key`. No
  deterministic kernel has been profiled yet, so every `profile first` row is a candidate,
  not a planned port.

<details>
<summary>Full inventory (274 functions)</summary>

| Function | File | Exported | Tests | Module | Rust candidate | Parity-critical |
|---|---|---|---|---|---|---|
| `complex_phenotypes` | grammar_complex.R | yes | direct | grammar | no — RNG | — |
| `additive` | grammar_layers.R | yes | direct | grammar | no — RNG | — |
| `.orthogonal_d_series` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `dominance` | grammar_layers.R | yes | direct | grammar | no — RNG | — |
| `epistasis` | grammar_layers.R | yes | direct | grammar | no — RNG | — |
| `vqtl` | grammar_layers.R | yes | direct | grammar | no — RNG | — |
| `.check_sim` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.resolve_n_qtn` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.last_layer_of_type` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.one_call_hint` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.resolve_prop` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.add_layer` | grammar_layers.R | no | indirect | grammar | no — RNG | — |
| `.draw_layer` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.dom_hetless` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.dom_partial_hetless` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.orthogonal_hetless_d` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.epi_hetless_d` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.rep_qtn` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.resolve_qtn_arg` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.resolve_epi_qtn` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.resolve_interaction_type` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.cite_vqtl` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.apply_phase` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `.require_markers` | grammar_layers.R | no | indirect | grammar | no — glue | — |
| `plot.phenotype_sim` | grammar_plot.R | S3 | indirect | grammar | no — glue | — |
| `.var_budget_matrix` | grammar_plot.R | no | direct | grammar | no — glue | — |
| `.plot_variance` | grammar_plot.R | no | direct | grammar | no — glue | — |
| `.plot_hist` | grammar_plot.R | no | indirect | grammar | no — glue | — |
| `.plot_effects` | grammar_plot.R | no | indirect | grammar | no — glue | — |
| `.plot_cor` | grammar_plot.R | no | indirect | grammar | no — glue | — |
| `.realize_phenotype` | grammar_realize.R | no | indirect | grammar | no — RNG | — |
| `.genetic_matrix` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.avg_effect` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.breeding_value_matrix` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.layer_scaled_effects` | grammar_realize.R | no | indirect | grammar | profile first | — |
| `.genetic_value_matrix` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.transcriptome_matrix` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.tx_raw` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.component_raw` | grammar_realize.R | no | direct | grammar | profile first | — |
| `.epi_unit_column` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.apply_vqtl` | grammar_realize.R | no | indirect | grammar | no — RNG | — |
| `.layer_qtn_effect` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.validate_rep` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.variance_budget` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.mediation_budget` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.orthogonal_var_split` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.seeded_residual` | grammar_realize.R | no | indirect | grammar | no — RNG | — |
| `.Random.seed_safe` | grammar_realize.R | no | direct | grammar | no — glue | — |
| `.restore_seed` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.realized_h2` | grammar_realize.R | no | direct | grammar | no — glue | — |
| `.check_h2_complete` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `.trait_mean` | grammar_realize.R | no | indirect | grammar | no — glue | — |
| `simulate_phenotype` | grammar_simulate_phenotype.R | yes | direct | grammar | no — RNG | — |
| `.check_arch_args` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.validate_count` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.validate_flag` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.validate_seed` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.validate_proportion` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.is_one_call` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.build_one_call` | grammar_simulate_phenotype.R | no | indirect | grammar | no — RNG | — |
| `.normalize_geno` | grammar_simulate_phenotype.R | no | direct | grammar | no — glue | — |
| `.select_individuals` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.expression_foundation` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.geno_cols` | grammar_simulate_phenotype.R | no | direct | grammar | no — glue | — |
| `.marker_maf_ref` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.layer_seed` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `print.phenotype_sim` | grammar_simulate_phenotype.R | S3 | indirect | grammar | no — glue | — |
| `.total_genetic_prop` | grammar_simulate_phenotype.R | no | direct | grammar | no — glue | — |
| `.total_variance_prop` | grammar_simulate_phenotype.R | no | direct | grammar | no — glue | — |
| `.expand_prop` | grammar_simulate_phenotype.R | no | indirect | grammar | no — glue | — |
| `.draw_qtn` | arch_independent.R | no | indirect | arch | no — RNG | — |
| `.draw_qtn_distinct_chr` | arch_independent.R | no | indirect | arch | no — RNG | — |
| `.draw_qtn_pairs` | arch_independent.R | no | indirect | arch | no — RNG | — |
| `.candidate_markers` | arch_independent.R | no | indirect | arch | profile first | — |
| `.type_occurrence` | arch_independent.R | no | indirect | arch | no — glue | — |
| `.draw_qtn_ld` | arch_ld.R | no | indirect | arch | no — RNG | — |
| `.pleio_draw` | effects_pleioarch.R | no | indirect | effects | no — RNG | — |
| `.pleio_check_zero_var` | effects_pleioarch.R | no | indirect | effects | no — glue | — |
| `.pleio_partition` | effects_pleioarch.R | no | indirect | effects | no — glue | — |
| `.pleio_single_unit_consequence` | effects_pleioarch.R | no | direct | effects | no — glue | — |
| `.pleio_nonadditive_draw` | effects_pleioarch.R | no | indirect | effects | no — RNG | — |
| `.pleio_units` | effects_pleioarch.R | no | indirect | effects | no — RNG | — |
| `.pleio_unit_effects` | effects_pleioarch.R | no | direct | effects | no — RNG | — |
| `.pleio_total_cor_check` | effects_pleioarch.R | no | direct | effects | no — glue | — |
| `.cite_pleioarch` | effects_pleioarch.R | no | indirect | effects | no — glue | — |
| `.cite_main` | effects_pleioarch.R | no | indirect | effects | no — glue | — |
| `.pleio_pi_vector` | effects_pleioarch.R | no | indirect | effects | no — glue | — |
| `.check_pleio_feasible` | effects_pleioarch.R | no | direct | effects | no — glue | — |
| `.draw_mvnorm` | effects_pleioarch.R | no | indirect | effects | no — RNG | — |
| `.draw_univariate` | effects_pleioarch.R | no | indirect | effects | no — RNG | — |
| `.pleio_cor_matrix` | effects_pleioarch.R | no | indirect | effects | profile first | — |
| `.effect_series` | effects_series.R | no | direct | effects | no — glue | — |
| `.draw_residual` | effects_series.R | no | indirect | effects | no — RNG | — |
| `as_numeric` | io_as_numeric.R | yes | direct | io | no — glue | — |
| `.het_to_letters` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `.call_to_letters` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `.gds_needed` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `detect_format` | io_detect_format.R | no | direct | io | no — glue | v1.3.0 |
| `.sniff_text_signature` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `.detect_df_format` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `parse_hapmap_chars_to_raw` | io_detect_format.R | no | direct | io | profile first | v1.3.0 |
| `compute_flip` | io_detect_format.R | no | direct | io | profile first | v1.3.0 |
| `assemble_output` | io_detect_format.R | no | indirect | io | no — glue | v1.3.0 |
| `format_conversion` | io_format_conversion.R | no | direct | io | no — glue | v1.3.0 |
| `.apply_coding` | io_read_formats.R | no | indirect | io | done (`numericalize_core`) | v1.3.0 |
| `handle_hapmap` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_numeric` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_table` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `.read_gds_to_raw` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_vcf` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_gds` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `.plink_cm` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_bed` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_ped` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `handle_finalreport` | io_read_formats.R | no | indirect | io | no — glue | v1.3.0 |
| `phenotypes_long` | io_write.R | yes | direct | io | no — glue | — |
| `phenotypes_wide` | io_write.R | yes | direct | io | no — glue | — |
| `write_phenotypes` | io_write.R | yes | direct | io | no — glue | — |
| `genetic_values` | io_write.R | yes | direct | io | no — glue | — |
| `mediation_split` | io_write.R | yes | direct | io | no — glue | — |
| `.qtn_var` | io_write.R | no | indirect | io | no — glue | — |
| `.tx_qtn_var` | io_write.R | no | indirect | io | no — glue | — |
| `qtn_table` | io_write.R | yes | direct | io | no — glue | — |
| `filter_geno` | qc_filter_geno.R | yes | direct | qc | no — glue | PLINK |
| `.ld_spec` | qc_filter_geno.R | no | indirect | qc | no — glue | PLINK |
| `.ld_prune` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.ld_sweep` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.plink_hap_rsq` | qc_filter_geno.R | no | direct | qc | profile first | PLINK |
| `.plink_em_hethet` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.plink_calc_lnlike` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.plink_cubic_roots` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.plink_calc_lnlike_quantile` | qc_filter_geno.R | no | indirect | qc | profile first | PLINK |
| `.plink_blocks_classify` | qc_filter_geno.R | no | indirect | qc | no — glue | PLINK |
| `.plink_blocks_chrom` | qc_filter_geno.R | no | direct | qc | profile first | PLINK |
| `.cite_ld` | qc_filter_geno.R | no | indirect | qc | no — glue | — |
| `.gabriel_blocks` | qc_ld_methods.R | no | indirect | qc | profile first | PLINK |
| `breed_composition` | cross_breed.R | yes | direct | cross | profile first | — |
| `heterosis` | cross_breed.R | yes | direct | cross | profile first | — |
| `.check_breeds` | cross_breed.R | no | indirect | cross | no — glue | — |
| `crossbreed` | cross_breed.R | yes | direct | cross | no — RNG | — |
| `synthetic_map` | cross_map.R | yes | direct | cross | no — glue | — |
| `.recycle_by_chr` | cross_map.R | no | indirect | cross | no — glue | — |
| `mate` | cross_mate.R | yes | direct | cross | no — RNG | — |
| `.check_plan` | cross_mate.R | no | indirect | cross | no — glue | — |
| `mating_design` | cross_mate.R | yes | direct | cross | no — RNG | — |
| `.cite_isqg` | cross_mating.R | no | indirect | cross | no — glue | — |
| `.draw_meiosis` | cross_mating.R | no | direct | cross | no — RNG | isqg |
| `.mate` | cross_mating.R | no | indirect | cross | done (`mate_haplotypes_core`) | — |
| `cross` | cross_mating.R | yes | direct | cross | no — RNG | — |
| `selfcross` | cross_mating.R | yes | direct | cross | no — RNG | — |
| `double_haploid` | cross_mating.R | yes | direct | cross | no — RNG | — |
| `.stable_key` | cross_pedigree.R | no | direct | cross | done (`stable_hash_core`) | — |
| `.founder_pedigree` | cross_pedigree.R | no | indirect | cross | no — glue | — |
| `.ensure_pedigree` | cross_pedigree.R | no | indirect | cross | no — glue | — |
| `.pedigree_union` | cross_pedigree.R | no | indirect | cross | no — glue | — |
| `.pedigree_ancestors` | cross_pedigree.R | no | indirect | cross | profile first | — |
| `.mating_pedigree` | cross_pedigree.R | no | indirect | cross | no — glue | — |
| `.pedigree_relabel` | cross_pedigree.R | no | indirect | cross | no — glue | — |
| `parentage` | cross_pedigree.R | yes | direct | cross | no — glue | — |
| `families` | cross_pedigree.R | yes | direct | cross | no — glue | — |
| `as_population` | cross_population.R | yes | direct | cross | no — glue | — |
| `.new_population` | cross_population.R | no | indirect | cross | no — glue | — |
| `.check_map` | cross_population.R | no | indirect | cross | no — glue | — |
| `n_individuals` | cross_population.R | yes | direct | cross | no — glue | — |
| `[.Population` | cross_population.R | S3 | indirect | cross | no — glue | — |
| `dosages` | cross_population.R | yes | direct | cross | profile first | — |
| `.resolve_geno_qtn` | cross_population.R | no | indirect | cross | no — glue | — |
| `additive_value` | cross_population.R | yes | direct | cross | profile first | — |
| `genotypic_value` | cross_population.R | yes | direct | cross | profile first | — |
| `phenotype_value` | cross_population.R | yes | direct | cross | no — RNG | — |
| `print.Population` | cross_population.R | S3 | indirect | cross | no — glue | — |
| `.check_population` | cross_population.R | no | indirect | cross | no — glue | — |
| `.check_single` | cross_population.R | no | indirect | cross | no — glue | — |
| `a_matrix` | select_blup.R | yes | direct | select | profile first | — |
| `predict_ebv` | select_blup.R | yes | direct | select | profile first | — |
| `.check_relationship` | select_blup.R | no | indirect | select | no — glue | — |
| `.blup_variances` | select_blup.R | no | indirect | select | no — glue | — |
| `.gblup_marker_effects` | select_blup.R | no | indirect | select | profile first | — |
| `prediction_accuracy` | select_blup.R | yes | direct | select | no — glue | — |
| `selection_methods` | select_blup.R | yes | direct | select | no — glue | — |
| `combining_ability` | select_combining.R | yes | direct | select | no — RNG | — |
| `.check_distinct` | select_combining.R | no | indirect | select | no — glue | — |
| `.check_ad` | select_combining.R | no | indirect | select | no — glue | — |
| `.expected_cross_means` | select_combining.R | no | indirect | select | profile first | — |
| `.simulate_cross_means` | select_combining.R | no | indirect | select | no — RNG | — |
| `.decompose_ca` | select_combining.R | no | indirect | select | profile first | — |
| `print.combining_ability` | select_combining.R | S3 | indirect | select | no — glue | — |
| `template_effects` | select_combining.R | yes | direct | select | no — glue | — |
| `select_ind` | select_ind.R | yes | direct | select | no — RNG | — |
| `.resolve_keep` | select_ind.R | no | indirect | select | no — glue | — |
| `.criterion_values` | select_ind.R | no | indirect | select | no — glue | — |
| `.quadratic_index_score` | select_ind.R | no | direct | select | no — glue | — |
| `.index_score` | select_ind.R | no | direct | select | no — glue | — |
| `.index_weights` | select_ind.R | no | direct | select | no — glue | — |
| `.combined_score` | select_ind.R | no | direct | select | no — glue | — |
| `.sel_top` | select_ind.R | no | indirect | select | no — glue | — |
| `.sel_within_family` | select_ind.R | no | indirect | select | no — glue | — |
| `.sel_among_family` | select_ind.R | no | indirect | select | no — glue | — |
| `.selection_result` | select_ind.R | no | indirect | select | no — glue | — |
| `.select_culling` | select_ind.R | no | indirect | select | no — glue | — |
| `mabc_select` | select_mabc.R | yes | direct | select | no — RNG | — |
| `recurrent_parent_recovery` | select_mabc.R | yes | direct | select | profile first | — |
| `.mabc_founders` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_markers` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_require_informative` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_check_flanks` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_parse_interval` | select_mabc.R | no | direct | select | no — glue | — |
| `.mabc_interval` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_recovery` | select_mabc.R | no | indirect | select | profile first | — |
| `.mabc_check_weights` | select_mabc.R | no | indirect | select | no — glue | — |
| `.mabc_weights` | select_mabc.R | no | direct | select | no — glue | — |
| `marker_select` | select_marker.R | yes | direct | select | no — RNG | — |
| `g_matrix` | select_ocs.R | yes | direct | select | profile first | — |
| `optimum_contribution` | select_ocs.R | yes | direct | select | profile first | — |
| `print.ocs` | select_ocs.R | S3 | indirect | select | no — glue | — |
| `sample_parents` | select_ocs.R | yes | direct | select | no — RNG | — |
| `.validate_coancestry_matrix` | select_ocs.R | no | indirect | select | no — glue | — |
| `.dosage_from` | select_ocs.R | no | indirect | select | profile first | — |
| `.frank_wolfe` | select_ocs.R | no | indirect | select | profile first | — |
| `.tune_lambda` | select_ocs.R | no | indirect | select | no — glue | — |
| `progeny_test` | select_progeny.R | yes | direct | select | no — RNG | — |
| `c.Population` | select_schemes.R | S3 | indirect | select | no — glue | — |
| `single_seed_descent` | select_schemes.R | yes | direct | select | no — RNG | — |
| `bulk` | select_schemes.R | yes | direct | select | no — RNG | — |
| `pedigree` | select_schemes.R | yes | direct | select | no — RNG | — |
| `recurrent_selection` | select_schemes.R | yes | direct | select | no — RNG | — |
| `.check_tandem` | select_schemes.R | no | indirect | select | no — glue | — |
| `.as_founder_pop` | select_schemes.R | no | indirect | select | no — glue | — |
| `.check_phenotyper` | select_schemes.R | no | indirect | select | no — glue | — |
| `.self_each` | select_schemes.R | no | indirect | select | no — RNG | — |
| `.intermate` | select_schemes.R | no | indirect | select | no — RNG | — |
| `.relabel` | select_schemes.R | no | direct | select | no — glue | — |
| `cross_usefulness` | select_usefulness.R | yes | direct | select | no — RNG | — |
| `.intensity_from_p` | select_usefulness.R | no | indirect | select | no — glue | — |
| `.additive_model` | select_usefulness.R | no | indirect | select | no — glue | — |
| `.additive_gv` | select_usefulness.R | no | indirect | select | profile first | — |
| `.resolve_pairs` | select_usefulness.R | no | indirect | select | no — glue | — |
| `.make_family` | select_usefulness.R | no | indirect | select | no — RNG | — |
| `observe_counts` | transcriptome_counts.R | yes | direct | transcriptome | no — RNG | — |
| `.tx_count_par` | transcriptome_counts.R | no | indirect | transcriptome | no — glue | — |
| `.attach_expression` | transcriptome_layer.R | no | indirect | transcriptome | no — RNG | — |
| `transcriptome` | transcriptome_layer.R | yes | direct | transcriptome | no — RNG | — |
| `.tx_cor_design_warning` | transcriptome_layer.R | no | indirect | transcriptome | no — glue | — |
| `.tx_grm` | transcriptome_mimic.R | no | direct | transcriptome | profile first | — |
| `.greml_h2` | transcriptome_mimic.R | no | direct | transcriptome | profile first | — |
| `.tx_estimate_factors` | transcriptome_mimic.R | no | indirect | transcriptome | profile first | — |
| `.tx_estimate_kappa` | transcriptome_mimic.R | no | direct | transcriptome | no — glue | — |
| `.tx_mimic_calibrate` | transcriptome_mimic.R | no | indirect | transcriptome | no — glue | — |
| `simulate_transcriptome` | transcriptome_simulate.R | yes | direct | transcriptome | no — RNG | — |
| `predict.transcriptome_sim` | transcriptome_simulate.R | S3 | indirect | transcriptome | no — RNG | — |
| `.tx_pergene` | transcriptome_simulate.R | no | indirect | transcriptome | no — RNG | — |
| `.tx_synthetic_coords` | transcriptome_simulate.R | no | indirect | transcriptome | no — RNG | — |
| `.tx_check_annotation` | transcriptome_simulate.R | no | indirect | transcriptome | no — glue | — |
| `print.transcriptome_sim` | transcriptome_simulate.R | S3 | indirect | transcriptome | no — glue | — |
| `base_line_multi_traits` | legacy_Base_line_multi_traits.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `base_line_single_trait` | legacy_Base_line_single_trait.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `genotypes` | legacy_Genotypes.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `phenotypes` | legacy_Phenotypes.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `qtn_linkage` | legacy_QTN_linkage.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `qtn_partially_pleiotropic` | legacy_QTN_partially_pleiotropic.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `qtn_pleiotropic` | legacy_QTN_pleiotropic.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `check_in` | legacy_check_in.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `constraint` | legacy_constraint.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `create_phenotypes` | legacy_create_phenotypes.R | yes | direct | legacy | no — frozen legacy | v1.3.0 |
| `genetic_effect` | legacy_genetic_effect.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `make_pd` | legacy_make_pd.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `qtn_from_user` | legacy_qtn_from_user.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `vQTL` | legacy_vQTL.R | no | indirect | legacy | no — frozen legacy | v1.3.0 |
| `numericalize_core` | extendr-wrappers.R | no | direct | rust (generated) | is Rust | v1.3.0 |
| `meiosis_core` | extendr-wrappers.R | no | direct | rust (generated) | is Rust | isqg |
| `mate_haplotypes_core` | extendr-wrappers.R | no | direct | rust (generated) | is Rust | isqg |
| `gamete_masks_core` | extendr-wrappers.R | no | direct | rust (generated) | is Rust | isqg |
| `stable_hash_core` | extendr-wrappers.R | no | direct | rust (generated) | is Rust | — |
| `.onAttach` | zzz.R | no | no | package | no — glue | — |

</details>

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
