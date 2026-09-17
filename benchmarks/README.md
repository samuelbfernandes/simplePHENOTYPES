# Transcriptome simulation benchmark suite

Runnable benchmarks for the transcriptome-simulation features of
**simplePHENOTYPES** (`simulate_transcriptome()`, `predict()`,
`observe_counts()`, and the expression-mediated `transcriptome()` phenotype
layer with `mediation_split()`).

Each script is **standalone**: it loads the package from the repository root with
`devtools::load_all(".")`, runs a small self-contained analysis on the bundled
maize panel (`SNP55K_maize282_maf04`, 10,650 markers x 280 individuals), prints a
result table to stdout, and writes a CSV plus a PNG under `benchmarks/output/`.
Sizes are kept modest so every script finishes in well under a minute, and all
random draws are seeded for reproducibility. Plotting is wrapped so a missing or
broken graphics device downgrades to a warning instead of crashing the run.

## How to run

From the package root (an rextendr package — the Rust side must be compiled,
which `devtools::load_all()` handles):

```sh
Rscript benchmarks/01_h2_calibration.R
Rscript benchmarks/02_eqtl_recovery.R
Rscript benchmarks/03_coexpression_fp_control.R
Rscript benchmarks/04_mediation_recovery.R
Rscript benchmarks/05_twas_power.R
```

Or all at once:

```sh
for f in benchmarks/0*_*.R; do Rscript "$f"; done
```

Outputs (CSV + PNG) land in `benchmarks/output/` (created on first run).
`benchmarks/_common.R` holds shared helpers (package loading, the output
directory, a crash-safe `bench_png()`, and a `bench_dosage()` genotype reader);
it is sourced by each script and is not run directly.

## What each script does

| # | Script | Question | Key output |
|---|--------|----------|------------|
| 01 | `01_h2_calibration.R` | Does the generator realize the requested per-gene expression **h2** and **cis fraction**? | Bias/RMSE tables + target-vs-realized scatter |
| 02 | `02_eqtl_recovery.R` | Are the planted **cis-eQTL** recoverable by a naive marginal scan? | Rank/detection-rate of the true cis-eQTL |
| 03 | `03_coexpression_fp_control.R` | **(headline)** Can co-expression exist with **no** genetic basis, and would a naive "genetic co-expression" test false-positive? | Within vs between-module correlation; naive-test FPR; block heatmap |
| 04 | `04_mediation_recovery.R` | Does `mediation_split()` recover the intended genetic-mediated share, and does derived (not real) expression enter H2? | Intended-vs-realized split; derived vs real H2 |
| 05 | `05_twas_power.R` | Does **TWAS** power rise with the transcriptome `prop`, and how do cis- vs trans-driven genes differ? | Power-vs-prop; cis-predicted cis vs trans association |

### 01 — h2 & cis-fraction calibration
Sweeps target h2 (with `cis_fraction = 0.5`) and target cis fraction (with
`h2 = 0.7`) over 200 genes each on the maize panel, and compares the **realized**
values in `$genes` (`h2_realized`, `cis_fraction_realized`) to the targets. Reports
per-target bias and RMSE and saves a target-vs-realized scatter. (Observed: per-gene
h2 bias ≈ +0.0003, RMSE ≈ 0.024; the realized cis fraction of genes that carry a
cis-eQTL matches the target essentially exactly, because the cis/trans budget is a
standardized decomposition.)

### 02 — cis-eQTL recovery
Simulates a cis-heavy transcriptome (`h2 = 0.8`, `cis_fraction = 0.95`), then for
40 genes runs a marginal association scan over a candidate set of {the gene's true
cis-eQTL from `$cis_eqtl`} + 200 random decoy markers, and reports where the true
cis-eQTL rank. (Observed: the true cis-eQTL is the top hit for ~100% of genes, vs a
~2.5% top-5 chance under random ranking.)

### 03 — co-expression without a genetic basis (headline)
Generates a **genotype-free** transcriptome (`geno = NULL`, every gene `h2 = 0`)
with real co-expression modules, shows within-module correlations far exceed
between-module (~14x), and runs a naive "genetic co-expression" test that flags a
module as genetically co-regulated when its within-module correlation is
significantly higher than background. Every module is flagged — a **100%
false-positive rate**, because the truth is `h2 = 0`. A genotype-driven foil
(`h2 > 0`) shows co-expression looks the same with or without a genetic basis, so
co-expression alone cannot distinguish them. This is a ground truth no
genotype→trait simulator provides.

### 04 — mediation recovery & the derived/real H2 asymmetry
Builds a derived expression-mediated phenotype and sweeps the per-gene expression
heritability `h2*`. Shows `mediation_split()`'s `genetic_mediated` share tracks the
intended `prop * h2*` (observed corr ≈ 0.997), and that feeding the **same
expression matrix** as a *real* observed source (`expression =`) yields H2 ≈ 0 and
`NULL` mediation — because only genome-**derived** expression is credited to
heritability (its genetic-mediated part appears in `genetic_values()`).

### 05 — TWAS power
For a derived transcriptome-mediated phenotype, (A) an observed-expression TWAS
(correlate each gene's expression with the phenotype) shows detection power rising
with `prop`; (B) a **cis-predicted** TWAS (correlate each gene's cis-eQTL-imputed
expression with the phenotype — the single-gene cis-TWAS setting) recovers
cis-driven causal genes but is blind to purely trans-driven ones (mean |r| = 0 for
purely-trans genes, since they have no cis component to impute). Note the
observed-expression TWAS also shows a rising non-causal "FPR": non-causal genes that
share a co-expression module with causal genes genuinely correlate with the
phenotype — a real confound the simulator lets you study.

## External-tool comparisons (planned — NOT installed here)

These benchmarks are self-contained and depend only on simplePHENOTYPES. The
comparisons below situate our transcriptome layer against existing simulators and
estimands. **They are intentionally not run or installed here** — each needs a
separate R/Bioconductor or Python environment. This section names each tool, the
metric it would be compared on, and why our tool differs, so the comparison can be
staged deliberately (e.g. for the Bioinformatics manuscript benchmark) rather than
pulled in as a heavy dependency.

- **PhenotypeSimulator** (R / CRAN). Genotype → multi-trait phenotype simulator
  with genetic, shared/independent noise, and correlation structure. *Metric:*
  realized trait heritability and cross-trait genetic correlation. *Why we differ:*
  it has **no gene-expression layer** — phenotypes are drawn directly from markers
  and latent structure. Our `transcriptome()` layer inserts an explicit
  expression-mediated component with a genetic/environmental split
  (`mediation_split()`); the comparison would show our expression layer as the
  extension of the same variance-partition idea (benchmarks 01, 04).

- **GWASBrewer** (R). Simulates GWAS **summary statistics** from a specified
  trait/variant network (direct and mediated effects, LD). *Metric:* recovery of
  the specified effect/network structure from simulated sumstats. *Why we differ:*
  GWASBrewer produces marker→trait sumstats over a trait network; we produce
  **individual-level** genome→expression→phenotype data with a per-gene eQTL truth
  table (`$cis_eqtl`, `$factor_eqtl`) and an expression **mediator**, enabling
  individual-level TWAS/mediation benchmarks (benchmarks 02, 05) rather than
  sumstat-level ones.

- **MESC** (Python; "mediated expression score regression"). Estimates the fraction
  of trait heritability **mediated by assayed gene expression** from GWAS + eQTL
  sumstats. *Metric:* estimated expression-mediated h2. *Why we differ:* MESC is an
  *estimator* of the mediated-h2 estimand; we are the *generator* that defines its
  ground truth. `mediation_split()`'s `genetic_mediated` share is exactly the
  quantity MESC estimates, so benchmark 04 provides a labelled dataset to validate
  an MESC run against a known mediated fraction.

- **twas_sim** (Python). Simulates a **single gene's** cis-window genotypes,
  expression, and a cis-TWAS association. *Metric:* TWAS power / type-I error for a
  cis-imputed gene. *Why we differ:* twas_sim is single-gene and cis-only; our
  generator is transcriptome-wide with cis **and** trans (latent-factor) eQTL plus
  co-expression modules. Benchmark 05's cis-predicted arm is precisely the
  twas_sim cis-only limit — and it demonstrates what twas_sim cannot model: purely
  trans-driven causal genes that a cis-TWAS is structurally blind to.

- **SPsimSeq** / **scDesign3** (R / Bioconductor). Simulate realistic bulk/
  single-cell RNA-seq **counts** with gene-gene correlation learned from a reference
  dataset. *Metric:* fidelity of simulated count distributions and co-expression to
  a reference. *Why we differ:* they reproduce co-expression but carry **no genetic
  truth** — there is no eQTL, no heritability, no genome to trace expression to. Our
  `observe_counts()` gives a comparable NB count layer, but on top of a transcriptome
  with a known genetic architecture; benchmark 03 (h2 = 0 co-expression) is the clean
  contrast: identical-looking co-expression, but ours is labelled as genetic or not.

## Reproducibility notes

- Seeds are fixed inside each script; re-running reproduces the tables and plots.
- A harmless `OMP: Warning #179` line may appear on stderr from a compiled
  dependency's OpenMP setup under a restricted temp directory; it does not affect
  results.
- `benchmarks/output/` is regenerated on each run and can be deleted safely.
