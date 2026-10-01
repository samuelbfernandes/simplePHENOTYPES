# Transcriptome simulation benchmark suite

Runnable benchmarks for the transcriptome-simulation features of
**simplePHENOTYPES** (`simulate_transcriptome()`, `predict()`,
`observe_counts()`, and the expression-mediated `transcriptome()` phenotype
layer with `mediation_split()`).

Each script is **standalone**: it loads the package from the repository root with
`devtools::load_all(".")`, runs a small self-contained analysis on the bundled
maize panel (`SNP55K_maize282_maf04`, 10,650 markers x 280 individuals), prints a
result table to stdout, and writes a CSV plus a PNG to an output directory (by
default a session temporary directory; see below -- nothing is written into the
repository unless you ask for it).
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

Outputs (CSV + PNG) land in `<tempdir()>/simplePHENOTYPES-benchmarks` by default
(the path is printed at the start of each run). To keep them, name a directory
with the `SP_BENCH_OUT` environment variable, e.g.
`SP_BENCH_OUT=benchmarks/output Rscript benchmarks/01_h2_calibration.R`
(`benchmarks/output/` is gitignored). `benchmarks/_common.R` holds shared helpers
(package loading, the output directory, a crash-safe `bench_png()`, and a
`bench_dosage()` genotype reader); it is sourced by each script and is not run
directly.

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
values in `$genes` to the targets. Two of the reported columns are informative and
two are not: the bounded `h2_allocated = Var(G)/(Var(G)+Var(R))` is a variance
allocation (not a heritability) that equals the target by construction on the
reference population (shown only as a sanity column), so the h2 panel uses the
realized heritability `h2_realized = Var(G)/Var(P)` (also stored as
`h2_var_ratio`); the cis panel uses
`cis_fraction_realized = v_cis/Var(G)` (covariance included), which scatters around
the target. Both panels have tolerances and the script stops if one is exceeded, so
neither can pass vacuously. (Observed: realized-h2 bias ≈ +0.0003, RMSE ≈ 0.024;
cis-fraction bias ≈ +0.002, RMSE ≈ 0.040, i.e. the realized cis share deviates from
the target by the cis-trans covariance.)

### 02 — cis-eQTL recovery
Simulates a cis-heavy transcriptome (`h2 = 0.8`, `cis_fraction = 0.95`), then for
two sets of 40 genes -- a **random** draw from all genes with a cis-eQTL (the
headline) and, for comparison, the 40 with the largest |cis effect| (the easiest
subset) -- runs a marginal association scan over a candidate set of {the gene's
true cis-eQTL from `$cis_eqtl`} + 200 random decoy markers, and reports where the
true cis-eQTL rank. The chance comparator is exact per gene,
`1 - choose(D, K)/choose(D + m, K)` for `m` true markers among `D` decoys and
top-`K` (about 4% on average here, since 47% of genes have more than one true
marker; the one-marker value `K/(D+1)` = 2.5% understates it). (Observed: the true
cis-eQTL is the top hit for ~100% of both sets.)

### 03 — co-expression without a genetic basis (headline)
Generates a **genotype-free** transcriptome (`geno = NULL`, every gene `h2 = 0`)
with real co-expression modules, shows within-module correlations far exceed
between-module (~14x), and runs a naive "genetic co-expression" test that flags a
module as genetically co-regulated when its within-module correlation is
significantly higher than background. Every module is flagged — a **100%
naive-inference error rate**, measured against the stated truth `h2 = 0` (the
Wilcoxon test itself is not a genetic test; it is presented as the naive,
wrong inference). A genotype-driven foil
(`h2 > 0`) shows co-expression looks the same with or without a genetic basis, so
co-expression alone cannot distinguish them. This is a ground truth no
genotype→trait simulator provides.

### 04 — mediation recovery & the derived/real H2 asymmetry
Builds a derived expression-mediated phenotype and sweeps the per-gene expression
heritability `h2*`. The rough expectation `prop * mean(h2)` holds only for
independent, equal-weight causal genes; the exact share is
`prop * Var(Σ w_g (G_g - mean)/s_Eg) / Var(Σ w_g (E_g - mean)/s_Eg) / V_P`
(no independence assumption), which the script computes from the expression
matrices independently of `mediation_split()` and checks against it (agreement to
~1e-16); `prop * h2*` is reported as the approximation (max deviation ≈ 0.03 with
20 causal genes). No correlation over the three sweep points is reported. Feeding the **same
expression matrix** as a *real* observed source (`expression =`) yields H2 ≈ 0 and
`NULL` mediation — because only genome-**derived** expression is credited to
heritability (its genetic-mediated part appears in `genetic_values()`). The
realized total share differs from `prop` because `V_P` is not exactly 1.

### 05 — TWAS power
For a derived transcriptome-mediated phenotype, (A) an observed-expression TWAS
(correlate each gene's expression with the phenotype) shows detection power rising
with `prop`, averaged over 20 phenotype replications per `prop` with the Monte
Carlo standard error reported (not a single realization); (B) a **cis-predicted**
TWAS (correlate each gene's cis-eQTL-imputed expression with the phenotype — the
single-gene cis-TWAS setting) recovers cis-driven causal genes; its power for
purely trans-driven ones is 0 **by construction** (they have no cis component to
impute), an illustration of the cis-only limit rather than an empirical finding.
Non-causal genes are split into *structural nulls* (modules with no causal gene:
no shared module and no direct causal-gene status, i.e. no designed path to the
phenotype; this is a statement about the simulated design, not a proof of zero
association, because marker linkage disequilibrium (LD) between eQTL can still
induce weak association, so their Bonferroni rate is the false-positive rate for
genes with no designed link and is expected to be small, not guaranteed zero) and
*module-linked* genes (they share a co-expression module with a causal gene and
genuinely correlate with the phenotype — a real confound the simulator lets you
study, with a rate that rises with `prop`).

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
- Output files are regenerated on each run and can be deleted safely (they live in a temporary directory unless `SP_BENCH_OUT` is set).
