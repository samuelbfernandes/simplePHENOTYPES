# DECISION-022 (DRAFT) — transcriptome simulation & expression-mediated phenotypes

> **DRAFT for review.** When approved, fold into `docs/DECISIONS.md` as DECISION-022
> and add the table row at the bottom. Companion spec: `SPEC-transcriptome.md`.

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

- **Scope now vs deferred.** Build the normalized generator + the transcriptome and
  genome→transcriptome phenotype bases first. **Defer** (explicit non-scope until
  the core variance semantics are stable): the RNA-seq **count** observation layer
  (`observe_counts()`, NB), directed regulatory networks, tissue specificity, and
  epistatic expression.

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
**Date:** DRAFT 2026-09-16 (pending review)

---

### Table row (append to the DECISIONS.md summary table on approval)
| 022 | transcriptome simulation: hybrid latent-factor eQTL (cis in a physical-bp window + trans via Q≪T latent regulatory factors + non-genetic co-expression modules + gene noise), normalized Gaussian scale; per-gene `h2_g`, cis fraction `omega_g`, residual module fraction `kappa_g` as non-competing knobs; **joint** cis/trans genetic scaling with a reported `cis_trans_cov` budget row (cf. `add_dom_cov`, DECISION-020); **fixed reference calibration** (DECISION-020/021); no annotation/reference data required (named `generic_bulk` profile; synthetic physical-bp coords; optional `mimic=` calibrates the expression generator only, GREML per-gene h²; ships an example annotation + expression dataset for SNP55K); phenotype bases inferred from inputs via a `transcriptome()` layer — markers-only (default) / `expression=` real transcriptome alone / `geno`+`transcriptome=` derived G→E→Y (`y = Z(δ+Bs)+Rs+η`, mediated+direct) / `geno`+`expression=` both real (`y = Zδ+Ẽs+η`, no simulated mediation, realized G–E cov reported) — all with genome-path rigor; counts/GRN/tissue/epistasis deferred; RNG in R (DECISION-006), no Rust yet | DRAFT (2026-09-16) |
