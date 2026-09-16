# SPEC-transcriptome.md — transcriptome simulation & expression-mediated phenotypes

> **DRAFT for review.** Companion to `SPEC.md`; extends the v2 grammar with (a) a
> generator of genetically controlled gene expression, and (b) phenotype
> simulation from a transcriptome. **Implementation status:** the parametric
> generator `simulate_transcriptome()` (§3), the `transcriptome()` phenotype
> layer for the genome-present modes (see the "v1 layer status" note in §3), and
> the **derived mediation split with covariance reporting** (`mediation_split()`;
> the genetic-mediated part of derived expression now counts toward realized `H²`)
> are implemented and tested, as are the **genotype-free mode 2**
> (`simulate_phenotype(expression = ...)` with no `geno`), `qtn_table()` gene rows,
> and **cross-population reuse** of a fixed architecture
> (`predict.transcriptome_sim()`). Still designed here but **not yet built**:
> `mimic` mode (§5), the counts layer (§7), and a genotype-free **generator**
> (`simulate_transcriptome(geno = NULL)`, a purely non-genetic transcriptome --
> distinct from the phenotype-side mode 2 above, which is done). Converged
> Claude + Codex design (see `project_transcriptome_simulation_design` memory and
> `DECISION-022-transcriptome-DRAFT.md`).

## 1. Goals and non-goals

### Goals
- Simulate a **genes × individuals expression matrix** that is genetically
  controlled (cis + trans eQTL) and realistically co-expressed, **without
  requiring reference expression data or a gene annotation**.
- Let a user **provide expression they want to mimic**; calibrate the generator
  to it (moments + low-rank co-expression, and heritability if genotypes are
  paired) and regenerate synthetic expression.
- Simulate a **phenotype from a transcriptome** — user-supplied or genome-derived
  — with the **same rigor as the genome path** (variance budget, realized-h²
  reporting, effect tables, seed reproducibility, input validation).
- One phenotype interface, three **bases**: genome; transcriptome (user); genome
  → transcriptome → phenotype (mediated + direct).

### Non-goals (this version)
- Mechanistic / directed gene-regulatory networks (knockouts, feedback, time
  courses). The trans structure is a phenomenological latent-factor model.
- Tissue specificity, alternative splicing, isoform-level expression.
- Epistasis in the expression model.
- RNA-seq **count** realism is an *optional* later observation layer, not the core.

## 2. Core object

`transcriptome_sim` — returned by `simulate_transcriptome()` (v1 return value):
- `expression` — genes × individuals (normalized ~Gaussian scale).
- `genetic_expression` — the genetic component `G` only (genes × individuals).
- `genes` — per-gene data frame: `gene_id`, `chr`, `tss`, `module`,
  `coordinate_source` (`"supplied"`/`"synthetic"`), `h2_target`, `h2_realized`,
  `cis_fraction_target`, `cis_fraction_realized`, `n_cis`, `trans_scale`.
- `cis_eqtl` — data frame `gene_id, snp, chr, pos, effect` (effective coefficient
  on centered dosage; the cis truth), or `NULL` if no gene has a cis-eQTL.
- `factor_eqtl` — data frame `factor, snp, chr, pos, hub_effect` (raw hub effects;
  the trans truth), or `NULL` if there is no genetic component.
- `loadings` — data frame `gene_id, factor, loading` (module membership; loading 1
  in v1).
- `reference` — list: `marker_mean` (per-marker reference means for centering),
  `ids`, `n_factors`, `kappa`, `cis_window`.
- `var_budget` — data frame `gene_id, v_cis, v_trans, cis_trans_cov, gr_cov`
  (`v_cis + v_trans + cis_trans_cov = Var(G)`; `gr_cov = 2 Cov(G, R)` is the
  finite-sample genetic-residual covariance, ~0 in expectation).
- `profile`, `seed`, `n_genes`, `n_ind`.

The effective trans coefficient on centered dosage `Z_k` for a hub of gene `g`'s
module is `genes$trans_scale[g] * factor_eqtl$hub_effect` (loading 1), so the whole
genetic component is reconstructable from the truth tables.

The truth tables make this a **known-ground-truth benchmark generator** for
eQTL / TWAS / mediation methods.

## 3. Function reference

> **Implementation status.** This section specifies the full design. **v1
> implements the parametric generator only** — the signature and behaviour marked
> *(v1)* below. The arguments marked *(planned)* — `mimic`, `counts`, and
> `geno = NULL` (purely non-genetic expression) — are designed here but **not yet
> in the function**; v1 requires `geno` and needs at least 3 individuals.

### `simulate_transcriptome()`
```r
# v1 signature (implemented):
simulate_transcriptome(
  geno,                           # (v1) Population / Population-backed phenotype_sim / data.frame / -1/0/1 matrix; REQUIRED
  n_genes        = 1000,
  annotation     = NULL,          # (v1) optional data.frame(gene_id, chr, tss); else synthetic coords
  cis_window     = 1e6,           # (v1) physical bp each side of the TSS
  h2             = "beta",        # (v1) per-gene total expression h2: "beta" profile, a scalar, or a length-n_genes vector
  cis_fraction   = 0.25,          # (v1) omega: cis share of MARGINAL genetic variance
  n_factors      = NULL,          # (v1) Q; default max(1, min(50, max(5, ceil(T/100)), n-2))
  residual_module_fraction = 0.15,# (v1) kappa
  profile        = "generic_bulk",# (v1)
  seed           = NULL           # (v1)
) -> transcriptome_sim
# planned additions: mimic = NULL (§5), counts = FALSE (§7), geno = NULL (non-genetic)
```
- **Parametric (default, v1):** parameters drawn from `profile`; synthetic gene
  coords when `annotation` is `NULL`; cis defined on **physical bp** (never cM).
  For **synthetic** coords, chromosomes are chosen weighted by eligible
  (MAF ≥ 0.05) marker count and each gene's TSS is snapped onto an eligible marker
  when a uniform placement finds none, so every synthetic gene has a cis marker. A
  bare genotype **matrix** (no chromosome info) is treated as a single chromosome
  (pos = column index). With a **supplied `annotation`**, a gene whose window
  contains no eligible marker realizes `cis_fraction_realized = 0` (all genetic
  variance from trans) — reported, not an error. **Trans hubs are drawn outside the
  cis windows of every gene in their module**, so factor-mediated variance is
  genuinely trans (distant). The genetic component and the residual are drawn
  **independently** (the residual keeps the shared-module structure, so genes
  co-express through `kappa`); realized `h2 = Var(G)/Var(E)` tracks the target up to
  a reported finite-sample `Cov(G, R)` (`var_budget$gr_cov`). Needs ≥ 3
  individuals.
- **`geno = NULL` (planned):** expression with no genetic component — only shared
  non-genetic modules + noise.
- **`mimic = E` (planned):** calibrate to `E` (§5), then generate.

### Phenotype bases (extend `simulate_phenotype()`)
Two new inputs set the expression source; the basis follows from which inputs are
present (no explicit `basis =` needed). `expression =` is **observed/real**
expression (a genes × individuals matrix, treated as data); `transcriptome =` is
**genome-derived** expression — either a `transcriptome_sim` from
`simulate_transcriptome()`, or `TRUE` to derive one internally from `geno` with the
default profile. At most one expression source (error if both are given). Genome
layers (`additive()`, …) require `geno`; the `transcriptome()` layer requires an
expression source. **Markers-only is the default** when neither is given.

| Inputs | Basis | Phenotype |
|---|---|---|
| `geno` (default) | markers only | `Y = Zδ + η` (current, unchanged) |
| `expression = E` (no `geno`) | transcriptome alone | `Y = Ẽs + η` (E is real data) |
| `geno` + `transcriptome = TRUE`/`<transcriptome_sim>` | genome → transcriptome → phenotype (derived) | `Y = Z(δ + Bs) + Rs + η` (mediated + direct) |
| `geno` + `expression = E` (both real) | genome + real transcriptome | `Y = Zδ + Ẽs + η` (two real predictors; realized G–E covariance reported, not asserted) |

```r
# markers only (default; current behavior)
simulate_phenotype(geno) |> additive(prop = 0.5)

# transcriptome alone (real expression, no genome)
simulate_phenotype(expression = E) |> transcriptome(prop = 0.5, n_genes = 20)

# genome -> derived transcriptome -> phenotype (mediation)
simulate_phenotype(geno, transcriptome = TRUE) |>       # or transcriptome = simulate_transcriptome(geno)
  transcriptome(prop = 0.3, n_genes = 20) |>            # genetically-mediated + env-mediated expression effect
  additive(prop = 0.2)                                  # optional DIRECT genome effect

# both real datasets: real markers AND observed expression
simulate_phenotype(geno, expression = E) |>
  additive(prop = 0.3) |>                               # direct marker effect
  transcriptome(prop = 0.4, n_genes = 20)               # observed-expression effect (no simulated mediation)
```

### The `transcriptome()` layer
> **v1 layer status (implemented):** the `transcriptome()` layer scores an
> attached expression source (`expression=` real, or `transcriptome=` derived) as
> a sparse linear function of z-scored expression, scaled to `prop`. It works both
> with a genome (`simulate_phenotype(geno, ...)`) and **without one** -- the
> genotype-free mode 2, `simulate_phenotype(expression = ...)` with no `geno`,
> where individuals come from the expression columns and only `transcriptome()`
> layers are valid. It is reported as a
> **distinct expression-mediated variance category, scaled by `prop` outside the
> marker `h2` budget**. **Mediation split (implemented):** for a *derived*
> transcriptome the component now decomposes into a genetic-mediated part `Tx_g`
> (traced to the genome through expression, using the identity
> `z_total = z_genetic + z_env` from `E = G + R` standardized by the same `sd(E)`)
> and an env-mediated part `Tx_e = Tx - Tx_g`. `Tx_g` **is** genetic and enters
> `genetic_values()` and realized `H²`; `mediation_split()` reports the realized
> genetic / env / covariance (`2·Cov(Tx_g,Tx_e)/V_P`) shares, which sum to the
> realized expression-mediated share (machine-precision closure, Codex-verified).
> A *real* `expression=` source asserts no genetic content, so `Tx_g = 0` and the
> whole component stays out of `H²`. `qtn_table()` gene rows, the genotype-free
> `simulate_phenotype(expression=)` mode 2, and cross-population reuse
> (`predict.transcriptome_sim()`) are **implemented**; **remaining follow-ups** are
> `mimic` calibration, the counts layer, and the genotype-free *generator*
> `simulate_transcriptome(geno = NULL)` (§3). v1 reports **marginal**
> variance shares per layer; when a transcriptome predictor is strongly
> (anti-)correlated with a marker layer (e.g. an expression gene equal to a causal
> marker's dosage -- a pathological input), their finite-sample marker-to-`Tx_g`
> covariance is included in the realized `H²` numerator (`Var(Zδ + Tx_g)`) but not
> attributed to any single marginal budget row, so a reported realized `H²` can
> still fall slightly outside `[0,1]` under such pathological inputs.

> **Orientation convention.** In every phenotype equation below, `E`, `R`, and `Z`
> are written in the standard **individuals × features** design-matrix orientation,
> so `E s`, `R s`, `Z δ` are matrix–vector products giving one value per individual.
> The `transcriptome_sim` stores `expression` / `genetic_expression` transposed
> (genes × individuals), so `E s` is computed as `t(tx$expression) %*% s`.

A first-class grammar layer whose predictors are **continuous standardized
expression** (not dosages), held to the same guarantees as `additive()`:
- `prop` — target share of phenotypic variance (marginal variance of the expression
  predictor `E s`), exactly like `additive(prop=)`.
- `n_genes` / `genes` — number, or explicit set, of causal genes (sparse).
- `slopes` — optional explicit slopes; else drawn `N(0,1)` and jointly rescaled by
  one common multiplier to hit `prop` on the reference (preserving relative slopes).
- Scores on the sim's expression source (real `expression=` or derived
  `transcriptome=`); reports a causal-gene effect table (the `qtn_table()` analog),
  realized variance share, and per-gene contributions.
- **Derived vs real semantics.** With a *derived* transcriptome the expression
  effect carries a genetically-mediated part (`Bs`) and an env-mediated part (`Rs`),
  and `additive()` supplies the **direct** genome effect (`δ`). `prop` sets the total
  expression-mediated variance share (marginal `Var(E s)`); its genetic-mediated vs
  env-mediated split, and hence the **overall `H²` is emergent and
  reported**, not separately forced (you cannot both pin `Var(E s) = prop` and pin
  the genetic-h² contribution `Var(Z(δ+Bs))` independently — they share `s`). To
  target overall `H²` instead, scale the total genetic score `Z(δ+Bs)` to it and let
  the expression-mediated share emerge; the SPEC picks the `prop`-on-`Var(E s)`
  convention and reports the rest. With a *real* `expression=` matrix there is **no
  simulated mediation** — the marker and expression layers are separate
  variance-budget components, and
  their realized covariance (expression is biologically downstream of the genome in
  real data) is reported, not forced to zero.

## 4. Model (normalized Gaussian scale)

`Z` = dosage centered on **reference** allele frequencies (do not assume −1/0/1 is
mean-zero). For gene `g`:
```
E_g   = mu_g + G_g + R_g
G_g   = sqrt(h2_g) * normalize_ref( sqrt(omega_g)*cis~_g + sqrt(1-omega_g)*trans~_g )
  cis_g   = sum_{j in window(g)} beta_gj  * Z_j
  trans_g = sum_q Lambda_gq * f_q ;  f_q = sum_{j in hub(q)} gamma_qj * Z_j
R_g   = sqrt(1-h2_g) * normalize_ref( sqrt(kappa_g)*mod~_g + sqrt(1-kappa_g)*eps~_g )
  mod_g = sum_q Lambda_gq * u_q ;  u_q ~ N(0,1)   # NON-genetic module factor (same loadings)
  eps_g ~ N(0, .)                                  # gene-specific noise
```
Key semantics:
- **Joint cis/trans scaling.** cis and trans covary (LD/structure), so
  `V_cis + V_trans != V_G`. Scale the *combined* genetic score to `h2_g`; report
  `V_G = V_cis + V_trans + 2*Cov(cis,trans)` with the covariance as its **own
  budget row** (the same pattern as the orthogonal model's `add_dom_cov` row,
  DECISION-020).
- **Co-expression decoupled from heritability.** Non-genetic module factors give
  co-expression even at `h2_g = 0`. `omega_g` controls the genetic cis/trans
  split; `kappa_g` controls residual co-expression; `h2_g` controls total genetic
  variance — three non-competing knobs.
- **Fixed reference calibration.** `mu_g`, allele freqs, and all sd's are frozen
  on the reference (founder) population and never re-estimated on
  descendant/crossed/selected populations — the fixed-scale principle of
  `additive_value()`/`phenotype_value()` (DECISION-020/021).
- **Realized, not asserted.** All reported h²/shares are computed from realized
  values (components not orthogonal at freq ≠ 0.5), per SPEC §2 V2.

Phenotype from expression (individuals × features orientation, per the layer
convention above — `E`/`R`/`Z` are individuals × genes/markers here):
`y = Z*delta + E*s + eta = Z*(delta + B*s) + R*s + eta` (B = the genetic effect
matrix of expression). `B*s` = mediated-genetic, `delta` = direct-genetic,
`R*s` = env-mediated. Do not force each share independently; scale
the total genetic score `Z*(delta+B*s)` to the phenotype's target h² and report
the mediated/direct split + covariance.

## 5. Mimic mode (calibrate to user expression)
`mimic = E_user` (genes × individuals). Estimate, on `E_user`:
- per-gene mean `mu_g` and total variance `V_g`;
- low-rank co-expression: a truncated SVD of the standardized matrix → factor
  scores and loadings `Lambda` (retain `n_factors` components), residual noise;
- if paired genotypes are supplied, per-gene `h2_g` via a **GREML-style** estimator
  (REML on a genomic relationship matrix, per gene), used to calibrate the target
  `h2_g` distribution — not to fit individual eQTL effects.
Then **generate** new synthetic expression matching those parameters (for the same
or new individuals). Mimic calibrates the **expression generator only**; phenotype
slopes are still generated de novo (the user does not want slopes fit to data).
Annotation is still optional.

## 6. Reproducibility
All stochastic draws in R (DECISION-006); one `seed` reproduces every basis.
Seed-threading extends the grammar's `(seed, layer_index, layer_type)` rule to the
`expression()` layer, and `simulate_transcriptome()` takes its own `seed`.

## 7. Deferred (explicit non-scope, revisit when the core is stable)
- Counts observation layer `observe_counts(tx, ...)`: `Y_gi ~ NB(mu, phi)`,
  `log mu = log L_i + alpha_g + sigma_g E_gi` (count-h² != latent-h²).
- Directed regulatory networks; tissue specificity; epistatic expression.

## 8. Testing requirements (mirror SPEC §8)
1. Realized per-gene `h2` and `cis_fraction` track targets within tolerance across
   the profile.
2. cis/trans covariance budget **closes**: `V_cis + V_trans + 2Cov = V_G`.
3. Co-expression: within-module > between-module; **non-zero at `h2 = 0`**.
4. Simulated cis-/trans-eQTL sit at their truth loci (recoverable by a mapping
   scan on the simulated data).
5. Fixed reference calibration: scoring a subset/descendant does **not** rescale
   (values comparable across generations).
6. Transcriptome-basis phenotype: realized `prop` ≈ target; effect table correct;
   equals the genome path's rigor.
7. Mediated phenotype: `y` decomposition `Z(delta+Bs) + Rs + eta` holds numerically.
8. Mimic mode reproduces the moments and leading factors of a held-out matrix.
9. Seed reproducibility across all three bases.
10. All theory claims pass an independent `dev/dual.sh` review before commit.

## 9. Resolved design decisions (2026-09-16)
1. **Basis API** — inferred from which inputs are supplied (`geno` / `expression` /
   `transcriptome`); no explicit `basis =` argument. Markers-only is the default.
2. **Layer name** — `transcriptome()`.
3. **Expression source semantics** — `expression =` is real/observed data
   (transcriptome alone when no `geno`; genome + real transcriptome when both);
   `transcriptome =` (`TRUE` or a `transcriptome_sim`) derives expression from the
   genome (G→E→Y, with simulated mediation). `simulate_transcriptome(geno = NULL)`
   (purely non-genetic expression) is planned, not in v1. See the §3 mode table.
4. **Mimic heritability estimator** — GREML-style (REML on a GRM, per gene) to
   calibrate the `h2_g` distribution.
5. **Example data** — ship an example dataset (a synthetic gene annotation, and an
   example expression matrix, for the bundled SNP55K panel) so all bases run
   out-of-the-box.

### Still open (later)
6. Whether `transcriptome()` should allow **non-linear** expression→phenotype
   effects (thresholds, expression×expression) later, or stay linear for v1.
   (Draft: linear v1.)
