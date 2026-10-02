# simplePHENOTYPES (development version)

## New features
* `filter_geno()` accepts a `Population` (founders or crossing progeny), so the
  MAF, heterozygosity and LD filters cover every genotype input
  `simulate_phenotype()` takes; the map, ids and pedigree are kept and only the
  failing markers are dropped.
* `dominance()` gains `effect =` (a geometric base, an explicit series of length
  `n_qtn`, or a per-trait list), matching `additive()` and v1 `dom_effect`, and
  `epistasis(effect =)` now also takes a per-trait list (one base or series per
  trait, like `additive()`); both are rejected under multi-trait
  `architecture = "pleiotropy"`, whose correlated draw sets the effects.
* Every `Population` now records its pedigree (`parentage()`, `families()`), kept
  through subsetting and pooling; `as_population()` gains `pool =`.
* `mating_design()` (random, factorial, nested, diallel, half-diallel) and `mate()`
  run mating plans across one or several populations.
* `combining_ability()`: GCA / SCA / testcross merit of candidates against
  testers (topcross, factorial, diallel), as exact expected cross means on a
  frozen architecture or from simulated progeny; `template_effects()` exports a
  simulation's realized per-locus `a` and `d`; `phenotype_value(d =)` scores the
  total genotypic value with a broad-sense `h2`.
* `progeny_test()`: half-sib progeny testing of parents on random mates.
* `select_ind(method = "culling")`: independent culling levels (simultaneous or
  sequential); tandem selection through a `trait` vector in `pedigree()` and
  `recurrent_selection()`.
* `marker_select()`: marker-assisted selection and staged gene pyramiding.
* `predict_ebv()`: known-variance BLUP breeding values (GBLUP or pedigree),
  with `a_matrix()` (numerator relationship matrix from the recorded pedigree)
  and `prediction_accuracy()`; `selection_methods()` lists the selection
  operators for front ends.
* `predict_ebv()` multi-trait BLUP: an individuals x traits `pheno` matrix (`NA`
  for missing records) with known genetic (`var_a`) and residual (`var_e`)
  covariance matrices predicts every trait for every individual, including traits
  an individual was not recorded on (Henderson & Quaas 1976).
* `write_phenotypes(file_type = "json")` writes the long or wide table as JSON (one
  object per row, UTF-8, values round-trip exactly); needs the suggested package
  jsonlite.
* `crossbreed()` (two-way, backcross, three-way, terminal, rotational),
  `breed_composition()` and `heterosis()`.
* `mabc_select()` and `recurrent_parent_recovery()`: marker-assisted
  backcross selection (foreground filter, flanking-marker recombinant
  selection, background recovery of the recurrent-parent genome; Frisch &
  Melchinger 2001, 2005).
* `additive(effect = list(...))`: one effect specification per trait in a
  single layer, e.g. to re-score per-trait effects frozen from an earlier
  (pleiotropic) simulation on fixed `qtn`.
* Under `architecture = "pleiotropy"`, `cor` now controls `dominance()` and
  `epistasis()` too (DECISION-023). Previously both layers gave every trait
  one identical effect series, so on the loci the traits shared the
  non-additive components correlated near +1 regardless of `cor` (e.g. a
  target of 0 gave 1.0 for pleiotropic epistasis).

## Behaviour changes
* Under `architecture = "ld"`, `epistasis()`, `dominance(same_as_add = FALSE)`
  and a second `additive()` layer now error: they could make one marker causal
  for both traits.
* Documentation now states that `cor` is a target: the realized correlation
  converges to it only as QTNs and individuals grow, for causal loci in
  approximate linkage equilibrium and without major QTNs.

## Audit fixes (independent dual-model audit, 2026-09)

Fixes for the defects confirmed by a two-model theory and implementation audit
(reports in the maintainers' `.tmp/audit-2026-09-29/`). Items that change seeded
output or reject previously accepted input are marked **(behaviour)**.

**Simulation grammar**
* **(behaviour)** Sub-seeds now use a position-sensitive hash of the draw label
  (DECISION-033). Replications 12/21, 13/31, ... and traits 12/21 used to share a
  sub-seed and returned identical QTNs, effects and residuals. Every seeded
  `simulate_phenotype()` result changes. The rule is (seed, layer type, occurrence
  of that type); inserting a same-type layer shifts later same-type layers.
* `additive()` + `dominance()` on shared loci report the realized Var(A), Var(D)
  and 2Cov(A, D) (`$ad_report`, and a note in `print()`). The realized-H2 gap is
  a coding-dependent term, not finite-sample noise; `additive(orthogonal = TRUE,
  ...)` is the recommended model (DECISION-034).
* **(behaviour)** `simulate_phenotype()` needs at least three individuals.
  `complex_phenotypes()` needs h2-complete inputs, is terminal (adding a layer
  errors) and no longer carries model-1 state (`mediation_split()` is `NULL`).
  `epistasis(qtn = list(...))` is one element per trait. Markers heterozygous in
  every individual are no longer drawn as QTNs.
* Under `architecture = "ld"`: the strongest-partner rule is documented, with an
  opt-in `partner = "random"`; r2 = 0 / 1 partners are never used. Geometric
  series that overflow or underflow are rejected.

**Genotype input**
* **(behaviour)** HapMap: the `alleles` column is validated against the observed
  calls; a mismatch warns and the orientation comes from the calls (previously the
  two homozygote classes could collapse to one code).
* **(behaviour)** All readers share one raw-allele contract (first-listed allele =
  allele 1). On tied (MAF = 0.5) markers a VCF read from a path and from a data
  frame now gives identical dosages; the tie rule is documented in `?as_numeric`.
  PLINK BED `allele` labels are `A1/A2` in file order.
* FinalReport literal `NA` alleles are no-calls (were a called homozygote). Calls
  are matched case-insensitively and the default `hets` is the full 18-code set.
  Duplicated marker or sample IDs are an error; `filter_geno()` validates its
  thresholds.

**Selection, relationship and prediction**
* **(behaviour)** `select_ind(trait = c(1, 2))` is an error (it used to be
  recycled into an alternating-trait criterion). Index scores match phenotypes to
  individuals by id. `bulk()` draws its seeds at random from the pooled progeny
  (seeded results change). `recurrent_selection()` errors when `n_parents` is at
  least the number of plants. Scheme wrappers restore the caller's RNG.
* `g_matrix()` / `optimum_contribution()` work on a subset `phenotype_sim`.
  `optimum_contribution()` needs a named `G`, warns for an unattainable target and
  centres merit internally, so constant offsets and very large merit spreads are
  handled. **(behaviour)** `sample_parents()` now
  allocates parent slots by largest remainder (`method = "allocate"`, default), so
  each individual's count stays within one slot of `n * c_i` (marginal count error
  is controlled); this approximates but does not preserve the group coancestry
  0.5 c'Gc, whose error depends on `n` and `G` (e.g. `n = 1`, `G = I`,
  `c = (0.5, 0.5)`: optimum 0.25, realized 0.5). In the F2 example the target
  coancestry is 0.0521 and the allocation gives 0.0518 at `n = 30` versus a
  multinomial mean of 0.0672; this is problem dependent, not a general guarantee
  (counterexample `c = (0.6, 0.4)`, `n = 1`, `G = diag(1, 3/7)`: allocation error
  0.286 versus multinomial 0.171); the previous independent weighted draw is
  `method = "multinomial"` (same results for the same seed). `cross_usefulness()`
  warns when a `"dh"`/`"selfcross"` family is built from heterozygous parents.
* `a_matrix()` gains `founder_f` for inbred-line founders. Several `predict_ebv()`,
  `prediction_accuracy()`, `combining_ability()` and `progeny_test()` inputs that
  were silently accepted are now errors. `template_effects()` errors instead of
  silently omitting a transcriptome layer.

**Transcriptome**
* **(behaviour)** `genes$h2_realized` is the realized heritability Var(G)/Var(P)
  computed from the realized genetic values and expression (it includes 2Cov(G, R)
  and is not bounded by 1 at small n); `h2_var_ratio` is kept as an identical alias.
  The bounded allocation Var(G)/(Var(G)+Var(R)), which equals the target on the
  reference panel by construction, is now `h2_allocated` (it is not a heritability).
  The mimic GREML identifiability guard now tests the spectrum of K after
  projecting out the intercept (a warning is issued and h2 is set to 0). The
  marginal-epistasis-share doc now gives epsilon/(epsilon+(1-epsilon)/s_ct^2).
  `cis_fraction_realized` is the realized
  share v_cis/Var(G) (it used to restate the target); new `epistasis_realized`.
  `mimic` estimates kappa without double counting genetic trans structure.
  `observe_counts()` returns integer counts.

**Crossing and Rust core**
* A malformed call to the Rust kernel is now an ordinary R error instead of aborting
  the R session (macOS gcc builds). Minimum supported Rust is 1.71.
* **(behaviour)** Chromosome order (and the seeded random stream) no longer depends
  on the storage type of `chr` or the locale; only character/factor `chr` that used
  to sort differently (e.g. `"1"`, `"10"`, `"2"`) changes.
* **(behaviour)** Seeded crossing (`cross()`, `selfcross()`, `double_haploid()`,
  `mate()`, `crossbreed()`) restores the caller's RNG state.
* `as_population()` warns when the map looks like Morgans; `cross()` warns when
  two panels list a marker's alleles in opposite order; `heterosis()` checks breed
  membership and states the Hardy-Weinberg condition of its retention fractions.

**Frozen v1 `create_phenotypes()` (bad inputs are rejected; valid output is unchanged)**
* **(behaviour)** Errors are re-signalled (the function used to print the message
  and return `NULL`) and the caller's RNG kind and state are restored.
* **(behaviour)** Rejected with an informative message (previously a cryptic error,
  `NULL`, `NA` phenotypes or silently wrong numbers): `QTN_list` with `ntraits = 1`;
  single-trait `"DE"`/`"ADE"`; `model = "D"` with indirect LD; `ld_max >= 1`; `h2`
  outside [0, 1] or below 0.05 with `rep > 1`; non-positive-definite `cor` (it is no
  longer repaired by the eigenvalue clamp); effect vectors of the wrong length;
  vQTL with several traits or `h2 = 0`; 0/1-only numeric genotypes.
* LD pairs that share a marker between traits, span chromosomes or fall outside
  `[ld_min, ld_max]` stop with an "LD contract" error. LD diagnostic files carry the
  right trait labels, and `Epistatic_QTNs.txt` lists each trait's own effects.
* The `create_phenotypes()` help now states the residual-seed formula the code
  actually uses and its consequences.

**Tooling**
* `evals/run.sh` never modifies working-tree sources; the commit-message and CI
  attribution guard share one pattern file (`.githooks/ai-patterns`) with a
  self-test (`dev/test-attribution-guard.sh`); the vignette uses the real
  `write_phenotypes()` signature; benchmarks write to a temporary directory.

### Review round 2 (Codex review of the audit fixes)

An independent Codex review of the fixes above found further defects; they are
fixed here. Items that change output or reject previously accepted input are
marked **(behaviour)**.

**Simulation grammar**
* `vqtl(same_as_add = TRUE)` after a `pleiotropy` additive layer whose traits retain
  different numbers of loci (e.g. `pi = c(1, 0.5)`) no longer errors "non-conformable
  arguments"; a reused dominance layer keeps the shared loci when a trait has
  `pi = 0`; `print()` shows the retained per-trait QTN counts.
* `$ad_report` gains `var_cA`, `var_cD` (aggregate component variances) so
  `realized = var_cA + var_cD + cov2_comp` closes for any number of additive/dominance
  layers; the printed note states `realized - requested/V_P` (not
  `realized - requested`) and that the cross term is one-signed across loci only when
  counted-allele frequencies lie on one side of 0.5.
* The orthogonal-model `d` guard names both causes (no heterozygotes, or heterozygous
  in every individual).
* Documentation: the layer sub-seed is collision-resistant over ordinary ranges, not
  injective (31-bit); duplicate `chr`/`pos` markers are accepted.

**Selection, relationship and prediction**
* `optimum_contribution()` is now invariant to a constant added to all merits (merit
  is centred before tuning), and a `target_coancestry` above the unconstrained
  optimum by any floating-point-resolvable margin (e.g. 0.5000005 vs 0.5) warns and
  uses `lambda = 0`; `select_ind()` validates a scalar `trait` for every non-culling
  method.
* `a_matrix(founder_f = )` accepts pedigree keys; a display id shared by founders in
  different pools is an ambiguity error.
* `prediction_accuracy()` rejects duplicated names when only one of `ebv`/`truth` is
  named.
* Supplied-`K` symmetrization no longer flushes subnormal entries to zero.

**Transcriptome**
* See the `h2_realized` / `h2_allocated` / `h2_var_ratio` entry under Transcriptome
  above, which is corrected in this round.

**Genotype input and crossing**
* `as_numeric()` now records the allele coded `+1` per marker in the
  `"counted_allele"` attribute of its result (dosages unchanged). `as_population()`
  keeps it and `cross()`/`c.Population()` stop when two populations count different
  alleles at a marker, which the `allele` label alone could not reveal (an all-`AA`
  and an all-`GG` panel converted separately both encode `+1`). The record is not
  written to text files.
* `as_numeric()` on a VCF file path now applies the same complete-diploid, biallelic
  rule as an in-memory VCF: haploid, partially missing and multiallelic calls become
  missing with a counted warning (previously SNPRelate kept them, and the valid
  diploids at that marker could be coded with the wrong sign).
* Numeric-format input (file or data frame) is normalized to the schema every reader
  emits (`snp`/`allele`/`chr` character, `pos` integer, `cm` double); a table with no
  positions no longer round-trips with a logical `pos`.
* A 9- or 10-column HapMap prefix is no longer detected as HapMap.
* `.same_map()` is symmetric (crossing A x B and B x A agree on map identity), and
  chromosome labels that tie numerically (`"1"`, `"01"`) are ordered by label, so a
  seeded mating no longer depends on row order.
* `heterosis()` documentation: the one-half retention of F1 heterosis holds per
  locus under exact criteria. A backcross to breed A retains 1/2 iff A is in
  Hardy-Weinberg proportions at the locus; an F2 retains 1/2 iff the two
  Hardy-Weinberg deviations sum to zero (e.g. deviations -0.12 and +0.12 with
  neither breed in HWE). A single fixed inbred line qualifies; a mixture of
  inbred lines is not sufficient.

**Frozen v1 `create_phenotypes()`**
* `create_phenotypes(seed = NULL)` under a non-Mersenne caller generator (e.g.
  L'Ecuyer-CMRG) now advances the caller's random-number stream, so two successive
  calls give different results; explicit seeds still restore the caller's RNG kind
  and state exactly.
* `create_phenotypes(architecture = "LD")` rejects seeds above about
  `.Machine$integer.max / 10` up front (the marker search derives retry seeds
  `seed * s + ...`, s <= 10); other architectures are unaffected.
* The indirect-LD contract check now also verifies the LD magnitude of every selected
  pair against `[ld_min, ld_max]` and against the reported LD (previously latent; no
  valid public output changed: 15/30 successes for model "A" seeds 1-30 are
  identical).
* Corrected help: `h2 = 0.05` is rejected when `rep > 1` (accepted range is
  `h2 > 0.05`); the seed-collision rule of the fully pleiotropic QTN draw is
  `seed >= 2 * rep`; `cor` semantics (trait 1 unchanged only for a unit diagonal); NA
  dosages error; duplicated `chr_pos` documented as accepted by the pleiotropic and
  partially pleiotropic architectures.
* Removed a false "none of the dominance QTNs has a heterozygous individual" warning
  for a locus that is heterozygous in every individual.

### Review round 3 (Codex re-review of round 2)

A second independent review closed the following; items that change output or
reject previously accepted input are marked **(behaviour)**.

* `optimum_contribution()`: centring the merit no longer overflows for huge finite
  merit ranges, and the band that detects a `target_coancestry` above the
  unconstrained optimum is relative to the target and optimum coancestries (round 4
  removes the `max(1, .)` floor that remained in this band).
* Supplied-`K` symmetrization uses the correctly rounded mean for near-symmetric
  subnormal entries.
* Transcriptome `genes$h2_realized`, `h2_var_ratio`, `h2_allocated`,
  `cis_fraction_realized` and `epistasis_realized` are scale-free in
  `simulate_transcriptome()` and `predict()` (the absolute 1e-12 cutoff is removed;
  an exactly zero denominator still reports 0).
* The orientation label fallback (cross-pool allele guard) is case-insensitive.
* An 11-column HapMap-like object without sample columns is no longer detected as
  HapMap.
* Whole-number double genotype columns are converted to integer.
* The v1 seed-overflow error message states the inclusive magnitude bound
  `abs(seed) <= N` when the accepted seeds include 0 (the ordinary case); see round 4
  for calls whose accepted range is not centred at 0.
* Documentation corrections: the `sample_parents()` coancestry comparison is an
  empirical statement about the F2 example only (problem dependent); multi-generation
  scheme accuracy under per-generation re-standardization stays near sqrt(h2) for additive-only
  architectures and declines with dominance; `heterosis()` retention criteria are
  exact (backcross needs only the recurrent breed in HWE; F2 needs the two HWE
  deviations to sum to zero); `Var(g) = prop_A + prop_D + 2Cov(c_A, c_D)` is for one
  additive and one dominance layer (general: `Var(c_A) + Var(c_D) + 2Cov(c_A, c_D)`);
  `h2_realized = h2/(1 + gr_cov)` holds only when `Var(G) + Var(R) = 1` (general
  form `Var(G)/(Var(G) + Var(R) + gr_cov)`, not bounded by 1); the marginal
  epistasis-share formula holds for the nondegenerate blend (the exact-cancellation
  fallback gives the target epsilon).

### Review round 4 (Codex re-review of round 3)

* **(behaviour, message only)** The v1 seed-overflow error no longer claims an
  inclusive magnitude bound when none applies. When the residual seed
  `(seed + rep) * round(10 * h2)` makes the accepted seeds an interval not centred at
  0 (for example `rep = 429496730`, `h2 = 0.5`, where seed 0 is rejected but
  -429496730 is accepted), the message states the accepted integers as `[lo, hi]`
  (or that no seed is accepted) and asks to reduce `rep` / `n_qtn`. The ordinary
  message `abs(seed) <= N` is unchanged and valid V1 outputs are bit-identical.
* `optimum_contribution()`: the band that detects a `target_coancestry` above the
  unconstrained optimum is now purely relative (16 machine epsilons times the larger
  of the target and the optimum coancestry), with no absolute floor, so tiny-scale
  `G` matrices are judged on their own scale.
* `simulate_transcriptome(mimic = )`: the per-gene rescale always hits the requested
  per-gene variance whenever the realized unit-scale variance is finite and
  positive (no absolute cutoff); a warning is issued when that realized variance is
  tiny (< 1e-12) and the rescale is ill-conditioned (it amplifies rounding noise).
  Only an exactly zero or non-finite realized variance keeps the unscaled fallback.
* Documentation: the `select_ind(on = "pheno")` response
  `R = i * Cov(A, P) / sigma_P` is stated as the linear-regression prediction (exact
  only if `E[A | P]` is linear, e.g. joint normality), reducing to `i * h2 * sigma_P`
  only when `Cov(A, P - A) = 0` (all non-additive parts, including epistasis under
  linkage disequilibrium, uncorrelated with the breeding value); `docs/DECISIONS.md`,
  `docs/THEORY_REVIEW.md` and `docs/ROADMAP.md` mirror the wording.

# simplePHENOTYPES 2.0.0

Version 2.0 is the release line that introduces the v2 simulation grammar and
the multi-generation / selection engine, alongside the frozen v1
`create_phenotypes()` (unchanged). It is the version `breedingDesigner` depends
on (`Imports: simplePHENOTYPES (>= 2.0)`); the exported engine surface is the
backend contract in `docs/BACKEND_CONTRACT.md`.

## Major changes
* New phenotype grammar: `simulate_phenotype()` with composable `additive()`,
  `dominance()`, `epistasis()` and `vqtl()` layers, `complex_phenotypes()`, and
  long/wide exporters. Orthogonal genotypic model via
  `additive(orthogonal = TRUE, a =, d =)`.
* Multi-generation genetics (the crossing schemes): `as_population()`,
  `cross()`, `selfcross()`, `double_haploid()`, `synthetic_map()`, `dosages()`
  — isqg-parity meiosis on a Rust core.
* Selection engine: `select_ind()` and the named schemes
  `single_seed_descent()`, `bulk()`, `pedigree()`, `recurrent_selection()`.
* Modern methods: `g_matrix()` (VanRaden), `optimum_contribution()` +
  `sample_parents()` (Meuwissen OCS), `cross_usefulness()`.
* Fixed-scale cross-generation accessors `additive_value()` and
  `phenotype_value()`.
* `filter_geno()` LD pruning (`indep_pairwise`, `indep`, `indep_pairphase`,
  Gabriel `blocks`) is byte-exact to PLINK 1.9.

## Notes
* `create_phenotypes()` is retained unchanged as frozen legacy (bugfix-only).

# simplePHENOTYPES 1.4.0
## Major changes
Implemented vQTL simulation
removed "Selected" from QTN output file name.

## Minor changes
Fixed bug that changed the working directory after the simulation
Fixed bug that make it stop when running h2 as a matrix
Fixed a bug that stopped the search for QTNs that fit ld_min and ld_max

# simplePHENOTYPES 1.3.1
## Minor changes
Fixed reference paper information

# simplePHENOTYPES 1.3.0
## Major changes
Implemented the parameter "ld_max" (replacing "ld") and "ld_min".
## Minor changes
Fixed bug that changed the working directory after the simulation

# simplePHENOTYPES 1.2.16
## Major changes
Included the parameter 'mean', so traits can be simulated with the desired mean (intercept) value.
Included QTN_list option for the LD architecture. 
set default seed generator as RNGversion('3.5.1') to ensure reproducibility.
## Minor changes
Renamed some output QTN info files to make it standard across different architectures

# simplePHENOTYPES 1.2.15
## Major changes
Included the parameter QTN_list = list(add = NULL, dom = NULL, epi = NULL) to give the user the possibility to select the specific markers to be used as QTNs.
## Minor changes
Included 'Master Seed' in the log file to facilitate reproducibility. Now it only saves individual seed numbers when verbose = TRUE (default).

# simplePHENOTYPES 1.2.14
## Minor changes
check if 'out_geno' is either 'numeric', 'plink' or 'gds'
replaced the dependence lqmm and uses the function make_pd() to make cor matrix positive definite

# simplePHENOTYPES 1.2.13
## Minor changes
Fixed bug when reading multiple files using geno_path

# simplePHENOTYPES 1.2.12
## Minor changes
Fix bug that stopped simplePHENOTYPES when using geno_obj and architecture = "LD"

# simplePHENOTYPES 1.2.11
## Minor changes
Fix bug in the simulation of single trait using multiple h2 values 

# simplePHENOTYPES 1.2.10
## Minor changes
Fix bug that made the direct LD option stop running


# simplePHENOTYPES 1.2.9
## Minor changes
Fixed bug that also removed the cause of LD when remove_QTN = TRUE with architecture = "LD"
Fixed bug when reading multiple files using geno_path

# simplePHENOTYPES 1.2.8
## Minor changes
Set all additive parameters to NULL when model is dominance or epistasis.

# simplePHENOTYPES 1.2.7
## Minor changes
Fixed bug when more than 9 traits were simulated under the "partially" architecture
Fixed bug when saving file name with very large name due to a large number of traits

# simplePHENOTYPES 1.2.6
## Minor changes
Fixed bug in the QTN MAF calculation on the LD architecture
Fixed bug when importing VCF and exporting BED files (implemented by chr_prefix)


# simplePHENOTYPES 1.2.4
## Major changes
**Input**
1. Implemented options for input format as VCF, plink bed/ped files, GDS.
1. Changed dosage (numeric format) information from 0, 1, and 2 to -1 (aa), 0 (Aa) and 1 (AA).
1. Implemented a new type of spurious pleiotropy, direct LD (type\_of\_ld = "indirect").
1. Included the option for assigning a residual correlation among traits.
1. Implemented a constrain option to select only heterozygote or only homozygote QTNs.
1. Included the warning\_file\_saver option to skip asking if the user wants to save one genotype file for each rep when vary\_QTN = FALSE.

**Output**
1. Included a new output file with the summary linkage disequilibrium information on the selected spurious pleiotropy QTNs.
1. Included MAF in the outputted QTN information file.
1. Calculates the proportion of phenotypic variation explained by each QTN (QTN\_variance = TRUE).
1. Includes the option to remove QTNs from the genotype file (remove_QTN = TRUE).
1. Renamed <Taxa> by <Trait> in Tassel output format.


## Minor changes

Fixed bug that didn't recognize geno\_obj as HapMap.
Fixed bug when simulating dominance will all SNPs being homozygotes.
Fixed bug when reaching the end of the file while looking for SNPs in LD.
Fixed bug in importing geno\_file from other directories.
Fixed bug in selecting QTNs when marker data < 6 SNPs.
Included file removal when simulation does not complete.
Renamed file outputted as numeric.
Renamed constrain option.
Implemented check for biallelic markers.
Changed the QTN file name.
Implemented an interactive question before Check to remove QTNs with vary\_QTN = T.
Included check.names as FALSE in all data.frames.
Check if geno\_file and geno\_path are NULL.
Incorrect output name.

