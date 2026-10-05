# BACKEND_CONTRACT.md — the simplePHENOTYPES engine surface

> simplePHENOTYPES is the **single computational backend** for breeding-program
> simulation. `breedingDesigner` (repo `fernandes-lab/breeding_designer`) is the
> **frontend + orchestration** and depends on this package
> (`Imports: simplePHENOTYPES (>= 2.0)`, DECISION-018). This file is the contract
> between them.

## The rule

There is **exactly one implementation** of every genetics, phenotype, and
selection computation, and it lives here, exported. Consumers (breedingDesigner,
and any other) **call only the exported functions listed below** — never a
`simplePHENOTYPES:::` internal, and never a reimplementation of engine genetics
in the consumer. If a consumer needs something the engine does not yet export,
the fix is to export it here (and bump the version), not to duplicate it
downstream.

Duplication that this contract exists to prevent (and retire):
`breedingDesigner`'s `.bv_index()` re-derives the fixed-scale additive value; it
must delegate to `additive_value()` unconditionally once `>= 2.0` is the pinned
minimum (see that repo's roadmap).

## Versioning (SemVer, as a dependency boundary)

- **Patch / minor** (`2.0.x`, `2.x.0`): additive — new exported functions or new
  optional arguments. Existing signatures and documented behavior are preserved.
  Safe for a consumer pinned at `>= 2.0`.
- **Major** (`3.0.0`): a breaking change to any function in this contract — a
  removed export, a renamed/removed argument, or a documented-behavior change.
  Consumers bump their pin deliberately.

`2.0.0` is the first version to declare this contract; it is where the crossing
schemes and the selection engine are implemented. `breedingDesigner` pins
`Imports: simplePHENOTYPES (>= 2.0)`.

The machine-checkable form of this list is
`tests/testthat/test-backend-contract.R`; that test fails if any contract
function stops being exported, which surfaces a break here before it reaches a
consumer.

## The contract surface

### Populations & crossing (multi-generation genetics)
`as_population`, `population_from_haplotypes`, `founders_coalescent`, `haplotypes`, `cross`, `selfcross`,
`double_haploid`, `dosages`, `n_individuals`, `synthetic_map`, and the `Population`
methods `[`, `c`, `print`. `population_from_haplotypes(cis, trans, map, ids, pool,
individuals_in_rows)` and its inverse `haplotypes(pop)` are the public known-phase
constructor and accessor (SPEC-0020 item 5, DECISION-039): 0/1 matrices, markers x
individuals by default, 1 = the counted (+1) allele, dosage = `cis + trans - 1`;
`map$counted` is kept only if the supplied map has it. A consumer must call these
instead of replacing the `cis`/`trans` slots of an `as_population()` result.
`founders_coalescent(n_ind, n_chr, seg_sites, inbred, species, split, theta, rho,
history, morgans, bp, pool, seed)` (SPEC-0020 item 9, DECISION-050) simulates founders with
historical LD natively (SMC', the model `runMacs()` runs) with the `runMacs()` GENERIC / MAIZE /
WHEAT / CATTLE presets; unlike `runMacs()` it is reproducible under `set.seed()` / `seed`.
Appended optional arguments (SPEC-0020 items 3, 6, 8): `interference = NULL` on `cross`,
`selfcross`, `double_haploid`, `mate`, `crossbreed` and, as the last formal, on `single_seed_descent`,
`bulk`, `pedigree`, `recurrent_selection`, `cross_usefulness`, `combining_ability` (simulated only;
an error with `method = "expected"`) and `progeny_test`; `reps = 1` on `simulate_phenotype`
and `complex_phenotypes`; `n_per_family` on `select_ind`; and a `source` attribute on
the `sample_parents()` result. Each is appended with a default that keeps today's
output and random stream.
Pedigree and mating plans (2.0.0.9001): `parentage`, `families`, `mating_design`,
`mate` (DECISION-024/025).
Crossbreeding (2.0.0.9001): `crossbreed`, `breed_composition`, `heterosis`
(DECISION-031).

### Phenotype grammar
`simulate_phenotype`, the layers `additive` (whose `effect` also takes a
per-trait list, added in 2.0.0.9000), `dominance`, `epistasis`, `vqtl`,
`complex_phenotypes`, `genetic_values`, `qtn_table`, the exporters
`phenotypes_long`, `phenotypes_wide`, `write_phenotypes`, and
`plot.phenotype_sim`.

### Selection engine
`select_ind` and the named schemes `single_seed_descent`, `bulk`, `pedigree`,
`recurrent_selection`. Added in 2.0.0.9001: `select_ind(method = "culling")` and
tandem selection (a `trait` vector on `pedigree` / `recurrent_selection`,
DECISION-028); `combining_ability`, `template_effects` (DECISION-026);
`progeny_test` (DECISION-027); `marker_select` (DECISION-029); `predict_ebv`,
`a_matrix`, `prediction_accuracy`, and the operator manifest `selection_methods`
(DECISION-030). Added in 2.0.0.9002: multi-trait BLUP through `predict_ebv` (an
individuals x traits `pheno` matrix with `var_a` / `var_e` covariance matrices,
DECISION-032).

### Modern methods
`g_matrix` (VanRaden), `optimum_contribution` + `sample_parents` (Meuwissen
OCS), `cross_usefulness`, and `print.ocs`; marker-assisted backcrossing
`mabc_select` + `recurrent_parent_recovery` (added in 2.0.0.9000).

### Fixed-scale cross-generation accessors
`additive_value`, `genotypic_value`, `phenotype_value` (which gains `d =` in
2.0.0.9001, DECISION-026, and the appended G x E arguments `gxe`, `gxe_intercept`,
`env`, `var_env` with the slope accessor `gxe_value`, DECISION-049) — the fixed-scale scorers
a downstream recurrent driver needs so a selection response is visible across
generations (DECISION-020 / DECISION-021). `genotypic_value()` is each
individual's own **per se** additive-plus-dominance total genotypic value
`G = A + D` (`a_j * dosage + d_j * (dosage == 0)`). A reciprocal-recurrent /
hybrid program uses it by scoring the **realized testcross / hybrid progeny** (so
dominance drives the progeny mean); it is a per se genotypic value, not a
parental SCA/GCA estimate, and not the transmissible breeding value
(use `select_ind(on = "bv")` for that).

### Genotype ingestion / QC
`as_numeric`, `filter_geno`.

## Crossing-layer preconditions and error contract

- **Errors, never crashes.** Every input the Rust kernel cannot honour (strand
  length or alphabet against the marker layout, a meiosis-event budget different
  from `n_prog * events_per_progeny`, non-finite or out-of-range chiasmata, flips
  other than 0/1, an unknown `design`, `code_as`, `model` or `impute`, a short or
  missing `flip`, raw dosages outside 0/1/2/NA) is an ordinary R error, raised
  before any output exists. A Rust panic would abort the R process on toolchains
  whose unwinder cannot cross R's frames (e.g. a gcc-linked macOS build), so the
  kernel entry points return `Result` (extendr feature `result_list`) and
  `R/extendr-wrappers.R` re-raises the error; regenerating that file with
  `rextendr::document()` drops the unwrapping and is caught by
  `tests/testthat/test-audit-rust-crossing.R`.
- **Allele orientation.** The -1/0/1 dosages are relative to the allele
  `as_numeric()` coded `+1` (by default the most frequent allele of that data
  set). `as_numeric()` records that allele per marker in the `"counted_allele"`
  attribute of its result (absent under `model = "Dom"`; `NA` where unknown);
  `as_population()` keeps it as `map$counted`. `cross()`, `c.Population()` and the
  selection schemes stop when both populations carry it and count different alleles
  at any marker. `as_numeric(counted_column = TRUE)` also writes the record as a character
  column `counted` immediately after `cm` (in the data frame and the numeric text file;
  default off, output unchanged). `as_population()` reads it (when present it is the record
  used; the `"counted_allele"` attribute, an R-object convenience that `[` does not realign,
  is then ignored), so files and row-subsetted panels keep the strong per-marker check.
  Without the record on both sides (numeric files written without the column, other
  software) the numeric `allele` label is compared (warning on opposite order, error on
  disjoint alleles). A supplied `counted` column is validated like `map$counted`; a numeric
  sixth column named `counted` is an individual. Crossing separately
  converted panels is still unsafe unless they are converted jointly or with the
  same `ref_allele` (`as_numeric(method = "reference", ref_allele = )`). A supplied
  `map$counted` is validated (character; one non-empty allele symbol per marker, `NA` =
  unknown; `""` and numeric columns are rejected). The map
  identity judgement (`.same_map()`) uses a symmetric relative tolerance
  (`1e-8 * max(1, |x|, |y|)`).
- **Map identity.** Two populations share a map when marker names, chromosome
  labels (compared as text) and `pos`/`cm` (element-wise relative tolerance 1e-8, symmetric in the two maps) agree; the
  `allele` column is not part of it. The same judgement is used for crossing,
  breed lists and pooling.
- **Chromosome order and the seeded stream.** Chromosomes are processed in a
  canonical, locale-independent order (numeric labels first, in numeric order;
  then other labels by prefix in byte order and trailing number), so integer and
  text `chr` give the same seeded progeny. Labels that tie on every canonical key
  (`"1"`, `"01"`) are ordered by the label in byte order, so the order is total. Text labels that used to sort as
  `"1", "10", "2"` are now `1, 2, 10`: seeded output changes only for such maps.
- **Length of a chromosome** for the Poisson crossover count (`interference = NULL`) is its *last* map
  position in Morgans (isqg convention), not its span; `cm` must be in
  centiMorgans (a map that looks like Morgans draws a warning).
- **Crossover interference.** `cross`, `selfcross`, `double_haploid`, `mate` and
  `crossbreed` take a trailing `interference = NULL` (appended, default unchanged); so do
  `single_seed_descent`, `bulk`, `pedigree`, `recurrent_selection`, `cross_usefulness`,
  `combining_ability` (simulated only) and `progeny_test`, which forward it to every
  meiosis they draw (DECISION-041). `NULL` is
  the Poisson model and the isqg random stream, bit-identical to versions without the
  argument; `list(nu =, p =)` (1 <= nu <= 1e6, p in [0, 1]) is the two-pathway gamma model
  (`?cross`, DECISION-041), drawn in R, consuming its own stream. The kernel is unchanged:
  it receives sorted chiasma positions in [0, L] under the same `counts`/`flips` contract.
  When `interference` is `NULL`, the option `simplePHENOTYPES.interference` (a `list(nu =, p =)`) is
  used if set (explicit argument, then option, then Poisson); unset leaves the Poisson/isqg stream
  bit-identical. `NULL` cannot switch a set option off for one call.
- **Batched meiosis.** The internal `mate_many_core()` (integer strands in and out, a
  mating table, one shared event stream consumed in mating order; DECISION-040) is what the
  crossing functions use; a batch equals running its matings sequentially, bit for bit.
  `mate_haplotypes_core()`, `meiosis_core()` and `gamete_masks_core()` are unchanged. The
  signature manifest in `tests/testthat/test-backend-contract.R` lists the appended
  `interference`; `cross`, `selfcross`, `double_haploid`, `mate`, `crossbreed`,
  `single_seed_descent`, `bulk`, `pedigree`, `recurrent_selection`, `cross_usefulness`,
  `combining_ability` and `progeny_test` keep their earlier formals in order (`interference`
  is the last).
- **Frozen signatures.** `create_phenotypes()` is untouched (DECISION-008). Every other
  change in this round is an appended optional argument or a new export, so a consumer
  pinned at `>= 2.0` is unaffected.
- **RNG.** `cross`, `selfcross`, `double_haploid`, `mate`, `mating_design(random)`
  and `crossbreed` restore the caller's RNG state when `seed =` is given.
- **Build profiles.** The development build (`DEBUG` set at install) is the debug
  profile and keeps `debug_assert!` live; the CRAN build is release. Validation
  is explicit and profile-independent, so both behave identically on invalid
  input. The declared minimum Rust version is 1.71 (that of the pinned
  `extendr-api 0.9.0`); `DESCRIPTION` must state the same.

## Not part of the contract

- `create_phenotypes()` — frozen v1 legacy (DECISION-008), bugfix-only. Consumers
  should build on the grammar (`simulate_phenotype()`), not this.
- Anything not exported (`:::`) — internal, may change without a version bump.
