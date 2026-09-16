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
`as_population`, `cross`, `selfcross`, `double_haploid`, `dosages`,
`n_individuals`, `synthetic_map`, and the `Population` methods `[`, `c`, `print`.

### Phenotype grammar
`simulate_phenotype`, the layers `additive`, `dominance`, `epistasis`, `vqtl`,
`complex_phenotypes`, `genetic_values`, `qtn_table`, the exporters
`phenotypes_long`, `phenotypes_wide`, `write_phenotypes`, and
`plot.phenotype_sim`.

### Selection engine
`select_ind` and the named schemes `single_seed_descent`, `bulk`, `pedigree`,
`recurrent_selection`.

### Modern methods
`g_matrix` (VanRaden), `optimum_contribution` + `sample_parents` (Meuwissen
OCS), `cross_usefulness`, and `print.ocs`.

### Fixed-scale cross-generation accessors
`additive_value`, `genotypic_value`, `phenotype_value` — the fixed-scale scorers
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

## Not part of the contract

- `create_phenotypes()` — frozen v1 legacy (DECISION-008), bugfix-only. Consumers
  should build on the grammar (`simulate_phenotype()`), not this.
- Anything not exported (`:::`) — internal, may change without a version bump.
