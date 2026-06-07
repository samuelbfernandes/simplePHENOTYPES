# ARCHITECTURE.md — simplePHENOTYPES Redesign
> Status: DRAFT — to be completed by Claude Code auditing R/ source
> Working directory: ~/Library/CloudStorage/OneDrive-UniversityofArkansas/UARK/collaboration/Software/simplePHENOTYPES

---

## 1. Project Background

simplePHENOTYPES (Fernandes & Lipka, BMC Bioinformatics 2020) simulates pleiotropic,
linked, and epistatic phenotypes from real genomic marker data. It supports additive,
dominance, and epistatic genetic architectures across multiple traits with controlled
degrees of pleiotropy (full, partial, spurious).

**Current problems driving this redesign:**
- Monolithic structure: one large script split into many tightly coupled functions
- Utility functions (format converters, etc.) are embedded but useful independently
- Single-generation only: no support for simulating crosses across generations
- No expression-based simulation
- Hard to maintain and extend
- Bugs present in current R source (see docs/BUGS.md)

---

## 2. Current Function Inventory
> TODO: Claude Code should complete this section by reading @R/
> Format: function_name | category | exported? | has_tests? | Rust_candidate?

### Known categories (to be completed):
- **Format converters**: VCF → numeric, HapMap, other marker formats → (-1, 0, 1) encoding
- **QTL/marker selection**: selecting causal loci from marker data
- **Effect size assignment**: additive, dominance, epistatic effect sampling
- **Phenotype generation**: computing genetic values + residuals + heritability scaling
- **Pleiotropy control**: full, partial, spurious pleiotropy logic
- **I/O utilities**: reading/writing phenotype and genotype files
- **Simulation orchestration**: main `create_phenotypes()` wrapper

---

## 3. Proposed Package Decomposition

### Option A: Two-package split (RECOMMENDED starting point)
```
genoUtils/          ← new standalone CRAN package
  - Format converters (VCF, HapMap, PLINK, numeric)
  - Marker data validation
  - No phenotype simulation logic

simplePHENOTYPES/   ← refactored, depends on genoUtils
  - Single-generation phenotype simulation (current capability, cleaned up)
  - Multi-generation simulation via isqg integration (new)
  - Expression-based simulation (new)
  - Imports genoUtils for format handling
```

### Option B: Single package with namespace separation
Keep everything in one CRAN package but organize into clear internal modules.
**Drawback:** doesn't solve the "utility functions useful to non-phenotype users" problem.

### Option C: Three packages
```
genoUtils/          ← format converters
simplePHENOTYPES/   ← single-generation phenotype sim
simBreed/           ← multi-generation breeding simulation
```
**Drawback:** more maintenance surface; premature separation before multi-gen is designed.

> DECISION NEEDED: See docs/DECISIONS.md DECISION-001

---

## 4. isqg Integration Design

isqg (Peixoto et al., G3 2019) is a binary framework for in silico quantitative genetics.
It implements:
- Binary strand-based genome encoding
- Recombination simulation via bitwise operations (implemented in C++)
- Mating/crossing schemes (F1, F2, DH, RIL, etc.)
- Custom recombination models via user-supplied C++ extensions

**Critical constraint: isqg is in the CRAN Archive (no longer actively maintained).**

### Integration options:
1. **Depend on archived version** — fragile; CRAN may reject packages depending on
   archived packages; no future bug fixes
2. **Fork isqg and maintain internally** — more work but controllable
3. **Absorb isqg C++ core into our Rust port** — cleanest long-term; bitwise operations
   map naturally to Rust; we own the code; no archived dependency

> DECISION NEEDED: See docs/DECISIONS.md DECISION-002

### Multi-generation simulation interface (draft):
```
Single generation (current):
  real_marker_data → select_QTL → assign_effects → generate_phenotype

Multi-generation (proposed):
  founder_genomes → define_cross_scheme → simulate_generations(n) →
  [each generation: recombination + selection + phenotype_generation]
```

---

## 5. Expression-Based Simulation Design

Current: simulation uses marker/SNP data as input.
New: simulation from gene expression data as input.

**Key differences to spec out:**
- Expression values are continuous, not discrete (-1, 0, 1)
- "Causal loci" become "causal transcripts/genes"
- Effect sizes operate on expression levels, not allele dosages
- Correlation structure between expression traits differs from LD in markers

> TODO: Requires separate SPEC section — leave as open question until marker-based
> redesign is complete.

---

## 6. Rust Port Candidates

Priority order (easiest → most complex):

| Priority | Function category | Rationale |
|----------|-------------------|-----------|
| 1 | Format converters (VCF → numeric, etc.) | Pure computation, no R idioms, high reuse |
| 2 | Effect size sampling | Numerical, statistical distributions |
| 3 | QTL/marker selection | Array operations, good Rust fit |
| 4 | Phenotype generation loop | Depends on above; ports after they're stable |
| 5 | Recombination simulation (isqg port) | Complex but critical for multi-gen |

**R-specific behaviors requiring explicit Rust handling:**
- `NA` values in marker data (no native equivalent in Rust)
- R's vectorization semantics (recycling rules)
- `set.seed()` reproducibility — must document Rust RNG seeding strategy
- S3/S4 class interfaces if any are used in existing API

---

## 7. Python API Design Notes

- Python package via PyO3 + maturin (see PROJECT_CONTEXT.md)
- Public API should mirror R function names where possible for cross-referencing
- NumPy arrays as natural equivalent to R numeric vectors
- pandas DataFrame as equivalent to R data.frame output

---

## 8. CRAN Constraints (applies to all packages in this ecosystem)

- Rust deps must be vendored: `rextendr::vendor_pkgs()` before each submission
- `DESCRIPTION` SystemRequirements: `Cargo (Rust's package manager), rustc`
- Bundle size limit: 5MB — format converter crate may be large; watch this
- All packages must build offline
- No nightly Rust features
- If depending on isqg: CRAN may flag archived dependencies → prefer Option 3 above

---

## 9. Open Design Questions

These must be resolved in docs/DECISIONS.md before writing SPEC.md:

1. **Package decomposition**: Option A, B, or C? (Section 3)
2. **isqg strategy**: depend, fork, or absorb? (Section 4)
3. **Package name for utilities**: `genoUtils`? `markerTools`? `genomicIO`?
4. **Backward compatibility**: does the refactored simplePHENOTYPES maintain the current
   `create_phenotypes()` API exactly, or is this a breaking v2?
5. **Multi-generation scope**: is this v2 of simplePHENOTYPES or a new package?
6. **Expression simulation**: v2 feature or deferred to v3?
7. **Minimum R version**: current is 4.1 — keep or raise?

---

## 10. Bootstrap Prompt for Claude Code

Once you've answered the open questions above, use this in Claude Code plan mode:

```
Read @docs/ARCHITECTURE.md and @docs/DECISIONS.md.
Complete Section 2 (Function Inventory) by auditing @R/ — list every exported
function, its category, whether it has tests, and whether it's a Rust port candidate.
Then update Section 3 to reflect the decisions made. Do not write any code.
Show me the completed sections for review before saving.
```
