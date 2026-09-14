# ROADMAP.md — simplePHENOTYPES

> Triage of the accumulated wish list, organized by topic and by whether it
> should land before this version is pushed. Status reflects the `restructure`
> branch as of 2026-09-09.
>
> Companion docs: `SPEC.md` (API contracts), `DECISIONS.md` (locked rationale),
> `NEXT_STEPS.md` (rolling pointer), `BUGS.md`.

Legend: `[x]` done · `[ ]` open · **Blocker** = do before pushing ·
**Later** = worth doing, not now · **Drop** = not worth building ·
**Obsolete** = the v2 architecture removed the need.

---

## 1. Already done

Items from the list that the v2 work has closed.

- [x] **Null / purely residual trait.** `simulate_phenotype(geno)` with no
      layers is exactly this: h² = 0, pure noise. It is the foundation of the
      grammar rather than a special option.
- [x] **All effects for each trait in one table with the QTNs.**
      `qtn_table()` returns `trait`, `layer`, `set`, `snp`, `chr`, `pos`,
      `maf`, `effect` in a single frame, across every layer and trait.
- [x] **Date/time in the log.** Superseded rather than implemented: nothing is
      written to disk, and the object carries `seed`, `n_traits`,
      `var_budget`, QTNs and effects. Provenance belongs to the script, not to
      a log file the package writes.
- [x] **Unit tests.** 271 tests: statistical validation of the grammar, exact
      isqg parity, the v1.3.0 frozen-legacy regression guard, format
      conversion, crossing genetics.
- [x] **Simulate based on a genetic map.** `synthetic_map()` builds one;
      `cross()` / `selfcross()` / `double_haploid()` use it for meiosis.
- [x] **`create_complex_phenotype()` — combine several architectures.**
      Shipped as `complex_phenotypes()`.
- [x] **Filter non-polymorphic SNPs.** `constraint()` / `maf_above` in the
      legacy engine; the grammar's `.candidate_markers()` respects MAF.
- [x] **Numericalizing HapMap with numeric columns** (bug). Handled: numeric
      HapMap input is detected and passed through.
- [x] **Multi-file output overwriting** (bug). Every run now lands in its own
      folder (`simplePHENOTYPES_output`, then `(1)`, `(2)`, …), so re-running
      never overwrites.
- [x] **Relative `home_dir`.** Works; output resolves under whatever path is
      given.
- [x] **`mean` / intercept per trait, and dropping the manual
      `pheno$Trait_1 <- pheno$Trait_1 + gv$gv_mu[1]` lines.** Done in the
      legacy engine (`mean = c(1, 1)`). **Still open for the grammar** — see
      §3.

---

## 2. Blockers — do before pushing

Small, and each one is either a correctness problem or something that will
embarrass the release.

- [x] **Finalize the PleioArch citation string.** Now reads
      `Prado et al. (in preparation).` (was a typo'd, unbalanced placeholder).
- [x] **Add the citation to the GitHub landing page.** README + `inst/CITATION`
      (modern `bibentry`) now list both papers; `citation()` shows the PleioArch
      entry alongside the 2020 paper.
- [x] **`vary_qtn`** now redraws an independent QTN set (loci + effects) per
      replication; `same_as_add` layers follow the additive loci per rep. (Was
      accepted and silently ignored.)
- [~] **Print the correlation when several h² values are simulated** (bug, v1 engine). Could not reproduce with current code -- passes on multi-h2 matrix + cor. Needs the exact failing call to confirm.
- [x] **Multi-file read messaging.** Each format already prints its own header; fixed the `message(files, sep=)` bug (message ignores `sep`) so filenames are newline-separated.
- [~] **Saving from `geno_path` (multiple files)** (bug). Could not reproduce -- multi-file HapMap read + save to disk works. Needs the exact failing case.

---

## 3. Worth doing soon — after the push, before the paper

Design work is understood; each is contained.

- [ ] **Validate `filter_geno()` LD pruning against PLINK 1.9 — next, and load-bearing.**
      `filter_geno()` now implements pairwise composite-r^2 pruning, phased
      (EM-haplotype) r^2 pruning, VIF pruning (PLINK `--indep`) and Gabriel-block
      tag selection. Because it will be used constantly, it must be **verified
      against PLINK 1.9 on the same dataset** so the kept/removed marker sets (and
      the block boundaries) match, not just look plausible — build a fixture that
      runs PLINK `--indep-pairwise` / `--indep` / `--blocks` on
      `SNP55K_maize282_maf04` and diffs the marker lists. Pressure-test edge cases
      (windows in variants vs kb, monomorphic/low-MAF loci, chromosome ends, the
      double-het EM at extreme allele frequencies). **If the R implementation is
      too slow** at whole-genome scale (the pairphase EM and the O(m^2)-within-window
      Gabriel classification are the suspects), port the deterministic inner loops
      to Rust — this fits the Rust boundary (DECISION-006): the r^2 / D' / EM
      computation and the block scan are pure and deterministic, so they may move
      to Rust while nothing about QTN sampling does. Keep the R version as the
      parity reference for the Rust port, mirroring the isqg approach.

- [x] **`mean` per trait in the grammar.** `simulate_phenotype(mean=)`, applied
      at the phenotype level; genetic values stay centered.
- [x] **Supply QTNs for one layer and randomize the rest.** Every layer takes
      `qtn =` (names or indices; epistasis takes a matrix); the others stay
      random, and fixed loci are exempt from `vary_qtn`.
- [x] **Proportion of variance explained by particular QTNs.** `qtn_table()`
      gains a `var_explained` column (additive/dominance; NA for epistasis,
      where variance is not per-locus).
- [x] **Subset individuals.** `simulate_phenotype(individuals=)` restricts the
      set (MAF recomputed on the subset; no genotype copy).
- [x] **Read `.gz` files.** `detect_format()` sees through a `.gz`/`.bz2`
      suffix; `fread` decompresses transparently.
- [x] **Diagnostic plots.** `plot.phenotype_sim()` (base graphics, no new dependency): variance partition, phenotype histogram, QTN effects, and the two-trait genetic-value scatter (or target-vs-realized h2 for one trait).
- [x] **Simulate alleles in repulsion.** `additive(phase = "repulsion")`
      alternates effect signs. Realized variance is still pinned to `prop`; the
      sign structure and cross-locus covariance are what change.
- [x] **Independent traits with QTNs on different chromosomes.** `distinct_chr = TRUE` partitions chromosomes among traits.

### Variance-decomposition fidelity (from the QG review)

The grammar currently uses a naive dosage/indicator coding (additive = -1/0/1
dosage, dominance = heterozygote indicator, epistasis = uncentered product), not
Fisher's orthogonal decomposition. Realized H2 tracks the sum of props closely
but not exactly when dominance/epistasis share loci with the additive layer at
non-0.5 allele frequencies. Improvements, roughly in order of effort:

- [x] **Centered additive-by-additive epistasis.** Each locus's design term is
      centered before the product, removing the first-order leakage into the
      additive main effects (`.component_raw()`, `R/grammar_realize.R`).
- [x] **Additive-by-dominance (a x d) and dominance-by-dominance (d x d)
      epistasis.** `epistasis(interaction_type=)` takes `"a"`/`"d"` per
      interacting position (`c("a","a")` default, `c("a","d")`, `c("d","d")`),
      each term centered. Note `"d"` terms are near-degenerate on a fully inbred
      panel (documented). The full-orthogonality caveat below still applies.
- [x] **Orthogonal genotypic model / average effects (breeding value /
      dominance deviation).** Shipped as `additive(orthogonal = TRUE, a =, d =)`
      (DECISION-020). The layer builds each locus's genotypic value from a
      per-locus additive effect `a` and dominance deviation `d`
      (`-a`/`+d`/`+a` for gene content 0/1/2) and scales the whole value to
      `prop`; the additive (breeding-value) component uses the average effect
      alpha_j = a_j + d_j(1 - 2 p_j) and the dominance deviation is the realized
      residual g - A (Cov(A, D) = 0 in expectation under random mating → per-locus
      HWE; LD is fine, nonrandom multilocus association is not; nonzero for
      structured / finite samples), and the
      additive/dominance variance shares **emerge** from (a, d, p) rather than a
      freely chosen `prop` — reported as the realized fractions Var(A)/Var(g) and
      Var(D)/Var(g) on separate variance-budget rows, plus an `add_dom_cov` row
      2Cov(A,D)/Var(g) so the three close to `prop`. The **degree of dominance**
      `d/abs(a)` (magnitude 1 = complete, > 1 = overdominance; `abs` because `a`
      may be negative) is now meaningful because the whole genotypic value is
      scaled together (not per-component). The breeding value ([select_ind()] `on = "bv"`, the OCS
      default, and the index methods) picks up the dominance-induced average
      effect automatically. Incompatible with `vary_qtn`, a separate `dominance()`
      layer, and the correlation-controlling `"pleiotropy"`/`"ld"` architectures.

---

## 4. Long tasks — worth waiting for

Each is a project. All are defensible; none should hold up this release.

- [ ] **Polyploids (0/1 per homologous chromosome).** The deepest change on the
      list. It touches the genome representation, meiosis (multivalent pairing,
      double reduction), dosage coding, and every variance formula that assumes
      diploidy. The Rust `Bits` representation generalizes cleanly to more than
      two strands, which is the good news; the quantitative genetics does not
      generalize as cleanly, and that is the real work. Treat as v3.
- [ ] **Categorical / threshold traits.** `y = 1 if l < γ₁`, etc. Conceptually
      straightforward — simulate the liability, then cut it — and a genuinely
      common need (disease status, ordinal scores). The work is in the API:
      thresholds by prevalence or by cut points, and what h² means on the
      observed versus liability scale. Good first "long" task.
- [ ] **eQTL simulation.** Already deferred to v3 by DECISION-005, on the
      grounds that transcript-level input has different semantics.
- [ ] **Haplotype-based simulation.** Now much closer than it was: the crossing
      code already carries phased haplotypes, so effects assigned to haplotypes
      rather than alleles is a real possibility. Needs a design.
- [ ] **GxE / GxB as an extra random term added to the genetic values.** The
      note in the list is right that adding it to the genetic values avoids
      touching the phenotype function. Fits the grammar as another layer.
- [ ] **Different means per subpopulation (PCA).** Requires a population-
      structure concept the package does not have. Pairs naturally with GxE.
- [ ] **Correlated allelic effects for high-throughput phenotyping.** This is
      the pleiotropy engine applied to many traits with a structured
      correlation matrix — which now works for any `n_traits`. Mostly a
      question of specifying the correlation structure (AR1 over time points,
      say) rather than new machinery.
- [ ] **Co-heritability (rG·h1·h2).** Small computation, but needs a decision
      on where it is reported.
- [ ] **Power calculation.** Scope creep risk: this is a study-design tool, not
      a simulator. Worth doing only as a vignette showing how to use the
      simulator for power, not as a function.
- [ ] **Drag-and-drop pipeline designer** (the "Shiny app", scoped up). A
      Lucidchart-style node canvas for building a breeding cycle — drag `cross`,
      `selfcross`, `double_haploid`, a **selection** node, `simulate_phenotype`,
      wire them, then execute. See §8 for the decided sequence and the
      architecture (one JS/React-Flow frontend, a DAG-JSON contract, thin R and
      Python executors). Keep the app itself out of the CRAN tarball (§5).
- [ ] **Python package** (PyPI + conda-forge/Bioconda). Reimplements the grammar
      in Python over the *shared* Rust core (PyO3/maturin); the stochastic R code
      cannot be shared. See §8.
- [ ] **GWAS wrapper (GAPIT, GEMMA).** See §5 — recommend a companion package.
- [ ] **Environtyping-based simulation.** Needs a defined data model first.

---

## 5. Not worth implementing here

Not rejections of the idea — rejections of putting them *in this package*.

- **GWAS wrappers (GAPIT / GEMMA).** simplePHENOTYPES simulates; it should not
  take on the analysis side. Wrapping GAPIT means inheriting its dependency
  tree and its breakage. A companion package, or a vignette showing the
  handoff, gives users the same thing without the maintenance burden. The
  package already writes GEMMA-format output, which is the right level of
  support.
- **GDS → VCF conversion (SeqArray).** Format conversion in that direction is
  SeqArray's job and it does it well. `as_numeric()` exists to get data *in*,
  not to become a format hub. Point at SeqArray in the docs.
- **Shiny app inside the package.** A Shiny dependency in a CRAN simulation
  package costs every user install time and CRAN-check surface for something
  most will never launch. Ship it as a separate repo/app that depends on
  simplePHENOTYPES. The drag-and-drop designer (§8a) is the modern form of this
  item — same conclusion: the app lives outside the CRAN tarball.
- **Power calculation as a function.** As above: better as a vignette.

---

## 6. Obsolete under the new architecture

- **"Option to simulate a null trait directly."** The grammar makes this the
  default state rather than an option (§1).
- **"Include date/time in the log file."** There is no log file. Nothing is
  written unless the user asks; provenance lives in the returned object and in
  the user's script. Re-adding a log would reintroduce exactly the
  write-side-effects the v2 design removed.
- **"Fix h² and vary the QTNs" (check with Kaio).** Half-obsolete: fixing h²
  is now how the grammar works (`h2` is a budget, `prop` optional). What
  remains is varying QTNs across replications — which is the `vary_qtn` item
  in §2, and should be tracked there rather than here.
- **Manual mean-shifting lines** (`pheno$Trait_1 <- ... + gv$gv_mu[1]`). Gone:
  removed in the legacy engine, and the grammar never had them.

---

## 7. Cross-cutting: what actually gates CRAN

Not from the wish list, but these outrank most of it.

- [ ] **Prado et al. (PleioArch) must be published first.** The CRAN submission is
      deliberately held until the correlation-control paper is out, so the
      `.cite_pleioarch()` notice can point at a real citation. **Until then the
      package is distributed on GitHub only** (`remotes::install_github()`), which
      also sidesteps the Bioconductor-dependency friction of a CRAN install.
- [x] Vendored Rust dependencies (`vendor.tar.xz`, 548 KB).
- [x] Bundle under the 5 MB cap.
- [ ] `devtools::check_win_devel()` — never run.
- [ ] Maintainer-address change explained in `cran-comments.md` (written;
      needs to survive an actual submission, since the old address is
      unreachable).

---

## 8. Post-grammar deliverables: designer and Python

Two large forward-looking pieces. The **maintainer's chosen sequence**
(2026-09-09): push the current grammar to GitHub → build the designer, shipping
it *within version 2* → update CRAN once the designer is ready → Python last.

Recorded so the plan is not lost, together with the engineering risks flagged
against it. These are open questions to resolve, not settled decisions.

### 8a. Drag-and-drop pipeline designer

- **Entry point:** a package function **`design_breeding_program()`** (exported
  from simplePHENOTYPES) launches the designer. The function lives in the
  package so the feature is discoverable from `library(simplePHENOTYPES)`, but it
  may *serve or open a GitHub-hosted app bundle* rather than carry the whole
  frontend in the tarball — keeping the React/Shiny weight out of the package
  itself (per §5). So the app can be shipped and iterated on GitHub while the
  launcher ships with the package.
- **Architecture (recommended):** build the node canvas **once** as a JS app
  ([React Flow](https://reactflow.dev) is the de-facto library for this), whose
  only output is a **DAG serialized to JSON**. Two *thin* backends execute that
  JSON: R (Shiny or `plumber`) and, later, Python (FastAPI). Do **not**
  reimplement the canvas twice (once in Shiny reactivity, once in Django) — that
  is the trap.
- **The real design work is the DAG-JSON schema**, the contract shared by the
  frontend and both executors. Define it *first*; the canvas is mechanical once
  the schema is stable. The schema should be runnable headless (paste JSON →
  execute) before any UI exists.
- **Designer core — BUILT (2026-09-11, DECISION-017).** The DAG-JSON contract and
  its headless R executor exist ahead of any UI, in `R/designer.R`, tested
  (`test-designer.R`):
  - `run_design(design, seed, data_objects)` — validates, topologically sorts, and
    executes a design (an R list, or a JSON string/path via `jsonlite` in
    Suggests). Returns a `design_result` keyed by node id, with `order` and
    `terminal` attributes and a `print` method.
  - `validate_design()` — unique ids, known node types, resolvable/required
    inputs, acyclic (errors name the first problem or the cycle).
  - `design_breeding_program()` — the exported launcher: runs a supplied design
    headless, or (no argument) points to the schema and the companion canvas.
  - **Node types** (each = one tested package function): `founders`
    (`as_population`), `cross`/`self`/`dh` (crossing core), `ssd`/`bulk`
    (schemes), `phenotype` (grammar model), `select` (`select_ind`), `ocs`
    (`optimum_contribution` + `sample_parents`), `usefulness`
    (`cross_usefulness`), `pedigree`/`recurrent` (selecting schemes; the genetic
    model is a `model` sub-spec and the additive loci are **frozen** once so every
    generation is scored on the same causal loci). Reproducible from one `seed`.
  - So "paste JSON → run" works today with no server; the designer adds
    orchestration, not new genetics.
- **Designer UI — Phase 2 (working prototype built 2026-09-11).** A node canvas
  (palette grouped by bucket, drag-to-wire ports, per-node inspector, live R /
  Python / DAG-JSON with Copy, client-side validate + Run) is built as a
  self-contained web page (the JS codegen mirrors `design_script()`). It opens on a
  worked program and lets a user assemble a design and copy runnable code today.
  Per §5 and the CRAN 5 MB / bundled-JS risk it stays a **companion app**, not in
  the tarball. Live in-browser **Run** needs the R engine, which a sandboxed page
  cannot reach; `inst/designer/plumber_api.R` exposes `validate_design` /
  `run_design` / `design_script` over localhost so a locally served canvas can
  execute against the machine's own data. Remaining Phase-2 work: harden the canvas
  (React Flow if it earns its weight), wire Run to the plumber runner, host the
  bundle, and have `design_breeding_program()` open/serve it. Everything the UI
  produces round-trips through `validate_design()` + `run_design()`.
- **Selection backend primitive — BUILT (2026-09-10, DECISION-015).** The
  selection engine and generation-advance functions now exist as plain, headless,
  tested R code, ahead of any UI:
  - `select_ind()` — truncation on `on = "pheno" | "gv" | "bv" | <vector/function>`
    (the custom vector/function hook is the genomic/phenomic-selection extension
    point for the single-score methods; the multi-trait `index`/`quadratic_index`
    methods score on true breeding values and ignore `on`), with methods mass /
    within-family / among-family / combined (Lush index) /
    Smith–Hazel `index` / quadratic_index / random; intensity by count /
    proportion / standardized
    `i`; both directions; returns a crossable `Population` subset carrying the
    realized selection differential `S` and intensity `i`.
  - Named scheme wrappers `single_seed_descent()`, `bulk()`, `pedigree()`,
    `recurrent_selection()` compose `select_ind()` with the crossing primitives;
    `c.Population()` pools progeny. These map to the designer's "selection bucket"
    nodes (SSD / bulk / pedigree / recurrent), with personalized nodes wiring the
    primitives directly.
  - Efficiency: ranking + stochastic orchestration in R, meiosis in Rust
    (DECISION-006); Populations are reference-based (no genotype copies).
  Design decisions are logged in DECISION-015 and mirrored in the designer
  manuscript (`breeding_designer/manuscript.tex`).
- **Modern methods — BUILT (2026-09-11, DECISION-016).** `g_matrix()` (VanRaden
  2008 genomic relationship), `optimum_contribution()` (Meuwissen 1997 OCS, with a
  dependency-free Frank–Wolfe optimizer) + `sample_parents()`, `cross_usefulness()`
  (Zhong & Jannink 2007; Lehermeier 2017 — cross ranking by µ + iσ, families
  simulated), and the quadratic genomic selection index form of Cerón-Rojas et al.
  (2026, *Nat Commun* 17:1991) as `select_ind(method = "quadratic_index")`. The
  latter is a *simulation of the index's behaviour* on the known-truth breeding
  values (the Fisher average-effect projection), not the paper's GEBV estimator:
  the package fits no genomic-prediction model, and the multi-trait index methods
  (`index`/`quadratic_index`) score on all traits' true breeding values, so they
  do not take an external per-trait GEBV. Single-score genomic/phenomic selection
  (e.g. `method = "mass"` on externally computed GEBVs) needs no engine code — it
  enters via the `on` custom-criterion hook; the index methods do not use `on`.
- **Modern methods — NEXT.**
  - **BQP relatedness-minimizing multi-trait index** (Montesinos-López et al. 2026,
    *Plant Methods*, doi:10.1186/s13007-025-01484-4): a discrete (binary) multi-trait
    index that also minimizes genetic relatedness among the selected set, solved by
    binary quadratic programming. Overlaps `optimum_contribution()` (same gain↔
    diversity goal, discrete vs. continuous); open decision: exact MIQP solver
    dependency vs. a dependency-free greedy+swap heuristic on the existing
    `g_matrix()`.
  - **PopVar-style cross selection from real data** (Mohammadi, Tiede & Smith 2015,
    *Crop Sci* 55:2068): the real-data counterpart of `cross_usefulness()` —
    predict biparental progeny mean/variance/correlated response from **estimated
    marker effects** (from a training set) rather than simulation with known
    effects. Needs a marker-effect estimation step (or a supplied effect vector).
- **THEORETICAL-CORRECTNESS SCRUTINY GATE (load-bearing, before the designer/
  selection engine ships).** The *entire* designer and selection-engine
  implementation must pass a deep review of the **theoretical correctness** of every
  method before release — not just "the tests pass," but that each estimator/optimizer
  matches its published definition. Concretely: (a) selection response tracks the
  breeder's equation R = i·h²·σ_P across intensities and heritabilities; (b) the
  combined-index and Smith–Hazel weights reproduce worked textbook examples;
  (c) OCS contributions match a reference solver (e.g. `optiSel`) on the same G and
  merit, and realized ΔF matches the constraint; (d) `g_matrix()` matches an
  independent VanRaden implementation; (e) `cross_usefulness()` µ/σ match analytic
  expectations and an AlphaSimR-simulated benchmark; (f) QGSI reproduces the
  Cerón-Rojas et al. reference behavior. This validation gates the manuscript and the
  CRAN update; log outcomes in DECISIONS.md and the designer manuscript.
- [x] **Fixed-scale additive value accessor (for cross-generation scoring).**
      Shipped as exported `additive_value(x, qtn, effect)` (`R/cross_population.R`):
      returns `sum_j dosage_ij * effect_j` over a given frozen architecture on the
      -1/0/1 dosage scale, with **no per-population centring or rescaling**, so mean
      additive value is comparable across generations (a selection response shows).
      This is the counterpart to `genetic_values()` (which re-centres/re-scales each
      layer to `prop` on the scored population and therefore cannot show a
      cross-generational trend). The breeding designer's `program_metrics()` can now
      delegate to it instead of re-deriving `breedingDesigner:::.bv_index` from
      `.freeze_loci()` loci/effects and `dosages()` (removing the reimplementation
      the review flagged, constitution #6 / ADR-0004). Accepts a `Population`, a
      Population-backed `phenotype_sim`, or a dosage matrix; `qtn` by name or index.
- [ ] **Fixed-scale PHENOTYPE accessor (for cross-generation selection).** The
      counterpart to `additive_value()` for *phenotypic* selection: a phenotype
      whose genetic component is a FIXED function of genotype (frozen loci/effects)
      plus a residual on a fixed variance, WITHOUT the per-population rescaling that
      `simulate_phenotype()` applies (it holds the genetic layer at `prop` on the
      scored population). Needed so a downstream recurrent driver can run
      `on = "pheno"` selection on a trait whose heritability actually declines as
      variance is exhausted. Symptom in the breeding designer (`program_metrics()`
      recurrent DAG execution, review O1): genomic `on = "gv"` selection is correct
      (rescaling is monotone within a population, so the ranking is unchanged), but
      `on = "pheno"` selection keeps its genetic share at `prop` every cycle, so
      selection accuracy does not decay and late-cycle response is optimistic. A
      `phenotype_value(x, qtn, effect, h2, ...)` (fixed genetic scale + residual
      draw) — or an option on the existing scorer to skip the per-population
      rescale — would let the designer drive faithful phenotypic selection.
- [ ] **Bump the package version when `additive_value()` (and the fixed-scale
      phenotype accessor) ship.** The accessor was added without a version bump
      (still `1.4.0.9001`), so a downstream package cannot pin
      `Imports: simplePHENOTYPES (>= <ver>)` to require it — the breeding designer
      therefore keeps a runtime fallback that re-derives the additive value
      (`breedingDesigner:::.bv_index`), the very reimplementation the review flags
      (O4). Bump the dev/release version on the next tag so downstream can pin the
      minimum and delete the fallback.
- **Open risk vs. the chosen sequence:**
  - Shipping the designer *inside* the CRAN tarball conflicts with §5 ("Shiny
    app inside the package") and the 5 MB budget (vendored Rust is already
    548 KB; a React bundle adds more; CRAN scrutinizes bundled minified JS).
    A **companion repo/app**, or a `run_designer()` launcher that serves a
    GitHub-hosted bundle, avoids this — and then the designer need not gate any
    CRAN release.
  - Gating the CRAN *update* on the designer keeps the finished grammar (and the
    unreachable-maintainer-address fix, §7) unreleased for the whole build.
    Recommended alternative: release the grammar as **2.0 now**, ship the
    designer as **2.1 / companion**. Maintainer has chosen to ship within v2;
    revisit if the designer slips.

### 8c. Digital-twin platform (north-star destination)

**Destination, not a near-term deliverable.** Today simplePHENOTYPES + the designer
is a **breeding-program simulator** (generic/parameterized, effects *specified*), the
same category as AlphaSimR — **not** a digital twin. A digital twin models a
*specific real* breeding program and is kept *synchronized* with its data. Naming the
gradient keeps the claims honest:

| Stage | What it is | Gap from here |
|-------|-----------|----------------|
| **Simulation** | generic/parameterized in-silico model | — (this is what exists) |
| **Digital shadow** | initialized from a *specific* program's real genotypes + fitted architecture; virtual mirrors real (one-way) | needs estimation-from-data |
| **Digital twin** | kept in sync each cycle with observed data; recommendations feed back to the real program | needs assimilation + validation + loop |

**What it takes to reach the destination (ordered):**
1. **Estimated, not specified, architecture.** Fit marker effects, variance
   components, `h2`, and genetic correlations from the program's own historical
   phenotypes (GBLUP/GWAS). This is the roadmap's PopVar-style / GS-predictor item
   (§8a modern-methods NEXT) — the first concrete step, and it turns the simulator
   into a **digital shadow**.
2. **Data assimilation loop (the defining feature).** Each season, ingest what was
   actually genotyped / planted / measured and update the twin's state, effects, and
   accuracy so it tracks reality. The biggest missing piece; separates a twin from a
   simulation.
3. **Environment / G×E.** Explicit multi-environment / weather-covariate model so
   simulated phenotypes match realized field means (today: `h2`-based residuals only).
4. **Calibration & validation against reality.** Backtest against the program's
   *observed* gain over past cycles; report uncertainty, not point estimates —
   upgrading the designer's *projected* diagnostics to *validated* forecasts.
5. **Decision feedback loop.** Optimize interventions on the twin (OCS / usefulness /
   designer), apply to the real program, recalibrate from the realized outcome.
6. **Platform engineering.** Each twin = a persistent, versioned object tied to a
   program id, with a datastore, update API, provenance, and live connectors (LIMS,
   field data, genotyping); the designer's pipeline tabs + `plumber` runner are the
   seed.

**Positioning rule:** call the current product a **breeding-program simulator /
in-silico designer**. "Digital twin of a breeding population" becomes accurate only
after steps 1–5. The genetics substrate (real genotypes, meiosis, variance
partition, selection, OCS/usefulness) is already strong; the missing work is
**estimation-from-data + synchronization**, not the genetics. (Mirrored in the
designer manuscript Discussion.)

### 8b. Python package

- **Registry:** PyPI is the direct CRAN analogue, but has **no gatekept review**
  (anyone can publish). For this audience, **conda-forge / Bioconda** is the
  channel that actually drives adoption. Target both.
- **What is shared vs. rewritten:** the deterministic Rust core (numericalize,
  meiosis) is shared via PyO3/maturin; the **stochastic grammar layer is R and
  must be reimplemented** in Python. Every grammar change is then a 2× edit.
- **Timing:** after the grammar/spec is stable (ideally after the first CRAN
  release), so the port is not chasing a moving target. Maintainer has placed it
  **after** the designer.
