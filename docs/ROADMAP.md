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
- [ ] **Fisher-orthogonal average effects (breeding value / dominance
      deviation).** The theoretically correct additive/dominance partition:
      derive the average effect alpha_j = a_j + d_j(1 - 2 p_j) per locus, build
      the additive (breeding-value) component from alpha and the dominance
      deviation as the orthogonal residual, so Cov(additive, dominance) = 0 and
      realized H2 = sum of components exactly at linkage equilibrium. This is a
      **re-parameterization, not a drop-in**: Vd becomes *derived* from
      (a, d, p) rather than a freely chosen `prop`, so it needs a new mode
      (e.g. a `genotypic_model()` builder, or `additive(orthogonal = TRUE)`
      taking per-locus a and d/degree, with Va/Vd/H2 emerging). This is also
      where a genuine **degree of dominance** would live -- the reason `degree`
      was removed from the current variance-partition grammar is that a scalar
      is washed out by per-component scaling; under the orthogonal genotypic
      model it becomes meaningful (d/a per locus). Roughly a day plus tests and
      docs; a real feature, but well understood.

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
- **Missing backend primitive — selection.** The maintainer's own example lists
  a "selection type" node, but there is **no selection / truncation /
  generation-advance function in the package today** (only the crossing
  functions). The designer therefore needs *new simulation code*, not just a UI.
  Build `select()` / generation-advance as plain, headless, tested R functions
  **before** the UI, and settle their API before freezing the grammar for CRAN.
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
