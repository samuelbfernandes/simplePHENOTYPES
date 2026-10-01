## PRIORITY — requests from breedingDesigner SPEC-0020 (2026-09-30)

The maintainer's decision (SPEC-0020 R19): the full-size Bančič cross-engine acceptance run
in breedingDesigner waits on item 1. Evidence: `breedingDesigner/docs/validation/
bancic-timing-2026-09-30.md` (one full-size DH cycle: 10,000 DH × 14,000 markers).

- [x] **1. `simulate_phenotype()` deparses the whole population to name it.**
      `R/grammar_simulate_phenotype.R:195`: `geno_name <- deparse(substitute(geno))`. When
      the caller passes the object inline (`do.call(simulate_phenotype, list(geno = pop, …))`,
      as breedingDesigner's `.build_model()` does), `deparse()` renders the entire
      10,000 × 14,000 `Population`: 47 s of a 63 s cycle (Rprof, 72 % of the simplePHENOTYPES
      arm; 60 % of the AlphaSimR arm, which runs the same grammar). Fix:
      `deparse(substitute(geno), nlines = 1)` (or `deparse1(..., collapse = "")` with a width
      cap), plus a test that a large inline `geno` costs milliseconds. Blocking.
      *Done 2026-09-30 (branch `feat/spec0020-engine-requests`): `.geno_label()` bounds the name (symbols and short calls unchanged; a 2000 x 14000 inline Population: 9.2 s -> 0 s); `as_numeric()` uses it too.*
- [x] **2. `double_haploid()` / `cross()` / `selfcross()` per-call cost.** 100 calls of
      `double_haploid(parent, n = 100)` on 14,000 markers take 12.9 s per cycle (AlphaSimR
      `makeDH` 0.3 s): each call serialises both parental haplotypes to bit strings and
      parses the progeny strings back (`cross_mating.R`, `mate_haplotypes_core`). Consider a
      vectorised multi-parent entry point (all parents in one Rust call) or raw-vector I/O.
      *Done 2026-09-30 (DECISION-040): `mate()` runs all plan rows in one Rust call (`mate_many_core`, integer I/O); 100 x `double_haploid(n = 100)` on 14,000 markers 14.1 s -> 2.3 s (6.1x), `mate()` DH plan 6.9x; draw-for-draw identical to the sequential path. `.stable_key()` in `R/cross_pedigree.R` could still be vectorised (~2 ms of ~23 ms per call).*
- [x] 3. Crossover interference model (gamma / count–location with interference; AlphaSimR
      `v`, `p`) as an option of the meiosis core — today Poisson only. Needed for
      like-for-like comparisons with AlphaSimR's default.
      *Done 2026-09-30 (DECISION-041): `interference = NULL` on `cross`/`selfcross`/`double_haploid`/`mate`/`crossbreed`; `list(nu, p)` = two-pathway gamma model, expected chiasma count per Morgan unchanged. Codex round 5: theory passed; `nu` limited to [1, 1e6] (2026-09-30); forwarded by every function that runs meiosis (see Block 3C).*
- [ ] 4. Additive-by-environment (G×E) trait layer (AlphaSimR `addTraitAG` semantics) — needed
      to reproduce Bančič Program 4.
      *Not implemented: waiting for the maintainer's instructions.*
- [x] 5. Public phased-haplotype constructor (`as_population(haplotypes = …)` or
      `population_from_haplotypes()`) — breedingDesigner currently builds the AlphaSimR
      view by replacing `cis`/`trans` after `as_population()` (documented shim).
      *Done 2026-09-30 (DECISION-039): `population_from_haplotypes()` / `haplotypes()` (markers x individuals; 1 = counted allele, dosage = cis + trans - 1; round trip identical).*
- [x] 6. Entry-mean replication in the grammar (`reps`: residual variance `var_e/reps`,
      AlphaSimR `setPheno(varE, reps)` semantics).
      *Done 2026-09-30 (DECISION-038): `reps` in `simulate_phenotype()` / `complex_phenotypes()`; residual variance `var_e/reps`, `h2` stays single-record; awaiting Codex theory review.*
- [x] 7. `sample_parents()`: keep the source id / slot of each drawn parent as an attribute
      (breedingDesigner recovers it by key today).
      *Done 2026-09-30: `attr(x, "source")` (slot, id, index, name) on `sample_parents()` output; draws unchanged.*
- [x] 8. `select_ind(method = "within_family")`: a per-family count (`n_per_family`) for
      unequal families (today `n` is a total apportioned proportionally).
      *Done 2026-09-30: `select_ind(method = "within_family", n_per_family = )` (one count or a vector named by family; undersized family is an error).*
- [ ] 9. Coalescent founders (MaCS-like) natively, so a run does not need AlphaSimR for
      historical-LD founders.
      *Not implemented: waiting for the maintainer's instructions.*

<!-- AUDIT:BEGIN -->
## Audit findings — 2026-09-17 (from dev/audit-all.sh, run on v2.0.0)

Transcripts: `dev/.audit/transcripts/audit-*-20260917-141443.md`. All 279 test
contexts pass after the fixes below (`devtools::test()` green).

**Fixed (verified with executed R):**
- [x] grammar V2 — `complex_phenotypes()` used `.genetic_matrix` not `.genetic_value_matrix`, dropping transcriptome-mediated genetic value (false "no genetic variation" error).
- [x] grammar V4 — dominance/`d`-epistasis on a *partially* hetless set: error only when the whole set is dead, else warn (was silently inert).
- [x] grammar V1 / Q2 "keep + clarify docs" — transcriptome `prop` is intentionally outside the marker h2 budget; `h2` param docs clarified.
- [x] effects-arch O1 — single shared pleiotropic QTN can't realize an intermediate `cor` (always ±1): now warns.
- [x] effects-arch O3 — `cor = 0` against a zero-variance trait is undefined (NA), not 0: now warns.
- [x] crossing O1 — `.as_founder_pop()` ignored `ind_idx`, so schemes reintroduced excluded individuals: now subsets.
- [x] crossing O2 — `recurrent_selection()` doc said "mating chain" but code samples independent random pairs (DECISION-015): doc corrected.
- [x] rust-core O1 — `as_numeric()` treated an in-memory character *matrix* as a file path: now distinguishes by `dim`.
- [x] io-formats O2 — multiallelic calls silently collapsed to biallelic (table + VCF): now set to NA with a message/warning. (Also fixed a latent `.het_to_letters()` crash on homozygote letters.)
- [x] io-formats O3 — empty VCF metadata column classified as a sample (`all(grepl(...,character(0)))` vacuous TRUE): now requires a real GT value.
- [x] io-formats O4 — variant-window LD pruning was a no-op on all-NA-chr data (`which(chr==NA)`): now groups NA-chr markers.
- [x] io-formats V2 — transcriptome `var_explained` was nominal component variance, not a fraction of realized V_P like marker rows: now divides by `var_p`.

**Not bugs (UNVERIFIABLE citation-page checks; implementations verified vs source/PLINK):**
grammar C3, effects-arch C3, crossing C3, io-formats O5.

**Design decision (Q1 = "extend control to non-additive") — DONE 2026-09-25:**
- [x] **Extend genetic-correlation control (`cor` / pleiotropy / LD architecture) to dominance and epistasis layers** — DECISION-023, `docs/SPEC-nonadditive-correlation.md`. Per-component PleioArch covariance; every component targets `cor` (realized correlation converges as units and individuals grow, given approximate linkage equilibrium), and the *total* targets `cor` for proportional (e.g. scalar) per-trait `prop`, else a warning gives its attenuated large-sample target (18/18 acceptance cells; additive/independent/ld bit-identical; independent review findings over five rounds fixed). Closes grammar P1/P4/X1 and effects-arch O2/X1. Tests: `test-nonadditive-cor.R` (fail on old code, pass now).
<!-- AUDIT:END -->


# TODO.md — simplePHENOTYPES v2 Development Checklist

> Track progress here. Update status as you go.
> See docs/NEXT_STEPS.md for rationale; docs/SPEC.md for grammar contracts.
> Start every Claude Code session: `@CLAUDE.md @docs/ARCHITECTURE.md @docs/SPEC.md`

---

## Block 1 — One-time setup
*Do once. ~30 min total. None of this touches the package code.*
*Status checked 2026-09-27 against the repo; the maintainer marked the machine-only items
(VS Code extensions, `GITHUB_TOKEN`) done on 2026-09-28.*

- [x] Save all docs from planning session to project:
  - [x] `CLAUDE.md` → project root (gitignored)
  - [x] `.gitignore` → project root
  - [x] `.Rbuildignore` → project root
  - [x] `docs/ARCHITECTURE.md`
  - [x] `docs/DECISIONS.md`
  - [x] `docs/SPEC.md`
  - [x] `docs/NEXT_STEPS.md`
  - [x] `docs/PROJECT_CONTEXT.md` *(not in the repo)*

- [x] Disable Claude AI attribution globally *(also enforced by `.githooks/commit-msg`
  and the `no-ai-attribution` CI check)*:
  ```bash
  mkdir -p ~/.claude
  # paste settings.json → ~/.claude/settings.json
  ```

- [x] Install VS Code extensions:
  `rust-analyzer`, `R` (Posit), `Claude Code`, `Error Lens`, `Even Better TOML`, `CodeLLDB`

- [x] Add MCP servers (terminal, one-time globally):
  ```bash
  claude mcp add github -- npx -y @modelcontextprotocol/server-github
  claude mcp add filesystem -- npx -y @modelcontextprotocol/server-filesystem ~/projects
  claude mcp add context7 -- npx -y @upstash/context7-mcp
  claude mcp add memory -- npx -y @modelcontextprotocol/server-memory
  claude mcp add git -- npx -y @modelcontextprotocol/server-git
  claude mcp add fetch -- npx -y @modelcontextprotocol/server-fetch
  claude mcp list    # verify
  ```

- [x] Export GitHub token (add to `~/.zshrc` or `~/.bashrc`):
  ```bash
  export GITHUB_TOKEN=ghp_your_existing_token
  ```

- [x] Install spec workflow slash commands (project root, one-time):
  ```bash
  npx @pimzino/claude-code-spec-workflow
  ```

---

## Block 2 — Before writing a single line of new code

- [x] **Capture v1.3.0 reference outputs** ← most important step
  - In a clean R session with simplePHENOTYPES 1.3.0 installed, run each README
    vignette example with `big_add_QTN_effect` removed.
  - Save the selected QTN names AND phenotype values as RDS files.
  - Store under `inst/extdata/v1_3_0_reference/`.
  - These files gate every future release; they must exist before any new code.
  - Claude Code prompt (new session, plan mode):
    > `@CLAUDE.md @docs/SPEC.md`
    > "Write `inst/extdata/v1_3_0_reference/capture_references.R` that installs
    > simplePHENOTYPES 1.3.0, runs each README example with big_add_QTN_effect removed,
    > and saves QTN names and phenotype matrix as RDS per example. Use
    > SNP55K_maize282_maf04. Plan first."

- [x] Run the capture script and verify RDS files exist:
  ```r
  list.files("inst/extdata/v1_3_0_reference/")
  ```

- [x] Commit captured references:
  ```bash
  git add inst/extdata/v1_3_0_reference/
  git commit -m "test: add v1.3.0 frozen reference outputs (big-QTN removed)"
  ```

---

## Block 2B — Build the testing foundation

*The package has no unit tests. Start here before touching grammar or Rust.*
*`as_numeric()` is the right first test: deterministic, already partially written in R,*
*input and output are both known, and it will be the first Rust port.*

### Step 1 — Set up testthat infrastructure
- [x] Add testthat to the package (if `tests/testthat/` does not yet exist):
  ```r
  usethis::use_testthat()
  ```
- [x] Verify the test harness runs (zero tests is fine):
  ```r
  devtools::test()
  ```

### Step 2 — Write the first test: `as_numeric()`
*Deterministic input → deterministic output. This is TDD: write the test first,*
*then the Rust port must pass it. No ambiguity about what "correct" means.*

- [x] Create `tests/testthat/test-as-numeric.R`.
  Include at minimum:
  - A small, hand-checkable example (e.g., 5 markers × 4 individuals; hand-verify
    the -1/0/1 encoding before writing the assertion).
  - The same conversion run through v1.3.0 (or the stored RDS) to assert parity.
  - NA / missing-marker handling.
  - A larger smoke test on `SNP55K_maize282_maf04` asserting output dimensions and
    value range (all values in {-1, 0, 1}).
  - Claude Code prompt:
    > `@CLAUDE.md @R/io_as_numeric.R`
    > "Write tests/testthat/test-as-numeric.R. Include: small hand-verifiable example
    > with expected -1/0/1 matrix written explicitly; NA handling; parity assertion
    > against inst/extdata/v1_3_0_reference/ for the same input. Plan first."

- [x] Make sure the test passes against the current R implementation:
  ```r
  testthat::test_file("tests/testthat/test-as-numeric.R")
  ```

- [x] Commit:
  ```bash
  git commit -m "test: add as_numeric() unit tests — first test in package"
  ```

### Step 3 — Write the v1.3.0 parity test harness
*Now that the reference RDS files exist and testthat is running, wire them up.*

- [x] Create `tests/testthat/test-v130-parity.R`:
  - For each RDS in `inst/extdata/v1_3_0_reference/`:
    load, configure the equivalent new-grammar call, assert identical QTNs and
    identical phenotype values.
  - These will **fail** until the grammar is implemented (Block 3) — that is expected.
    Mark them with `testthat::skip("grammar not yet implemented")` initially so
    `devtools::test()` stays green.
  - Remove the skips one by one as each grammar function is completed.

- [x] Commit:
  ```bash
  git commit -m "test: add v1.3.0 parity harness (skipped until grammar implemented)"
  ```

### Step 4 — Port `as_numeric()` to Rust
*The test already exists. Port, then remove the R fallback and confirm the test still*
*passes. This proves the pattern for every future surgical Rust port.*

- [x] Set up the rextendr scaffold (if `src/rust/` does not yet exist):
  ```r
  rextendr::use_extendr()
  ```
- [x] Implement `as_numeric()` in `src/rust/src/numeric.rs`.
  Claude Code prompt (new session, plan mode):
  > `@CLAUDE.md @docs/ARCHITECTURE.md @R/io_as_numeric.R
  >  @tests/testthat/test-as-numeric.R`
  > "Port R/io_as_numeric.R to src/rust/src/numeric.rs. The Rust output must pass the
  > existing testthat tests exactly. R passes the genotype matrix; Rust returns -1/0/1
  > integers. No RNG — fully deterministic. Plan first."

- [x] Confirm Rust passes the existing test without modification:
  ```r
  testthat::test_file("tests/testthat/test-as-numeric.R")
  ```
- [x] Run full check:
  ```r
  devtools::test()
  rcmdcheck::rcmdcheck()
  ```
- [x] Commit:
  ```bash
  git commit -m "port: as_numeric() to Rust — tests pass"
  ```

### Step 5 — Add a few more deterministic tests before grammar work
*Build a minimal safety net. These don't need to be exhaustive — just enough to*
*catch regressions during the restructure.*

- [x] `test-format-conversion.R` — test HapMap / VCF input parsing (deterministic).
- [x] `test-shim.R` — smoke-test that `create_phenotypes()` still runs without error
  on the SNP55K dataset (not checking values yet — just that it doesn't crash).
- [x] Commit:
  ```bash
  git commit -m "test: add format-conversion and shim smoke tests"
  ```

---

## Block 3 — Grammar implementation
*Only start here after Block 2B is complete and all existing tests are green.*

- [x] **Function inventory** — done 2026-09-28: 274 functions (dead `.causal_loci()`
  removed), table in ARCHITECTURE.md §5.1. Original prompt:
  > `@CLAUDE.md @docs/ARCHITECTURE.md @R/`
  > "Audit every function in R/. For each: name | exported? | has_tests? | new module
  > | Rust candidate | parity-critical. Write table to ARCHITECTURE.md §5. No code."

- [x] Restructure `R/` into modules — implemented as **filename prefixes** in a flat `R/`
  (`grammar_*`, `arch_*`, `effects_*`, `io_*`), because R does not source `R/` subdirs
  (see ARCHITECTURE.md §4).

- [x] Implement grammar functions in R:
  - [x] `simulate_phenotype()` foundation (+ one-call detection folded in)
  - [x] `additive()` layer
  - [x] `dominance()` layer (same_as_add default)
  - [x] `epistasis()` layer (2-way default)
  - [x] `vqtl()` layer (same_as_add default)
  - [x] `complex_phenotypes()`
  - [x] ~~`sim_phenotypes()`~~ dropped — one-call folded into `simulate_phenotype()`
        (DECISION; create_phenotypes() also remains for v1-style one-call).
  - [x] PleioArch pleiotropy engine (`cor` targeted for any number of traits; the
        Cholesky fallback was removed by DECISION-013; see DECISION-023 for what is and
        is not exact).

- [x] Mark `create_phenotypes()` as frozen legacy (superseded `@description` note; no
  lifecycle dependency to avoid an unused-import NOTE). Does NOT delegate (DECISION-008).
- [x] `test-v130-parity.R` repointed at `create_phenotypes()` (regression guard); added
  `test-grammar.R` statistical/structural tests (variance identity, seed invariance,
  realized `cor`, QTN structure, complex weighting) per DECISION-009.

- [x] Capture isqg reference outputs into `inst/extdata/isqg_v1_outputs/` (6 fixtures,
  each storing the drawn randomness alongside isqg's output).
- [x] Port isqg meiosis / cross / DH to Rust (`src/rust/src/meiosis.rs`, `genome.rs`),
  gated by **exact** bit-parity (DECISION-012) in `tests/testthat/test-isqg-parity.R`.
  R draws all randomness; the Rust core never calls an RNG.

- [x] Genetic map: `synthetic_map()` (exported) plus a centromere-suppressed synthetic
  map written into `SNP55K_maize282_maf04$cm`, which was previously logical all-NA.
  Regenerated by `data-raw/make_SNP55K_map.R`.
- [x] Write multi-generation functions: `cross()`, `selfcross()`, `double_haploid()`,
  plus `as_population()` / `dosages()` / `n_individuals()` / `[` / `print` and the
  `.normalize_geno()` hook so `simulate_phenotype()` accepts a `Population` (SPEC §4.1).
- [x] Vignette demonstrating every user-facing function
  (`vignettes/simplePHENOTYPES-v2.Rmd`).

---

## Block 3B — Engine work requested by breedingDesigner

breedingDesigner (BD, `../breeding_designer`) runs every simulation through this package
and must not reimplement genetics, so each item below is blocked on the engine. Refs point
into the BD repo. Methods/theory items need an independent theory review
(`docs/THEORY_REVIEW.md`) before they ship.

**Now — blocking current BD users**
- [x] **Publish the `restructure` engine where users install from.** *Done: merged to
  `master` (all 16 exports below present) and tagged `v2.0.0`; the `restructure`
  branch is deleted. BD-side follow-up: switch `Remotes:` to
  `samuelbfernandes/simplePHENOTYPES@v2.0.0` (or `master`) and pin `(>= 2.0.0)`.* A plain
  `remotes::install_github("samuelbfernandes/simplePHENOTYPES")` (or the released
  version) lacks operators BD dispatches, and runs fail with `'select_ind' is not an
  exported object`. BD needs all of: `as_population`, `cross`, `selfcross`,
  `double_haploid`, `single_seed_descent`, `bulk`, `simulate_phenotype`, `select_ind`,
  `optimum_contribution`, `cross_usefulness`, `recurrent_selection`, `sample_parents`,
  `pedigree`, `n_individuals`, `dosages`, `synthetic_map` (BD `R/run_app.R`
  `engine_required_exports()`; BD pins `simplePHENOTYPES (>= 1.4.0.9002)`,
  `Remotes: samuelbfernandes/simplePHENOTYPES@restructure`). Tag a version when merged.

**Next — methods BD has specced and is waiting on**
- [x] **Marker-assisted backcross: `mabc_select()` + `recurrent_parent_recovery()`**
  *Done (PR #8, `R/select_mabc.R`, `tests/testthat/test-mabc.R`; unselected recovery
  0.5000 / 0.7513 / 0.8719 / 0.9351 at F1–BC3; theory review PASS).*
  (full design in `docs/ROADMAP.md` → "Modern methods — NEXT"; BD
  `docs/reviews/backcross-theory-review.md`, `specs/SPEC-0012`). BD now ships the
  unselected structural backcross and a lineage-level recovery test (BD
  `tests/testthat/test-backcross-recovery.R`: synthetic informative-marker fixture,
  60 lineages, predeclared tolerance, selfing negative control) that can be reused
  to validate the unselected baseline here.
- [x] **Fixed-scale TOTAL genotypic value (additive + dominance) accessor** — *already
  provided by `genotypic_value(x, qtn, a, d)` (`G = A + D`, fixed scale; in the backend
  contract).* Sibling of
  `additive_value()`, with no per-population re-centring/re-scaling. Gates
  heterosis-capturing reciprocal recurrent selection (RRS-1b); until then BD runs an
  additive-only RRS and says so (BD `specs/SPEC-0011`, `R/rrs.R`). Not yet in
  `docs/ROADMAP.md`.
- [x] **Per-trait `effect` list in `additive()`.** *Done (PR #8): `additive(effect =
  list(...))`; re-scoring a frozen pleiotropic template reproduces its values.* Frozen multi-trait architectures with
  different per-trait effects currently need one trait-masked layer per trait in BD
  (BD `specs/SPEC-0002`); a per-trait effect list makes it a single layer.

**Specced in `docs/SPEC-block3b.md` — done (branch `feat/block3b-engine`; every batch
independently theory-reviewed, AGREE)**
Source: BD `docs/BREEDING_METHODS_CATALOG.md` ("Engine:" notes).
- [x] Pedigree in `Population` + mating plans: `parentage()`, `families()`,
  `mating_design()`, `mate()` (DECISION-024/025; F1/F2).
- [x] Combining-ability scorer (GCA / SCA, testcross merit):
  `combining_ability()` (expected + simulated), `template_effects()`,
  `phenotype_value(d =)` (DECISION-026).
- [x] Marker-based selection for MAS / gene pyramiding, and a marker index for MARS:
  `marker_select()`; `additive_value()` is the MARS index (DECISION-029).
- [x] BLUP / EBV prediction verb: `predict_ebv()` (GBLUP / pedigree, known variances),
  `a_matrix()`, `prediction_accuracy()`, `selection_methods()` (DECISION-030).
- [x] Multi-trait BLUP: `predict_ebv()` with an individuals x traits `pheno` matrix and
  `var_a` / `var_e` covariance matrices (DECISION-032).
- [ ] Single-step (pedigree + genomic, `H`) BLUP — still deferred by maintainer
  decision D14 (`docs/SPEC-block3b.md`).
- [x] Progeny-mean scorer (progeny testing): `progeny_test()` (DECISION-027); families
  in a design via `mating_design()` / `mate()` + `families()`.
- [x] Multi-trait sequential rules: `select_ind(method = "culling")` and tandem
  selection (a `trait` vector on `pedigree()` / `recurrent_selection()`, DECISION-028).
- [x] Multi-population mating for crossbreeding: `crossbreed()` (two-way / backcross /
  three-way / terminal / rotational), `breed_composition()`, `heterosis()`
  (DECISION-031); `mate()` executes BD SPEC-0007 plans.

**Still open**
- [ ] Polyploid (tetrasomic) model — BD's `Wheat_div` is allotetraploid coded diploid
  per subgenome as an approximation (already `docs/ROADMAP.md` → "Polyploids").
  Maintainer decisions D22/D23: deferred to v3; the concrete autotetraploid driver is
  **potato**. (Disomic allotetraploids such as wheat are already modelled correctly
  by the diploid-per-subgenome coding.)
- [ ] Python package — BD's canvas already generates Python for the planned
  `simplephenotypes` API (already `docs/ROADMAP.md` §8b).

---

## Block 3C — Post-audit follow-ups (independent dual-model audit, 2026-09/10)

Status 2026-10-01: audit and review rounds 1-4 are on master (NEWS "Audit fixes",
DECISION-033 to 037). The SPEC-0020 engine requests (items 1, 2, 3, 5, 6, 7, 8 above, the
interference propagation, and their review fixes) are in **PR #13**
(`feat/spec0020-engine-requests`, DECISION-038 to 041): full test suite green (63 files, 5761
expectations), local `R CMD check` 0 errors / 0 warnings, Codex re-reviews fixed.

**Housekeeping (do first)**
- [ ] Merge PR #13 once CI is green (the session's auto-fix watcher may not survive archiving;
  re-check `gh pr checks 13`). After merging, the main checkout's uncommitted `TODO.md` edits are
  already contained in this file (copied 2026-09-30): discard them (`git checkout TODO.md`)
  before `git pull`.
- [ ] Items **4** (G x E trait layer) and **9** (coalescent founders) in the PRIORITY list are
  NOT implemented: waiting for the maintainer's instructions.
- [ ] Evidence is gitignored and lives only in the audit worktree
  (`.claude/worktrees/kind-shamir-c7d1df/.tmp/`): copy `.tmp/audit-2026-09-29/` (reports,
  equation-to-code PDF) and `.tmp/codex-review*/` to a permanent folder before the worktree is
  deleted.
- [ ] Remove the leftover agent worktree/branch `worktree-agent-abf14fefa4d4acf18`.

**Reviews still owed (the other model must review genetics changes, AGENTS.md)**
- [ ] Codex re-review of the round-4 changes (script to write, adapt
  `.tmp/codex-review-round3.sh`): R4-1 (`.tune_lambda()` purely relative above-optimum band,
  `R/select_ocs.R`) and R4-5 (`.tx_mimic_scale()` always rescales to the requested per-gene
  variance, counted warning when ill-conditioned, `R/transcriptome_simulate.R`).
- [ ] Codex review of the last fixes of PR #13 (validation and wording only, 24 tests): the
  overwrite warning for default-named `as_numeric()` output files, `.check_counted()` rejecting
  matrices, the DECISION-038 index row.
- [ ] Codex review of round-3 items never sent to Codex: 11-column HapMap guard (R3-6),
  integer genotype schema (R3-7), heterosis retention text (R3-12), `h2_*` documentation
  (R3-14), case-insensitive orientation labels (R3-5).
- [ ] Add the interference rubric item (M4) to `docs/THEORY_REVIEW.md` after the review (no
  exact text was proposed yet); verify the cited McPeek & Speed (1995) and Housworth & Stahl
  (2003) references in `?cross` (author/year/journal only, unverified).
- [ ] Pre-PR gate: a full `R CMD check` **with vignettes** (needs pandoc). The `v1-to-v2`
  vignette failure escaped the test suite and local no-vignette checks; consider a test that
  evaluates every vignette's code.

**Open follow-ups from the new features**
- [ ] `.stable_key()` in `R/cross_pedigree.R` could be vectorised (~2 ms of ~23 ms per call).
- [ ] `cross_usefulness()` `"dh"` / `"selfcross"`: the interference dispersion is not tested
  (forwarding only, via a call counter).
- [ ] Optional: a scheme-level interference default (an option or an `as_population()` argument)
  instead of passing `interference =` to every function.
- [ ] The default-name overwrite warning of `as_numeric()` fires before conversion, so a failed
  conversion still warns.
- [ ] `reps` with a derived transcriptome layer is conditional on a fixed transcriptome
  covariate (DECISION-038); revisit if a per-record transcriptome environment is wanted.

**Known gaps left open on purpose**
- [ ] `counted_allele` lives on the R object only: it is not written to numeric text files or
  kept through row-subsetting, so those cases fall back to the `allele` label check
  (DECISION-036). Needs an output-contract change to persist it.
- [ ] testthat edition 3 is deferred (9 tests fail under it).
- [ ] V1 direct LD with dominance meets the LD contract for only a minority of seeds (3 of 29
  for `model = "D"`); the contract error is intentional, a better search is not implemented.
- [ ] Two fixed-table PLINK tests not written: `.plink_calc_lnlike`, `.plink_blocks_classify`.
- [ ] Not all 328 proposed tests of the audit were adopted (see
  `.tmp/audit-2026-09-29/reconciliation/*` section 6).
- [ ] Unverifiable citation pages (e.g. PRED-F3 Ceron-Rojas, AUX-F20 CRAN baseline version,
  Meuwissen 1997 / Baik et al. 2005 equation pages): confirm against the sources before the
  Python port quotes them.

**Before the Python re-creation**
- [ ] Use `.tmp/audit-2026-09-29/simplePHENOTYPES_equation_code_map.pdf` (equation to code map)
  and `V1_AUDIT_REPORT.md` / `V2_AUDIT_REPORT.md` as the porting checklist; regenerate the PDF
  line numbers after the fixes (they refer to the pre-fix tree).

---

## Block 4 — CRAN submission

- [x] All parity tests passing (no `skip()`s remaining in `test-v130-parity.R`).
- [ ] `devtools::check_win_devel()` passes.
- [x] `rextendr::vendor_pkgs()` run; bundle < 5MB *(`src/rust/vendor.tar.xz` 532 KB;
  package tarball 2.8 MB, 2026-09-29)*.
- [x] `cran-comments.md` written explaining v2 changes.
- [ ] `devtools::release()`.

---

## Quick Reference — starting a Claude Code session

```
@CLAUDE.md @docs/ARCHITECTURE.md @docs/DECISIONS.md @docs/SPEC.md @docs/NEXT_STEPS.md
```
One conversation per feature. Plan mode for all new work.
