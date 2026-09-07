# TODO.md — simplePHENOTYPES v2 Development Checklist

> Track progress here. Update status as you go.
> See docs/NEXT_STEPS.md for rationale; docs/SPEC.md for grammar contracts.
> Start every Claude Code session: `@CLAUDE.md @docs/ARCHITECTURE.md @docs/SPEC.md`

---

## Block 1 — One-time setup
*Do once. ~30 min total. None of this touches the package code.*

- [ ] Save all docs from planning session to project:
  - [ ] `CLAUDE.md` → project root (gitignored)
  - [ ] `.gitignore` → project root
  - [ ] `.Rbuildignore` → project root
  - [ ] `docs/ARCHITECTURE.md`
  - [ ] `docs/DECISIONS.md`
  - [ ] `docs/SPEC.md`
  - [ ] `docs/NEXT_STEPS.md`
  - [ ] `docs/PROJECT_CONTEXT.md`

- [ ] Disable Claude AI attribution globally:
  ```bash
  mkdir -p ~/.claude
  # paste settings.json → ~/.claude/settings.json
  ```

- [ ] Install VS Code extensions:
  `rust-analyzer`, `R` (Posit), `Claude Code`, `Error Lens`, `Even Better TOML`, `CodeLLDB`

- [ ] Add MCP servers (terminal, one-time globally):
  ```bash
  claude mcp add github -- npx -y @modelcontextprotocol/server-github
  claude mcp add filesystem -- npx -y @modelcontextprotocol/server-filesystem ~/projects
  claude mcp add context7 -- npx -y @upstash/context7-mcp
  claude mcp add memory -- npx -y @modelcontextprotocol/server-memory
  claude mcp add git -- npx -y @modelcontextprotocol/server-git
  claude mcp add fetch -- npx -y @modelcontextprotocol/server-fetch
  claude mcp list    # verify
  ```

- [ ] Export GitHub token (add to `~/.zshrc` or `~/.bashrc`):
  ```bash
  export GITHUB_TOKEN=ghp_your_existing_token
  ```

- [ ] Install spec workflow slash commands (project root, one-time):
  ```bash
  npx @pimzino/claude-code-spec-workflow
  ```

---

## Block 2 — Before writing a single line of new code

- [ ] **Capture v1.3.0 reference outputs** ← most important step
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

- [ ] Run the capture script and verify RDS files exist:
  ```r
  list.files("inst/extdata/v1_3_0_reference/")
  ```

- [ ] Commit captured references:
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
- [ ] Add testthat to the package (if `tests/testthat/` does not yet exist):
  ```r
  usethis::use_testthat()
  ```
- [ ] Verify the test harness runs (zero tests is fine):
  ```r
  devtools::test()
  ```

### Step 2 — Write the first test: `as_numeric()`
*Deterministic input → deterministic output. This is TDD: write the test first,*
*then the Rust port must pass it. No ambiguity about what "correct" means.*

- [ ] Create `tests/testthat/test-as-numeric.R`.
  Include at minimum:
  - A small, hand-checkable example (e.g., 5 markers × 4 individuals; hand-verify
    the -1/0/1 encoding before writing the assertion).
  - The same conversion run through v1.3.0 (or the stored RDS) to assert parity.
  - NA / missing-marker handling.
  - A larger smoke test on `SNP55K_maize282_maf04` asserting output dimensions and
    value range (all values in {-1, 0, 1}).
  - Claude Code prompt:
    > `@CLAUDE.md @R/as_numeric.R`
    > "Write tests/testthat/test-as-numeric.R. Include: small hand-verifiable example
    > with expected -1/0/1 matrix written explicitly; NA handling; parity assertion
    > against inst/extdata/v1_3_0_reference/ for the same input. Plan first."

- [ ] Make sure the test passes against the current R implementation:
  ```r
  testthat::test_file("tests/testthat/test-as-numeric.R")
  ```

- [ ] Commit:
  ```bash
  git commit -m "test: add as_numeric() unit tests — first test in package"
  ```

### Step 3 — Write the v1.3.0 parity test harness
*Now that the reference RDS files exist and testthat is running, wire them up.*

- [ ] Create `tests/testthat/test-v130-parity.R`:
  - For each RDS in `inst/extdata/v1_3_0_reference/`:
    load, configure the equivalent new-grammar call, assert identical QTNs and
    identical phenotype values.
  - These will **fail** until the grammar is implemented (Block 3) — that is expected.
    Mark them with `testthat::skip("grammar not yet implemented")` initially so
    `devtools::test()` stays green.
  - Remove the skips one by one as each grammar function is completed.

- [ ] Commit:
  ```bash
  git commit -m "test: add v1.3.0 parity harness (skipped until grammar implemented)"
  ```

### Step 4 — Port `as_numeric()` to Rust
*The test already exists. Port, then remove the R fallback and confirm the test still*
*passes. This proves the pattern for every future surgical Rust port.*

- [ ] Set up the rextendr scaffold (if `src/rust/` does not yet exist):
  ```r
  rextendr::use_extendr()
  ```
- [ ] Implement `as_numeric()` in `src/rust/src/numeric.rs`.
  Claude Code prompt (new session, plan mode):
  > `@CLAUDE.md @docs/ARCHITECTURE.md @R/as_numeric.R
  >  @tests/testthat/test-as-numeric.R`
  > "Port R/as_numeric.R to src/rust/src/numeric.rs. The Rust output must pass the
  > existing testthat tests exactly. R passes the genotype matrix; Rust returns -1/0/1
  > integers. No RNG — fully deterministic. Plan first."

- [ ] Confirm Rust passes the existing test without modification:
  ```r
  testthat::test_file("tests/testthat/test-as-numeric.R")
  ```
- [ ] Run full check:
  ```r
  devtools::test()
  rcmdcheck::rcmdcheck()
  ```
- [ ] Commit:
  ```bash
  git commit -m "port: as_numeric() to Rust — tests pass"
  ```

### Step 5 — Add a few more deterministic tests before grammar work
*Build a minimal safety net. These don't need to be exhaustive — just enough to*
*catch regressions during the restructure.*

- [ ] `test-format-conversion.R` — test HapMap / VCF input parsing (deterministic).
- [ ] `test-shim.R` — smoke-test that `create_phenotypes()` still runs without error
  on the SNP55K dataset (not checking values yet — just that it doesn't crash).
- [ ] Commit:
  ```bash
  git commit -m "test: add format-conversion and shim smoke tests"
  ```

---

## Block 3 — Grammar implementation
*Only start here after Block 2B is complete and all existing tests are green.*

- [ ] **Function inventory** — Claude Code, new session:
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
  - [x] PleioArch pleiotropy engine (`cor`, exact for 2 traits; Cholesky fallback >2).

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

## Block 4 — CRAN submission

- [ ] All parity tests passing (no `skip()`s remaining in `test-v130-parity.R`).
- [ ] `devtools::check_win_devel()` passes.
- [ ] `rextendr::vendor_pkgs()` run; bundle < 5MB.
- [ ] `cran-comments.md` written explaining v2 changes.
- [ ] `devtools::release()`.

---

## Quick Reference — starting a Claude Code session

```
@CLAUDE.md @docs/ARCHITECTURE.md @docs/DECISIONS.md @docs/SPEC.md @docs/NEXT_STEPS.md
```
One conversation per feature. Plan mode for all new work.
