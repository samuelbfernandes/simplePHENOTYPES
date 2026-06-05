# CLAUDE.md — R Package → Rust Migration Project

## Project Overview
Updating and migrating an R package to Rust for performance-critical functions, with bindings
for both R (CRAN) and Python. The R source is the behavioral ground truth.

## Restructure goals (complete before porting)

The `restructure` branch refactors the R source before any Rust work begins:

- **Split `create_phenotypes.R`** into focused, single-responsibility functions. The original
  function grew too large; the goal is one function per logical step (input handling, QTN
  selection, genetic effects, phenotype assembly, output).
- **Expose internal utilities as proper functions.** Example: numericalization logic was buried
  inside `create_phenotypes()` and is being extracted into `as_numeric()` (not exported to
  users) so it can be tested and later ported independently.
- **New internal files created during restructure**:
  - `R/as_numeric.R` — numeric coercion (internal)
  - `R/check_in.R` — pre-flight parameter validation (internal)
  - `R/format_conversion.R` — format conversion utilities (internal)
  - `R/handling_input_formats.R` — genotype input dispatch (internal)
  - `R/table_to_numeric.R` — character table → numeric matrix (internal)
  - `R/vQTL.R` — variance QTL simulation (exported)
- The R API seen by users (`create_phenotypes()` signature) must remain backward-compatible.

## Architecture
my-package/
├── CLAUDE.md              ← you are here (auto-loaded every session)
├── .claudeignore
├── Cargo.toml             ← Cargo workspace root
├── docs/
│   ├── SPEC.md            ← what to build (living document)
│   ├── BUGS.md            ← bug triage log
│   └── DECISIONS.md       ← architectural decisions + rationale
├── core/                  ← pure Rust logic (no R or Python deps)
├── r-pkg/                 ← extendr/rextendr R bindings
│   └── R/                 ← original R source (READ ONLY during port)
└── py-pkg/                ← PyO3 + maturin Python bindings

## Language & Toolchain
- **Core logic**: Rust (stable, no nightly features)
- **R bindings**: extendr + rextendr (CRAN-compatible)
- **Python bindings**: PyO3 + maturin
- **Rust edition**: 2021
- **R minimum version**: 4.1

## Build & Test Commands

### Rust core
```bash
cargo build -p core
cargo test -p core
cargo clippy -- -D warnings
cargo fmt
cargo doc --open
```

### R package
```bash
Rscript -e "devtools::document()"
Rscript -e "devtools::test()"
Rscript -e "rcmdcheck::rcmdcheck()"
Rscript -e "rextendr::vendor_pkgs()"   # run before every CRAN submission
```

### Python package
```bash
maturin develop --release
pytest
maturin build --release
```

## Workflow: Bugs First, Port Second

### The non-negotiable sequence

```
Phase 0 → Bootstrap workspace (CLAUDE.md, Cargo workspace, scaffolds)
Phase 1 → Fix all bugs in R source  ← DO THIS BEFORE ANY PORTING
Phase 2 → Write SPEC.md             ← plan before code, always
Phase 3 → Port function by function ← one new conversation per function
Phase 4 → New features              ← Rust-first, never R-first
```

**Always follow this order — never port broken code:**
1. Fix bugs in R source (R/ is ground truth)
2. Add/complete R tests (pin correct behavior)
3. Port to Rust core (tests prove Rust matches corrected R)
4. Add new features directly in Rust core
### Bug fix loop (Phase 1)
```
/bug-create → /bug-analyze → /bug-fix → /bug-verify
→ Rscript -e "devtools::test()"
→ commit: "fix: <description> [R source]"
```

### Bug fix loop (per bug)
/bug-create → /bug-analyze → /bug-fix → /bug-verify → Propose fix in plan mode →
Review diff → Accept → Run devtools::test() → Commit
Commit format: `fix: <description> [R source]`

### Port loop (per function)
Plan only (no code) → Review plan → Implement →
cargo clippy → cargo test → devtools::test() → Commit
Commit format: `port: <function_name>() to Rust core + R/Python bindings`

## Spec-Driven Development Phases

| Step | Document | Goal/Contents |
|------|----------|----------|
| 0 | CLAUDE.md | Bootstrap workspace |
| 1 | `docs/SPEC.md` | What to build: public API, R-specific behaviors (NA handling, recycling, vectorization), architecture per layer, CRAN constraints, open questions |
| 2 | Plan review | Claude proposes implementation plan in plan mode; human approves |
| 3 | Code | Implementation follows the approved plan |


**Always open a new conversation per function during Phase 3.**
**Always use plan mode before touching any file.**

## Rust Conventions
- Error handling: `thiserror` in `core/`, `anyhow` in binaries
- No `unwrap()` in library code — use `?` propagation
- No `unsafe` without a `// SAFETY: <reason>` comment
- Prefer `&str` over `String` in function arguments
- Document lifetime annotations with an inline comment
- No nightly-only features (CRAN builds on stable Rust)

## Binding Conventions
- R bindings: `#[extendr]` macro + `extendr_module!`
- Python bindings: `#[pyfunction]` / `#[pyclass]` macros
- **Both binding layers import from `core::` only — zero business logic
  in r-pkg/ or py-pkg/**
- Match R function signatures exactly in the public API
- NA handling must be explicit in both binding layers

## Testing Rules
- Every ported function needs:
  - A unit test in `core/` 
  - An integration test comparing Rust output vs. R reference output
- Run both `cargo test` and `devtools::test()` before every commit
- Never commit a port without passing tests in both languages

## CRAN Constraints
- Rust dependencies must be vendored: `rextendr::vendor_pkgs()`
- Package must compile offline
- DESCRIPTION `SystemRequirements` must list: `Cargo (Rust's package manager), rustc`
- Built package must stay under 5MB — note any exceptions in cran-comments.md
- No nightly Rust features

## VS Code Context Tips
- Reference files with `@filename` (fuzzy match) or type `@` for picker
- Reference specific lines: select code → Alt+K (Win) / Option+K (Mac)
- Reference terminal output: `@terminal:name`
- Use **plan mode** for all new work; **default mode** for bug fixes
- Open a **new conversation** for each function port (keeps context clean)
- Auto-loaded files: CLAUDE.md (this file), any file in docs/ you @-mention

## MCP Servers (CLI only — not available in VS Code extension)
```bash
claude mcp add github       # browse issues, R package history
claude mcp add filesystem   # explicit cross-project file access
```

**Self-review prompt for SPEC.md:**
> "Review this spec from three perspectives: (1) a Rust developer — flag anything underspecified for the core crate; (2) a CRAN maintainer — flag anything that violates CRAN policies; (3) a Python user — flag API decisions that would feel unidiomatic in Python."

## Avoid
- Don't modify R/ source files during porting (it's the reference)
- Don't put logic in binding layers — core/ only
- Don't use unwrap() in library code
- Don't use nightly Rust features
- Don't skip the plan review step
- Don't port a function that still has a known bug
- Don't add large crate dependencies (watch the 5MB CRAN vendor limit)
- Don't run /bug-fix without first running /bug-analyze

## Claude Code Configuration

### VS Code extension usage
- Use **plan mode** for all new work and ports (Claude proposes, you approve)
- Use **default mode** for bug fixes (review each diff before accepting)
- Use **auto-accept** only for formatting/linting passes
- Open a **new conversation** for each function port — keeps context window small