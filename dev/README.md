# dev/ — two-model development pipeline

Claude Code + Codex, used so an **independent** model catches implementation bugs and —
above all — **errors in the genetic theory** before they land.

## One-time setup (per clone)

```bash
./dev/setup.sh
```

Sets `core.hooksPath` to the committed `.githooks/`, verifies the AI-attribution guard is
active, and checks that `claude` and `codex` are on PATH. If `codex` is missing:
`npm i -g @openai/codex` (the VS Code extension alone is not enough — the pipeline drives
the CLI).

## Daily use

```bash
# Independent genetic-theory review of what you're about to commit:
dev/dual.sh review --staged
dev/dual.sh review R/select_ind.R R/select_ocs.R      # or specific files

# Let one model implement and the other review+critique, iterating to green:
dev/dual.sh loop "add rrBLUP GEBV as an `on` criterion, per docs/THEORY_REVIEW.md S5"

# Just the objective gate:
dev/dual.sh check
```

## Installed-package gate (`dev/test-installed.sh`)

Run before opening a PR. `devtools::test()` loads the sources, so it cannot see problems
that only appear in an installed package (R CMD check / CI): tests that read `R/`, `docs/`
or `benchmarks/`, undeclared optional packages, files written next to the sources.

    bash dev/test-installed.sh                 # build, install, run every test file
    FILTER='grammar|cross' bash dev/test-installed.sh   # only matching test files
    KEEP=1 bash dev/test-installed.sh          # keep the temporary directory

It builds the working tree with `R CMD build` (no vignettes), installs it with
`R CMD INSTALL --no-docs` into a temporary library (compiles the Rust core: several
minutes), copies `tests/` from the built tarball to a temporary directory, and runs
`testthat::test_dir(package = "simplePHENOTYPES", load_package = "installed")` with
`NOT_CRAN=true`. Everything lives under `mktemp -d` in `${TMPDIR:-/tmp}`. Exit status is
non-zero on any failed test or error (skips are fine; the source-reading tests and
`test-vignettes.R` skip there by design, so also run `devtools::test()`).

## Skeptical debate to consensus (`dev/debate.sh`)

Stronger than a one-shot review: a **SKEPTIC** must be *convinced* and a **DEFENDER** must
fix or rebut each objection **with an executed run** (rhetoric loses to a failing/passing
`Rscript`). It proceeds only when the skeptic explicitly agrees **and** `devtools::test()`
is green; on deadlock it **escalates to you** and never auto-proceeds or commits.

```bash
dev/debate.sh R/select_ind.R              # debate the correctness of a file
dev/debate.sh --staged               # debate the staged diff
dev/debate.sh --task "add rrBLUP GEBV criterion per docs/THEORY_REVIEW.md M1"
```
Default: `claude` defends, `codex` is skeptic (swap with `DEFENDER=codex SKEPTIC=claude`).
Keep them **different models** — cross-model independence is what makes agreement mean
something; two of the same model can share a blind spot and confirm each other. The full
dialogue is saved to a transcript file (path printed at start and end).

> Why it works when `dual.sh` doesn't quite: it maintains a shared transcript (the
> stateless CLIs get memory), gives the defender a real rebuttal turn, parses a structured
> JSON verdict (`{"verdict":"AGREE|BLOCK",...}`) from both sides, caps rounds with a
> human-escalation referee, and grounds every claim in executed R rather than persuasion.

By default `claude` implements and `codex` reviews. Swap with
`IMPL=codex REVIEWER=claude dev/dual.sh ...`. **The implementer never reviews its own
genetics change** — that independence is the whole point.

The reviewer applies `docs/THEORY_REVIEW.md` (per-item PASS / FAIL / UNVERIFIABLE) and must
end with a fenced JSON verdict, `{"verdict":"AGREE|BLOCK","open":[...],"confidence":0-1,
"summary":"..."}`. The loop stops only when `devtools::test()` is green **and** the
reviewer's verdict is `AGREE`. It never commits for you — you inspect the diff and commit.

> **The reviewer is read-only by contract, not by sandbox.** `codex` runs with
> `-s workspace-write` (R needs temp files to produce executed evidence), so "do not edit
> source" is enforced by the prompt only; the `claude` reviewer runs in
> `--permission-mode plan`. After any review run `git status` / `git diff`, or use
> `ISOLATE=1` so a stray edit lands in a throwaway worktree.

## Isolation, provenance, and structured verdicts

- **Isolation** — `ISOLATE=1 dev/debate.sh --task "…"` runs the agent in a throwaway git
  worktree (created in `$TMPDIR`, on branch `agent/<ts>`); your working tree is never
  touched until you merge that branch. Recommended for any run that edits code. The
  scripts default `TMPDIR` to `<repo>/.tmp` (the macOS per-user temp is unwritable in some
  terminals), which on a OneDrive-synced checkout is inside the synced tree — set
  `PIPELINE_TMPDIR` to a directory outside it to avoid the churn. Provenance is always
  written to the main checkout's `dev/.audit/`, not into the worktree.
- **Provenance** — every review/debate/eval run appends ONE JSON record (models + versions,
  rubric hash, verdict, commit, transcript path) to `dev/.audit/log.jsonl` and copies the
  transcript to `dev/.audit/transcripts/`. `dev/.audit/` is gitignored (local, not shipped).
- **Structured verdicts** — reviewer/skeptic replies end in a fenced JSON
  `{"verdict":"AGREE|BLOCK","open":[…],"confidence":…}`, parsed by `dev/lib/parse_verdict.py`
  — the harness never scrapes prose to decide agreement.

## Evals — is the reviewer any good?

`evals/` meta-tests the reviewer against seeded theory bugs (see `evals/README.md`):
```bash
evals/run.sh --check    # no model: golden set still applies (also runs in CI)
evals/run.sh --eval     # score recall/precision of the reviewer on 5 seeded bugs
                        # (bugs are seeded in a throwaway git worktree, never your tree)
```

## Guarantees

- **No AI co-author, ever.** `.githooks/commit-msg` rejects any trailer
  (`Co-Authored-By`, `Assisted-by`, `Generated-by`, `Reviewed-by`, `Signed-off-by`, ...)
  naming an assistant, plus "generated/made/written with|by <assistant>" lines and 🤖 —
  whichever model wrote the message. The patterns live in one file, `.githooks/ai-patterns`,
  sourced by both the hook and the CI workflow (`attribution-guard.yml`) so they cannot
  drift; `bash dev/test-attribution-guard.sh` checks 40+ sample messages. Human co-authors
  (including GitHub `noreply` addresses) are allowed. Bypassable locally with
  `--no-verify`; the CI guard is the backstop.
- **Both models read the same rules.** `AGENTS.md` (committed) is canonical; `CLAUDE.md`
  (gitignored) just points to it, and Codex reads `AGENTS.md` natively. No drift.

## Adjusting CLI flags

If your `codex`/`claude` version uses different non-interactive flags, edit the
`run_claude_*` / `run_codex_*` functions at the top of `dual.sh` — they are the only
place invocation flags live.
