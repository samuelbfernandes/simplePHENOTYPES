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
dev/dual.sh review R/select.R R/ocs.R      # or specific files

# Let one model implement and the other review+critique, iterating to green:
dev/dual.sh loop "add rrBLUP GEBV as an `on` criterion, per docs/THEORY_REVIEW.md S5"

# Just the objective gate:
dev/dual.sh check
```

## Skeptical debate to consensus (`dev/debate.sh`)

Stronger than a one-shot review: a **SKEPTIC** must be *convinced* and a **DEFENDER** must
fix or rebut each objection **with an executed run** (rhetoric loses to a failing/passing
`Rscript`). It proceeds only when the skeptic explicitly agrees **and** `devtools::test()`
is green; on deadlock it **escalates to you** and never auto-proceeds or commits.

```bash
dev/debate.sh R/select.R              # debate the correctness of a file
dev/debate.sh --staged               # debate the staged diff
dev/debate.sh --task "add rrBLUP GEBV criterion per docs/THEORY_REVIEW.md M1"
```
Default: `claude` defends, `codex` is skeptic (swap with `DEFENDER=codex SKEPTIC=claude`).
Keep them **different models** — cross-model independence is what makes agreement mean
something; two of the same model can share a blind spot and confirm each other. The full
dialogue is saved to a transcript file (path printed at start and end).

> Why it works when `dual.sh` doesn't quite: it maintains a shared transcript (the
> stateless CLIs get memory), gives the defender a real rebuttal turn, parses a structured
> `VERDICT: AGREE|BLOCK` from both sides, caps rounds with a human-escalation referee, and
> grounds every claim in executed R rather than persuasion.

By default `claude` implements and `codex` reviews. Swap with
`IMPL=codex REVIEWER=claude dev/dual.sh ...`. **The implementer never reviews its own
genetics change** — that independence is the whole point.

The reviewer applies `docs/THEORY_REVIEW.md` and must end with `THEORY: PASS|FAIL`. The
loop stops only when `devtools::test()` is green **and** the reviewer returns PASS. It
never commits for you — you inspect the diff and commit.

## Guarantees

- **No AI co-author, ever.** `.githooks/commit-msg` rejects any `Co-Authored-By` naming
  an assistant, plus "generated with"/🤖 lines — whichever model wrote the message.
- **Both models read the same rules.** `AGENTS.md` (committed) is canonical; `CLAUDE.md`
  (gitignored) just points to it, and Codex reads `AGENTS.md` natively. No drift.

## Adjusting CLI flags

If your `codex`/`claude` version uses different non-interactive flags, edit the
`run_claude_*` / `run_codex_*` functions at the top of `dual.sh` — they are the only
place invocation flags live.
