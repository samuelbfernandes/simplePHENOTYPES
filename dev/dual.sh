#!/usr/bin/env bash
# dual.sh — two-model development driver (Claude Code + Codex).
#
# The point: one model implements, the OTHER independently reviews — focused on
# genetic-theory correctness (docs/THEORY_REVIEW.md), not just passing tests.
#
#   dev/dual.sh review [--staged | <paths...>]     # standalone independent theory review
#   dev/dual.sh loop "<task>" [max_iters]          # implement -> review -> fix, until green
#   dev/dual.sh check                              # just run the objective gate (tests)
#
# Env knobs:
#   IMPL=claude|codex     implementer  (default: claude)
#   REVIEWER=codex|claude reviewer     (default: the model that is NOT the implementer)
#   MAX_ITERS=4           loop cap
#
# Flags for each CLI are centralized in run_claude()/run_codex() below — if your codex
# version uses different flags, fix them there only.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

RUBRIC="docs/THEORY_REVIEW.md"
IMPL="${IMPL:-claude}"
REVIEWER="${REVIEWER:-$([ "$IMPL" = claude ] && echo codex || echo claude)}"
MAX_ITERS="${MAX_ITERS:-4}"

have() { command -v "$1" >/dev/null 2>&1; }

# --- model invocations (adjust flags here if your CLI versions differ) ------------
# Claude Code headless: -p prints the reply. acceptEdits lets it write during 'loop';
# review runs in plan mode (read-only). Scope tools tightly for safety.
run_claude_impl()   { claude -p "$1" --permission-mode acceptEdits \
                        --allowedTools "Edit,Write,Read,Grep,Glob,Bash(Rscript:*),Bash(cd src/rust*)"; }
run_claude_review() { claude -p "$1" --permission-mode plan \
                        --allowedTools "Read,Grep,Glob,Bash(Rscript:*),Bash(git diff:*)"; }
# Codex headless: exec = non-interactive. full-auto writes in the workspace; read-only
# for review. (OpenAI Codex CLI: `codex exec`, `--full-auto`, `-s read-only`.)
# stdin from /dev/null so codex never blocks waiting on it inside a script.
run_codex_impl()    { codex exec -s workspace-write --skip-git-repo-check "$1" </dev/null; }
run_codex_review()  { codex exec -s read-only     --skip-git-repo-check "$1" </dev/null; }

impl()   { case "$IMPL"     in claude) run_claude_impl   "$1";; codex) run_codex_impl   "$1";; esac; }
review() { case "$REVIEWER" in claude) run_claude_review "$1";; codex) run_codex_review "$1";; esac; }

require_model() {
  if ! have "$1"; then
    echo "✖ '$1' CLI not found on PATH." >&2
    [ "$1" = codex ] && echo "  Install: npm i -g @openai/codex  (you have it in VS Code; the pipeline needs the CLI)." >&2
    exit 127
  fi
}

# --- objective gate --------------------------------------------------------------
gate() {
  echo "→ objective gate: devtools::test()"
  Rscript -e "options(crayon.enabled=FALSE); devtools::test(stop_on_failure=TRUE)"
}

# --- standalone theory review ----------------------------------------------------
cmd_review() {
  require_model "$REVIEWER"
  local scope diff
  if [ "${1:-}" = "--staged" ]; then
    diff="$(git diff --staged)"; scope="staged changes"
  elif [ $# -gt 0 ]; then
    diff="$(git diff -- "$@"; echo; echo '--- current contents ---'; cat "$@" 2>/dev/null)"; scope="$*"
  else
    diff="$(git diff)"; scope="unstaged working tree"
  fi
  [ -z "$diff" ] && { echo "Nothing to review in: $scope"; exit 0; }

  local prompt
  prompt="You are the INDEPENDENT reviewer (model: $REVIEWER). Review the change below for
GENETIC-THEORY correctness and implementation bugs. Apply the rubric in $RUBRIC exactly,
including its required output format and final 'THEORY: PASS|FAIL' line. Do NOT edit any
file. Try to falsify each claim with a concrete numerical counterexample. Do not invent
citations or page numbers — mark UNVERIFIABLE if you cannot confirm.

Scope: $scope

<change>
$diff
</change>"
  echo "→ independent theory review by: $REVIEWER  (scope: $scope)"
  review "$prompt"
}

# --- implement -> review -> fix loop --------------------------------------------
cmd_loop() {
  local task="${1:?usage: dev/dual.sh loop \"<task>\" [max_iters]}"
  MAX_ITERS="${2:-$MAX_ITERS}"
  require_model "$IMPL"; require_model "$REVIEWER"
  echo "→ implementer=$IMPL  reviewer=$REVIEWER  max_iters=$MAX_ITERS"

  impl "Implement the following in this package, following AGENTS.md. Do not commit.
Task: $task"

  local i
  for i in $(seq 1 "$MAX_ITERS"); do
    echo "=== iteration $i: gate + independent review ==="
    if ! gate; then
      impl "devtools::test() failed. Fix the failures without touching unrelated code, then stop."
      continue
    fi
    local verdict
    verdict="$(cmd_review 2>&1)"; echo "$verdict"
    if printf '%s' "$verdict" | grep -qE 'THEORY:[[:space:]]*PASS'; then
      echo "✓ tests green AND theory review PASS. Review the diff, then commit yourself."
      return 0
    fi
    impl "The independent theory reviewer ($REVIEWER) reported issues below. Address each
FAIL against its primary source; do not touch unrelated code; then stop.
$verdict"
  done
  echo "⚠ reached max_iters ($MAX_ITERS) without a clean PASS — inspect manually."; return 1
}

case "${1:-}" in
  review) shift; cmd_review "$@";;
  loop)   shift; cmd_loop   "$@";;
  check)  gate;;
  *) sed -n '2,12p' "$0"; exit 2;;
esac
