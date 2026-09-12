#!/usr/bin/env bash
# debate.sh — adversarial debate-to-consensus between two models (Claude + Codex).
#
# A SKEPTIC (reviewer) must be *convinced* the implementation is correct; a DEFENDER
# (implementer) must fix or rebut each objection with EXECUTED EVIDENCE. The change
# proceeds only when the skeptic explicitly agrees AND the objective gate is green.
# On deadlock it escalates to you — it never auto-proceeds on disagreement, and it
# never commits.
#
#   dev/debate.sh <paths...>        # debate the correctness of these files
#   dev/debate.sh --staged         # debate the staged diff
#   dev/debate.sh --task "<task>"  # DEFENDER implements first, then debate
#
# Env: DEFENDER=claude|codex (default claude), SKEPTIC=codex|claude (default the other),
#      MAX_ROUNDS=5.
#
# WHY it's more than dual.sh: (1) a shared transcript gives the stateless CLIs memory
# each round; (2) the defender gets a rebuttal turn, not just "fix"; (3) both emit a
# machine-parseable verdict; (4) a max-rounds+escalate referee; (5) claims are settled by
# runnable R, not rhetoric — two models agreeing is not proof unless it's grounded.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

RUBRIC="docs/THEORY_REVIEW.md"
DEFENDER="${DEFENDER:-claude}"
SKEPTIC="${SKEPTIC:-$([ "$DEFENDER" = claude ] && echo codex || echo claude)}"
MAX_ROUNDS="${MAX_ROUNDS:-5}"
TRANSCRIPT="${TMPDIR:-/tmp}/debate-$(date +%Y%m%d-%H%M%S).md"

have() { command -v "$1" >/dev/null 2>&1; }
require() { have "$1" || { echo "✖ '$1' not on PATH." >&2; [ "$1" = codex ] && echo "  npm i -g @openai/codex, or: brew install codex" >&2; exit 127; }; }
[ "$DEFENDER" = "$SKEPTIC" ] && echo "⚠ defender and skeptic are the SAME model — cross-model independence is the whole point." >&2

# --- model calls. Defender may edit; skeptic is read-only. Flags live only here. ----
call_claude_rw() { claude -p "$1" --permission-mode acceptEdits \
                     --allowedTools "Edit,Write,Read,Grep,Glob,Bash(Rscript:*),Bash(cd src/rust*)"; }
call_claude_ro() { claude -p "$1" --permission-mode plan \
                     --allowedTools "Read,Grep,Glob,Bash(Rscript:*),Bash(git diff:*)"; }
call_codex_rw()  { codex exec -s workspace-write --skip-git-repo-check "$1" </dev/null; }
call_codex_ro()  { codex exec -s read-only       --skip-git-repo-check "$1" </dev/null; }
defender() { case "$DEFENDER" in claude) call_claude_rw "$1";; codex) call_codex_rw "$1";; esac; }
skeptic()  { case "$SKEPTIC"  in claude) call_claude_ro "$1";; codex) call_codex_ro "$1";; esac; }

log() { printf '\n\n---\n## %s\n\n%s\n' "$1" "$2" >> "$TRANSCRIPT"; echo "── $1 ──"; }

# --- scope ------------------------------------------------------------------------
TASK=""
if [ "${1:-}" = "--task" ]; then TASK="${2:?task text}"; SCOPE="task: $TASK";
elif [ "${1:-}" = "--staged" ]; then SCOPE="staged diff"; CHANGE="$(git diff --staged)";
elif [ $# -gt 0 ]; then SCOPE="files: $*"; CHANGE="$(git diff -- "$@" 2>/dev/null; echo; echo '--- contents ---'; cat "$@" 2>/dev/null)";
else SCOPE="unstaged working tree"; CHANGE="$(git diff)"; fi

require "$DEFENDER"; require "$SKEPTIC"
echo "defender=$DEFENDER  skeptic=$SKEPTIC  max_rounds=$MAX_ROUNDS"
echo "transcript: $TRANSCRIPT"
{ echo "# Debate — $SCOPE"; echo "defender=$DEFENDER skeptic=$SKEPTIC  $(date)"; } > "$TRANSCRIPT"

gate() { echo "→ gate: devtools::test()"; Rscript -e "options(crayon.enabled=FALSE); devtools::test(stop_on_failure=TRUE)"; }

RULES='RULES OF THE DEBATE:
- Ground every claim in EXECUTED evidence. An objection stands only if you can exhibit
  concrete inputs that produce a wrong result (run R with Rscript to show it). A rebuttal
  stands only if you can exhibit a run showing correct behavior. Rhetoric loses to a
  failing/passing run.
- Do not invent citations or page numbers; mark UNVERIFIABLE if you cannot confirm.
- Be specific: cite file:line and the equation as implemented.'

if [ -n "$TASK" ]; then
  log "DEFENDER implements (round 0)" "$(defender "Implement in this package, following AGENTS.md. Do not commit. Task: $TASK
Then briefly describe what you changed and why it is correct. $RULES")"
  CHANGE="$(git diff)"
fi

# --- SKEPTIC opens ----------------------------------------------------------------
OPEN="$(skeptic "You are the SKEPTIC (model: $SKEPTIC), reviewing a change for GENETIC-THEORY
correctness. Apply the rubric $RUBRIC. Try hard to FALSIFY it. List each objection as
O1, O2, … with: rubric id, file:line, the flaw, and a concrete falsifying case (run R to
demonstrate if you can). $RULES

End your message with two lines EXACTLY:
VERDICT: BLOCK        (or 'VERDICT: AGREE' only if you genuinely found nothing)
OPEN: O1,O2,...       (the ids still unresolved; empty if AGREE)

Scope: $SCOPE
<change>
$CHANGE
</change>")"
log "SKEPTIC opening ($SKEPTIC)" "$OPEN"

verdict_of() { printf '%s' "$1" | grep -oiE 'VERDICT:[[:space:]]*(AGREE|BLOCK)' | tail -1 | grep -oiE 'AGREE|BLOCK' | tr a-z A-Z; }

if [ "$(verdict_of "$OPEN")" = AGREE ]; then
  echo "Skeptic found nothing on opening. Running gate…"
else
  for r in $(seq 1 "$MAX_ROUNDS"); do
    echo "=== round $r ==="
    REB="$(defender "You are the DEFENDER (model: $DEFENDER). The skeptic raised objections below.
For EACH open item: either FIX the code (edit files), or REBUT it by exhibiting an executed
run that proves correctness, or CONCEDE. Do not touch unrelated code; do not commit. $RULES

End with one line per item, e.g. 'STANCE O1: FIXED', 'STANCE O2: REBUT', 'STANCE O3: CONCEDE',
then a one-paragraph summary.

<skeptic>
$OPEN
</skeptic>")"
    log "DEFENDER round $r ($DEFENDER)" "$REB"
    CHANGE="$(git diff)"

    OPEN="$(skeptic "You are the SKEPTIC ($SKEPTIC). The defender responded below and may have edited
the code. Re-examine EACH item. Concede only when an executed run convinces you; otherwise
PRESS with new evidence. Do not soften just to reach agreement. $RULES

End with two lines EXACTLY:
VERDICT: AGREE|BLOCK
OPEN: <ids still unresolved, empty if AGREE>

<defender>
$REB
</defender>
<current-change>
$CHANGE
</change>")"
    log "SKEPTIC round $r ($SKEPTIC)" "$OPEN"
    [ "$(verdict_of "$OPEN")" = AGREE ] && break
  done
fi

echo
if [ "$(verdict_of "$OPEN")" = AGREE ]; then
  if gate; then
    echo "✓ CONSENSUS: skeptic agrees AND devtools::test() is green."
    echo "  Review the diff and the transcript, then commit yourself: $TRANSCRIPT"
    exit 0
  else
    echo "✖ Skeptic agreed but the OBJECTIVE GATE FAILED. Do not proceed — tests must pass."
    echo "  Transcript: $TRANSCRIPT"; exit 1
  fi
else
  echo "⚠ NO CONSENSUS after $MAX_ROUNDS rounds — ESCALATED TO YOU."
  echo "  Unresolved: $(printf '%s' "$OPEN" | grep -iE '^OPEN:' | tail -1)"
  echo "  Read the full dialogue: $TRANSCRIPT"; exit 2
fi
