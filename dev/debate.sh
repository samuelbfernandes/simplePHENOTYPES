#!/usr/bin/env bash
# debate.sh — adversarial debate-to-consensus between two models (Claude + Codex).
#
# A SKEPTIC (reviewer) must be *convinced* the implementation is correct; a DEFENDER
# (implementer) must fix or rebut each objection with EXECUTED EVIDENCE. The change
# proceeds only when the skeptic emits a structured verdict AGREE AND the objective gate
# is green. On deadlock it escalates to you; it never auto-proceeds and never commits.
#
#   dev/debate.sh <paths...>        # debate the correctness of these files
#   dev/debate.sh --staged         # debate the staged diff
#   dev/debate.sh --task "<task>"  # DEFENDER implements first, then debate
#
# Env: DEFENDER=claude|codex (default claude), SKEPTIC=codex|claude (default the other),
#      MAX_ROUNDS=5, ISOLATE=1 (run in a throwaway git worktree; recommended for --task).
#
# What makes this more than a one-shot review: a shared transcript gives the stateless
# CLIs memory each round; the defender gets a rebuttal turn; BOTH sides emit a
# machine-parseable JSON verdict; a max-rounds+escalate referee; claims settled by
# runnable R, not rhetoric; and every run is recorded to dev/.audit/.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Writable repo-local temp before any git/R call (see dual.sh for why).
export TMPDIR="${PIPELINE_TMPDIR:-$DIR/../.tmp}"; mkdir -p "$TMPDIR"
cd "$(git rev-parse --show-toplevel)"
. "$DIR/lib/verdict.sh"; . "$DIR/lib/audit.sh"; . "$DIR/lib/worktree.sh"

RUBRIC="docs/THEORY_REVIEW.md"
DEFENDER="${DEFENDER:-claude}"
SKEPTIC="${SKEPTIC:-$([ "$DEFENDER" = claude ] && echo codex || echo claude)}"
MAX_ROUNDS="${MAX_ROUNDS:-5}"
ISOLATE="${ISOLATE:-0}"
TRANSCRIPT="${TMPDIR:-/tmp}/debate-$(date +%Y%m%d-%H%M%S).md"

have() { command -v "$1" >/dev/null 2>&1; }
require() { have "$1" || { echo "✖ '$1' not on PATH." >&2; [ "$1" = codex ] && echo "  brew install codex  (or npm i -g @openai/codex)" >&2; exit 127; }; }
[ "$DEFENDER" = "$SKEPTIC" ] && echo "⚠ defender and skeptic are the SAME model — cross-model independence is the point." >&2

# --- model calls. Defender may edit; skeptic is read-only. Flags live only here. ----
call_claude_rw() { claude -p "$1" --permission-mode acceptEdits \
                     --allowedTools "Edit,Write,Read,Grep,Glob,Bash(Rscript:*),Bash(cd src/rust*)"; }
call_claude_ro() { claude -p "$1" --permission-mode plan \
                     --allowedTools "Read,Grep,Glob,Bash(Rscript:*),Bash(git diff:*)"; }
call_codex_rw()  { codex exec -s workspace-write --skip-git-repo-check "$1" </dev/null; }
# skeptic is workspace-write so it can RUN R for evidence (read-only blocks R's temp);
# the skeptic prompt forbids editing source.
call_codex_ro()  { codex exec -s workspace-write --skip-git-repo-check "$1" </dev/null; }
defender() { case "$DEFENDER" in claude) call_claude_rw "$1";; codex) call_codex_rw "$1";; esac; }
skeptic()  { case "$SKEPTIC"  in claude) call_claude_ro "$1";; codex) call_codex_ro "$1";; esac; }

log() { printf '\n\n---\n## %s\n\n%s\n' "$1" "$2" >> "$TRANSCRIPT"; echo "── $1 ──"; }
gate() { echo "→ gate: devtools::test()"; Rscript -e "options(crayon.enabled=FALSE); devtools::test(stop_on_failure=TRUE)"; }

# --- scope -----------------------------------------------------------------------
TASK=""
if [ "${1:-}" = "--task" ]; then TASK="${2:?task text}"; SCOPE="task: $TASK";
elif [ "${1:-}" = "--staged" ]; then SCOPE="staged diff"; CHANGE="$(git diff --staged)";
elif [ $# -gt 0 ]; then SCOPE="files: $*"; CHANGE="$(git diff -- "$@" 2>/dev/null; echo; echo '--- contents ---'; cat "$@" 2>/dev/null)";
else SCOPE="unstaged working tree"; CHANGE="$(git diff)"; fi

require "$DEFENDER"; require "$SKEPTIC"

# --- optional isolation ----------------------------------------------------------
WT=""; BR=""
if [ "$ISOLATE" = 1 ]; then
  IFS='|' read -r WT BR <<< "$(wt_create)"
  echo "→ isolated worktree: $WT  (branch $BR)"
  cd "$WT"
  trap 'echo "worktree kept for review: $WT (branch $BR); remove with: git worktree remove --force \"$WT\" && git branch -D \"$BR\""' EXIT
fi

echo "defender=$DEFENDER  skeptic=$SKEPTIC  max_rounds=$MAX_ROUNDS  isolate=$ISOLATE"
echo "transcript: $TRANSCRIPT"
{ echo "# Debate — $SCOPE"; echo "defender=$DEFENDER skeptic=$SKEPTIC isolate=$ISOLATE  $(date)"; } > "$TRANSCRIPT"

RULES='RULES OF THE DEBATE:
- Ground every claim in EXECUTED evidence. An objection stands only if you can exhibit
  concrete inputs producing a wrong result (run R with Rscript to show it). A rebuttal
  stands only if you can exhibit a run showing correct behavior. Rhetoric loses to a run.
- Do not invent citations or page numbers; mark UNVERIFIABLE if you cannot confirm.
- Be specific: cite file:line and the equation as implemented.'

if [ -n "$TASK" ]; then
  log "DEFENDER implements (round 0)" "$(defender "Implement in this package, following AGENTS.md. Do not commit. Task: $TASK
Then briefly describe what you changed and why it is correct. $RULES")"
  CHANGE="$(git diff)"
fi

# --- SKEPTIC opens ---------------------------------------------------------------
OPEN="$(skeptic "You are the SKEPTIC (model: $SKEPTIC), reviewing a change for GENETIC-THEORY
correctness. You REVIEW ONLY — do not modify, create, or delete any source file; you may
run Rscript and create only temporary files to gather evidence. Apply the rubric $RUBRIC.
Try hard to FALSIFY it. List each objection as O1, O2, … with rubric id, file:line, the
flaw, and a concrete falsifying case. $RULES

$(verdict_instruction)

Scope: $SCOPE
<change>
$CHANGE
</change>")"
log "SKEPTIC opening ($SKEPTIC)" "$OPEN"

finish() {  # <verdict>
  local v="$1" saved; saved="$(audit_save_transcript "$TRANSCRIPT")"
  audit_record "debate" "$SCOPE" "$v" "$DEFENDER" "$SKEPTIC" "$saved"
  echo "audit: $(audit_dir)/log.jsonl   transcript: $saved"
}

if [ "$(parse_verdict "$OPEN")" != AGREE ]; then
  for r in $(seq 1 "$MAX_ROUNDS"); do
    echo "=== round $r ==="
    REB="$(defender "You are the DEFENDER (model: $DEFENDER). The skeptic raised objections below.
For EACH open item: FIX the code, or REBUT by exhibiting an executed run that proves
correctness, or CONCEDE. Do not touch unrelated code; do not commit. $RULES
End with one line per item ('STANCE O1: FIXED|REBUT|CONCEDE') and a short summary.
<skeptic>
$OPEN
</skeptic>")"
    log "DEFENDER round $r ($DEFENDER)" "$REB"
    CHANGE="$(git diff)"
    OPEN="$(skeptic "You are the SKEPTIC ($SKEPTIC). You REVIEW ONLY — do not modify any file; you
may run Rscript and create only temp files. The defender responded below and may have
edited the code. Re-examine EACH item. Concede only when an executed run convinces you;
else PRESS with new evidence. Do not soften just to agree. $RULES

$(verdict_instruction)
<defender>
$REB
</defender>
<current-change>
$CHANGE
</current-change>")"
    log "SKEPTIC round $r ($SKEPTIC)" "$OPEN"
    [ "$(parse_verdict "$OPEN")" = AGREE ] && break
  done
fi

echo
if [ "$(parse_verdict "$OPEN")" = AGREE ]; then
  if gate; then
    echo "✓ CONSENSUS: skeptic AGREE and devtools::test() green. Review the diff + transcript, then commit."
    finish AGREE
    [ "$ISOLATE" = 1 ] && echo "changes are on branch $BR (worktree $WT) — merge it when satisfied."
    exit 0
  else
    echo "✖ Skeptic agreed but the OBJECTIVE GATE FAILED. Do not proceed."; finish BLOCK; exit 1
  fi
else
  echo "⚠ NO CONSENSUS after $MAX_ROUNDS rounds — ESCALATED TO YOU."
  echo "  unresolved(JSON): $(verdict_json "$OPEN")"; finish BLOCK; exit 2
fi
