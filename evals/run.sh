#!/usr/bin/env bash
# run.sh — evaluate the review harness against a golden set of SEEDED THEORY BUGS.
#
# You cannot trust a reviewer you have never tested. Each entry in mutations.json injects
# one known genetic-theory error into the real source; a good reviewer must BLOCK it
# (recall), and must NOT fabricate a BLOCK on clean code (precision). This meta-tests the
# reviewer/debate loop itself.
#
#   evals/run.sh --check          # (default) no model: verify every mutation applies &
#                                 # reverts cleanly. Cheap; run in CI.
#   evals/run.sh --eval           # run the reviewer on each seeded bug + a clean control,
#                                 # score recall/precision. Calls a model (costs tokens).
#
# Env for --eval: REVIEWER=codex|claude (default codex), keep worktrees with KEEP=1.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEV="$DIR/../dev"
. "$DEV/lib/verdict.sh"; . "$DEV/lib/audit.sh"; . "$DEV/lib/worktree.sh"
MUT="$DIR/mutations.json"
REVIEWER="${REVIEWER:-codex}"
RUBRIC="docs/THEORY_REVIEW.md"

n() { jq 'length' "$MUT"; }
field() { jq -r ".[$1].$2" "$MUT"; }

# apply/revert an exact-substring mutation with a python helper (no regex surprises).
apply_mut() { # <file> <find> <replace>
  python3 - "$1" "$2" "$3" <<'PY'
import sys
f,find,repl=sys.argv[1],sys.argv[2],sys.argv[3]
s=open(f).read()
if find not in s:
    sys.stderr.write("FIND-NOT-PRESENT\n"); sys.exit(3)
open(f,"w").write(s.replace(find,repl,1)); print("OK")
PY
}

# ---------- --check : golden set is valid, no model -------------------------------
cmd_check() {
  local total; total="$(n)"; local ok=0 present=0 skipped=0
  [ "$total" -eq 0 ] && { echo "no mutations defined in $MUT (empty golden set) — OK"; return 0; }
  echo "→ validating $total seeded mutations apply & revert cleanly (no model)"
  for i in $(seq 0 $((total-1))); do
    local id file find repl bak
    id="$(field "$i" id)"; file="$(field "$i" file)"
    find="$(field "$i" find)"; repl="$(field "$i" replace)"
    if [ ! -f "$file" ]; then echo "  ○ $id: $file not present yet — skipped"; skipped=$((skipped+1)); continue; fi
    present=$((present+1)); bak="$(mktemp)"; cp "$file" "$bak"
    if apply_mut "$file" "$find" "$repl" >/dev/null 2>&1 && ! diff -q "$file" "$bak" >/dev/null; then
      echo "  ✓ $id ($file) — applies"
      ok=$((ok+1))
    else
      echo "  ✖ $id ($file) — find string not present (source moved? update mutations.json)"
    fi
    cp "$bak" "$file"; rm -f "$bak"          # always restore exact bytes
  done
  echo "check: $ok/$present present mutations valid ($skipped skipped, file not committed yet)"
  [ "$ok" -eq "$present" ]                    # green when every PRESENT target is valid
}

# ---------- --eval : does the reviewer catch them? --------------------------------
# Review the FULL CONTENTS of one file. We review the working-tree file directly (not a
# git diff) because the target sources may be untracked WIP — a git diff cannot show an
# untracked file, and a worktree cut from HEAD would not contain it at all.
review_file() { # <file>
  local body prompt; body="$(cat "$1")"
  prompt="You are an INDEPENDENT reviewer (model: $REVIEWER) checking ONE source file for
GENETIC-THEORY correctness. It may or may not contain a seeded error. Apply $RUBRIC and
try to FALSIFY each formula against its primary source. Do not edit files.
$(verdict_instruction)
File: $1
<file>
$body
</file>"
  case "$REVIEWER" in
    codex)  codex exec -s read-only --skip-git-repo-check "$prompt" </dev/null;;
    claude) claude -p "$prompt" --permission-mode plan --allowedTools "Read,Grep,Glob";;
  esac
}

cmd_eval() {
  command -v "$REVIEWER" >/dev/null || { echo "✖ reviewer '$REVIEWER' not on PATH"; exit 127; }
  local total; total="$(n)"; local caught=0 localized=0 present=0
  [ "$total" -eq 0 ] && { echo "no mutations defined in $MUT — nothing to eval"; return 0; }
  echo "→ eval: reviewer=$REVIEWER over $total seeded bugs (+1 clean control)"; echo

  # precision control: a real target file, UNMUTATED, should NOT be blocked.
  local ctrl clean clean_v="n/a"; ctrl="$(field 0 file)"
  if [ -f "$ctrl" ]; then clean="$(review_file "$ctrl")"; clean_v="$(parse_verdict "$clean")"; fi
  echo "clean control ($ctrl) → verdict=$clean_v (want AGREE)"; echo

  for i in $(seq 0 $((total-1))); do
    local id file find repl rubric bak out v
    id="$(field "$i" id)"; file="$(field "$i" file)"; find="$(field "$i" find)"
    repl="$(field "$i" replace)"; rubric="$(field "$i" rubric)"
    if [ ! -f "$file" ]; then echo "  ○ $id: $file not present — skipped"; continue; fi
    present=$((present+1))
    bak="$(mktemp)"; cp "$file" "$bak"
    trap 'cp "$bak" "$file" 2>/dev/null; rm -f "$bak"' INT TERM   # restore on interrupt
    if ! apply_mut "$file" "$find" "$repl" >/dev/null 2>&1; then
      echo "  ✖ $id: could not seed (find string not present)"; cp "$bak" "$file"; rm -f "$bak"; trap - INT TERM; continue
    fi
    out="$(review_file "$file")"; v="$(parse_verdict "$out")"
    cp "$bak" "$file"; rm -f "$bak"; trap - INT TERM              # restore immediately
    local hit=no; [ "$v" = BLOCK ] && { caught=$((caught+1)); hit=yes; }
    local loc=no; printf '%s' "$out" | grep -qw "$rubric" && { localized=$((localized+1)); loc=yes; }
    printf '  %s %-26s verdict=%-7s (want BLOCK)  rubric-cited(%s)=%s\n' \
      "$([ "$hit" = yes ] && echo ✓ || echo ✖)" "$id" "$v" "$rubric" "$loc"
    audit_record "eval:$id" "$rubric" "$v" "-" "$REVIEWER" "-" >/dev/null 2>&1 || true
  done

  echo; echo "SCORECARD (reviewer=$REVIEWER)"
  echo "  recall   (bugs blocked):     $caught/$present present"
  echo "  localized(right rubric id):  $localized/$present"
  echo "  precision(clean not blocked): $([ "$clean_v" = AGREE ] && echo 'PASS (AGREE)' || echo "SUSPECT ($clean_v)")"
  [ "$caught" -eq "$present" ] && [ "$clean_v" = AGREE ]
}

case "${1:---check}" in
  --check) cmd_check;;
  --eval)  cmd_eval;;
  *) sed -n '2,15p' "$0"; exit 2;;
esac
