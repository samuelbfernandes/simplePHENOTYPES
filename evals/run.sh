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
review_diff() { # <scope-label>  (reads the worktree's diff via git); echoes model text
  local prompt diff; diff="$(git diff)"
  prompt="You are an INDEPENDENT reviewer (model: $REVIEWER) checking a change for
GENETIC-THEORY correctness. Apply $RUBRIC. Try to falsify. Do not edit files.
$(verdict_instruction)
Scope: $1
<change>
$diff
</change>"
  case "$REVIEWER" in
    codex)  codex exec -s read-only --skip-git-repo-check "$prompt" </dev/null;;
    claude) claude -p "$prompt" --permission-mode plan --allowedTools "Read,Grep,Glob,Bash(git diff:*)";;
  esac
}

cmd_eval() {
  command -v "$REVIEWER" >/dev/null || { echo "✖ reviewer '$REVIEWER' not on PATH"; exit 127; }
  local total; total="$(n)"; local caught=0; local localized=0
  [ "$total" -eq 0 ] && { echo "no mutations defined in $MUT — nothing to eval"; return 0; }
  echo "→ eval: reviewer=$REVIEWER over $total seeded bugs + 1 clean control"; echo

  # precision control: clean code should NOT be blocked
  local wt br; IFS='|' read -r wt br <<< "$(wt_create)"
  ( cd "$wt"; echo "# eval clean control" >> .eval_touch 2>/dev/null || true )
  local clean; clean="$(cd "$wt" && review_diff "CLEAN CONTROL (no seeded bug)")"
  local clean_v; clean_v="$(parse_verdict "$clean")"
  wt_cleanup "$wt" "$br"
  echo "clean control → verdict=$clean_v (want AGREE)"; echo

  for i in $(seq 0 $((total-1))); do
    local id file find repl rubric
    id="$(field "$i" id)"; file="$(field "$i" file)"; find="$(field "$i" find)"
    repl="$(field "$i" replace)"; rubric="$(field "$i" rubric)"
    IFS='|' read -r wt br <<< "$(wt_create)"
    ( cd "$wt" && apply_mut "$file" "$find" "$repl" >/dev/null ) || { echo "  ✖ $id: could not seed"; wt_cleanup "$wt" "$br"; continue; }
    local out v; out="$(cd "$wt" && review_diff "seeded bug $id ($rubric)")"; v="$(parse_verdict "$out")"
    local hit=no; [ "$v" = BLOCK ] && { caught=$((caught+1)); hit=yes; }
    local loc=no; printf '%s' "$out" | grep -qw "$rubric" && { localized=$((localized+1)); loc=yes; }
    printf '  %s %-26s verdict=%-7s (want BLOCK)  rubric-cited(%s)=%s\n' \
      "$([ "$hit" = yes ] && echo ✓ || echo ✖)" "$id" "$v" "$rubric" "$loc"
    audit_record "eval:$id" "$rubric" "$v" "-" "$REVIEWER" "-" >/dev/null 2>&1 || true
    KEEP="${KEEP:-0}"; wt_cleanup "$wt" "$br" "$KEEP"
  done

  echo; echo "SCORECARD (reviewer=$REVIEWER)"
  echo "  recall   (bugs blocked):     $caught/$total"
  echo "  localized(right rubric id):  $localized/$total"
  echo "  precision(clean not blocked): $([ "$clean_v" = AGREE ] && echo 'PASS (AGREE)' || echo "SUSPECT ($clean_v)")"
  [ "$caught" -eq "$total" ] && [ "$clean_v" = AGREE ]
}

case "${1:---check}" in
  --check) cmd_check;;
  --eval)  cmd_eval;;
  *) sed -n '2,15p' "$0"; exit 2;;
esac
