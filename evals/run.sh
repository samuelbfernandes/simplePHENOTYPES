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
# Env for --eval: REVIEWER=codex|claude (default codex), keep the worktree with KEEP=1.
#   (REVIEWER_CMD=<cmd> substitutes any command that takes the prompt as $1 — test hook.)
# Neither mode writes to a source file in the working tree; --eval seeds bugs only inside a
# throwaway git worktree and verifies the targets byte-for-byte on exit.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Writable repo-local temp before any git/R call (see dev/dual.sh for why).
export TMPDIR="${PIPELINE_TMPDIR:-$DIR/../.tmp}"; mkdir -p "$TMPDIR"
cd "$(git rev-parse --show-toplevel)"
DEV="$DIR/../dev"
. "$DEV/lib/verdict.sh"; . "$DEV/lib/audit.sh"; . "$DEV/lib/worktree.sh"
MUT="$DIR/mutations.json"
REVIEWER="${REVIEWER:-codex}"
RUBRIC="docs/THEORY_REVIEW.md"

n() { jq 'length' "$MUT"; }
field() { jq -r ".[$1].$2" "$MUT"; }

# Apply an exact-substring mutation to a COPY (never in place): reads <src>, writes the
# mutated text to <dst>. Exit 3 when the find string is absent. Python, not sed, so there
# are no regex surprises. Source files in the working tree are never written by this script.
apply_mut() { # <src> <dst> <find> <replace>
  python3 - "$1" "$2" "$3" "$4" <<'PY'
import sys
src,dst,find,repl=sys.argv[1:5]
s=open(src).read()
if find not in s:
    sys.stderr.write("FIND-NOT-PRESENT\n"); sys.exit(3)
open(dst,"w").write(s.replace(find,repl,1)); print("OK")
PY
}

mk_tmp() { mktemp "${TMPDIR:-/tmp}/$1.XXXXXX"; }   # bare mktemp ignores TMPDIR on macOS

# ---------- --check : golden set is valid, no model -------------------------------
# Mutations are applied to a temporary COPY of each target, so --check never writes to a
# source file (an interrupt cannot leave a seeded bug behind).
cmd_check() {
  local total; total="$(n)"; local ok=0 present=0 skipped=0
  [ "$total" -eq 0 ] && { echo "no mutations defined in $MUT (empty golden set) — OK"; return 0; }
  echo "→ validating $total seeded mutations apply & differ from the source (no model)"
  local out; out="$(mk_tmp evals-check)"
  # shellcheck disable=SC2064
  trap "rm -f '$out'" EXIT
  for i in $(seq 0 $((total-1))); do
    local id file find repl
    id="$(field "$i" id)"; file="$(field "$i" file)"
    find="$(field "$i" find)"; repl="$(field "$i" replace)"
    if [ ! -f "$file" ]; then echo "  ○ $id: $file not present yet — skipped"; skipped=$((skipped+1)); continue; fi
    present=$((present+1))
    if apply_mut "$file" "$out" "$find" "$repl" >/dev/null 2>&1 && ! cmp -s "$file" "$out"; then
      echo "  ✓ $id ($file) — applies"
      ok=$((ok+1))
    else
      echo "  ✖ $id ($file) — find string not present (source moved? update mutations.json)"
    fi
  done
  echo "check: $ok/$present present mutations valid ($skipped skipped, file not committed yet)"
  [ "$ok" -eq "$present" ]                    # green when every PRESENT target is valid
}

# ---------- --eval : does the reviewer catch them? --------------------------------
# Every seeded bug lives ONLY inside a throwaway `git worktree` under $TMPDIR (the same
# mechanism as dev/debate.sh ISOLATE=1). The working tree is never written: the target
# files are copied INTO the worktree first (so untracked / uncommitted WIP is what gets
# reviewed), the mutation is applied there, and the reviewer runs with the worktree as its
# cwd. A snapshot of every target is taken up front and verified byte-for-byte on EXIT as a
# belt-and-braces guard; a mismatch is restored from the snapshot and reported loudly.
review_file() { # <file>  (reads the file relative to the current directory)
  local body prompt; body="$(cat "$1")"
  prompt="You are an INDEPENDENT reviewer (model: $REVIEWER) checking ONE source file for
GENETIC-THEORY correctness. It may or may not contain a seeded error. Apply $RUBRIC and
try to FALSIFY each formula against its primary source. Do not edit files.
$(verdict_instruction)
File: $1
<file>
$body
</file>"
  if [ -n "${REVIEWER_CMD:-}" ]; then   # test hook: any command taking the prompt as $1
    "$REVIEWER_CMD" "$prompt" </dev/null
    return
  fi
  case "$REVIEWER" in
    codex)  codex exec -s workspace-write --skip-git-repo-check "$prompt" </dev/null;;
    claude) claude -p "$prompt" --permission-mode plan --allowedTools "Read,Grep,Glob,Bash(Rscript:*)";;
    *) echo "unknown REVIEWER '$REVIEWER' (codex|claude)" >&2; return 127;;
  esac
}

EVAL_ROOT=""; EVAL_SNAP=""; EVAL_WT=""; EVAL_BR=""; EVAL_FILES=""
eval_cleanup() {  # EXIT trap for --eval: verify sources untouched, then remove the worktree
  local rc=$? f snapf
  trap - EXIT INT TERM
  if [ -n "$EVAL_SNAP" ] && [ -d "$EVAL_SNAP" ]; then
    for f in $EVAL_FILES; do
      snapf="$EVAL_SNAP/$(printf '%s' "$f" | tr '/' '_')"
      if [ -f "$snapf" ] && ! cmp -s "$snapf" "$EVAL_ROOT/$f"; then
        echo "✖ INTEGRITY: $f differs from its pre-eval snapshot — restoring exact bytes" >&2
        cp "$snapf" "$EVAL_ROOT/$f" && cmp -s "$snapf" "$EVAL_ROOT/$f" \
          && echo "  restored $f" >&2 || echo "  ✖ RESTORE FAILED for $f — recover with git checkout / your editor history" >&2
      fi
    done
    rm -rf "$EVAL_SNAP"
  fi
  if [ -n "$EVAL_WT" ]; then
    if [ "${KEEP:-0}" = 1 ]; then echo "worktree kept: $EVAL_WT (branch $EVAL_BR)" >&2
    else wt_cleanup "$EVAL_WT" "$EVAL_BR"; fi
  fi
  exit "$rc"
}

cmd_eval() {
  if [ -z "${REVIEWER_CMD:-}" ]; then
    command -v "$REVIEWER" >/dev/null || { echo "✖ reviewer '$REVIEWER' not on PATH"; exit 127; }
  fi
  local total; total="$(n)"; local caught=0 localized=0 present=0 errors=0
  [ "$total" -eq 0 ] && { echo "no mutations defined in $MUT — nothing to eval"; return 0; }
  echo "→ eval: reviewer=$REVIEWER over $total seeded bugs (+1 clean control)"; echo

  # Isolation + integrity guard, armed BEFORE anything is created.
  EVAL_ROOT="$(pwd)"; EVAL_SNAP="$(mktemp -d "${TMPDIR:-/tmp}/evals-snap.XXXXXX")"
  EVAL_FILES="$(jq -r '[.[].file]|unique|.[]' "$MUT" | tr '\n' ' ')"
  trap eval_cleanup EXIT
  trap 'exit 130' INT TERM
  local f
  for f in $EVAL_FILES; do
    [ -f "$f" ] && cp "$f" "$EVAL_SNAP/$(printf '%s' "$f" | tr '/' '_')"
  done
  IFS='|' read -r EVAL_WT EVAL_BR <<< "$(wt_create)"
  echo "→ throwaway worktree: $EVAL_WT (branch $EVAL_BR); the working tree is never modified"
  for f in $EVAL_FILES; do   # WIP contents (even untracked) are what gets reviewed
    if [ -f "$f" ]; then mkdir -p "$EVAL_WT/$(dirname "$f")"; cp "$f" "$EVAL_WT/$f"; fi
  done

  # precision control: a real target file, UNMUTATED, should NOT be blocked.
  local ctrl clean clean_v="n/a" rc=0; ctrl="$(field 0 file)"
  if [ -f "$ctrl" ]; then
    clean="$(cd "$EVAL_WT" && review_file "$ctrl")" || rc=$?
    if [ "$rc" -ne 0 ]; then clean_v="ERROR(rc=$rc)"; errors=$((errors+1)); else clean_v="$(parse_verdict "$clean")"; fi
  fi
  echo "clean control ($ctrl) → verdict=$clean_v (want AGREE)"; echo

  for i in $(seq 0 $((total-1))); do
    local id file find repl rubric out v
    id="$(field "$i" id)"; file="$(field "$i" file)"; find="$(field "$i" find)"
    repl="$(field "$i" replace)"; rubric="$(field "$i" rubric)"
    if [ ! -f "$file" ]; then echo "  ○ $id: $file not present — skipped"; continue; fi
    present=$((present+1))
    # start from the pristine WIP copy every time, then seed the bug in the worktree only
    if ! apply_mut "$EVAL_SNAP/$(printf '%s' "$file" | tr '/' '_')" "$EVAL_WT/$file" "$find" "$repl" >/dev/null 2>&1; then
      echo "  ✖ $id: could not seed (find string not present)"; continue
    fi
    rc=0; out="$(cd "$EVAL_WT" && review_file "$file")" || rc=$?
    if [ "$rc" -ne 0 ]; then v="ERROR"; errors=$((errors+1)); else v="$(parse_verdict "$out")"; fi
    local hit=no; [ "$v" = BLOCK ] && { caught=$((caught+1)); hit=yes; }
    local loc=no; printf '%s' "${out:-}" | grep -qw "$rubric" && { localized=$((localized+1)); loc=yes; }
    printf '  %s %-26s verdict=%-7s (want BLOCK)  rubric-cited(%s)=%s\n' \
      "$([ "$hit" = yes ] && echo ✓ || echo ✖)" "$id" "$v" "$rubric" "$loc"
    audit_record "eval:$id" "$rubric" "$v" "-" "$REVIEWER" "-" >/dev/null 2>&1 || true
  done

  echo; echo "SCORECARD (reviewer=$REVIEWER)"
  echo "  recall   (bugs blocked):     $caught/$present present"
  echo "  localized(right rubric id):  $localized/$present"
  echo "  reviewer errors (non-zero exit): $errors"
  echo "  precision(clean not blocked): $([ "$clean_v" = AGREE ] && echo 'PASS (AGREE)' || echo "SUSPECT ($clean_v)")"
  [ "$errors" -eq 0 ] && [ "$caught" -eq "$present" ] && [ "$clean_v" = AGREE ]
}

case "${1:---check}" in
  --check) cmd_check;;
  --eval)  cmd_eval;;
  *) sed -n '2,15p' "$0"; exit 2;;
esac
