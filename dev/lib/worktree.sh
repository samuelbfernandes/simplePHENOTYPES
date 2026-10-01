#!/usr/bin/env bash
# worktree.sh — run an agent in an isolated git worktree so its edits never touch your
# working tree until you choose to bring them in. Sourced, not executed.
#
# Worktrees are created under $TMPDIR. The pipeline scripts default TMPDIR to <repo>/.tmp
# (gitignored) because the macOS per-user temp is unwritable in some terminals; on a
# OneDrive-synced checkout that means the throwaway worktree is inside the synced tree —
# set PIPELINE_TMPDIR to a directory outside it to avoid the churn. A run lands on branch
# agent/<ts>; nothing merges automatically — you review the branch and merge or delete it.

# wt_create -> prints "<worktree_path>|<branch>"
wt_create() {
  local ts br wt
  ts="$(date +%Y%m%d-%H%M%S)-$$"
  br="agent/$ts"
  wt="${TMPDIR:-/tmp}/agent-wt-$ts"
  git worktree add -q -b "$br" "$wt" HEAD >&2
  echo "$wt|$br"
}

# wt_cleanup <worktree_path> <branch> [keep]
# keep=1 leaves the branch (edits worth reviewing); default removes both.
wt_cleanup() {
  git worktree remove --force "$1" 2>/dev/null || true
  [ "${3:-}" = "1" ] || git branch -D "$2" 2>/dev/null || true
}
