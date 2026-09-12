#!/usr/bin/env bash
# worktree.sh — run an agent in an isolated git worktree so its edits never touch your
# working tree until you choose to bring them in. Sourced, not executed.
#
# Worktrees are created OUTSIDE the repo (in $TMPDIR) so a OneDrive-synced project does
# not churn on throwaway checkouts. A run lands on branch agent/<ts>; nothing merges
# automatically — you review the branch and merge or delete it.

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
