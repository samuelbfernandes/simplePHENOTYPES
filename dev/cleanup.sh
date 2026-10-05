#!/usr/bin/env bash
# Remove regenerable, git-ignored build and render artifacts from a checkout.
#
#   dev/cleanup.sh [checkout]          # dry run: list what would be removed
#   dev/cleanup.sh [checkout] --apply  # remove it
#
# Only ignored files that are re-created by rendering, checking or testing are
# touched. Tracked files, the compiled package (src/), .tmp/ (it holds audit
# evidence), .claude/, context/ and graphify-out/ are never removed.
set -euo pipefail

root="${1:-.}"
apply=false
[[ "${2:-}" == "--apply" || "${1:-}" == "--apply" ]] && apply=true
[[ "${1:-}" == "--apply" ]] && root="."
cd "$root"
git rev-parse --is-inside-work-tree >/dev/null

targets=(
  .quarto
  docs/ARCHITECTURE_DIAGRAM.html
  docs/ARCHITECTURE_DIAGRAM_files
  tests/testthat/_problems
  tests/testthat/testthat-problems.rds
  tests/testthat/Rplots.pdf
  .Rhistory
)
# vignette purl / render output (vignettes/<name>.R and .html next to <name>.Rmd)
for rmd in vignettes/*.Rmd; do
  stem="${rmd%.Rmd}"
  targets+=("$stem.R" "$stem.html")
done

found=0
for t in "${targets[@]}"; do
  [[ -e "$t" ]] || continue
  # never touch a tracked path
  if git ls-files --error-unmatch "$t" >/dev/null 2>&1; then
    echo "skip (tracked): $t"
    continue
  fi
  found=1
  if $apply; then
    rm -rf -- "$t"
    echo "removed: $t"
  else
    echo "would remove: $t"
  fi
done
[[ $found -eq 1 ]] || echo "nothing to clean"
$apply || echo "(dry run; add --apply to remove)"
