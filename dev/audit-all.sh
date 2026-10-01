#!/usr/bin/env bash
# audit-all.sh — audit the whole package in cohesive groups (one reviewer call each),
# most theory-critical first, and record any issues at the TOP of TODO.md.
#
#   dev/audit-all.sh            # live v2 code: high + medium tiers        [default]
#   dev/audit-all.sh high       # only the theory-critical groups
#   dev/audit-all.sh all        # also the frozen v1 legacy engine (big, low ROI)
#   dev/audit-all.sh legacy     # only the legacy groups
#
# One `dual.sh review` (independent reviewer = codex) per group — fewer calls than
# per-file, coherent context. Groups follow the R/ filename prefixes (ARCHITECTURE.md
# §4) and cover every live v2 file. Designer files are OMITTED (moved to BD).
# Each group's transcript is saved under dev/.audit/transcripts/.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TMPDIR="${PIPELINE_TMPDIR:-$DIR/../.tmp}"; mkdir -p "$TMPDIR"
cd "$(git rev-parse --show-toplevel)"
. "$DIR/lib/verdict.sh"; . "$DIR/lib/audit.sh"
TS="$(date +%Y%m%d-%H%M%S)"; audit_init; ADIR="$(audit_dir)/transcripts"

# tier|label|paths (most theory-critical first)
# NB: not named GROUPS — that is a reserved bash array (the caller's group IDs).
AUDIT_GROUPS=(
"high|grammar|R/grammar_simulate_phenotype.R R/grammar_layers.R R/grammar_realize.R R/grammar_complex.R R/grammar_plot.R"
"high|effects-arch|R/effects_pleioarch.R R/effects_series.R R/arch_independent.R R/arch_ld.R R/qc_ld_methods.R"
"high|crossing-schemes|R/cross_population.R R/cross_mating.R R/cross_map.R R/cross_pedigree.R R/cross_mate.R R/cross_breed.R R/select_schemes.R"
"high|selection|R/select_ind.R R/select_ocs.R R/select_usefulness.R R/select_marker.R R/select_mabc.R"
"high|prediction|R/select_blup.R R/select_combining.R R/select_progeny.R"
"high|transcriptome|R/transcriptome_simulate.R R/transcriptome_layer.R R/transcriptome_counts.R R/transcriptome_mimic.R"
"high|rust-core|src/rust/src/meiosis.rs src/rust/src/genome.rs src/rust/src/numeric.rs src/rust/src/hash.rs src/rust/src/lib.rs R/extendr-wrappers.R R/io_as_numeric.R"
"med|io-formats|R/io_read_formats.R R/io_detect_format.R R/io_format_conversion.R R/io_write.R R/qc_filter_geno.R"
"legacy|legacy-core|R/legacy_create_phenotypes.R R/legacy_Phenotypes.R R/legacy_check_in.R R/legacy_constraint.R"
"legacy|legacy-qtn|R/legacy_QTN_linkage.R R/legacy_QTN_pleiotropic.R R/legacy_QTN_partially_pleiotropic.R R/legacy_qtn_from_user.R R/legacy_genetic_effect.R R/legacy_vQTL.R"
"legacy|legacy-baseline|R/legacy_Base_line_multi_traits.R R/legacy_Base_line_single_trait.R R/legacy_Genotypes.R R/legacy_make_pd.R"
)

DRY=0
if [ "${1:-}" = "--list" ]; then DRY=1; shift; fi   # print resolved groups, no model calls
case "${1:-live}" in
  high)    TIERS="high";;
  live|"") TIERS="high med";;
  all)     TIERS="high med legacy";;
  legacy)  TIERS="legacy";;
  *) echo "usage: dev/audit-all.sh [--list] [high|live|all|legacy]"; exit 2;;
esac
[ "$DRY" = 1 ] && TIERS="high med legacy"   # --list shows every group
echo "→ audit-all tiers: $TIERS   reviewer: ${REVIEWER:-codex}$([ "$DRY" = 1 ] && echo '   (DRY LIST)')"

results=()   # "label|verdict|summary|transcript"
for g in "${AUDIT_GROUPS[@]}"; do
  tier="${g%%|*}"; rest="${g#*|}"; label="${rest%%|*}"; paths="${rest#*|}"
  case " $TIERS " in *" $tier "*) : ;; *) continue;; esac
  existing=(); for p in $paths; do [ -e "$p" ] && existing+=("$p"); done
  [ ${#existing[@]} -eq 0 ] && { echo "○ $label: no files present — skip"; continue; }
  if [ "$DRY" = 1 ]; then printf "  %-18s %2d files: %s\n" "$label" "${#existing[@]}" "${existing[*]}"; continue; fi
  echo "=== auditing: $label (${#existing[@]} files) ==="
  tfile="$ADIR/audit-$label-$TS.md"
  # dual.sh writes the ONE transcript (to $tfile) and the ONE audit record for this group.
  out="$(REVIEW_TRANSCRIPT="$tfile" REVIEW_LABEL="audit:$label" dev/dual.sh review "${existing[@]}" 2>>"$tfile.err")" || true
  printf '%s\n' "$out"
  v="$(parse_verdict "$out")"
  s="$(verdict_json "$out" | jq -r '.summary // ""' 2>/dev/null || true)"
  o="$(verdict_json "$out" | jq -r '(.open // []) | join(",")' 2>/dev/null || true)"
  [ -n "$o" ] && s="$s (open: $o)"
  s="${s//|/ }"                       # keep '|' out of the summary so field-split is safe
  results+=("$label|$v|${s}|$tfile")
  echo "  -> $label: $v"
done

# --- write findings to the TOP of TODO.md (idempotent, marker-delimited) ----------
if [ ${#results[@]} -gt 0 ]; then
  [ -f TODO.md ] || printf '# TODO\n\n' > TODO.md    # create it if the repo has none
  block="<!-- AUDIT:BEGIN -->
## Audit findings — $(date +%F) (TOP PRIORITY — from dev/audit-all.sh)
"
  issues=0
  for r in "${results[@]}"; do
    label="${r%%|*}"; r1="${r#*|}"; v="${r1%%|*}"; r2="${r1#*|}"; s="${r2%|*}"; t="${r2##*|}"
    if [ "$v" = AGREE ]; then
      block+="- [x] AUDIT ${label}: AGREE — no issues  (\`${t}\`)
"
    else
      issues=$((issues+1))
      block+="- [ ] **AUDIT ${label}: ${v}** — ${s:-see transcript}  (\`${t}\`)
"
    fi
  done
  block+="<!-- AUDIT:END -->
"
  tmp="$(mktemp "${TMPDIR:-/tmp}/audit-todo.XXXXXX")"   # bare mktemp ignores TMPDIR on macOS
  { printf '%s\n' "$block"; sed '/<!-- AUDIT:BEGIN -->/,/<!-- AUDIT:END -->/d' TODO.md; } > "$tmp" && mv "$tmp" TODO.md
  echo "→ TODO.md updated: $issues group(s) need attention (review & commit TODO.md yourself)"
fi

echo; echo "SUMMARY"
if [ ${#results[@]} -gt 0 ]; then
  for r in "${results[@]}"; do
    l="${r%%|*}"; rr="${r#*|}"; v="${rr%%|*}"
    printf "  %-18s %s\n" "$l" "$v"
  done
else
  echo "  (nothing audited)"
fi
