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
# per-file, coherent context. Selection (select/ocs/usefulness) is OMITTED (audit it
# separately with dev/dual.sh review). Designer files are OMITTED (migrating to BD).
# Each group's transcript is saved under dev/.audit/transcripts/.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export TMPDIR="${PIPELINE_TMPDIR:-$DIR/../.tmp}"; mkdir -p "$TMPDIR"
cd "$(git rev-parse --show-toplevel)"
. "$DIR/lib/verdict.sh"; . "$DIR/lib/audit.sh"
TS="$(date +%Y%m%d-%H%M%S)"; audit_init; ADIR="$(audit_dir)/transcripts"

# tier|label|paths (most theory-critical first)
GROUPS=(
"high|grammar|R/grammar_simulate_phenotype.R R/grammar_layers.R R/grammar_realize.R R/grammar_complex.R"
"high|effects-arch|R/effects_pleioarch.R R/effects_series.R R/arch_independent.R R/arch_ld.R R/ld_methods.R"
"high|crossing-schemes|R/cross_population.R R/cross_mating.R R/cross_map.R R/schemes.R"
"high|rust-core|src/rust/src/meiosis.rs src/rust/src/genome.rs src/rust/src/numeric.rs src/rust/src/lib.rs R/extendr-wrappers.R R/as_numeric.R R/numericalization.R R/table_to_numeric.R"
"med|io-formats|R/handling_input_formats.R R/file_loader.R R/detect_format.R R/format_conversion.R R/io_write.R R/filter_geno.R"
"legacy|legacy-core|R/legacy_create_phenotypes.R R/legacy_Phenotypes.R R/legacy_check_in.R R/legacy_constraint.R"
"legacy|legacy-qtn|R/legacy_QTN_linkage.R R/legacy_QTN_pleiotropic.R R/legacy_QTN_partially_pleiotropic.R R/legacy_qtn_from_user.R R/legacy_genetic_effect.R R/legacy_vQTL.R"
"legacy|legacy-baseline|R/legacy_Base_line_multi_traits.R R/legacy_Base_line_single_trait.R R/legacy_Genotypes.R R/legacy_make_pd.R"
)

case "${1:-live}" in
  high)    TIERS="high";;
  live|"") TIERS="high med";;
  all)     TIERS="high med legacy";;
  legacy)  TIERS="legacy";;
  *) echo "usage: dev/audit-all.sh [high|live|all|legacy]"; exit 2;;
esac
echo "→ audit-all tiers: $TIERS   reviewer: ${REVIEWER:-codex}"

results=()   # "label|verdict|summary|transcript"
for g in "${GROUPS[@]}"; do
  IFS='|' read -r tier label paths <<< "$g"
  case " $TIERS " in *" $tier "*) : ;; *) continue;; esac
  existing=(); for p in $paths; do [ -e "$p" ] && existing+=("$p"); done
  [ ${#existing[@]} -eq 0 ] && { echo "○ $label: no files present — skip"; continue; }
  echo "=== auditing: $label (${#existing[@]} files) ==="
  tfile="$ADIR/audit-$label-$TS.md"
  out="$(dev/dual.sh review "${existing[@]}" 2>>"$tfile.err" | tee "$tfile")" || true
  v="$(parse_verdict "$out")"
  s="$(verdict_json "$out" | jq -r '.summary // ""' 2>/dev/null || true)"
  o="$(verdict_json "$out" | jq -r '(.open // []) | join(",")' 2>/dev/null || true)"
  [ -n "$o" ] && s="$s (open: $o)"
  results+=("$label|$v|$s|$tfile")
  echo "  -> $label: $v"
done

# --- write findings to the TOP of TODO.md (idempotent, marker-delimited) ----------
if [ ${#results[@]} -gt 0 ] && [ -f TODO.md ]; then
  block="<!-- AUDIT:BEGIN -->
## Audit findings — $(date +%F) (TOP PRIORITY — from dev/audit-all.sh)
"
  issues=0
  for r in "${results[@]}"; do
    IFS='|' read -r label v s t <<< "$r"
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
  tmp="$(mktemp)"
  { printf '%s\n' "$block"; sed '/<!-- AUDIT:BEGIN -->/,/<!-- AUDIT:END -->/d' TODO.md; } > "$tmp" && mv "$tmp" TODO.md
  echo "→ TODO.md updated: $issues group(s) need attention (review & commit TODO.md yourself)"
fi

echo; echo "SUMMARY"
for r in "${results[@]}"; do IFS='|' read -r l v s t <<< "$r"; printf "  %-18s %s\n" "$l" "$v"; done
