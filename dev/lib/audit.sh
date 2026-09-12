#!/usr/bin/env bash
# audit.sh — provenance for every model-assisted run. Sourced, not executed.
# Records who/what/when so a result is reproducible and attributable to a run, not a
# vibe. The audit dir is gitignored (local provenance, not shipped).

audit_dir()  { echo "$(git rev-parse --show-toplevel)/dev/.audit"; }
audit_init() { mkdir -p "$(audit_dir)/transcripts"; }

# Version of the rubric actually used (git blob hash of THEORY_REVIEW.md), so a verdict
# can be tied to the exact criteria that produced it.
audit_rubric_hash() {
  git hash-object docs/THEORY_REVIEW.md 2>/dev/null | cut -c1-12 || echo "none"
}

# audit_record <kind> <scope> <verdict> <defender> <skeptic> <transcript_path>
audit_record() {
  audit_init
  jq -cn \
    --arg ts        "$(date -u +%FT%TZ)" \
    --arg kind      "${1:-}" \
    --arg scope     "${2:-}" \
    --arg verdict   "${3:-}" \
    --arg defender  "${4:-}" \
    --arg skeptic   "${5:-}" \
    --arg transcript "${6:-}" \
    --arg rubric    "$(audit_rubric_hash)" \
    --arg claude    "$(command -v claude >/dev/null && claude --version 2>/dev/null || echo NA)" \
    --arg codex     "$(command -v codex  >/dev/null && codex  --version 2>/dev/null || echo NA)" \
    --arg commit    "$(git rev-parse --short HEAD 2>/dev/null || echo NA)" \
    --arg branch    "$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo NA)" \
    '{ts:$ts,kind:$kind,scope:$scope,verdict:$verdict,
      defender:$defender,skeptic:$skeptic,transcript:$transcript,
      rubric_hash:$rubric,models:{claude:$claude,codex:$codex},
      commit:$commit,branch:$branch}' \
    >> "$(audit_dir)/log.jsonl"
}

# audit_save_transcript <src> -> echoes the saved path (under dev/.audit/transcripts/)
audit_save_transcript() {
  audit_init
  local dest; dest="$(audit_dir)/transcripts/$(basename "$1")"
  cp "$1" "$dest" 2>/dev/null && echo "$dest" || echo "$1"
}
