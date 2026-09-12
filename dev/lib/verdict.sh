#!/usr/bin/env bash
# verdict.sh — structured-verdict contract shared by the review/debate/eval harnesses.
# A parseable JSON verdict beats regex-scraping prose. Sourced, not executed.

LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The instruction appended to every reviewer/skeptic prompt.
verdict_instruction() {
  cat <<'EOF'
You MUST end your reply with a single fenced JSON object on its own — nothing after it:

```json
{"verdict":"AGREE|BLOCK","open":["O1","O2"],"confidence":0.0,"summary":"one line"}
```

- "verdict": AGREE only if you are convinced the implementation is correct; BLOCK otherwise.
- "open": ids of the objections still unresolved (empty when AGREE).
- "confidence": 0..1, your calibrated confidence in the verdict.
Ground the verdict in executed evidence, not rhetoric. Do not invent citations/pages.
EOF
}

# parse_verdict <text>  -> AGREE | BLOCK | UNKNOWN
parse_verdict() { printf '%s' "$1" | python3 "$LIB_DIR/parse_verdict.py"; }
# verdict_json <text>   -> normalized JSON (or {})
verdict_json()  { printf '%s' "$1" | python3 "$LIB_DIR/parse_verdict.py" --json; }
