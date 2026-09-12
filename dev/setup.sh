#!/usr/bin/env bash
# setup.sh — one-time bootstrap for the two-model dev pipeline (run per clone).
#   ./dev/setup.sh
# Idempotent. Sets the git hooks path, makes scripts executable, verifies tooling.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

echo "→ point git at committed hooks (.githooks/)"
git config core.hooksPath .githooks
chmod +x .githooks/* dev/*.sh 2>/dev/null || true

echo "→ verify the co-author guard is active"
if printf 'test\n\nCo-Authored-By: Claude <noreply@anthropic.com>\n' | .githooks/commit-msg /dev/stdin >/dev/null 2>&1; then
  echo "  ✖ hook did NOT block an AI trailer — check .githooks/commit-msg" >&2; exit 1
else
  echo "  ✓ AI co-author trailers are blocked"
fi

echo "→ check models"
if command -v claude >/dev/null 2>&1; then echo "  ✓ claude  $(claude --version 2>/dev/null)"; else echo "  ✖ claude not on PATH"; fi
if command -v codex  >/dev/null 2>&1; then echo "  ✓ codex   $(codex --version 2>/dev/null)"; else
  echo "  ⚠ codex not on PATH — install:  npm i -g @openai/codex"
  echo "    (you have Codex in VS Code; the pipeline drives the CLI. Review step needs it.)"
fi

echo "→ done. Usage:"
echo "    dev/dual.sh review --staged      # independent genetic-theory review before commit"
echo "    dev/dual.sh loop \"<task>\"        # implement (one model) -> review (the other) -> fix"
