#!/usr/bin/env python3
"""Extract a structured verdict from a model's stdout.

Models are asked to end with a JSON object {"verdict":"AGREE|BLOCK","open":[...],
"confidence":0-1,...}. Free-text prose is unreliable to scrape, so we parse the JSON.

Usage:
  parse_verdict.py            < model_output      # prints: AGREE | BLOCK | UNKNOWN
  parse_verdict.py --json     < model_output      # prints the normalized JSON (or {})
  parse_verdict.py --field open < model_output     # prints a field (open as CSV)
"""
import sys, json, re

def extract(text):
    best = None
    # flat JSON objects containing "verdict" (open:[...] has brackets, not braces)
    for m in re.finditer(r'\{[^{}]*"verdict"[^{}]*\}', text, re.S | re.I):
        try:
            best = json.loads(m.group(0))
        except Exception:
            pass
    if best is None:
        mm = re.findall(r'verdict"?\s*[:=]\s*"?(AGREE|BLOCK)', text, re.I)
        if mm:
            best = {"verdict": mm[-1].upper(), "open": []}
    return best or {}

def main():
    args = sys.argv[1:]
    obj = extract(sys.stdin.read())
    if "--json" in args:
        print(json.dumps(obj)); return
    if "--field" in args:
        f = args[args.index("--field") + 1]
        v = obj.get(f, "")
        print(",".join(map(str, v)) if isinstance(v, list) else str(v)); return
    v = str(obj.get("verdict", "")).upper()
    print(v if v in ("AGREE", "BLOCK") else "UNKNOWN")

if __name__ == "__main__":
    main()
