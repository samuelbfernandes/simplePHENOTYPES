# evals/ — meta-testing the review harness

You cannot trust a reviewer you have never tested. This is a **golden set of seeded
theory bugs**: each entry in `mutations.json` injects one known quantitative-genetics
error into the real source. A good reviewer must **catch** it (recall) and must **not**
fabricate objections on clean code (precision).

```bash
evals/run.sh --check   # (default) verify each mutation still applies & reverts. No model,
                       # no tokens — run this in CI so the golden set can't silently rot.
evals/run.sh --eval    # run the reviewer on each seeded bug + a clean control; score
                       # recall / localization / precision. Calls a model (costs tokens).
REVIEWER=claude evals/run.sh --eval    # swap the reviewer model
```

`--eval` seeds each bug inside a throwaway **git worktree** (never your working tree),
runs the reviewer, parses its **JSON verdict**, and records each result to `dev/.audit/`.

## The seeded bugs (all map to `docs/THEORY_REVIEW.md`)

| id | file | breaks | rubric |
|----|------|--------|--------|
| O1-vanraden-denominator | R/ocs.R | drops the 2 in VanRaden `2·Σp(1−p)` | O1 |
| O2-coancestry-half | R/ocs.R | drops the ½ in group coancestry `½c'Gc` | O2 |
| S1-intensity-divide-p | R/select.R | drops `/p` in `i(p)=φ(Φ⁻¹(1−p))/p` | S1 |
| S2-smith-hazel-swap | R/select.R | swaps P,G in `b=P⁻¹Ga` | S2 |
| U1-usefulness-intensity | R/usefulness.R | drops `i` in `U=μ+iσ` | U1 |

## Adding a bug
Append to `mutations.json`: a unique `find` substring from the source, the wrong
`replace`, the `rubric` id it should trip, and a `bug` note. Run `evals/run.sh --check`
to confirm it applies. Keep `find` distinctive so it matches exactly one line.

> If `--check` reports "find string not present", the source moved — update the `find`
> string. That failure in CI is the golden set doing its job.
