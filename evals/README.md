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

`--eval` seeds each bug inside a throwaway **git worktree** under `$TMPDIR` (never your
working tree: `run.sh` writes to no source file in either mode, and `--check` mutates
temporary copies). The current contents of each target file, including uncommitted or
untracked work, are copied into the worktree first. It then runs the reviewer there, parses
its **JSON verdict**, and records each result to `dev/.audit/`. Every target is snapshotted
first and verified byte-for-byte on exit (any mismatch is restored and reported); the
worktree is removed on exit (`KEEP=1` keeps it). A reviewer that exits non-zero counts as
an error (reported in the scorecard, non-zero exit), never a silent pass.

## The seeded bugs (all map to `docs/THEORY_REVIEW.md`)

| id | file | breaks | rubric |
|----|------|--------|--------|
| O1-vanraden-denominator | R/select_ocs.R | drops the 2 in VanRaden `2·Σp(1−p)` | O1 |
| O2-coancestry-half | R/select_ocs.R | drops the ½ in group coancestry `½c'Gc` | O2 |
| S1-intensity-divide-p | R/select_ind.R | drops `/p` in `i(p)=φ(Φ⁻¹(1−p))/p` | S1 |
| S2-smith-hazel-swap | R/select_ind.R | swaps P,G in `b=P⁻¹Ga` | S2 |
| U1-usefulness-intensity | R/select_usefulness.R | drops `i` in `U=μ+iσ` | U1 |

## Adding a bug
Append to `mutations.json`: a unique `find` substring from the source, the wrong
`replace`, the `rubric` id it should trip, and a `bug` note. Run `evals/run.sh --check`
to confirm it applies. Keep `find` distinctive so it matches exactly one line.

> If `--check` reports "find string not present", the source moved — update the `find`
> string. That failure in CI is the golden set doing its job.
