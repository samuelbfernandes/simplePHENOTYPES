# BUGS.md — simplePHENOTYPES v2

> Known bugs and their status. Use the `/bug-*` workflow
> (`/bug-create → /bug-analyze → /bug-fix → /bug-verify`, commit prefix `fix:`).
> Bug fixes are the priority gate before new-feature work (DECISION + TODO_newfeatures).

## Status legend
`[ ]` open · `[~]` in progress · `[x]` fixed · `[?]` needs investigation · `[-]` won't fix

## Open / known

- `[?]` **v1 `cor` argument is buggy** — Cholesky-correlated genetic effects are only
  accurate for trait 1; PVE is misreported after transformation. v2 does NOT carry `cor`
  forward; correlation is reimplemented as `rho_g` via PleioArch (DECISION-010). The buggy
  `cor` remains only inside the frozen `create_phenotypes()`; decide whether to fix in
  place or leave documented-as-is.
## Fixed

- `[x]` **LD direct QTN selection "diverges" from v1.3.0 — NOT a regression.** With
  `architecture = "LD", type_of_ld = "direct", seed = 200`, `create_phenotypes()` selects
  a different 3rd causal marker (`ss196453915`, chr 4) than the captured v1.3.0 reference
  (`ss196527787`, chr 3). Root cause: the reference was captured from CRAN 1.3.0 (tag at
  `01d6b0a`), which had a bug in the direct-LD `ld_min`/`ld_max` candidate-acceptance loop.
  That bug was fixed **after** 1.3.0 by commits `68a227e` and `b95529a` (both "fixed bug on
  ld_min/ld_max"), which the frozen-legacy engine already includes — hence the corrected
  3rd-marker pick. Indirect (untouched path) still matches 1.3.0 exactly, which is why only
  direct diverged. Resolution (DECISION-009): the `ld_direct.rds` reference is **re-blessed**
  from the current frozen-legacy output (`capture_references.R` 5d now captures it from the
  current package, not 1.3.0); `test-v130-parity.R` test 4 is un-skipped and passes.

- `[x]` **Test-artifact pollution** — stray `tests/testthat/file*.gds` files (written by
  `create_phenotypes()` GDS conversion into the testthat working dir) and untracked
  `tests/test_numeric.txt` removed. `tests/testthat/setup.R` now auto-cleans the `.gds`
  files via a `teardown_env()` deferral after each run; `.Rbuildignore` guards them from
  the CRAN bundle (`.gitignore` already excluded them from commits).

- `[x]` Non-biallelic SNPs now set to NA; clearer error for `method = "reference"` without
  `ref_allele` (commit 4069b47).
- `[x]` CRAN check issues resolved — 0 errors / 0 warnings / 0 notes (commit ac63155).

> Add new bugs above as `/bug-create` generates them.
