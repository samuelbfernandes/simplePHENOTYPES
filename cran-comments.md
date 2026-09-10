## Maintainer address change

The maintainer address changes in this submission, from `samuelf@illinois.edu`
to `fernandessb101@gmail.com`.

`samuelf@illinois.edu` was an institutional address at a former affiliation and
is no longer accessible, so the confirmation email CRAN normally sends to the
*old* address cannot be answered. `fernandessb101@gmail.com` is a personal
address that will remain reachable.

The package continues to be maintained by the same person (Samuel B. Fernandes),
under the same GitHub account that has held it since the first CRAN release
(https://github.com/samuelbfernandes/simplePHENOTYPES). I am happy to provide any
additional confirmation of identity you need — please contact me at the new
address.

## Test environments

* macOS (local), R 4.4.3, `R CMD check --as-cran`
* win-builder (devel and release)

## R CMD check results

No ERRORs. Locally, the only WARNING and NOTEs are environmental, not package
issues:

* WARNING "A complete check needs the 'checkbashisms' script" — the `checkbashisms`
  helper is not installed on the local macOS machine; there are no shell scripts
  in the package that would trip it.
* NOTE "checking HTML version of manual" — emitted by an older local `tidy` that
  does not recognise the HTML5 `<main>` element R now generates; not reproduced
  on CRAN's infrastructure.
* NOTE "checking for future file timestamps" — the local sandbox clock, unrelated
  to the package.

The expected NOTE on CRAN is:

* "New maintainer" — the maintainer address change described above.
