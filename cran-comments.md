## Submission

This is a major update from the CRAN version 1.3.0 to the 2.x line. It adds a
composable phenotype-simulation grammar (`simulate_phenotype()` with `additive()`,
`dominance()`, `epistasis()`, `vqtl()` layers), multi-generation crossing and
selection tools, and a small Rust core (via extendr) for genotype conversion and
meiosis. The v1 interface `create_phenotypes()` is kept unchanged, so existing
code keeps working. There are no reverse dependencies on CRAN.

### Rust code

* All Rust dependencies are vendored in `src/rust/vendor.tar.xz` (xz-compressed,
  from `cargo vendor`), so the package builds offline; nothing is downloaded at
  install time.
* `cargo build` runs with `-j 2 --offline` on CRAN (`tools/config.R`), and the
  `cargo`/`rustc` versions are reported during configure (`tools/msrv.R`).
* `SystemRequirements: Cargo (Rust's package manager), rustc >= 1.65.0, xz`.
* The authors and licences of the vendored crates are declared in the DESCRIPTION
  `Copyright` field, which points to `inst/COPYRIGHTS`.

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

* macOS Tahoe 26.7 (aarch64), R 4.6.1, `R CMD check --as-cran`
* win-builder (devel and release) — to be run before submission

## R CMD check results

0 ERRORs. Locally the remaining WARNING and NOTEs come from the check machine,
not the package:

* WARNING "A complete check needs the 'checkbashisms' script" — the helper is not
  installed on the local macOS machine.
* NOTE "checking HTML version of manual" — the local HTML Tidy is too old to
  validate R's generated HTML.

The expected NOTE on CRAN is:

* "New maintainer" — the maintainer address change described above.

URLs and DOIs were checked with `urlchecker::url_check()` and
`tools:::check_doi_db()`: all resolve.
