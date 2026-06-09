## Submission summary

DRDID 1.3.0 is a performance and robustness update.

* The 2x2 estimators are roughly 1.5x faster (via `fastglm`'s low-level entry point
  `fastglmPure` and `crossprod()`-based influence-function computations), and the
  weighted bootstrap benefits proportionally. Point estimates, standard errors, and
  influence functions are unchanged up to floating-point precision (~1e-14).

* The propensity-score Hessian is now checked for singularity before inversion
  (matching the existing outcome-regression check); `std_ipw_did_rc` gained the
  convergence / NA-coefficient guards the other estimators already had; and the
  repeated cross-section pre-processing now warns on dropped collinear covariates
  (matching the panel path). These only add clear errors/warnings on ill-conditioned
  inputs that previously produced silently-incorrect output.

See NEWS.md for the full list of changes.

## Test environments

* Local: macOS 26.5 (aarch64-apple-darwin23), R 4.6.0
* win-builder: R-devel and R-release   <!-- to run: devtools::check_win_devel(); check_win_release() -->
* R-hub v2 (GitHub Actions), incl. the compiled-code sanitizer/valgrind runs
  <!-- to run: rhub::rhub_check() -->
* macOS builder   <!-- to run: devtools::check_mac_release() -->

## R CMD check results

`R CMD check --as-cran` on the local environment above: 0 errors | 0 warnings | 0
notes. (The local run reports one NOTE, "Skipping checking HTML validation: 'tidy'
doesn't look like recent enough HTML Tidy", which is a property of the local
machine, not the package.)

## Reverse dependencies

There are 2 CRAN reverse dependencies, `did` and `ptetools`. Both were checked
against this version and pass with 0 errors / 0 warnings / 0 notes. The numeric
change introduced here is below the precision of any downstream expectation, and
the new errors/warnings only fire on ill-conditioned designs that the reverse
dependencies do not exercise.
