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
* win-builder: R-devel and R-release
* macOS builder (mac.R-project.org): R-release
* R-hub v2: linux, macos, macos-arm64, windows (R-devel), and the AddressSanitizer
  builds clang-asan and gcc-asan (the package contains compiled C++)

## R CMD check results

0 errors | 0 warnings | 0 notes, on every environment above (including both
AddressSanitizer builds).

## Reverse dependencies

There are 2 CRAN reverse dependencies, `did` and `ptetools`. Both were checked
against this version and pass with 0 errors / 0 warnings / 0 notes. The numeric
change introduced here is below the precision of any downstream expectation, and
the new errors/warnings only fire on ill-conditioned designs that the reverse
dependencies do not exercise.
