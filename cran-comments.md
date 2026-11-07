## Submission summary

DRDID 1.4.1 is a test-only compatibility update for the upcoming fastglm 0.1.3,
following a report from the fastglm maintainer (Jared Huling).

fastglm 0.1.3 reports aliased coefficients of rank-deficient designs as `NA`
instead of zero. DRDID still rejects these designs, but with a different error
message, so five regression-test expectations matching the old message would
fail. The tests now accept either diagnostic.

Compared with the DRDID 1.4.0 tarball on CRAN, only `DESCRIPTION`, `NEWS.md`,
and `tests/testthat/test_regression_lock.R` differ. The R code, C++ sources,
documentation, and vignettes are unchanged, as are estimator behavior and
results.

We apologize for the short interval since 1.4.0 (published 2026-10-22). This
update keeps DRDID's checks passing once fastglm 0.1.3 reaches CRAN.

## Test environments

* Local: macOS 26.5.1 (aarch64-apple-darwin23), R 4.6.0, 2026-11-07, checked
  separately against fastglm 0.1.2 (current CRAN) and fastglm 0.1.3 (upstream
  commit 889c7c65eabf58f0abef55fedd71dacae74c58da).
* GitHub Actions: Ubuntu 24.04.5, R-devel (2026-11-06 r90845), R 4.6.2,
  and R 4.5.3.
* GitHub Actions: macOS 26.7.2 (arm64), R 4.6.2.
* GitHub Actions: Windows Server 2022 (x64), R 4.6.2 UCRT.

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes with each fastglm
version. All 4,413 testthat assertions passed in both, with no warnings or
skips. With fastglm 0.1.3, the 1.4.0 test file reproduces the five failures.

The GitHub Actions checks used the release candidate's source commit
f772e1d9 with fastglm 0.1.2. Each platform reported Status: OK
(0 errors | 0 warnings | 0 notes), and all 4,413 testthat assertions passed
on each platform with no warnings or skips.

https://github.com/pedrohcgs/DRDID/actions/runs/42406215997

## Reverse dependencies

Both CRAN reverse dependencies were checked from their current CRAN sources
against DRDID 1.4.1 and fastglm 0.1.2 in an isolated library on the local
environment above, with `--as-cran --no-manual` and `NOT_CRAN=false`.

* `did` 2.5.1: 0 errors | 1 warning | 0 notes. Examples, tests, and rebuilt
  vignettes passed, with 1,715 test assertions and 23 upstream skips.
* `ptetools` 1.0.1: 0 errors | 1 warning | 0 notes. Examples and tests passed,
  with 8 assertions and 3 upstream skips.

In both cases the warning is solely CRAN incoming feasibility for the
already-published version and its older source build timestamp; it is
unrelated to DRDID compatibility. Since the DRDID package code is unchanged
from 1.4.0, reverse dependencies are unaffected.
