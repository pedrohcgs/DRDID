## Submission summary

DRDID 1.4.0 is a bug-fix update for inference, panel preprocessing, and
weighted estimation.

* Correct the sign of the pre-period comparison-group outcome-regression term in
  the influence functions of `drdid_rc()` and `drdid_rc1()`.

* Preserve individual pre/post pairing when panel outcomes, covariates, or
  sampling weights are missing, and preserve user column names during
  preprocessing in formula-based estimators.

* Use individual-level influence functions and analytical inference in
  `twfe_did_panel()`, reject incomplete or mismatched panel inputs, and preserve
  covariate dimensions when removing a supplied intercept in both TWFE estimators.

* Account for observation-weight normalization in `ipw_did_rc()` inference and
  make its weighted bootstrap refit the same statistic as the point estimator.

* Require numerically verified original IPT fits in improved DR estimators and
  every weighted-bootstrap replicate. Failed IPT fits now stop rather than
  silently using ordinary logit inference. Stabilize weighted IPT arithmetic,
  ignore zero-weight contributions, and consistently weight the fallback
  objective, gradient, and Hessian.

* Make the improved panel weighted bootstrap refit the point estimator's
  statistic with trimming and supplied designs without an intercept.

* Validate sampling weights, binary treatment coding, bootstrap repetition
  counts, and regression cell sizes. Handle exactly degenerate bootstrap draws
  explicitly and reject nonfinite or otherwise invalid bootstrap scales.

* Correct the minimum required version of BMisc to 1.4.8.

These corrections can change standard errors and confidence intervals. Panel
preprocessing and weighted IPT fallback corrections can also change point
estimates for affected inputs.

See NEWS.md for the full list of changes.

## Test environments

* Local: macOS 26.5.1 (aarch64-apple-darwin23), R 4.6.0, 2026-10-22.
* GitHub Actions: Ubuntu 24.04.5, R-devel (2026-10-21 r90732), R 4.6.1,
  and R 4.5.3.
* GitHub Actions: macOS 26.6.2 (arm64), R 4.6.1.
* GitHub Actions: Windows Server 2022 (x64), R 4.6.1 UCRT.

## R CMD check results

Local `R CMD check --as-cran`: 0 errors | 0 warnings | 0 notes.
All 4,413 testthat assertions passed. A separate spelling check found no errors.

The GitHub Actions checks used the release candidate's source commit
6dbc303cc16512c1578646836dbeb49737b5d4fc, before updating the version/date,
NEWS heading, and spelling wordlist. Each platform reported 0 errors and
0 warnings, with one spelling NOTE for technical terms subsequently added to
the wordlist. All 4,413 testthat assertions passed on each platform.

https://github.com/pedrohcgs/DRDID/actions/runs/39825405732

## Reverse dependencies

Both CRAN reverse dependencies were checked from their current CRAN sources
against DRDID 1.4.0 in an isolated library on the local environment above.

* `ptetools` 1.0.1: 0 errors | 0 warnings | 0 notes. Tests passed with
  8 assertions and 3 upstream skips.
* `did` 2.5.1: 0 errors | 1 warning | 0 notes with `--as-cran --no-manual`
  and `NOT_CRAN=false`. Examples, tests, and rebuilt vignettes passed, with
  1,715 test assertions and 23 upstream skips. The warning is solely CRAN
  incoming feasibility for the already-published `did` version and its older
  source build timestamp; it is unrelated to DRDID compatibility.
