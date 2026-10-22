# DRDID 1.4.0

  * Panel preprocessing now excludes whole units with missing covariates or
    sampling weights before checking time invariance. Retained units remain
    matched by ID, and observed time-varying covariates or weights still fail.

  * Formula-based estimators preserve user columns named `y`, `D`, `post`, or
    `w` when extracting outcomes, weights, and covariates. Renaming a variable
    no longer changes its meaning through internal preprocessing columns.

  * Improved DR estimators now require a numerically verified solution of the
    original exponential-IPT equations. Logit MLE is used only for initialization;
    failed IPT fits stop and suggest explicitly rerunning the traditional DR
    estimator. The secondary IPT candidate is checked against the original
    equations even when its continued objective reports convergence. The returned
    `ps.flag` identifies the verified fitting route (0 or 1); value 2 is removed.
    Fitting normalizes the common scale of weights to avoid scale-dependent
    optimizer stopping. Numerically ill-conditioned designs are rejected.
    Every weighted-bootstrap replicate must pass the same check, and failures
    stop with the replicate number without skipping or redrawing the fit.

  * The weighted bootstrap in `drdid_imp_panel()` now uses the point estimator's
    treated-weight denominator for its combined residual contributions. This
    makes each draw refit the same implemented statistic, including active
    trimming and supplied designs without an intercept. Weighted-bootstrap
    draws, standard errors, and intervals can change;
    point estimates, analytic inference, and multiplier inference are unchanged.

  * IPT criteria now evaluate weighted nonlinear contributions without first
    creating avoidable overflowing or underflowing intermediate exponentials.
    Weighted values, gradients, and Hessians retain representable terms and
    reject unresolved opposing overflow terms. The original criteria and
    continuation rule are preserved. Fitting reuses fixed-data checks and
    avoids unnecessary exceptional arithmetic while retaining these safeguards.

  * Zero-weight rows no longer enter IPT objective, gradient, or Hessian
    arithmetic. Extreme finite covariates in those rows cannot create NaNs
    through zero times an overflowing intermediate value. The original
    objective divisor, model dimension, and secondary continuation threshold
    are preserved.

  * The IPT propensity-score fallback used by the improved DR estimators now
    weights its full objective consistently with its gradient and Hessian.
    This corrects optimization with unequal observation or bootstrap weights
    when the primary `trust` fit does not converge.

  * Sampling weights are checked before normalization; negative, nonfinite,
    all-zero, or incorrectly sized low-level weights now give a clear error.
    Legitimate zero weights and existing missing-data exclusion are preserved.

  * The weighted bootstrap in `ipw_did_rc()` now refits the same statistic as
    the point estimator, including normalization of the combined observation
    and bootstrap weights and the full-sample pre/post time shares. Trimming
    affects comparison contributions, as in the point estimator. Weighted
    bootstrap draws, standard errors, and intervals can change; the original
    ATT, analytic inference, and multiplier bootstrap are unchanged.

  * `ipw_did_rc()` now includes the effect of sample-mean normalization of
    observation weights in its influence function and the resulting analytic
    and multiplier-bootstrap inference under iid sampling. The ATT point
    estimate is unchanged.

  * Bootstrap inference now returns a zero standard error and point interval when
    every centered bootstrap draw is exactly zero. Nonfinite draws or a zero IQR
    with nonzero draws stop with a clear error; no draws are removed or redrawn,
    and positive-scale inference retains its existing arithmetic.

  * `twfe_did_panel()` and `twfe_did_rc()` now preserve matrix dimensions when
    removing a supplied intercept. An intercept plus one covariate now fits
    the same regression as supplying that covariate as a numeric vector or a
    matrix with one column, instead of failing.
    Intercept-only inputs retain their existing no-covariate behavior.

  * `twfe_did_panel()` now accounts for dependence between each individual's
    pre/post observations in its analytical standard error. With covariates,
    the returned influence function now has one value per individual instead
    of one per period row. This matches individual-clustered sandwich inference
    without small-sample corrections and the other panel estimators' convention.
    The point estimate and both bootstrap procedures are preserved, up to
    floating-point precision. The no-covariate inference path is unchanged.
    Direct calls now reject missing or mismatched individual inputs before
    fitting, preventing row omission from breaking pre/post alignment.

  * Bootstrap estimation now rejects `nboot` values that are not finite scalar
    integers of at least 2 before fitting or drawing bootstrap weights. The
    low-level `NULL` default remains 999; `boot = FALSE` continues to ignore it.

  * Panel preprocessing for `drdid()`, `ordid()`, and `ipwdid()` now removes
    both observations of a unit when either outcome is missing,
    matches pre/post outcomes by individual ID, and checks group sizes after
    exclusions. Previously, dropping individual rows with missing outcomes could
    silently pair different individuals and produce incorrect point estimates
    and standard errors. Warnings now report the number of affected panel units.

  * `drdid_rc()` and `drdid_rc1()`: the estimation effect of the pre-period
    outcome regression for the comparison group now enters the influence function
    with the sign given in Appendix B of Sant'Anna and Zhao (2020). It previously
    entered with the opposite sign. Point estimates are unchanged. Analytic and
    multiplier-bootstrap standard errors change. The two versions agree
    asymptotically when the propensity score model is correctly specified, but
    not when only the outcome regressions are.

  * Treatment indicators must use 0 for comparison and 1 for treated units.
    Invalid coding is rejected explicitly rather than silently changing weights.

  * The internal GLM fitter now rejects designs with fewer rows than columns
    before entering `fastglm`. Such designs cannot have full column rank.
    Callers receive an ordinary R error when a regression cell is undersized;
    supported fits retain the same calculations.

# DRDID 1.3.0

This release substantially speeds up the 2x2 estimators and tightens the handling
of ill-conditioned designs. Point estimates, standard errors, and influence
functions are unchanged up to floating-point precision (~1e-14).

## Performance

  * The propensity-score and outcome regressions now use `fastglm`'s low-level
    entry point (`fastglmPure`) in place of the `fastglm()` wrapper, skipping the
    per-call input-coercion and family/deviance bookkeeping that was paid on every
    fit. This is applied throughout: the point estimators, every weighted-bootstrap
    helper (which refits on each of `nboot` iterations), and the IPT/calibration
    propensity initializer.

  * The influence-function computations use `crossprod()` and BLAS dot products in
    place of `colMeans(. * int.cov)` and `t(int.cov) %*% (W * int.cov)`.

  * Together these make the 2x2 estimators roughly 1.5x faster, with the weighted
    bootstrap benefiting proportionally.

## Robustness

  * The propensity-score Hessian is now checked for singularity with `rcond()`
    before inversion, matching the long-standing check on the outcome-regression
    design. A near-singular propensity-score design previously yielded a silently
    incorrect standard error; it now stops with an informative message.

  * `std_ipw_did_rc` now warns when the propensity-score estimation does not
    converge and stops when the estimated coefficients are `NA`, matching the other
    propensity-score estimators (these guards were previously missing, so a
    rank-deficient design returned a corrupted standard error silently).

  * The repeated cross-section pre-processing now warns when collinear covariates
    are dropped, matching the panel pre-processing (previously the repeated
    cross-section path dropped them without notice).

## Internal

  * Added a numeric regression-lock test suite that pins the ATT and standard error
    of every exported estimator, asserts the influence-function invariants, and
    checks that singular designs and negative weights are rejected.

# DRDID 1.2.3
  * Fix typo on non-stabilized IPW with trimming
  * Unify degree of freedom adjustments in analytical std errors
# DRDID 1.2.2
  * Add trimming argument to avoid severe overlap problems. The default is to trim the propensity score in the comparison group that is above 0.995.
  
# DRDID 1.2.1
  * Fix typo on returning influence functions for TWFE regressions.

# DRDID 1.2.0
  * Improve code to avoid redundant data checks
  
  * Use fastglm instead of parglm for improved speed

  
# DRDID 1.1.0
  * Restore solve as default to invert matrix, as it is faster than qr.solve for small matrices.
  
  * Improve error handling for non-invertible matrices.
  
  * Changing estimation methods for fastglm and parglm (in place of lm and glm).
  
  * Do not let the estimated propensity score be above 1 - 1e-6 (instead of 1 - e-16).

# DRDID 1.0.7
  * Speed up data processing using Rcpp
  
  * The weights are now enforced to be normalized and have mean 1 across all observations.

# DRDID 1.0.6
  * Use qr.solve as default (instead of solve)
  
  * Drop collinear variables in pre_process_drdid.R (useful in drdid command but not other commands)
  
  * Add compatibility with R 3.5
  
  * Improve invertibility of outcome regression design matrix
  
# DRDID 1.0.4
  * Fixed links
  
# DRDID 1.0.3
  * Add new flags for non-unique unit identifier
  
  * Better handle of factor variables as covariates

# DRDID 1.0.2
  * Fix issue with NA in covariates
  
# DRDID 1.0.1
  * Allows for treating covariates as factor and alike when computing DiD
  
  * Improve error and warning handling due to collinearity and convergence issues.
  
# DRDID 1.0.0
  * First official version of package, functions for computing a variety of difference-in-differences (DiD) estimators for the ATT. 
  
  * Documentation is improved compared to the devel version, including examples for every function now.
  
  * Created wrapper function `drdid`, `ordid` and `ipwdid` to implement doubly-robust, outcome regression and inverse probability weighted DID estimators.
  
  * Add dataset used in the empirical application of Sant'Anna and Zhao (2020).


