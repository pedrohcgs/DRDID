# DRDID 1.2.4
  * Speed: the 2x2 estimators are about 1.5x faster, with estimates, standard errors, and influence functions unchanged up to floating-point precision. The propensity-score and outcome regressions now use `fastglm`'s low-level entry point (`fastglmPure`), which skips the per-call wrapper overhead, and the influence-function computations use `crossprod()` in place of `colMeans()` / `t()`. The weighted bootstrap benefits proportionally since it refits on every iteration.

  * The propensity-score Hessian is now checked for singularity with `rcond()` before inversion, matching the long-standing check on the outcome-regression design. A near-singular propensity-score design previously yielded a silently incorrect standard error; it now stops with an informative message.

  * `std_ipw_did_rc` now warns when the propensity-score estimation does not converge and stops when the estimated coefficients are `NA`, matching the other propensity-score estimators (these guards were previously missing, so a rank-deficient design returned a corrupted standard error silently).

  * The repeated cross-section pre-processing now warns when collinear covariates are dropped, matching the panel pre-processing (previously the repeated cross-section path dropped them without notice).

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


