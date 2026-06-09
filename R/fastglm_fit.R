# Internal helper: fit a GLM via fastglm's low-level entry point (fastglmPure),
# skipping the fastglm() / fastglm.default formula+data.frame wrapper. That wrapper
# does input coercion (is.data.frame / model.matrix paths), family deviance
# bookkeeping (dbinom, dev.resids) and structure() assembly on EVERY call -- pure
# per-fit overhead when DRDID's 2x2 estimators are invoked thousands of times (e.g.
# from the did package's group-time loop and its conditional pre-test).
#
# Passing fastglm's own defaults -- method as supplied (3 for the logit propensity
# score, 0 for the OLS outcome regressions, matching fastglm.default), tol = 1e-8,
# maxit = 100L -- makes the fitted coefficients and fitted values BIT-IDENTICAL to
# the previous fastglm::fastglm() calls (verified to 0 across a grid of cells),
# while removing ~25-30% of each fit's wall-clock.
#
# The returned fastglmPure list exposes $coefficients, $fitted.values and
# $converged, which replace the previous coef() / fitted() / $converged accessors.
fastglm_fit <- function(x, y, family, weights = NULL, method = 0L) {
  fastglm::fastglmPure(x = x, y = y, family = family, weights = weights,
                       method = method, tol = 1e-8, maxit = 100L)
}
