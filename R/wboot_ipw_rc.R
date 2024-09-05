# Bootstrapped IPW Difference-in-Differences with Repeated Cross Section
# 2 periods and 2 groups

wboot_ipw_rc <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  # ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- suppressWarnings(parglm::parglm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  # Lambda estimation
  lambda.b <- mean(b.weights * post)
  # Compute IPW estimator
  att.b <- mean(b.weights * (D - ps.b * (1 - D)/(1 - ps.b)) * ((post - lambda.b)/(lambda.b * (1 - lambda.b))) * y) / mean(b.weights * D)
  #-----------------------------------------------------------------------------
  return(att.b)
}
