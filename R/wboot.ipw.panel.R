# Bootrstapped IPW Difference-in-Differences with Panel Data
# 2 periods and 2 groups

wboot.ipw.panel <- function(nn, n, deltaY, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  # Compute IPW estimator
  att.b <- mean(b.weights * (D - ps.b * (1 - D)/(1 - ps.b)) * deltaY) / mean(b.weights * D)

  #-----------------------------------------------------------------------------
  return(att.b)
}
