# Bootstrapped IPW Difference-in-Differences with Repeated Cross Section
# 2 periods and 2 groups

wboot_ipw_rc <- function(nn, n, y, post, D, int.cov, i.weights,
                         trim.level = 0.995){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  # ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- suppressWarnings(fastglm::fastglm(x = int.cov,
                                            y = D,
                                            family = stats::binomial(),
                                            weights =  b.weights,
                                            intercept = FALSE,
                                            method = 3)$fitted.values)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  trim.ps <- (ps.b < 1.01)
  trim.ps[D==0] <- (ps.b[D==0] < trim.level)
  # Lambda estimation
  lambda.b <- mean(trim.ps * b.weights * post)
  # Compute IPW estimator
  att.b <- mean(b.weights *trim.ps*
                  (D - ps.b * (1 - D)/(1 - ps.b)) * ((post - lambda.b)/(lambda.b * (1 - lambda.b))) * y) / mean(b.weights * D)
  #-----------------------------------------------------------------------------
  return(att.b)
}
