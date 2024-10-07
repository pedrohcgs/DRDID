# Bootrstapped standardized IPW Difference-in-Differences with Panel Data
# 2 periods and 2 groups

wboot.std.ipw.panel <- function(nn, n, deltaY, D, int.cov, i.weights){
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
  # Compute  standardized IPW estimator
  w.treat.b <- b.weights * D
  w.cont.b <- b.weights * (1 - D) * ps.b / (1 - ps.b)

  aipw.1 <- mean(w.treat.b * deltaY) / mean(w.treat.b)
  aipw.0 <- mean(w.cont.b * deltaY) / mean(w.cont.b)

  att.b <- aipw.1 - aipw.0
  #-----------------------------------------------------------------------------
  return(att.b)
}
