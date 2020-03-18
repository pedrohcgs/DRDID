# Bootrstapped standardized IPW Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_std_ipw_rc <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  # Compute  standardized IPW estimator
  w.treat.pre.b <- b.weights * D * (1 - post)
  w.treat.post.b <- b.weights * D * post
  w.cont.pre.b <- b.weights * ps.b * (1 - D) * (1 - post)/ (1 - ps.b)
  w.cont.post.b <- b.weights * ps.b * (1 - D) * post/ (1 - ps.b)

  ipw.1 <- mean(w.treat.post.b * y) / mean(w.treat.post.b) -
    mean(w.treat.pre.b * y) / mean(w.treat.pre.b)
  ipw.0 <- mean(w.cont.post.b * y) / mean(w.cont.post.b) -
    mean(w.cont.pre.b * y) / mean(w.cont.pre.b)

  att.b <- ipw.1 - ipw.0
  #-----------------------------------------------------------------------------
  return(att.b)
}
