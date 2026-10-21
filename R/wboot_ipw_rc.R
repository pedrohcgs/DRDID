# Bootstrapped IPW Difference-in-Differences with Repeated Cross Section
# 2 periods and 2 groups

wboot_ipw_rc <- function(nn, n, y, post, D, int.cov, i.weights,
                         trim.level = 0.995){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  b.weights <- b.weights / mean(b.weights)
  # Propensity score estimation
  # ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- suppressWarnings(fastglm_fit(int.cov, D, stats::binomial(), b.weights, method = 3)$fitted.values)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  trim.ps <- (ps.b < 1.01)
  trim.ps[D==0] <- (ps.b[D==0] < trim.level)
  # Match the full-sample treatment and time shares in ipw_did_rc.
  Pi.b <- mean(b.weights * D)
  lambda.b <- mean(b.weights * post)
  one.minus.lambda.b <- mean(b.weights * (1 - post))

  # Refit the same four components, trimming comparison observations only.
  w.treat.pre <- b.weights * D * (1 - post)
  w.treat.post <- b.weights * D * post
  w.cont.pre <- trim.ps * b.weights * ps.b * (1 - D) * (1 - post) / (1 - ps.b)
  w.cont.post <- trim.ps * b.weights * ps.b * (1 - D) * post / (1 - ps.b)

  att.b <- (mean(w.treat.post * y / (Pi.b * lambda.b)) -
              mean(w.treat.pre * y / (Pi.b * one.minus.lambda.b))) -
    (mean(w.cont.post * y / (Pi.b * lambda.b)) -
       mean(w.cont.pre * y / (Pi.b * one.minus.lambda.b)))
  #-----------------------------------------------------------------------------
  return(att.b)
}
