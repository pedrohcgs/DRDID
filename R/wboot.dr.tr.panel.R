# Bootrstapped "Traditional" Doubly Robust Difference-in-Differences with panel data
# 2 periods and 2 groups

wboot.dr.tr.panel <- function(nn, n, deltaY, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- as.vector(ps.b)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  reg.coeff.b <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
                                     subset = D==0,
                                     weights = b.weights))
  out.reg.b <- as.vector(tcrossprod(reg.coeff.b, int.cov))
  # Compute AIPW estimator
  att.b <- aipw.did.panel(deltaY, D, ps.b, out.reg.b, b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}
