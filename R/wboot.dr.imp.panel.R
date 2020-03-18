# Bootrstapped Improved Doubly Robust Difference-in-Differences with panel Data
# 2 periods and 2 groups

wboot.dr.imp.panel <- function(nn, n, deltaY, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  ps.b <- pscore.cal(D, int.cov, b.weights, n)$pscore
  ps.b <- pmin(ps.b, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.reg.b <- wols.br.panel(deltaY, D, int.cov, ps.b, i.weights = b.weights)$out.reg
  # Compute AIPW estimator
  att.b <- aipw.did.panel(deltaY, D, ps.b, out.reg.b, b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}
