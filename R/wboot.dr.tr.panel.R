# Bootrstapped "Traditional" Doubly Robust Difference-in-Differences with panel data
# 2 periods and 2 groups

wboot.dr.tr.panel <- function(nn, n, deltaY, D, int.cov, i.weights,
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
  ps.b <- as.vector(ps.b)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  trim.ps <- (ps.b < 1.01)
  trim.ps[D==0] <- (ps.b[D==0] < trim.level)

  #Compute the Outcome regression for the control group
  # reg.coeff.b <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
  #                                    subset = D==0,
  #                                    weights = b.weights))
  control_filter <- (D == 0)
  reg.coeff.b <- stats::coef(fastglm::fastglm(
                              x = int.cov[control_filter, , drop = FALSE],
                              y = deltaY[control_filter],
                              weights = b.weights[control_filter],
                              family = gaussian(link = "identity")
  ))
  out.reg.b <- as.vector(tcrossprod(reg.coeff.b, int.cov))
  # Compute AIPW estimator
  att.b <- aipw.did.panel(deltaY, D, ps.b, out.reg.b, b.weights, trim.ps)
  #-----------------------------------------------------------------------------
  return(att.b)
}
