# Bootrstapped "Traditional" Doubly Robust Difference-in-Differences
# 2 periods and 2 groups

wboot_drdid_rc <- function(nn, n, y, post, D, int.cov, i.weights,
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
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  # reg.coeff.pre.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                          subset = ((D==0) & (post==0)),
  #                                          weights = b.weights))
  control_pre <- (D == 0) & (post == 0)
  reg.coeff.pre.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[control_pre, , drop = FALSE],
                                  y = y[control_pre],
                                  weights = b.weights[control_pre],
                                  family = gaussian(link = "identity")
  ))
  out.y.cont.pre.b <-   as.vector(tcrossprod(reg.coeff.pre.b, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  # reg.coeff.post.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                           subset = ((D==0) & (post==1)),
  #                                           weights = b.weights))
  control_post <- (D == 0) & (post == 1)
  reg.coeff.post.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[control_post, , drop = FALSE],
                                  y = y[control_post],
                                  weights = b.weights[control_post],
                                  family = gaussian(link = "identity")
  ))
  out.y.cont.post.b <-   as.vector(tcrossprod(reg.coeff.post.b, int.cov))


  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  # reg.treat.coeff.pre.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                              subset = ((D==1) & (post==0)),
  #                                              weights = b.weights))
  treat_pre <- (D == 1) & (post == 0)
  reg.treat.coeff.pre.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[treat_pre, , drop = FALSE],
                                  y = y[treat_pre],
                                  weights = b.weights[treat_pre],
                                  family = gaussian(link = "identity")
  ))
  out.y.treat.pre.b <-   as.vector(tcrossprod(reg.treat.coeff.pre.b, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  # reg.treat.coeff.post.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                               subset = ((D==1) & (post==1)),
  #                                               weights = b.weights))
  treat_post <- (D == 1) & (post == 1)
  reg.treat.coeff.post.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[treat_post, , drop = FALSE],
                                  y = y[treat_post],
                                  weights = b.weights[treat_post],
                                  family = gaussian(link = "identity")
  ))
  out.y.treat.post.b <-   as.vector(tcrossprod(reg.treat.coeff.post.b, int.cov))



  # Compute AIPW estimator
  att.b <- aipw_did_rc(y, post, D, ps.b,
                       out.y.treat.post.b, out.y.treat.pre.b,
                       out.y.cont.post.b, out.y.cont.pre.b,
                       b.weights,
                       trim.ps)
  #-----------------------------------------------------------------------------
  return(att.b)
}
