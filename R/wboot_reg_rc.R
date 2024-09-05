# Bootrstapped Regression-based Robust Difference-in-Differences with Repeated Cross Section Data
# 2 periods and 2 groups

wboot_reg_rc <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Outcome regression for the control group
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


  out.reg.pre.b <- as.vector(tcrossprod(reg.coeff.pre.b, int.cov))
  out.reg.post.b <- as.vector(tcrossprod(reg.coeff.post.b, int.cov))
  # Compute OR estimator
  att.b <- mean( b.weights * D * post * y)/mean(b.weights * D * post) -
    mean( b.weights * D * (1 - post) * y)/mean(b.weights * D * (1 - post)) -
    mean(b.weights * D * (out.reg.post.b - out.reg.pre.b))/mean(b.weights * D)
  #-----------------------------------------------------------------------------
  return(att.b)
}
