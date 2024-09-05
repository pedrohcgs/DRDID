# Bootrstapped "Improved" Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_drdid_rc1 <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  # Propensity score estimation
  # ps.b <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- suppressWarnings(parglm::parglm(D ~ -1 + int.cov, family = "binomial", weights = b.weights)$fitted.values)
  ps.b <- as.vector(ps.b)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  # reg.coeff.pre.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                        subset = ((D==0) & (post==0)),
  #                                        weights = b.weights))
  control_pre <- (D == 0) & (post == 0)
  reg.coeff.pre.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[control_pre, , drop = FALSE],
                                  y = y[control_pre],
                                  weights = b.weights[control_pre],
                                  family = gaussian(link = "identity")
  ))
  out.y.pre.b <-   as.vector(tcrossprod(reg.coeff.pre.b, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  # reg.coeff.post.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                         subset = ((D==0) & (post==1)),
  #                                         weights = b.weights))
  control_post <- (D == 0) & (post == 1)
  reg.coeff.post.b <- stats::coef(fastglm::fastglm(
                                  x = int.cov[control_post, , drop = FALSE],
                                  y = y[control_post],
                                  weights = b.weights[control_post],
                                  family = gaussian(link = "identity")
  ))
  out.y.post.b <-   as.vector(tcrossprod(reg.coeff.post.b, int.cov))
  # Combine the ORs
  out.y.b <- post * out.y.post.b + (1 - post) * out.y.pre.b

  # Compute AIPW estimator
  att.b <- aipw_did_rc1(y, post, D, ps.b, out.y.b, b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}
