# Bootrstapped Improved and locally efficient Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_drdid_imp_rc <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Pscore using the pscore.cal
  ps.b <- pscore.cal(D, int.cov, i.weights = b.weights, n = n)
  ps.b <- as.vector(ps.b$pscore)
  ps.b <- pmin(ps.b, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.cont.pre.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = T, treat = F)
  out.y.cont.pre.b <-  as.vector(out.y.cont.pre.b$out.reg)
  out.y.cont.post.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = F, treat = F)
  out.y.cont.post.b <-  as.vector(out.y.cont.post.b$out.reg)

  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                 subset = ((D==1) & (post==0)),
                                                 weights = b.weights))
  out.y.treat.pre.b <-   as.vector(tcrossprod(reg.treat.coeff.pre.b, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post.b <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                  subset = ((D==1) & (post==1)),
                                                  weights = b.weights))
  out.y.treat.post.b <-   as.vector(tcrossprod(reg.treat.coeff.post.b, int.cov))

  # Compute AIPW estimator
  att.b <- aipw_did_rc(y, post, D, ps.b,
                       out.y.treat.post.b, out.y.treat.pre.b,
                       out.y.cont.post.b, out.y.cont.pre.b,
                       b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}
