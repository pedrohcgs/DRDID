# Bootrstapped Improved Doubly Robust Difference-in-Differences with Repeated Cross Section data
# 2 periods and 2 groups

wboot_drdid_imp_rc1 <- function(nn, n, y, post, D, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Pscore using the pscore.cal
  ps.b <- pscore.cal(D, int.cov, i.weights = b.weights, n = n)
  ps.b <- as.vector(ps.b$pscore)
  ps.b <- pmin(ps.b, 1 - 1e-6)
  #Compute the Outcome regression for the control group
  out.y.pre.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = TRUE, treat = FALSE)
  out.y.pre.b <-  as.vector(out.y.pre.b$out.reg)
  out.y.post.b <- wols_rc(y, post, D, int.cov, ps.b, b.weights, pre = FALSE, treat = FALSE)
  out.y.post.b <-  as.vector(out.y.post.b$out.reg)

  # Combine the ORs
  out.y.b <- post * out.y.post.b + (1 - post) * out.y.pre.b

  # Compute AIPW estimator
  att.b <- aipw_did_rc1(y, post, D, ps.b, out.y.b, b.weights)
  #-----------------------------------------------------------------------------
  return(att.b)
}
