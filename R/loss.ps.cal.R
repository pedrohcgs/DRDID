# Loss function for estimation of the calibrated PS, using trust Based on Tan (2019).
loss.ps.cal <- function(gam, D, int.cov, iw){
  #gam: argument (pscore parameters)
  #D: Treatment/Group indicator
  #int.cov: covariate matrix; n x k (with intercept!)
  # iw: sampling weights (useful for weighted bootstrap, too)

  n <- dim(int.cov)[1]
  k <- dim(int.cov)[2]

  if (!any(is.na(gam))) {
    # Zero-weight rows contribute nothing. Exclude them before indices,
    # exponentials and matrix products can create 0 * Inf or 0 * NaN.
    active <- iw != 0
    if (!any(active)) {
      return(list(value = 0, gradient = rep(0, k), hessian = matrix(0, k, k)))
    }
    if (!all(active)) {
      D <- D[active]
      int.cov <- int.cov[active, , drop = FALSE]
      iw <- iw[active]
    }
    sample.fraction <- length(iw) / n
    ps.ind <- as.vector(int.cov %*% gam)
    #exp.ps.ind <- exp(-ps.ind)
    exp.ps.ind <- exp(ps.ind)
    #val <- base::mean(ifelse(D, exp.ps.ind, ps.ind) * iw)
    val <- - base::mean(ifelse(D, ps.ind, -exp.ps.ind) * iw) * sample.fraction

    #grad <- apply(ifelse(D, -exp.ps.ind, 1) * iw * int.cov, 2, mean)
    grad <- - apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean) * sample.fraction

    #hess <- (t(int.cov) %*% (ifelse(D, exp.ps.ind, 0) * iw * int.cov))/n
    hess <- - (t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov))/n

  } else {
    val <- Inf
    grad <- rep(NA, k)
    hess <- matrix(NA, k, k)
  }
  list(value=val, gradient=grad, hessian=hess)
}
