# Loss for the original exponential-IPT criterion, using trust (Tan, 2019).
loss.ps.cal <- function(gam, D, int.cov, iw, .ipt = NULL) {
  k <- ncol(int.cov)
  if (any(is.na(gam))) {
    return(list(value = Inf, gradient = rep(NA, k), hessian = matrix(NA, k, k)))
  }

  # Exclude zero weights before indices or nonlinear arithmetic; keep n in
  # the mean criterion. Treated rows require only their linear contribution.
  if (is.null(.ipt)) .ipt <- ipt.loss.prepare(D, int.cov, iw)
  n <- .ipt$n
  D <- .ipt$D
  int.cov <- .ipt$X
  iw <- .ipt$iw
  if (!length(iw)) {
    return(list(value = 0, gradient = rep(0, k), hessian = matrix(0, k, k)))
  }
  sample.fraction <- .ipt$sample.fraction
  ps.ind <- as.vector(int.cov %*% gam)
  control <- .ipt$control
  exp.ps.ind <- numeric(length(iw))
  exp.ps.ind[control] <- exp(ps.ind[control])

  # Retain the ordinary calculation path when its terms are in range.
  val <- -mean(ifelse(D, ps.ind, -exp.ps.ind) * iw) * sample.fraction
  grad <- -apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean) * sample.fraction
  hess <- -(t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov)) / n
  if (all(is.finite(c(val, grad, hess))) &&
      ipt.loss.in.range(ps.ind, control, .ipt, n)) {
    return(list(value = val, gradient = grad, hessian = hess))
  }

  # Outside the safe box, check the multiplications themselves. A rounded
  # subnormal intermediate can precede an ordinary-range complete product.
  exp.safe <- !control | (is.finite(exp.ps.ind) &
                           exp.ps.ind >= .Machine$double.xmin)
  mass <- ipt.loss.product(iw, exp.ps.ind, safe.b = exp.safe)
  value <- ipt.loss.product(-iw, ps.ind)
  value$value[control] <- mass$value[control]
  value$safe[control] <- mass$safe[control]
  loading <- -iw
  loading[control] <- mass$value[control]
  loading.safe <- !control | mass$safe
  terms <- list(value = value$value, value.safe = value$safe,
                gradient = loading, gradient.safe = loading.safe,
                hessian = mass$value, hessian.safe = mass$safe,
                hessian.column.first = TRUE)
  if (all(is.finite(c(val, grad, hess))) && ipt.loss.products.in.range(terms, .ipt)) {
    return(list(value = val, gradient = grad, hessian = hess))
  }
  repair <- function(rows, j = NULL, l = NULL) {
    ipt.loss.recover(ps.ind, D, int.cov, iw, rows, j, l)
  }
  ipt.loss.aggregate(list(value = val, gradient = grad, hessian = hess),
                     terms, int.cov, repair, divisor = n)
}
