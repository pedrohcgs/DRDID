# Loss for the original exponential-IPT criterion, using trust (Tan, 2019).
loss.ps.cal <- function(gam, D, int.cov, iw) {
  n <- nrow(int.cov)
  k <- ncol(int.cov)
  if (any(is.na(gam))) {
    return(list(value = Inf, gradient = rep(NA, k), hessian = matrix(NA, k, k)))
  }

  # Exclude zero weights before indices or nonlinear arithmetic; keep n in
  # the mean criterion. Treated rows require only their linear contribution.
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
  control <- D == 0
  exp.ps.ind <- numeric(length(iw))
  exp.ps.ind[control] <- exp(ps.ind[control])

  # Retain the ordinary calculation path when its terms are in range.
  val <- -mean(ifelse(D, ps.ind, -exp.ps.ind) * iw) * sample.fraction
  grad <- -apply(ifelse(D, 1, -exp.ps.ind) * iw * int.cov, 2, mean) * sample.fraction
  hess <- -(t(int.cov) %*% (ifelse(D, 0, -exp.ps.ind) * iw * int.cov)) / n
  if (all(is.finite(c(val, grad, hess))) &&
      ipt.loss.in.range(ps.ind, control, int.cov, iw, n)) {
    return(list(value = val, gradient = grad, hessian = hess))
  }

  log.iw <- log(abs(iw))
  log.mass <- log.iw + ps.ind
  log.value <- ifelse(control, log.mass, log.iw + log(abs(ps.ind)))
  sign.value <- ifelse(control, sign(iw), -sign(iw) * sign(ps.ind))
  log.gradient <- ifelse(control, log.mass, log.iw)
  sign.gradient <- ifelse(control, sign(iw), -sign(iw))
  log.hessian <- ifelse(control, log.mass, -Inf)
  sign.hessian <- ifelse(control, sign(iw), 0)
  if (any(!is.finite(c(val, grad, hess))) ||
      any(!is.finite(exp.ps.ind[control]) |
            exp.ps.ind[control] < .Machine$double.xmin) ||
      ipt.loss.needs.log(log.value, sign.value, log.gradient, sign.gradient,
                         log.hessian, sign.hessian, int.cov)) {
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
    repair <- function(rows, j = NULL, l = NULL) {
      ipt.loss.recover(ps.ind, D, int.cov, iw, rows, j, l)
    }
    return(ipt.loss.aggregate(list(value = val, gradient = grad, hessian = hess),
                               terms, int.cov, repair, divisor = n))
  }
  list(value = val, gradient = grad, hessian = hess)
}
