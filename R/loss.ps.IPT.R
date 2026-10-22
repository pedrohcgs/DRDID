# Continued IPT criterion of Graham, Pinto and Egel (2012, 2016).
loss.ps.IPT <- function(gamma1, n, D, int.cov, iw) {
  # The declared n determines the original continuation threshold, including
  # zero-weight rows excluded from arithmetic.
  active <- iw != 0
  if (!all(active)) {
    D <- D[active]
    int.cov <- int.cov[active, , drop = FALSE]
    iw <- iw[active]
  }
  k <- ncol(int.cov)
  if (!length(iw)) {
    res <- 0
    attr(res, "gradient") <- matrix(0, k, 1)
    attr(res, "hessian") <- matrix(0, k, k)
    return(res)
  }
  v <- as.vector(int.cov %*% gamma1)
  control <- D == 0
  exponential <- control & v < log(n - 1)
  tail <- control & !exponential
  t <- v[tail] - log(n - 1)

  # At controls, -phi(v)-v equals exp(v) below the join and the
  # centered quadratic below in the tail. Evaluate these directly to avoid
  # cancellation of v or of the constant in phi'(v)+1.
  value <- -iw * v
  score <- -iw
  curvature <- numeric(length(iw))
  ev <- exp(v[exponential])
  value[exponential] <- iw[exponential] * ev
  score[exponential] <- iw[exponential] * ev
  curvature[exponential] <- iw[exponential] * ev
  value[tail] <- iw[tail] * (n - 1) * (1 + t + 0.5 * t^2)
  score[tail] <- iw[tail] * (n - 1) * (1 + t)
  curvature[tail] <- iw[tail] * (n - 1)
  res <- sum(value)
  gradient <- t(int.cov) %*% score
  hessian <- t(curvature * int.cov) %*% int.cov
  if (all(is.finite(c(res, gradient, hessian))) &&
      ipt.loss.in.range(v, control, int.cov, iw, n)) {
    attr(res, "gradient") <- gradient
    attr(res, "hessian") <- hessian
    return(res)
  }

  log.iw <- log(abs(iw))
  log.value <- log.iw + log(abs(v))
  sign.value <- -sign(iw) * sign(v)
  log.gradient <- log.iw
  sign.gradient <- -sign(iw)
  log.hessian <- rep(-Inf, length(iw))
  sign.hessian <- numeric(length(iw))
  log.value[exponential] <- log.gradient[exponential] <-
    log.hessian[exponential] <- log.iw[exponential] + v[exponential]
  sign.value[control] <- sign.gradient[control] <- sign.hessian[control] <- sign(iw[control])
  # log(1+t+t^2/2), without squaring a large t first.
  log.polynomial <- numeric(length(t))
  small <- t <= 1
  log.polynomial[small] <- log1p(t[small] + 0.5 * t[small]^2)
  log.polynomial[!small] <- 2 * log(t[!small]) +
    log(0.5 + 1 / t[!small] + (1 / t[!small])^2)
  log.value[tail] <- log.iw[tail] + log(n - 1) + log.polynomial
  log.gradient[tail] <- log.iw[tail] + log(n - 1) + log1p(t)
  log.hessian[tail] <- log.iw[tail] + log(n - 1)

  if (any(!is.finite(c(res, gradient, hessian))) ||
      any(!is.finite(ev) | ev < .Machine$double.xmin) ||
      ipt.loss.needs.log(log.value, sign.value, log.gradient, sign.gradient,
                         log.hessian, sign.hessian, int.cov)) {
    stable <- ipt.loss.aggregate(log.value, sign.value, log.gradient, sign.gradient,
                                  log.hessian, sign.hessian, int.cov)
    res <- stable$value
    gradient <- matrix(stable$gradient, k, 1)
    dimnames(gradient) <- list(colnames(int.cov), NULL)
    hessian <- stable$hessian
  }
  attr(res, "gradient") <- gradient
  attr(res, "hessian") <- hessian
  res
}
