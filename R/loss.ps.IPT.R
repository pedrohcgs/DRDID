# Continued IPT criterion of Graham, Pinto and Egel (2012, 2016).
loss.ps.IPT <- function(gamma1, n, D, int.cov, iw, .ipt = NULL) {
  # The declared n determines the original continuation threshold, including
  # zero-weight rows excluded from arithmetic.
  if (is.null(.ipt)) .ipt <- ipt.loss.prepare(D, int.cov, iw)
  D <- .ipt$D
  int.cov <- .ipt$X
  iw <- .ipt$iw
  k <- ncol(int.cov)
  if (!length(iw)) {
    res <- 0
    # Empty products retain the ordinary derivative column-name conventions.
    attr(res, "gradient") <- crossprod(int.cov, numeric(0))
    attr(res, "hessian") <- crossprod(int.cov)
    return(res)
  }
  v <- as.vector(int.cov %*% gamma1)
  control <- .ipt$control
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
      ipt.loss.in.range(v, control, .ipt, n)) {
    attr(res, "gradient") <- gradient
    attr(res, "hessian") <- hessian
    return(res)
  }

  # Complete products can be normal even after an intermediate lost range.
  # Build their reliability masks whenever the safe fast path was unavailable.
  value.safe <- ipt.loss.product(-iw, v)$safe
  gradient.safe <- hessian.safe <- rep(TRUE, length(iw))
  mass <- ipt.loss.product(iw[exponential], ev,
                           safe.b = is.finite(ev) & ev >= .Machine$double.xmin)
  value.safe[exponential] <- gradient.safe[exponential] <-
    hessian.safe[exponential] <- mass$safe
  coefficient <- ipt.loss.product(iw[tail], n - 1)
  polynomial <- 1 + t + 0.5 * t^2
  value.safe[tail] <- ipt.loss.product(coefficient$value, polynomial,
                                       coefficient$safe)$safe
  gradient.safe[tail] <- ipt.loss.product(coefficient$value, 1 + t,
                                          coefficient$safe)$safe
  hessian.safe[tail] <- coefficient$safe
  terms <- list(value = value, value.safe = value.safe,
                gradient = score, gradient.safe = gradient.safe,
                hessian = curvature, hessian.safe = hessian.safe)
  if (all(is.finite(c(res, gradient, hessian))) && ipt.loss.products.in.range(terms, .ipt)) {
    attr(res, "gradient") <- gradient
    attr(res, "hessian") <- hessian
    return(res)
  }
  repair <- function(rows, j = NULL, l = NULL) {
    ipt.loss.recover(v, D, int.cov, iw, rows, j, l, continued.n = n)
  }
  stable <- ipt.loss.aggregate(list(value = res, gradient = gradient, hessian = hessian),
                                terms, int.cov, repair)
  res <- stable$value
  gradient <- matrix(stable$gradient, k, 1)
  if (!is.null(colnames(int.cov))) dimnames(gradient) <- list(colnames(int.cov), NULL)
  hessian <- stable$hessian
  attr(res, "gradient") <- gradient
  attr(res, "hessian") <- hessian
  res
}
