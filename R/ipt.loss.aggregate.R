# Exceptional-range arithmetic for the two IPT criteria. Sum signed, scaled
# contributions without first forming a possibly overflowing weighted term.
# This does not change either criterion or certify arbitrary cancellation.
ipt.loss.in.range <- function(index, control, X, iw, n) {
  # A conservative ordinary-range box avoids rowwise logarithms on routine
  # optimizer iterations. Within these bounds, row contributions, X factors
  # and their sums are far from either floating-point range boundary.
  if (n < 2 || n > 1e50 || any(abs(iw) < 1e-50 | abs(iw) > 1e50)) return(FALSE)
  abs.X <- abs(X)
  if (any(abs.X > 1e50 | (abs.X != 0 & abs.X < 1e-50))) return(FALSE)
  if (any(abs(index[control]) > 100)) return(FALSE)
  linear.index <- abs(index[!control])
  !any(linear.index > 1e100 | (linear.index != 0 & linear.index < 1e-100))
}

ipt.log.sum <- function(logabs, direction, divisor = 1) {
  keep <- direction != 0 & logabs != -Inf
  if (anyNA(keep)) return(NaN)
  if (!any(keep)) return(0)
  logabs <- logabs[keep]
  direction <- direction[keep]
  infinite <- logabs == Inf
  if (any(infinite)) {
    if (length(unique(direction[infinite])) == 1) return(direction[infinite][1] * Inf)
    return(NaN)
  }
  largest <- max(logabs)
  scaled <- sum(direction * exp(logabs - largest))
  if (is.na(scaled)) return(NaN)
  if (scaled == 0) return(0)
  sign(scaled) * exp(largest + log(abs(scaled)) - log(divisor))
}

ipt.loss.needs.log <- function(log.value, sign.value, log.gradient,
                               sign.gradient, log.hessian, sign.hessian, X) {
  lower <- log(.Machine$double.xmin)
  upper <- log(.Machine$double.xmax)
  outside <- function(x, direction) {
    result <- direction != 0 & (x < lower | x > upper)
    anyNA(result) || any(result)
  }
  # Vectorize over rows; an optimizer can evaluate this with many observations.
  smallest <- rep(Inf, nrow(X))
  largest <- rep(-Inf, nrow(X))
  for (j in seq_len(ncol(X))) {
    log.column <- log(abs(X[, j]))
    smallest <- pmin(smallest, ifelse(X[, j] == 0, Inf, log.column))
    largest <- pmax(largest, log.column)
  }
  has.X <- largest != -Inf
  sign.gradient <- sign.gradient * has.X
  sign.hessian <- sign.hessian * has.X
  outside(log.value, sign.value) ||
    outside(log.gradient + smallest, sign.gradient) ||
    outside(log.gradient + largest, sign.gradient) ||
    outside(log.hessian + 2 * smallest, sign.hessian) ||
    outside(log.hessian + 2 * largest, sign.hessian)
}

ipt.loss.aggregate <- function(log.value, sign.value, log.gradient,
                               sign.gradient, log.hessian, sign.hessian,
                               X, divisor = 1) {
  k <- ncol(X)
  log.X <- log(abs(X))
  sign.X <- sign(X)
  gradient <- vapply(seq_len(k), function(j) {
    ipt.log.sum(log.gradient + log.X[, j], sign.gradient * sign.X[, j], divisor)
  }, numeric(1))
  hessian <- matrix(0, k, k)
  names(gradient) <- colnames(X)
  dimnames(hessian) <- list(colnames(X), colnames(X))
  for (j in seq_len(k)) {
    for (l in seq_len(j)) {
      hessian[j, l] <- ipt.log.sum(log.hessian + log.X[, j] + log.X[, l],
                                  sign.hessian * sign.X[, j] * sign.X[, l], divisor)
      hessian[l, j] <- hessian[j, l]
    }
  }
  list(value = ipt.log.sum(log.value, sign.value, divisor),
       gradient = gradient, hessian = hessian)
}
