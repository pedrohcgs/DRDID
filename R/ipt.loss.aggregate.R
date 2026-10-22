# Exceptional-range arithmetic for the two IPT criteria. Reliable ordinary
# products retain their mantissas; exceptional products retain binary scaling
# until the complete signed sum and original divisor have been applied.
ipt.loss.prepare <- function(D, X, iw) {
  n <- nrow(X)
  active <- iw != 0
  if (!all(active)) {
    D <- D[active]
    X <- X[active, , drop = FALSE]
    iw <- iw[active]
  }
  abs.X <- abs(X)
  X.range <- if (length(abs.X)) range(abs.X) else c(0, 0)
  if (isTRUE(X.range[1] == 0 && X.range[2] > 0)) {
    X.range[1] <- min(abs.X[abs.X != 0])
  }
  ordinary.range <- !any(abs(iw) < 1e-50 | abs(iw) > 1e50) &&
    all(is.finite(X.range)) && X.range[2] <= 1e50 &&
    (X.range[1] == 0 || X.range[1] >= 1e-50)
  list(D = D, X = X, iw = iw, n = n, control = D == 0,
       sample.fraction = length(iw) / n, ordinary.range = ordinary.range,
       X.range = X.range)
}

ipt.loss.in.range <- function(index, control, prepared, n) {
  if (n < 2 || n > 1e50 || !prepared$ordinary.range) return(FALSE)
  if (!length(index)) return(TRUE)
  bounds <- range(abs(index))
  # These global bounds imply both group-specific checks below. All-zero
  # indices are safe; mixed zero/nonzero or wider ranges use the exact checks.
  if (isTRUE(bounds[2] <= 100 && (bounds[1] >= 1e-100 || bounds[2] == 0))) return(TRUE)
  if (any(abs(index[control]) > 100)) return(FALSE)
  linear.index <- abs(index[!control])
  !any(linear.index > 1e100 | (linear.index != 0 & linear.index < 1e-100))
}

# Certify every nonzero derivative product with scalar magnitude bounds.
# Both Hessian multiplication orders are covered by checking each successive
# X factor. The margins keep a rounded bound away from the normal/finite
# boundaries; underflow or overflow while forming a bound simply rejects it.
ipt.loss.normal.products <- function(loading, X.range, degree) {
  if (X.range[2] == 0) return(TRUE)
  magnitude <- abs(loading[loading != 0])
  if (!length(magnitude)) return(TRUE)
  bounds <- range(magnitude)
  for (i in seq_len(degree)) {
    bounds <- bounds * X.range
    if (any(!is.finite(bounds)) || bounds[1] < 4 * .Machine$double.xmin ||
        bounds[2] > .Machine$double.xmax / 4) return(FALSE)
  }
  TRUE
}

ipt.loss.products.in.range <- function(terms, prepared) {
  all(terms$value.safe) && all(terms$gradient.safe) && all(terms$hessian.safe) &&
    ipt.loss.normal.products(terms$gradient, prepared$X.range, 1) &&
    ipt.loss.normal.products(terms$hessian, prepared$X.range, 2)
}

# A product is reliable only if no earlier multiplication lost range.
# An exact zero factor annihilates even an out-of-range nonlinear loading.
ipt.loss.product <- function(a, b, safe.a = TRUE, safe.b = TRUE) {
  zero <- (a == 0 & safe.a) | (b == 0 & safe.b)
  zero[is.na(zero)] <- FALSE
  value <- a * b
  value[zero] <- 0
  safe <- zero | (safe.a & safe.b & is.finite(value) &
                    abs(value) >= .Machine$double.xmin)
  safe[is.na(safe)] <- FALSE
  list(value = value, safe = safe)
}

# log2 selects an integer power only: division by that exact power retains
# the represented mantissa, including subnormal inputs and values near xmax.
ipt.loss.binary <- function(value) {
  value <- unname(value)
  exponent <- pmax(-1074, pmin(1023, floor(log2(abs(value)))))
  exponent[value == 0] <- 0
  cbind(mantissa = value / 2^exponent, exponent = exponent)
}

ipt.loss.binary.product <- function(factors, index = NULL) {
  if (any(factors == 0, na.rm = TRUE)) return(c(0, 0))
  if (any(!is.finite(factors))) return(c(NaN, 0))
  direction <- prod(sign(factors))
  if (!is.null(index)) {
    if (is.na(index)) return(c(NaN, 0))
    # Exponential products have at most q and two X factors. Beyond these
    # bounds no nonzero finite factors can rescue their floating-point range.
    if (index > 4096) return(c(direction * Inf, 0))
    if (index < -4096) return(c(0, 0))
    pieces <- 2^ceiling(log2(max(1, abs(index) / 500)))
    factors <- c(factors, rep(exp(index / pieces), pieces))
  }
  mantissa <- 1
  exponent <- 0
  for (factor in abs(factors)) {
    binary <- ipt.loss.binary(factor)
    mantissa <- mantissa * binary[1, 1]
    exponent <- exponent + binary[1, 2]
    shift <- floor(log2(mantissa))
    mantissa <- mantissa / 2^shift
    exponent <- exponent + shift
  }
  unname(c(direction * mantissa, exponent))
}

# Recover complete products, including X factors before rounding. The centered
# tail is factored as s*s*p, with s=max(1,t), without forming a large square.
ipt.loss.recover <- function(index, D, X, iw, rows, j = NULL, l = NULL,
                             continued.n = NULL) {
  result <- matrix(0, length(rows), 2)
  for (r in seq_along(rows)) {
    i <- rows[r]
    factors <- iw[i]
    if (!is.null(j)) factors <- c(factors, X[i, j])
    if (!is.null(l)) factors <- c(factors, X[i, l])
    if (D[i] != 0) {
      if (!is.null(l)) next
      factors[1] <- -factors[1]
      if (is.null(j)) factors <- c(factors, index[i])
      result[r, ] <- ipt.loss.binary.product(factors)
    } else if (is.null(continued.n) || index[i] < log(continued.n - 1)) {
      result[r, ] <- ipt.loss.binary.product(factors, index[i])
    } else {
      offset <- index[i] - log(continued.n - 1)
      factors <- c(factors, continued.n - 1)
      if (is.null(j)) {
        scale <- max(1, offset)
        polynomial <- if (offset <= 1) 1 + offset + 0.5 * offset^2 else
          0.5 + 1 / offset + (1 / offset)^2
        factors <- c(factors, scale, scale, polynomial)
      } else if (is.null(l)) {
        factors <- c(factors, 1 + offset)
      }
      result[r, ] <- ipt.loss.binary.product(factors)
    }
  }
  result
}

ipt.loss.binary.sum <- function(binary, divisor = 1) {
  if (anyNA(binary)) return(NaN)
  mantissa <- unname(binary[, 1])
  exponent <- unname(binary[, 2])
  infinite <- is.infinite(mantissa)
  if (any(infinite)) {
    # A cutoff infinity has unknown magnitude relative to other extended
    # overflowing terms. Opposite signs are unresolved, not a signed infinity.
    overflowing <- infinite | (exponent > 1023 & mantissa != 0)
    directions <- unique(sign(mantissa[overflowing]))
    if (length(directions) > 1) return(NaN)
    return(directions[1] * Inf)
  }
  keep <- mantissa != 0
  mantissa <- mantissa[keep]
  exponent <- exponent[keep]
  if (!length(mantissa)) return(0)
  # Sum nearby exponents first. If large terms cancel, their residual is
  # reinserted before the next range; tiny terms are not discarded up front.
  while (length(mantissa) > 1) {
    largest <- max(exponent)
    near <- exponent >= largest - 512
    if (sum(near) == 1) {
      mantissa <- mantissa[near]
      exponent <- exponent[near]
      break
    }
    scaled <- mantissa[near] * 2^(exponent[near] - largest)
    total <- sum(scaled[order(abs(scaled), decreasing = TRUE)])
    mantissa <- mantissa[!near]
    exponent <- exponent[!near]
    if (total != 0) {
      residual <- ipt.loss.binary(total)
      mantissa <- c(mantissa, residual[1, 1])
      exponent <- c(exponent, largest + residual[1, 2])
    }
    if (!length(mantissa)) return(0)
  }
  # Apply N before the final normal/subnormal range rounding, rather than
  # rounding separately divided ordinary and exceptional contributions.
  quotient <- ipt.loss.binary(mantissa / divisor)
  mantissa <- unname(quotient[1, 1])
  exponent <- exponent + quotient[1, 2]
  if (exponent > 1023) return(sign(mantissa) * Inf)
  if (exponent < -1075) return(0)
  if (exponent == -1075) return((mantissa / 2) * 2^-1074)
  unname(mantissa * 2^exponent)
}

ipt.loss.aggregate <- function(ordinary, terms, X, repair, divisor = 1) {
  component <- function(raw, safe, previous, j = NULL, l = NULL) {
    if (all(safe) && is.finite(previous)) return(previous)
    binary <- matrix(0, length(raw), 2)
    binary[safe, ] <- ipt.loss.binary(raw[safe])
    if (any(!safe)) binary[!safe, ] <- repair(which(!safe), j, l)
    ipt.loss.binary.sum(binary, divisor)
  }
  k <- ncol(X)
  gradient <- numeric(k)
  hessian <- matrix(0, k, k)
  names(gradient) <- colnames(X)
  if (!is.null(colnames(X))) dimnames(hessian) <- list(colnames(X), colnames(X))
  for (j in seq_len(k)) {
    product <- ipt.loss.product(terms$gradient, X[, j], terms$gradient.safe)
    gradient[j] <- component(product$value, product$safe, ordinary$gradient[j], j)
    for (l in seq_len(k)) {
      # The primary crossproduct loads its column first; the secondary
      # crossproduct loads its row first. Check the order actually evaluated.
      first.column <- if (isTRUE(terms$hessian.column.first)) l else j
      second.column <- if (isTRUE(terms$hessian.column.first)) j else l
      first <- ipt.loss.product(terms$hessian, X[, first.column], terms$hessian.safe)
      product <- ipt.loss.product(first$value, X[, second.column], first$safe)
      hessian[j, l] <- component(product$value, product$safe, ordinary$hessian[j, l], j, l)
    }
  }
  list(value = component(terms$value, terms$value.safe, ordinary$value),
       gradient = gradient, hessian = hessian)
}
