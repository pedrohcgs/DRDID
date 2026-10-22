# Check the ORIGINAL exponential-IPT score, including candidates produced by the
# continued secondary objective. This is a numerical stationarity check, not a
# proof of existence or a guarantee about subsequent trimming or outcome fits.
# The roundoff allowances are practical safeguards, not interval arithmetic.
ipt.valid <- function(gamma, D, int.cov, i.weights){
  k <- ncol(int.cov)
  if (!is.numeric(gamma) || length(gamma) != k || any(!is.finite(gamma)) ||
      any(!is.finite(i.weights)) || any(i.weights < 0)) return(FALSE)
  active <- i.weights > 0
  x <- int.cov[active, , drop = FALSE]
  d <- D[active]
  if (nrow(x) < k || k == 0 || any(!is.finite(x)) ||
      any(!d %in% c(0, 1)) || !any(d == 1) || sum(d == 0) < k) return(FALSE)

  # Equilibrate columns, then work in an orthonormal basis of the supplied span.
  # No intercept is added. Reject numerically rank-deficient designs rather than
  # squaring their condition number in the original crossproduct matrix.
  basis <- function(a){
    z <- svd(a, nu = min(dim(a)), nv = min(dim(a)))
    # A merely nonzero singular value is insufficient: an inaccurate retained
    # span can hide imbalance. Keep a precision margin for the 1e-6 score gate.
    relative.cutoff <- max(max(dim(a)) * .Machine$double.eps,
                           sqrt(.Machine$double.eps))
    if (length(z$d) < k || any(!is.finite(z$d)) ||
        min(z$d) <= relative.cutoff * max(z$d)) return(NULL)
    z
  }
  scales <- apply(abs(x), 2, max)
  if (any(scales == 0)) return(FALSE)
  scaled <- sweep(x, 2, scales, "/")
  design <- tryCatch(basis(scaled), error = function(e) NULL)
  if (is.null(design)) return(FALSE)
  z <- design$u

  # Bound uncertainty in the retained span using its reconstruction residual,
  # allowing for the column divisions and length-k reconstruction products.
  # A small span error can still matter when control curvature is weak below.
  eps <- .Machine$double.eps
  product.factor <- (k + 2) * eps / (1 - (k + 2) * eps)
  ud <- sweep(z, 2, design$d, "*")
  reconstructed <- tcrossprod(ud, design$v)
  absolute.product <- tcrossprod(abs(ud), abs(design$v))
  v.error <- norm(crossprod(design$v) - diag(k), "F") +
    product.factor * norm(crossprod(abs(design$v)), "F")
  if (!is.finite(v.error) || v.error >= 1) return(FALSE)
  span.error <- (norm(scaled - reconstructed, "F") +
                  product.factor * (norm(scaled, "F") + norm(absolute.product, "F"))) /
    (min(design$d) * sqrt(1 - v.error))
  eta <- as.vector(x %*% gamma)
  if (any(!is.finite(eta))) return(FALSE)

  # Evaluate only positive-weight controls' exponentials, on a common log scale.
  # This avoids 0 * Inf and overflow in exp(eta) when the effective mass is finite.
  log.weights <- log(i.weights[active])
  log.mass <- log.weights
  log.mass[d == 0] <- log.mass[d == 0] + eta[d == 0]
  if (any(!is.finite(log.mass))) return(FALSE)
  common.scale <- max(log.mass)
  shifted <- log.mass - common.scale
  mass <- exp(shifted)

  # Cancellation in X gamma can corrupt the exponential score even when its
  # span is accurately represented. Bound control index arithmetic before
  # exponentiating. Treated indices do not enter their score masses.
  absolute.index <- as.vector(abs(x[d == 0, , drop = FALSE]) %*% abs(gamma))
  if (any(!is.finite(absolute.index))) return(FALSE)
  log.error <- product.factor * (1 + abs(log.weights) +
                                  abs(log.mass) + abs(common.scale))
  log.error[d == 0] <- log.error[d == 0] +
    product.factor * (absolute.index / (1 - product.factor) +
                        abs(log.weights[d == 0]) + abs(eta[d == 0]))
  # The computed common scale is an arbitrary shared normalization, so its
  # own error need not be propagated. Include log/add/subtract/exp roundoff;
  # the normal minimum is a conservative absolute allowance for underflow.
  mass.error <- (mass * (expm1(log.error) + product.factor) +
                   .Machine$double.xmin * exp(log.error)) / (1 - product.factor)
  if (any(!is.finite(mass.error))) return(FALSE)
  score <- as.vector(crossprod(z, ifelse(d == 1, mass, -mass)))
  controls <- z[d == 0, , drop = FALSE] * sqrt(mass[d == 0])
  curvature <- tryCatch(basis(controls), error = function(e) NULL)
  if (is.null(curvature)) return(FALSE)

  # H = V diag(s^2) V'; solve by the singular values of sqrt(mass) * Z.
  # Both diagnostics are invariant to common weight scaling and nonsingular
  # model-basis changes in exact arithmetic (subject to numerical rank limits).
  projected <- as.vector(crossprod(curvature$v, score)) / curvature$d
  step <- as.vector(curvature$v %*% (projected / curvature$d))
  newton <- sqrt(sum(projected^2) / sum(mass))
  index.step <- max(abs(z %*% step))

  # Couple span accuracy to control curvature: separate condition cutoffs can
  # accept an inaccurate span whose error is amplified by the Newton solve.
  # Q=Z+E is a basis of the original model, with ||E|| bounded by span.error.
  # Weyl's bound controls its curvature; perturbing H b = S controls the change
  # in fitted-index steps. Include practical dot-product/SVD/solve allowances.
  dot.factor <- (nrow(z) + k + 2) * eps / (1 - (nrow(z) + k + 2) * eps)
  score.error <- span.error * sqrt(sum((mass + mass.error)^2)) +
    sqrt(sum(colSums(abs(z) * mass.error)^2)) +
    dot.factor * sqrt(sum(colSums(abs(z) * (mass + mass.error))^2))
  control.mass <- mass[d == 0]
  control.mass.error <- mass.error[d == 0]
  root.mass.error <- sqrt(control.mass.error)
  positive <- control.mass > 0
  root.mass.error[positive] <- control.mass.error[positive] /
    (sqrt(pmax(control.mass[positive] - control.mass.error[positive], 0)) +
       sqrt(control.mass[positive]))
  control.error <- sqrt(max(control.mass + control.mass.error)) * span.error +
    norm(z[d == 0, , drop = FALSE] * root.mass.error, "F") +
    dot.factor * norm(controls, "F")
  lower <- min(curvature$d) - control.error
  if (!is.finite(lower) || lower <= 0) return(FALSE)
  hessian.error <- 2 * max(curvature$d) * control.error + control.error^2
  step.norm <- sqrt(sum(step^2))
  solve.residual <- score - as.vector(crossprod(controls, controls %*% step))
  step.error <- (score.error + sqrt(sum(solve.residual^2)) +
                   hessian.error * step.norm) / lower^2
  index.error <- (sqrt(max(rowSums(z^2))) + span.error) * step.error +
    span.error * step.norm + product.factor * max(abs(z) %*% abs(step))
  # Requiring a small change in every active fitted index also prevents large
  # covariate shifts from hiding imbalance in a coordinate-relative score ratio.
  tolerance <- 1e-6
  is.finite(newton) && is.finite(index.step) && is.finite(index.error) &&
    newton <= tolerance && index.step + index.error <= tolerance
}
