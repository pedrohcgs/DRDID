# Check the ORIGINAL exponential-IPT score, including candidates produced by the
# continued secondary objective. This is a numerical stationarity check, not a
# proof of existence or a guarantee about subsequent trimming or outcome fits.
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
    if (length(z$d) < k || any(!is.finite(z$d)) ||
        min(z$d) <= max(dim(a)) * .Machine$double.eps * max(z$d)) return(NULL)
    z
  }
  scales <- apply(abs(x), 2, max)
  if (any(scales == 0)) return(FALSE)
  design <- tryCatch(basis(sweep(x, 2, scales, "/")), error = function(e) NULL)
  if (is.null(design)) return(FALSE)
  z <- design$u
  eta <- as.vector(x %*% gamma)
  if (any(!is.finite(eta))) return(FALSE)

  # Evaluate only positive-weight controls' exponentials, on a common log scale.
  # This avoids 0 * Inf and overflow in exp(eta) when the effective mass is finite.
  log.mass <- log(i.weights[active])
  log.mass[d == 0] <- log.mass[d == 0] + eta[d == 0]
  if (any(!is.finite(log.mass))) return(FALSE)
  mass <- exp(log.mass - max(log.mass))
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
  # Requiring a small change in every active fitted index also prevents large
  # covariate shifts from hiding imbalance in a coordinate-relative score ratio.
  tolerance <- 1e-6
  is.finite(newton) && is.finite(index.step) &&
    newton <= tolerance && index.step <= tolerance
}
