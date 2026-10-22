context("Weighted nonlinear IPT arithmetic")

expect_ipt_arithmetic_ratio <- function(actual, expected, tolerance = 1e-10) {
  expect_true(all(is.finite(actual)))
  expect_equal(as.vector(actual) / as.vector(expected),
               rep(1, length(expected)), tolerance = tolerance)
}

test_that("primary weighted exponentials survive intermediate overflow and underflow", {
  X <- matrix(1, 1, 1)
  overflow <- DRDID:::loss.ps.cal(720, 0, X, 2 * exp(-720))
  for (part in overflow) expect_ipt_arithmetic_ratio(part, 2)

  # Split the exponential into representable factors for an independent
  # reference, rather than forming exp(-1000) before applying its weight.
  expected <- (exp(-500) * 1e150)^2
  underflow <- DRDID:::loss.ps.cal(-1000, 0, X, 1e300)
  for (part in underflow) expect_ipt_arithmetic_ratio(part, expected)
  secondary <- DRDID:::loss.ps.IPT(-1000, 4, 0, X, 1e300)
  expect_ipt_arithmetic_ratio(as.numeric(secondary), expected)
  expect_ipt_arithmetic_ratio(attr(secondary, "gradient"), expected)
  expect_ipt_arithmetic_ratio(attr(secondary, "hessian"), expected)

  # A positive subnormal exponential has already lost precision before a large
  # weight restores an ordinary-range contribution. Rescue partial underflow,
  # not only an exponential that has rounded all the way to zero.
  expected <- (exp(-370) * 1e150)^2
  partial <- DRDID:::loss.ps.cal(-740, 0, X, 1e300)
  for (part in partial) expect_ipt_arithmetic_ratio(part, expected, 5e-13)
  partial <- DRDID:::loss.ps.IPT(-740, 4, 0, X, 1e300)
  expect_ipt_arithmetic_ratio(as.numeric(partial), expected, 5e-13)
  expect_ipt_arithmetic_ratio(attr(partial, "gradient"), expected, 5e-13)
  expect_ipt_arithmetic_ratio(attr(partial, "hessian"), expected, 5e-13)
})

test_that("ordinary no-intercept zero rows keep the ordinary calculation path", {
  X <- matrix(c(0, -1, 1, 2), ncol = 1)
  D <- c(0, 1, 0, 1)
  w <- c(1, 0.5, 2, 1)
  gamma <- 0.2
  # This private clone makes use of exceptional aggregation observable, without
  # changing the package namespace or relying on rounding differences.
  ordinary <- DRDID:::loss.ps.cal
  private <- new.env(parent = environment(ordinary))
  private$ipt.loss.aggregate <- function(...) stop("unnecessary fallback")
  private$ipt.loss.needs.log <- function(...) stop("unnecessary rowwise logarithms")
  environment(ordinary) <- private
  result <- ordinary(gamma, D, X, w)
  mass <- w * (1 - D) * exp(as.vector(X * gamma))
  expect_identical(result$value, -mean(ifelse(D, as.vector(X * gamma),
                                              -exp(as.vector(X * gamma))) * w))
  expect_identical(result$gradient, -apply(ifelse(D, 1, -exp(as.vector(X * gamma))) *
                                           w * X, 2, mean))
  expect_identical(result$hessian, crossprod(X, mass * X) / nrow(X))
  secondary <- DRDID:::loss.ps.IPT
  environment(secondary) <- private
  expect_equal(secondary(gamma, nrow(X), D, X, w),
               DRDID:::loss.ps.IPT(gamma, nrow(X), D, X, w))
})

test_that("representable derivatives survive an underflowing control mass", {
  X <- matrix(1e200, 1, 1)
  # eta=-1000 has an unrepresentable mass, but multiplying by X or X^2
  # before rounding gives representable gradient and curvature.
  expected_gradient <- exp(-500) * (exp(-500) * 1e200)
  expected_hessian <- (exp(-500) * 1e200)^2
  primary <- DRDID:::loss.ps.cal(-1e-197, 0, X, 1)
  secondary <- DRDID:::loss.ps.IPT(-1e-197, 4, 0, X, 1)
  expect_equal(primary$value, 0)
  expect_equal(as.numeric(secondary), 0)
  expect_ipt_arithmetic_ratio(primary$gradient, expected_gradient)
  expect_ipt_arithmetic_ratio(primary$hessian, expected_hessian)
  expect_ipt_arithmetic_ratio(attr(secondary, "gradient"), expected_gradient)
  expect_ipt_arithmetic_ratio(attr(secondary, "hessian"), expected_hessian)

  # Retain the original n and apply X factors before a tiny mass is rounded
  # away by division. The second observation has zero weight.
  tiny <- .Machine$double.xmin * .Machine$double.eps
  primary <- DRDID:::loss.ps.cal(log(tiny) / 2, c(0, 0), matrix(2, 2, 1), c(1, 0))
  expect_equal(as.numeric(primary$gradient) / tiny, 1)
  expect_equal(as.numeric(primary$hessian) / tiny, 2)
})

test_that("primary means remain finite when an unnormalized Hessian sum overflows", {
  X <- matrix(1, 6, 1)
  D <- rep(0, 6)
  w <- rep(4e307, 6)
  result <- DRDID:::loss.ps.cal(log(2), D, X, w)
  for (part in result) expect_ipt_arithmetic_ratio(part, 8e307)
  # The secondary criterion is a sum, not a mean: its actual value is outside
  # the representable range and must not be replaced by a finite surrogate.
  secondary <- DRDID:::loss.ps.IPT(log(2), 6, D, X, w)
  expect_equal(as.numeric(secondary), Inf)
  expect_true(all(attr(secondary, "gradient") == Inf))
  expect_true(all(attr(secondary, "hessian") == Inf))
})

test_that("continued controls avoid cancellation and overflowing tail polynomials", {
  result <- DRDID:::loss.ps.IPT(-50, 4, 0, matrix(1, 1, 1), 1e200)
  expected <- 1e200 * exp(-50)
  expect_ipt_arithmetic_ratio(as.numeric(result), expected)
  expect_ipt_arithmetic_ratio(attr(result, "gradient"), expected)
  expect_ipt_arithmetic_ratio(attr(result, "hessian"), expected)

  tail <- DRDID:::loss.ps.IPT(1e200, 4, 0, matrix(1, 1, 1), 1e-300)
  expect_ipt_arithmetic_ratio(as.numeric(tail), 1.5e100)
  expect_ipt_arithmetic_ratio(attr(tail, "gradient"), 3e-100)
  expect_ipt_arithmetic_ratio(attr(tail, "hessian"), 3e-300)
  # A treated observation has only a linear term, even at the same index.
  treated <- DRDID:::loss.ps.IPT(1e200, 4, 1, matrix(1, 1, 1), 1e-200)
  primary <- DRDID:::loss.ps.cal(1e200, 1, matrix(1, 1, 1), 1e-200)
  expect_ipt_arithmetic_ratio(as.numeric(treated), -1)
  expect_ipt_arithmetic_ratio(attr(treated, "gradient"), -1e-200)
  expect_equal(attr(treated, "hessian"), matrix(0, 1, 1))
  expect_ipt_arithmetic_ratio(primary$value, -1)
  expect_ipt_arithmetic_ratio(primary$gradient, -1e-200)
  expect_equal(primary$hessian, matrix(0, 1, 1))
})

test_that("continued control values preserve the original join and declared n", {
  for (n in c(2, 4, 12)) {
    for (offset in c(-0.01, 0, 0.01, 1)) {
      gamma <- log(n - 1) + offset
      result <- DRDID:::loss.ps.IPT(gamma, n, 0, matrix(1, 1, 1), 2.5)
      if (offset < 0) {
        value <- gradient <- hessian <- 2.5 * exp(gamma)
      } else {
        value <- 2.5 * (n - 1) * (1 + offset + offset^2 / 2)
        gradient <- 2.5 * (n - 1) * (1 + offset)
        hessian <- 2.5 * (n - 1)
      }
      expect_ipt_arithmetic_ratio(as.numeric(result), value)
      expect_ipt_arithmetic_ratio(attr(result, "gradient"), gradient)
      expect_ipt_arithmetic_ratio(attr(result, "hessian"), hessian)
    }
  }
})

test_that("extreme-range values and derivatives describe the same criteria", {
  # Finite differences of the scalar objective and of its supplied gradient
  # exercise the fallback, independently of its internal log aggregation.
  cases <- list(
    list(gamma = 720, h = 1e-4, n = 4, w = 2 * exp(-720), secondary = FALSE),
    list(gamma = -1000, h = 1e-4, n = 4, w = 1e300, secondary = FALSE),
    list(gamma = 1e200, h = 1e195, n = 4, w = 1e-300, secondary = TRUE)
  )
  for (case in cases) {
    evaluate <- function(gamma) {
      if (case$secondary) {
        result <- DRDID:::loss.ps.IPT(gamma, case$n, 0, matrix(1, 1, 1), case$w)
        list(value = as.numeric(result), gradient = attr(result, "gradient"),
             hessian = attr(result, "hessian"))
      } else {
        DRDID:::loss.ps.cal(gamma, 0, matrix(1, 1, 1), case$w)
      }
    }
    center <- evaluate(case$gamma)
    plus <- evaluate(case$gamma + case$h)
    minus <- evaluate(case$gamma - case$h)
    numerical_gradient <- (plus$value - minus$value) / (2 * case$h)
    numerical_hessian <- (plus$gradient - minus$gradient) / (2 * case$h)
    expect_ipt_arithmetic_ratio(numerical_gradient, center$gradient, 1e-7)
    expect_ipt_arithmetic_ratio(numerical_hessian, center$hessian, 1e-7)
  }
})

test_that("real optimizers can use the repaired positive-weight criteria", {
  # The original exponential root is 720. Treated eta=1440 must not be
  # exponentiated; tiny positive control weights offset exp(720).
  D <- rep(c(1, 0), each = 4)
  X <- matrix(ifelse(D == 1, 2, 1), ncol = 1)
  w <- ifelse(D == 1, 1, 2 * exp(-720))
  fit <- trust::trust(DRDID:::loss.ps.cal, parinit = 719, rinit = 1, rmax = 1000,
                      iterlim = 1000, D = D, int.cov = X, iw = w)
  expect_true(fit$converged)
  expect_equal(as.numeric(fit$argument), 720, tolerance = 1e-8)
  expect_lt(max(abs(fit$gradient)), 1e-8)

  # This root is below the continuation threshold. Directly evaluating its
  # control exponential avoids cancellation in both phi+v and phi'+1.
  D <- c(1, 0, 1, 0)
  X <- matrix(1, 4, 1)
  w <- ifelse(D == 1, 1e20 * exp(-50), 1e20)
  secondary <- stats::nlm(DRDID:::loss.ps.IPT, -49, n = 4, D = D, int.cov = X,
                          iw = w, gradtol = 1e-8, iterlim = 1000,
                          check.analyticals = FALSE)
  expect_true(secondary$code %in% c(1L, 2L))
  expect_equal(as.numeric(secondary$estimate), -50, tolerance = 1e-6)
  expect_lt(max(abs(secondary$gradient)), 1e-8)
})
