context("Weighted nonlinear IPT arithmetic")

expect_ipt_arithmetic_ratio <- function(actual, expected, tolerance = 1e-10) {
  expect_true(all(is.finite(actual)))
  expect_equal(as.vector(actual) / as.vector(expected),
               rep(1, length(expected)), tolerance = tolerance)
}

test_that("all-zero continued derivatives preserve named and unnamed attributes", {
  named <- matrix(c(1, 1, -1, 1), 2, 2,
                   dimnames = list(c("pre", "post"), c("intercept", "trend")))
  for (X in list(named, unname(named))) {
    result <- DRDID:::loss.ps.IPT(c(0.2, 0.4), 2, c(1, 0), X, c(0, 0))
    expected_gradient <- matrix(0, 2, 1)
    expected_hessian <- matrix(0, 2, 2)
    if (!is.null(colnames(X))) {
      dimnames(expected_gradient) <- list(colnames(X), NULL)
      dimnames(expected_hessian) <- list(colnames(X), colnames(X))
    }
    expect_identical(as.numeric(result), 0)
    expect_identical(attr(result, "gradient"), expected_gradient)
    expect_identical(attr(result, "hessian"), expected_hessian)
  }
})

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
  expect_null(dimnames(attr(secondary, "gradient")))

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
  private$ipt.loss.product <- function(...) stop("unnecessary reliability masks")
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

test_that("subnormal tail scores cannot erase a representable gradient contrast", {
  tiny <- 2^-1074
  X <- matrix(c(2^27, 1), ncol = 1)
  gamma <- 1 + 2^-28
  # The exact control score is tiny*(2^27+1.5). Rounding it before X
  # changes its complete gradient to the treated term, falsely yielding zero.
  result <- DRDID:::loss.ps.IPT(gamma, 2, c(0, 1), X,
                               c(tiny, 0x1.0000004p-1020))
  expect_identical(as.numeric(attr(result, "gradient")), -2^-1048)
})

test_that("continued Hessians recover subnormal intermediates in either column order", {
  tiny <- 2^-1074
  X <- rbind(c(2^27 + 0.5, 2^27), c(-1, 1))
  gamma <- c(0, 1)
  # The positive off-diagonal term is 2^-1020 + 2^-1048. The second
  # row subtracts 2^-1020. Both complete row terms are ordinary-range numbers.
  for (order in list(1:2, 2:1)) {
    result <- DRDID:::loss.ps.IPT(gamma[order], 2, c(0, 0), X[, order],
                                 c(tiny, 2^-1020))
    expect_identical(attr(result, "hessian")[1, 2], 2^-1048)
    expect_identical(attr(result, "hessian")[2, 1], 2^-1048)
  }

  # The same loss of precision also matters without cancellation. Applying
  # the large X factor first supplies an independent representable reference.
  X <- rbind(c(1e8 + 0.5, 1e100), c(1e8 + 0.5, 1e100))
  expected <- 2 * (tiny * X[1, 2]) * X[1, 1]
  result <- DRDID:::loss.ps.IPT(c(1e92, 0), 2, c(0, 0), X, rep(tiny, 2))
  expect_ipt_arithmetic_ratio(attr(result, "hessian")[1, 2], expected, 1e-12)
  expect_ipt_arithmetic_ratio(attr(result, "hessian")[2, 1], expected, 1e-12)
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
  # exercise the fallback, independently of its internal binary aggregation.
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

test_that("an unrelated subnormal row cannot erase finite ordinary value differences", {
  q2 <- 0x1.d945df4f8ee2ep+1009
  X <- diag(3)
  D <- c(0, 1, 0)
  primary_values <- secondary_values <- numeric(2)
  control <- sqrt(exp(700))
  offset <- control - log(2)
  secondary_reference <- 2 * (1 + offset + 0.5 * offset^2) - q2
  for (i in seq_along(c(-700, -740))) {
    third <- c(-700, -740)[i]
    primary <- DRDID:::loss.ps.cal(c(700, q2, third), D, X, rep(1, 3))
    secondary <- DRDID:::loss.ps.IPT(c(control, q2, third), 3, D, X, rep(1, 3))
    primary_values[i] <- primary$value
    secondary_values[i] <- as.numeric(secondary)
    expect_lt(primary$value, 0)
    expect_lt(as.numeric(secondary), 0)
    # The surviving differences are hundreds of ulps. Ordinary reductions
    # differ slightly near cancellation; zero is not an acceptable answer.
    expect_ipt_arithmetic_ratio(primary$value, (exp(700) - q2) / 3, 0.01)
    expect_ipt_arithmetic_ratio(as.numeric(secondary), secondary_reference, 0.01)
    expect_null(names(primary$value))
    expect_null(names(secondary))
  }
  expect_equal(primary_values[2] / primary_values[1], 1, tolerance = 0.01)
  expect_equal(secondary_values[2] / secondary_values[1], 1, tolerance = 0.01)
})

test_that("ordinary value and derivative contrasts retain the actual weight mantissas", {
  q <- c(0x1.d945df4f8ec8ep+1009, 0x1.d945df4f8ee2ep+1009, 1)
  for (scale in c(1, exp(-600))) {
    w <- scale * q
    # References use the actual represented scaled weights, rather than
    # assuming multiplication by a common floating-point scale is exact.
    difference <- w[1] - w[2]
    X <- matrix(c(0, 1, -740), ncol = 1)
    primary <- DRDID:::loss.ps.cal(1, c(0, 1, 0), X, w)
    secondary <- DRDID:::loss.ps.IPT(1, 3, c(0, 1, 0), X, w)
    expect_ipt_arithmetic_ratio(primary$value, difference / 3, 0.01)
    expect_ipt_arithmetic_ratio(as.numeric(secondary), difference, 0.01)

    X <- rbind(c(1, 1, 0), c(1, -1, 0), c(0, 0, 1))
    primary <- DRDID:::loss.ps.cal(c(0, 0, -740), rep(0, 3), X, w)
    secondary <- DRDID:::loss.ps.IPT(c(0, 0, -740), 3, rep(0, 3), X, w)
    expect_ipt_arithmetic_ratio(primary$gradient[2], difference / 3, 0.01)
    expect_ipt_arithmetic_ratio(primary$hessian[1, 2], difference / 3, 0.01)
    expect_ipt_arithmetic_ratio(attr(secondary, "gradient")[2], difference, 0.01)
    expect_ipt_arithmetic_ratio(attr(secondary, "hessian")[1, 2], difference, 0.01)
  }
})

test_that("two exceptional products retain differences in their weight mantissas", {
  q <- c(0x1.d945df4f8ec8ep+1009, 0x1.d945df4f8ee2ep+1009)
  X <- rbind(c(1, 1), c(1, -1))
  for (scale in c(1, exp(-600))) {
    w <- scale * q
    # Both control exponentials are subnormal. Compute their common factor
    # without first rounding exp(-740), and subtract represented weights first.
    difference <- (w[1] - w[2]) * exp(-370) * exp(-370)
    primary <- DRDID:::loss.ps.cal(c(-740, 0), c(0, 0), X, w)
    secondary <- DRDID:::loss.ps.IPT(c(-740, 0), 2, c(0, 0), X, w)
    expect_ipt_arithmetic_ratio(primary$gradient[2], difference / 2, 0.01)
    expect_ipt_arithmetic_ratio(primary$hessian[1, 2], difference / 2, 0.01)
    expect_ipt_arithmetic_ratio(attr(secondary, "gradient")[2], difference, 0.01)
    expect_ipt_arithmetic_ratio(attr(secondary, "hessian")[1, 2], difference, 0.01)
  }
})

test_that("recovered tiny terms combine with cancellation before the original divisor", {
  X <- diag(3)
  w <- rep(1e300, 3)
  expected <- exp(-500) * 1e300 * exp(-500)
  primary <- DRDID:::loss.ps.cal(c(0, 1, -1000), c(0, 1, 0), X, w)
  secondary <- DRDID:::loss.ps.IPT(c(0, 1, -1000), 3, c(0, 1, 0), X, w)
  expect_ipt_arithmetic_ratio(primary$value, expected / 3)
  expect_ipt_arithmetic_ratio(as.numeric(secondary), expected)

  tiny <- .Machine$double.xmin * .Machine$double.eps
  # Individually dividing tiny and 0.49*tiny by two would round both to zero.
  primary <- DRDID:::loss.ps.cal(c(0, log(0.49)), c(0, 0), diag(2), rep(tiny, 2))
  expect_identical(as.numeric(primary$value) / tiny, 1)
  expect_null(names(primary$value))
})

test_that("exponential Hessian reliability follows the actual multiplication order", {
  X <- matrix(c(1e300, 1e-100), nrow = 1)
  gamma <- c(log(1e-300) / 1e300, 0)
  eta <- as.numeric(X %*% gamma)
  expected <- (exp(eta / 2) * X[1, 1]) * (exp(eta / 2) * X[1, 2])
  # mass*X2 underflows before multiplication by X1, while the reverse order
  # is representable. Column permutations and both criteria must retain the
  # same curvature, despite their different crossproduct operand order.
  for (order in list(1:2, 2:1)) {
    primary <- DRDID:::loss.ps.cal(gamma[order], 0, X[, order, drop = FALSE], 1)
    secondary <- DRDID:::loss.ps.IPT(gamma[order], 4, 0, X[, order, drop = FALSE], 1)
    expect_ipt_arithmetic_ratio(primary$hessian[1, 2], expected)
    expect_ipt_arithmetic_ratio(primary$hessian[2, 1], expected)
    expect_ipt_arithmetic_ratio(attr(secondary, "hessian")[1, 2], expected)
    expect_ipt_arithmetic_ratio(attr(secondary, "hessian")[2, 1], expected)
  }
})

test_that("exceptional objectives remain nonflat inside a rounded-log bin", {
  q2 <- 0x1.d945df4f8ee2ep+1009
  gamma <- c(700, q2, -740)
  plus <- minus <- gamma
  plus[2] <- plus[2] + 16 * 2^(1009 - 52)
  minus[2] <- minus[2] - 16 * 2^(1009 - 52)
  expect_identical(log(plus[2]), 700)
  expect_identical(log(minus[2]), 700)
  evaluate <- function(g) DRDID:::loss.ps.cal(g, c(0, 1, 0), diag(3), rep(1, 3))
  numerical <- (evaluate(plus)$value - evaluate(minus)$value) / (plus[2] - minus[2])
  expect_gt(abs(numerical), 0.1)
  # Use the represented step and allow ordinary cancellation quantization.
  expect_equal(as.numeric(numerical) / evaluate(gamma)$gradient[2], 1, tolerance = 0.1)
})

test_that("far-out-of-range indices have bounded recovery and exact zero derivatives", {
  X <- matrix(c(0, 1), nrow = 1)
  positive <- DRDID:::loss.ps.cal(c(0, 1e308), 0, X, 1)
  negative <- DRDID:::loss.ps.cal(c(0, -1e308), 0, X, 1)
  expect_equal(positive$value, Inf)
  expect_equal(positive$gradient, c(0, Inf))
  expect_equal(positive$hessian, matrix(c(0, 0, 0, Inf), 2, 2))
  expect_equal(negative$value, 0)
  expect_equal(negative$gradient, c(0, 0))
  expect_equal(negative$hessian, matrix(0, 2, 2))
})

test_that("opposing extended overflow terms do not produce a false signed infinity", {
  hi <- 0x1.0000000000001p+12
  lo <- 0x1.fffffffffffffp+11
  # 2*exp(lo) dominates exp(hi), but the latter crosses the bounded
  # exponential-recovery cutoff. Refuse to assign a sign to unresolved
  # opposing overflow terms rather than claiming the smaller term dominates.
  for (orientation in c(1, -1)) {
    X <- rbind(c(hi, -orientation), c(lo, orientation))
    for (order in list(1:2, 2:1)) {
      result <- DRDID:::loss.ps.cal(c(1, 0), c(0, 0), X[order, ], c(1, 2)[order])
      expect_identical(result$value, Inf)
      expect_true(is.nan(result$gradient[2]))
      expect_true(all(is.nan(result$hessian[c(2, 3)])))
      expect_identical(result$hessian[c(1, 4)], c(Inf, Inf))
    }
  }
})

test_that("cutoff infinities retain known same-sign and bounded-opposition results", {
  for (direction in c(1, -1)) {
    # A finite mantissa with an extended exponent can represent a term beyond
    # double range. Same signs remain unambiguous; bounded ordinary terms
    # cannot reverse an exponential product beyond the far-out cutoff.
    same_sign <- rbind(c(direction * Inf, 0), c(direction, 5910))
    bounded_opposition <- rbind(c(direction * Inf, 0), c(-direction, 1023))
    expect_identical(DRDID:::ipt.loss.binary.sum(same_sign), direction * Inf)
    expect_identical(DRDID:::ipt.loss.binary.sum(bounded_opposition), direction * Inf)
  }
})
