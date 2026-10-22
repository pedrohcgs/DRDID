context("Zero-weight rows do not enter IPT objective arithmetic")

test_that("zero-weight extreme rows cannot corrupt either IPT criterion", {
  X <- cbind(1, rep(c(-1, 1), each = 4))
  D <- c(1, 0, 0, 0, 1, 1, 1, 0)
  gamma <- c(0, log(3))
  w <- rep(1, 8)
  primary <- DRDID:::loss.ps.cal(gamma, D, X, w)
  # Keep the same declared n in the continued objective: its continuation
  # threshold belongs to the original problem, including zero-weight rows.
  secondary <- DRDID:::loss.ps.IPT(gamma, 9, D, X, w)
  for (extreme in c(1e200, 1e308)) {
    for (group in c(0, 1)) {
      augmented_X <- rbind(X, c(1, extreme))
      augmented_D <- c(D, group)
      augmented_w <- c(w, 0)
      result <- DRDID:::loss.ps.cal(gamma, augmented_D, augmented_X, augmented_w)
      # The primary objective is a mean over the original nine rows.
      expect_equal(result$value, primary$value * 8 / 9, tolerance = 1e-12)
      expect_equal(result$gradient, primary$gradient * 8 / 9, tolerance = 1e-12)
      expect_equal(result$hessian, primary$hessian * 8 / 9, tolerance = 1e-12)
      continued <- DRDID:::loss.ps.IPT(gamma, 9, augmented_D, augmented_X, augmented_w)
      expect_equal(as.numeric(continued), as.numeric(secondary), tolerance = 1e-12)
      expect_equal(attr(continued, "gradient"), attr(secondary, "gradient"),
                   tolerance = 1e-12)
      expect_equal(attr(continued, "hessian"), attr(secondary, "hessian"),
                   tolerance = 1e-12)
    }
  }
})

test_that("an all-zero-weight criterion is zero even at extreme finite indices", {
  X <- cbind(1, c(-1e308, 1e308, 1e200, -1e200))
  D <- c(0, 1, 0, 1)
  gamma <- c(0, 2)
  primary <- DRDID:::loss.ps.cal(gamma, D, X, rep(0, 4))
  expect_identical(primary$value, 0)
  expect_equal(primary$gradient, c(0, 0))
  expect_equal(primary$hessian, matrix(0, 2, 2))
  secondary <- DRDID:::loss.ps.IPT(gamma, 4, D, X, rep(0, 4))
  expect_equal(as.numeric(secondary), 0)
  expect_equal(as.vector(attr(secondary, "gradient")), c(0, 0))
  expect_equal(attr(secondary, "hessian"), matrix(0, 2, 2))
})

test_that("excluding zero weights preserves model dimension and original continuation", {
  # Only two active rows, but the declared problem still has n=8. Index 1.5
  # is below log(7), whereas replacing n by the active row count changes it
  # to the polynomial branch and therefore changes the objective.
  n <- 8
  X <- matrix(c(1, 1, rep(1e308, 6)), ncol = 1)
  D <- c(1, 0, rep(0, 6))
  w <- c(2, 3, rep(0, 6))
  gamma <- 1.5
  value <- -2 * gamma + 3 * exp(gamma)
  gradient <- -2 + 3 * exp(gamma)
  hessian <- 3 * exp(gamma)
  primary <- DRDID:::loss.ps.cal(gamma, D, X, w)
  secondary <- DRDID:::loss.ps.IPT(gamma, n, D, X, w)
  expect_equal(primary$value, value / n, tolerance = 1e-12)
  expect_equal(primary$gradient, gradient / n, tolerance = 1e-12)
  expect_equal(primary$hessian, matrix(hessian / n, 1, 1), tolerance = 1e-12)
  expect_equal(as.numeric(secondary), value, tolerance = 1e-12)
  expect_equal(as.vector(attr(secondary, "gradient")), gradient, tolerance = 1e-12)
  expect_equal(attr(secondary, "hessian"), matrix(hessian, 1, 1), tolerance = 1e-12)
})
