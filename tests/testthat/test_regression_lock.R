# =============================================================================
# Numeric regression lock + robustness guards.
#
# The rest of the DRDID test suite asserts only structural EQUIVALENCE (ATT == ATT
# across equivalent calls), so a uniform numeric shift in the estimation engine
# (e.g. the propensity / OLS solver) would pass every existing test undetected.
# This file pins the actual ATT and standard-error values of the six core 2x2
# estimators on a fixed dataset, and exercises the singular-design / bad-weight
# guards. Tolerance 1e-6 catches any meaningful engine change while tolerating
# cross-platform BLAS round-off (~1e-12).
# =============================================================================

make_lock_data <- function() {
  set.seed(20240115)
  n <- 1000
  X <- cbind(1, rnorm(n), rbinom(n, 1, 0.5))
  D <- rbinom(n, 1, plogis(-0.2 + 0.5 * X[, 2] - 0.3 * X[, 3]))
  y0 <- 1 + X[, 2] + 0.5 * X[, 3] + rnorm(n)
  y1 <- y0 + 0.8 * D + rnorm(n)
  post <- rbinom(n, 1, 0.5)
  y <- 1 + X[, 2] + 0.5 * X[, 3] + 0.8 * D * post + rnorm(n)
  list(X = X, D = D, w = rep(1, n), y0 = y0, y1 = y1, post = post, y = y)
}

test_that("core estimators reproduce locked ATT / se values (numeric regression lock)", {
  d <- make_lock_data()
  tol <- 1e-6
  expect_equal(drdid_panel(d$y1, d$y0, d$D, d$X, d$w)$ATT,       0.6551964166, tolerance = tol)
  expect_equal(drdid_panel(d$y1, d$y0, d$D, d$X, d$w)$se,        0.0687789916, tolerance = tol)
  expect_equal(std_ipw_did_panel(d$y1, d$y0, d$D, d$X, d$w)$ATT, 0.6555031513, tolerance = tol)
  expect_equal(std_ipw_did_panel(d$y1, d$y0, d$D, d$X, d$w)$se,  0.0689834951, tolerance = tol)
  expect_equal(reg_did_panel(d$y1, d$y0, d$D, d$X, d$w)$ATT,     0.6550830747, tolerance = tol)
  expect_equal(reg_did_panel(d$y1, d$y0, d$D, d$X, d$w)$se,      0.0686698690, tolerance = tol)
  expect_equal(drdid_rc(d$y, d$post, d$D, d$X, d$w)$ATT,         0.8186719766, tolerance = tol)
  expect_equal(drdid_rc(d$y, d$post, d$D, d$X, d$w)$se,          0.1331289343, tolerance = tol)
  expect_equal(std_ipw_did_rc(d$y, d$post, d$D, d$X, d$w)$ATT,   0.6780615210, tolerance = tol)
  expect_equal(std_ipw_did_rc(d$y, d$post, d$D, d$X, d$w)$se,    0.1925907969, tolerance = tol)
  expect_equal(reg_did_rc(d$y, d$post, d$D, d$X, d$w)$ATT,       0.8694663521, tolerance = tol)
  expect_equal(reg_did_rc(d$y, d$post, d$D, d$X, d$w)$se,        0.1682145937, tolerance = tol)
})

test_that("influence-function analytic invariants hold (mean(IF)=0, se = sd(IF) scaling)", {
  d <- make_lock_data()
  r <- drdid_panel(d$y1, d$y0, d$D, d$X, d$w, inffunc = TRUE)
  n <- length(d$D)
  expect_equal(length(r$att.inf.func), n)
  expect_lt(abs(mean(r$att.inf.func)), 1e-8)
  expect_equal(stats::sd(r$att.inf.func) * sqrt(n - 1) / n, r$se, tolerance = 1e-10)
})

test_that("rank-deficient designs are rejected by every core estimator (incl. std_ipw_did_rc)", {
  set.seed(1)
  n <- 800; x <- rnorm(n); Xrd <- cbind(1, x, x)   # collinear -> rank deficient
  D <- rbinom(n, 1, plogis(x)); w <- rep(1, n)
  y0 <- x + rnorm(n); y1 <- y0 + D + rnorm(n)
  post <- rbinom(n, 1, 0.5); y <- x + D * post + rnorm(n)
  expect_error(drdid_panel(y1, y0, D, Xrd, w),       "singular")
  expect_error(std_ipw_did_panel(y1, y0, D, Xrd, w), "singular")
  expect_error(reg_did_panel(y1, y0, D, Xrd, w),     "singular")
  expect_error(drdid_rc(y, post, D, Xrd, w),         "singular")
  # std_ipw_did_rc previously SILENTLY returned a corrupted SE here; now it errors.
  expect_error(std_ipw_did_rc(y, post, D, Xrd, w),   "singular")
  expect_error(reg_did_rc(y, post, D, Xrd, w),       "singular")
})

test_that("negative weights are rejected", {
  d <- make_lock_data()
  bad_w <- ifelse(d$D == 1, -1, 1)
  expect_error(drdid_panel(d$y1, d$y0, d$D, d$X, bad_w), "non-negative")
  expect_error(drdid_rc(d$y, d$post, d$D, d$X, bad_w),   "non-negative")
})
