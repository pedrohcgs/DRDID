context("TWFE covariate dimensions survive removal of a supplied intercept")

make_twfe_shape_data <- function() {
  set.seed(2861)
  n <- 72
  D <- rep(c(0, 0, 1, 1), length.out = n)
  post <- rep(c(0, 1, 0, 1), length.out = n)
  x <- rnorm(n)
  z <- rnorm(n)
  y0 <- 0.7 + 0.3 * D + 0.5 * x - 0.3 * z + rnorm(n)
  dy <- 0.4 + 1.2 * D + 0.25 * x + 0.15 * z * D + rnorm(n)
  list(n = n, D = D, post = post, x = x, z = z, y0 = y0,
       y1 = y0 + dy, y = y0 + post * dy, w = exp(0.2 * x - 0.1 * z))
}

# A separate numeric-design WLS calculation avoids the production formula and
# covariate preprocessing. Panel rows retain a common weight and one paired IF.
twfe_shape_reference <- function(d, panel, X, w) {
  if (panel) {
    dd <- rep(d$D, 2)
    post <- rep(0:1, each = d$n)
    y <- c(d$y0, d$y1)
    XX <- if (is.null(X)) NULL else rbind(X, X)
    ww <- rep(w / mean(w), 2)
  } else {
    dd <- d$D
    post <- d$post
    y <- d$y
    XX <- X
    ww <- w / mean(w)
  }
  Z <- cbind(1, dd, post, XX, interaction = dd * post)
  fit <- stats::lm.wfit(Z, y, ww)
  H <- crossprod(Z, ww * Z) / nrow(Z)
  phi <- as.vector((ww * fit$residuals * Z) %*% solve(H)[, ncol(Z)])
  if (panel) phi <- (phi[seq_len(d$n)] + phi[d$n + seq_len(d$n)]) / 2
  se <- if (panel) sd(phi) * sqrt(d$n - 1) / d$n else sd(phi) / sqrt(d$n)
  list(ATT = unname(tail(fit$coefficients, 1)), se = se, IF = phi)
}

twfe_shape_call <- function(d, panel, X, w, ...) {
  if (panel) {
    twfe_did_panel(d$y1, d$y0, d$D, X, w, ...)
  } else {
    twfe_did_rc(d$y, d$post, d$D, X, w, ...)
  }
}

test_that("an intercept and one covariate agree with independent weighted WLS", {
  d <- make_twfe_shape_data()
  wz <- d$w
  wz[seq(5, d$n, by = 11)] <- 0
  for (panel in c(TRUE, FALSE)) {
    for (w in list(rep(1, d$n), d$w, wz)) {
      expected <- twfe_shape_reference(d, panel, matrix(d$x, ncol = 1), w)
      for (X in list(d$x, matrix(d$x, ncol = 1), cbind(1, d$x))) {
        fit <- twfe_shape_call(d, panel, X, w, inffunc = TRUE)
        expect_equal(unname(fit$ATT), expected$ATT, tolerance = 1e-11)
        expect_equal(fit$se, expected$se, tolerance = 1e-11)
        expect_equal(fit$att.inf.func, expected$IF, tolerance = 1e-10)
        expect_length(fit$att.inf.func, d$n)
        expect_equal(as.numeric(c(fit$lci, fit$uci)),
                     expected$ATT + c(-1, 1) * 1.96 * expected$se, tolerance = 1e-11)
        expect_null(fit$boots)
      }
    }
  }
})

test_that("one-covariate bootstrap draws use the same unit or row weights", {
  d <- make_twfe_shape_data()
  w <- d$w
  w[seq(5, d$n, by = 11)] <- 0
  B <- 19
  for (panel in c(TRUE, FALSE)) {
    expected <- twfe_shape_reference(d, panel, matrix(d$x, ncol = 1), w)
    for (method in c("multiplier", "weighted")) {
      set.seed(829)
      if (method == "multiplier") {
        p <- (1 + sqrt(5)) / (2 * sqrt(5))
        draws <- vapply(seq_len(B), function(b) {
          v <- ifelse(rbinom(d$n, 1, p) == 1, (1 - sqrt(5)) / 2, (1 + sqrt(5)) / 2)
          mean(v * expected$IF)
        }, numeric(1))
      } else {
        draws <- vapply(seq_len(B), function(b)
          twfe_shape_reference(d, panel, matrix(d$x, ncol = 1), w * rexp(d$n))$ATT,
          numeric(1))
      }
      expected_rng <- .Random.seed
      expected_se <- IQR(draws) / (qnorm(0.75) - qnorm(0.25))
      for (X in list(d$x, matrix(d$x, ncol = 1), cbind(1, d$x))) {
        set.seed(829)
        fit <- twfe_shape_call(d, panel, X, w, boot = TRUE, boot.type = method,
                                nboot = B, inffunc = TRUE)
        expect_equal(unname(fit$ATT), expected$ATT, tolerance = 1e-11)
        expect_equal(as.vector(fit$boots), draws, tolerance = 1e-11)
        expect_equal(fit$se, expected_se, tolerance = 1e-11)
        expect_equal(fit$att.inf.func, expected$IF, tolerance = 1e-10)
        expect_identical(.Random.seed, expected_rng)
      }
    }
  }
})

test_that("intercept-only inputs keep the no-covariate TWFE path", {
  d <- make_twfe_shape_data()
  w <- d$w
  w[seq(5, d$n, by = 11)] <- 0
  for (panel in c(TRUE, FALSE)) {
    expected <- twfe_shape_reference(d, panel, NULL, w)
    for (method in c("analytic", "multiplier", "weighted")) {
      set.seed(307)
      baseline <- twfe_shape_call(d, panel, NULL, w, boot = method != "analytic",
                                   boot.type = if (method == "analytic") "weighted" else method,
                                   nboot = 19, inffunc = TRUE)
      baseline_rng <- .Random.seed
      expect_equal(unname(baseline$ATT), expected$ATT, tolerance = 1e-11)
      expect_equal(as.vector(baseline$att.inf.func), expected$IF, tolerance = 1e-10)
      if (method == "analytic") expect_equal(baseline$se, expected$se, tolerance = 1e-11)
      for (X in list(rep(1, d$n), matrix(1, d$n, 1))) {
        set.seed(307)
        fit <- twfe_shape_call(d, panel, X, w, boot = method != "analytic",
                                boot.type = if (method == "analytic") "weighted" else method,
                                nboot = 19, inffunc = TRUE)
        expect_equal(fit, baseline, tolerance = 1e-11)
        expect_identical(.Random.seed, baseline_rng)
      }
    }
  }
})

test_that("two-covariate fits and optional IF output retain their behavior", {
  d <- make_twfe_shape_data()
  XX <- cbind(d$x, d$z)
  for (panel in c(TRUE, FALSE)) {
    expected <- twfe_shape_reference(d, panel, XX, d$w)
    for (X in list(XX, cbind(1, XX))) {
      fit <- twfe_shape_call(d, panel, X, d$w, inffunc = TRUE)
      expect_equal(unname(fit$ATT), expected$ATT, tolerance = 1e-11)
      expect_equal(fit$se, expected$se, tolerance = 1e-11)
      expect_equal(fit$att.inf.func, expected$IF, tolerance = 1e-10)
    }
    expect_null(twfe_shape_call(d, panel, cbind(1, d$x), d$w)$att.inf.func)
  }
})
