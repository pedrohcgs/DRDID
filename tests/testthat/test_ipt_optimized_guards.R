context("IPT optimization safety across evaluations")

test_that("repeated coefficient evaluations cannot retain a range certificate", {
  for (k in c(1L, 20L)) {
    X <- matrix(1, 2, k, dimnames = list(NULL, paste0("x", seq_len(k))))
    if (k > 1L) X[, -1L] <- 0
    D <- c(0, 1)
    w <- c(1e300, 1)
    # The first and last calls have normal exponentials; the middle call has
    # an underflowed exponential rescued by its positive control weight.
    for (eta in c(0, -1000, 0)) {
      gamma <- c(eta, rep(0, k - 1L))
      mass <- if (eta == 0) 1e300 else (exp(-500) * 1e150)^2
      primary <- DRDID:::loss.ps.cal(gamma, D, X, w)
      secondary <- DRDID:::loss.ps.IPT(gamma, 7, D, X, w)
      expect_equal(primary$hessian[1, 1] / (mass / 2), 1, tolerance = 1e-12)
      expect_equal(attr(secondary, "hessian")[1, 1] / mass, 1,
                   tolerance = 1e-12)
      expect_identical(names(primary$gradient), colnames(X))
      expect_identical(dimnames(attr(secondary, "gradient")),
                       list(colnames(X), NULL))
      if (k > 1L) {
        expect_true(all(primary$hessian[-1L, ] == 0))
        expect_true(all(attr(secondary, "hessian")[-1L, ] == 0))
      }
    }
  }
})

test_that("tiny positive weights and signed zero weights stay distinct", {
  tiny <- 2^-1074
  X <- matrix(c(2^500, 1, 1), ncol = 1,
              dimnames = list(NULL, "x"))
  D <- c(0, 1, 0)
  gamma <- 0
  expected <- tiny * 2^1000
  for (zero in c(0, -0)) {
    w <- c(tiny, zero, zero)
    primary <- DRDID:::loss.ps.cal(gamma, D, X, w)
    secondary <- DRDID:::loss.ps.IPT(gamma, 11, D, X, w)
    expect_identical(as.numeric(primary$hessian), expected / 3)
    expect_identical(as.numeric(attr(secondary, "hessian")), expected)
    # Excluded rows may contain arbitrarily large finite model values.
    changed <- X
    changed[2:3, 1] <- c(1e308, -1e308)
    expect_identical(DRDID:::loss.ps.cal(gamma, D, changed, w), primary)
    expect_identical(DRDID:::loss.ps.IPT(gamma, 11, D, changed, w), secondary)
  }
})

test_that("reused caller variables cannot contaminate independent objectives", {
  for (k in c(1L, 20L)) {
    X <- matrix(seq_len(6L * k) / (6L * k), 6, k,
                dimnames = list(NULL, paste0("x", seq_len(k))))
    D <- c(0, 1, 0, 1, 0, 1)
    w <- c(0.5, 1, 1.5, 2, 0, 1)
    gamma <- rep(0.1 / k, k)
    for (iteration in 1:3) {
      if (iteration == 2L) {
        X[, 1] <- -X[, 1]
        D <- 1 - D
        w <- rev(w)
      }
      if (iteration == 3L) {
        X <- X / 2
        w <- w * 3
      }
      eta <- as.vector(X %*% gamma)
      mass <- w * (1 - D) * exp(eta)
      score <- mass - w * D
      value <- sum(mass - w * D * eta)
      primary <- DRDID:::loss.ps.cal(gamma, D, X, w)
      secondary <- DRDID:::loss.ps.IPT(gamma, 17, D, X, w)
      expect_equal(primary$value, value / 6, tolerance = 1e-12)
      expect_equal(unname(primary$gradient), as.vector(crossprod(X, score)) / 6,
                   tolerance = 1e-12)
      expect_equal(primary$hessian, crossprod(X, mass * X) / 6,
                   tolerance = 1e-12)
      expect_equal(as.numeric(secondary), value, tolerance = 1e-12)
      expect_equal(attr(secondary, "gradient"), crossprod(X, score),
                   tolerance = 1e-12)
      expect_equal(attr(secondary, "hessian"), crossprod(X, mass * X),
                   tolerance = 1e-12)
    }
  }
})

test_that("continued objectives use declared n after zero-weight filtering", {
  X <- matrix(c(1, 1e308, -1e308), ncol = 1,
              dimnames = list(NULL, "x"))
  D <- c(0, 1, 0)
  w <- c(2.5, 0, -0)
  gamma <- 2
  for (n in c(2, 100, 2)) {
    eta <- gamma
    if (eta < log(n - 1)) {
      value <- score <- curvature <- 2.5 * exp(eta)
    } else {
      offset <- eta - log(n - 1)
      value <- 2.5 * (n - 1) * (1 + offset + offset^2 / 2)
      score <- 2.5 * (n - 1) * (1 + offset)
      curvature <- 2.5 * (n - 1)
    }
    result <- DRDID:::loss.ps.IPT(gamma, n, D, X, w)
    expect_equal(as.numeric(result), value, tolerance = 1e-12)
    expect_equal(as.numeric(attr(result, "gradient")), score, tolerance = 1e-12)
    expect_equal(as.numeric(attr(result, "hessian")), curvature, tolerance = 1e-12)
  }
})

test_that("both optimizer routes refresh coefficient-dependent arithmetic", {
  # Substitute optimizer call heads in a private clone, keeping the real loss
  # functions and the arguments supplied by pscore.cal. This probes successive
  # evaluations of the same prepared problem on both routes.
  fit <- DRDID:::pscore.cal
  replace <- function(node) {
    if (identical(node, quote(trust::trust))) return(as.name(".probe_primary"))
    if (identical(node, quote(stats::nlm))) return(as.name(".probe_secondary"))
    if (is.call(node)) node <- as.call(lapply(node, replace))
    node
  }
  body(fit) <- replace(body(fit))
  private <- new.env(parent = environment(fit))
  seen <- character()
  private$.probe_primary <- function(objfun, parinit, rinit, rmax, iterlim, ...) {
    for (eta in c(0, -1000, 0)) {
      result <- objfun(eta / 2^500, ...)
      expected <- if (eta == 0) 2^999 else (exp(-500) * 2^500)^2 / 2
      expect_equal(as.numeric(result$hessian) / expected, 1, tolerance = 1e-12)
    }
    seen <<- c(seen, "primary")
    list(argument = 1e4)
  }
  private$.probe_secondary <- function(f, p, iterlim, gradtol, check.analyticals, ...) {
    for (eta in c(0, -1000, 0)) {
      result <- f(eta / 2^500, ...)
      expected <- if (eta == 0) 2^1001 else 2 * (exp(-500) * 2^500)^2
      expect_equal(as.numeric(attr(result, "hessian")) / expected, 1,
                   tolerance = 1e-12)
    }
    seen <<- c(seen, "secondary")
    list(estimate = 0)
  }
  environment(fit) <- private
  result <- expect_warning(fit(c(0, 1, 0, 1), matrix(2^500, 4, 1), rep(1, 4), 4),
                           "Used the secondary IPT algorithm")
  expect_identical(seen, c("primary", "secondary"))
  expect_identical(result$pscore, rep(0.5, 4))
  expect_equal(result$flag, 1)
})
