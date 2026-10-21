context("Weighted IPT fallback objective and optimization")

ipt_objective_numerical_gradient <- function(gamma, n, D, X, w, h = 1e-5) {
  vapply(seq_along(gamma), function(j) {
    plus <- minus <- gamma
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    (as.numeric(DRDID:::loss.ps.IPT(plus, n, D, X, w)) -
       as.numeric(DRDID:::loss.ps.IPT(minus, n, D, X, w))) / (2 * h)
  }, numeric(1))
}

ipt_objective_numerical_hessian <- function(gamma, n, D, X, w, h = 1e-6) {
  vapply(seq_along(gamma), function(j) {
    plus <- minus <- gamma
    plus[j] <- plus[j] + h
    minus[j] <- minus[j] - h
    (as.vector(attr(DRDID:::loss.ps.IPT(plus, n, D, X, w), "gradient")) -
       as.vector(attr(DRDID:::loss.ps.IPT(minus, n, D, X, w), "gradient"))) / (2 * h)
  }, numeric(length(gamma)))
}

test_that("IPT objective derivatives agree in the exponential, tail and join regions", {
  X <- cbind(1, c(-1.5, -1, -0.5, 0, 0.4, 0.8, 1.2, 1.6))
  n <- nrow(X)
  D <- c(0, 1, 0, 0, 1, 0, 1, 0)
  w <- c(0.2, 0.5, 1, 1.6, 0.7, 2.1, 1.4, 0.9)
  zero_w <- w
  zero_w[c(2, 6)] <- 0
  weights <- list(rep(1, n), w, w / mean(w), zero_w, 3.7 * w)
  vstar <- log(n - 1)
  coefficients <- list(exponential = c(-0.3, 0.2),
                       tail = c(vstar + 0.8, 0.15),
                       join = c(vstar, 0.25))
  for (region in names(coefficients)) {
    gamma <- coefficients[[region]]
    for (iw in weights) {
      result <- DRDID:::loss.ps.IPT(gamma, n, D, X, iw)
      numerical_gradient <- ipt_objective_numerical_gradient(gamma, n, D, X, iw)
      numerical_hessian <- ipt_objective_numerical_hessian(gamma, n, D, X, iw)
      expect_lt(max(abs(as.vector(attr(result, "gradient")) - numerical_gradient)),
                 1e-6)
      expect_lt(max(abs(attr(result, "hessian") - numerical_hessian)), 1e-4)

      if (region == "exponential") {
        # Below the extrapolation threshold, the fallback must nest the
        # independently implemented primary weighted IPT criterion (sum vs mean).
        primary <- DRDID:::loss.ps.cal(gamma, D, X, iw)
        expect_equal(as.numeric(result), n * primary$value, tolerance = 1e-12)
        expect_equal(as.vector(attr(result, "gradient")),
                     n * primary$gradient, tolerance = 1e-12)
        expect_equal(attr(result, "hessian"), n * primary$hessian,
                     tolerance = 1e-12)
      }
    }
  }
})

test_that("the entire IPT fallback criterion is weighted", {
  X <- cbind(1, c(-1, -0.5, 0, 0.4, 1, 1.5))
  D <- c(0, 1, 0, 1, 0, 1)
  w <- c(0.5, 0, 1.2, 2, 0, 1)
  gamma <- c(0.1, 0.4)
  result <- DRDID:::loss.ps.IPT(gamma, nrow(X), D, X, w)
  scaled <- DRDID:::loss.ps.IPT(gamma, nrow(X), D, X, 4.3 * w)
  expect_equal(as.numeric(scaled), 4.3 * as.numeric(result), tolerance = 1e-12)
  changed <- X
  changed[w == 0, 2] <- c(3, -2)
  # Rows with zero weight cannot contribute either part of the criterion.
  expect_equal(as.numeric(DRDID:::loss.ps.IPT(gamma, nrow(X), D, changed, w)),
               as.numeric(result), tolerance = 1e-12)
})

make_forced_ipt_fallback <- function() {
  # Clone the production function and replace only its primary optimizer call.
  # Namespace-qualified calls cannot be masked by a normal lexical binding.
  # This private AST substitution avoids any trust/DRDID namespace mutation and
  # works without a dependency on a particular testthat mocking API version.
  forced <- DRDID:::pscore.cal
  replacements <- 0L
  replace_trust <- function(node) {
    if (identical(node, quote(trust::trust))) {
      replacements <<- replacements + 1L
      return(as.name(".ipt_forced_trust_failure"))
    }
    if (is.call(node)) node <- as.call(lapply(node, replace_trust))
    node
  }
  body(forced) <- replace_trust(body(forced))
  stopifnot(replacements == 1L)
  private <- new.env(parent = environment(forced))
  private$.ipt_forced_trust_failure <- function(..., parinit) {
    list(converged = FALSE, argument = parinit)
  }
  environment(forced) <- private
  forced
}

test_that("the actual nlm fallback solves weighted IPT after primary optimizer failure", {
  # Three covariate support points with a misspecified logit give distinct MLE
  # and IPT solutions. Both groups have overlap at every support point.
  X <- cbind(1, rep(c(-1, 0, 1), each = 40))
  D <- c(rep(c(1, 0), c(8, 32)), rep(c(1, 0), c(28, 12)),
         rep(c(1, 0), c(20, 20)))
  n <- nrow(X)
  w <- exp(0.4 * X[, 2])
  zero_w <- w
  zero_w[c(1, 41, 81)] <- 0
  original_body <- body(DRDID:::pscore.cal)
  original_trust <- trust::trust
  forced <- make_forced_ipt_fallback()
  for (iw in list(rep(1, n), w / mean(w), zero_w / mean(zero_w), 3.7 * w)) {
    # The fault injection exercises branch correctness, not the natural
    # probability that the primary optimizer fails on observed data.
    result <- expect_warning(forced(D, X, iw, n),
                              "trust algorithm did not converge")
    expect_equal(result$flag, 1)
    moment <- colMeans(iw * (D - (1 - D) * result$pscore / (1 - result$pscore)) * X)
    expect_lt(max(abs(moment)), 1e-6)
    expect_lt(max(stats::qlogis(result$pscore)), log(n - 1))

    # The untouched primary optimizer uses its separate weighted criterion.
    primary <- DRDID:::pscore.cal(D, X, iw, n)
    expect_equal(primary$flag, 0)
    expect_equal(result$pscore, primary$pscore, tolerance = 1e-6)
  }
  expect_identical(body(DRDID:::pscore.cal), original_body)
  expect_identical(trust::trust, original_trust)
})
