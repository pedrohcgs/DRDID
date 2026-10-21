context("Improved estimators require a verified original IPT fit")

required_ipt_error <- paste0(
  "IPT propensity-score fit did not converge. Consider explicitly using the ",
  "traditional DR estimator by rerunning `drdid(..., estMethod = \"trad\")`."
)

make_required_ipt_fixture <- function() {
  # Both groups have the same symmetric covariate distribution, giving the
  # exact, identified ordinary IPT root gamma=(0,0).
  list(D = rep(c(1, 1, 0, 0), 2),
       X = cbind(1, rep(c(-1, 1), 4)), w = rep(1, 8), gamma = c(0, 0))
}

make_required_ipt_routes <- function(primary, secondary, initial = NULL) {
  # A private clone keeps the real candidate validator and substitutes only
  # optimizer/initialization call heads. Nothing is assigned to a namespace.
  forced <- DRDID:::pscore.cal
  counts <- c(primary = 0L, secondary = 0L, initial = 0L)
  replace_calls <- function(node) {
    if (identical(node, quote(trust::trust))) {
      counts["primary"] <<- counts["primary"] + 1L
      return(as.name(".required_ipt_primary"))
    }
    if (identical(node, quote(stats::nlm))) {
      counts["secondary"] <<- counts["secondary"] + 1L
      return(as.name(".required_ipt_secondary"))
    }
    if (identical(node, as.name("fastglm_fit"))) {
      counts["initial"] <<- counts["initial"] + 1L
      return(as.name(".required_ipt_initial"))
    }
    if (is.call(node)) node <- as.call(lapply(node, replace_calls))
    node
  }
  body(forced) <- replace_calls(body(forced))
  stopifnot(all(counts == 1L))
  private <- new.env(parent = environment(forced))
  calls <- new.env(parent = emptyenv())
  calls$primary <- calls$secondary <- 0L
  private$.required_ipt_primary <- function(objfun, parinit, ...) {
    calls$primary <- calls$primary + 1L
    calls$primary_start <- parinit
    if (is.function(primary)) primary(objfun, parinit, ...) else primary
  }
  private$.required_ipt_secondary <- function(f, p, ...) {
    calls$secondary <- calls$secondary + 1L
    calls$secondary_start <- p
    if (is.function(secondary)) secondary(f, p, ...) else secondary
  }
  private$.required_ipt_initial <- if (is.null(initial)) DRDID:::fastglm_fit else initial
  environment(forced) <- private
  list(fit = forced, calls = calls)
}

test_that("ordinary IPT validation is invariant to equivalent model coordinates", {
  d <- make_required_ipt_fixture()
  expect_true(DRDID:::ipt.valid(d$gamma, d$D, d$X, d$w))
  transforms <- list(diag(c(1, 1e-8)), diag(c(1, 1e8)),
                      matrix(c(1, 0, 1e6, 1), 2),
                      matrix(c(1, 1, 1, -1), 2))
  for (A in transforms) {
    expect_true(DRDID:::ipt.valid(as.vector(solve(A, d$gamma)),
                                  d$D, d$X %*% A, d$w))
  }
  # Also preserve a nonzero root, so a shifted covariate requires the matching
  # intercept adjustment rather than leaving all coefficients equal to zero.
  nonzero_D <- c(1, 0, 0, 0, 1, 1, 1, 0)
  nonzero_X <- cbind(1, rep(c(-1, 1), each = 4))
  nonzero_gamma <- c(0, log(3))
  for (A in c(list(diag(2)), transforms)) {
    expect_true(DRDID:::ipt.valid(as.vector(solve(A, nonzero_gamma)),
                                  nonzero_D, nonzero_X %*% A, d$w))
  }
  for (scale in c(1e-12, 1e12)) {
    expect_true(DRDID:::ipt.valid(d$gamma, d$D, d$X, scale * d$w))
  }

  # Root's finite-overlap counterexample: the raw per-column relative score
  # becomes tiny after shifting x, although the same material imbalance remains.
  D <- c(1, 1, 1, 1, 0, 0, 0, 0)
  X <- cbind(1, c(1, 1, 0, -1, -1, -1, 0, 1))
  expect_equal(colMeans((2 * D - 1) * X), c(0, 0.25))
  expect_false(DRDID:::ipt.valid(c(0, 0), D, X, rep(1, 8)))
  X[, 2] <- X[, 2] + 1e6
  expect_false(DRDID:::ipt.valid(c(0, 0), D, X, rep(1, 8)))
})

test_that("IPT validation handles zero weights and rejects nonfinite or unidentified fits", {
  d <- make_required_ipt_fixture()
  expect_true(DRDID:::ipt.valid(d$gamma, c(d$D, 1),
                                rbind(d$X, c(1, 1e150)), c(d$w, 0)))
  for (gamma in list(c(NA, 0), c(Inf, 0), c(NaN, 0), c(1000, 0), 0)) {
    expect_false(DRDID:::ipt.valid(gamma, d$D, d$X, d$w))
  }
  expect_false(DRDID:::ipt.valid(c(0, 0, 0), d$D,
                                 cbind(d$X, d$X[, 2]), d$w))
})

test_that("a verified primary candidate is accepted regardless of optimizer status", {
  d <- make_required_ipt_fixture()
  for (converged in c(TRUE, FALSE)) {
    routes <- make_required_ipt_routes(
      list(converged = converged, argument = d$gamma),
      function(...) stop("secondary must not run"))
    result <- routes$fit(d$D, d$X, d$w, length(d$D))
    expect_equal(result$flag, 0)
    expect_equal(result$pscore, rep(0.5, length(d$D)))
    expect_equal(routes$calls$secondary, 0)
  }
})

test_that("invalid primary candidates and exceptions can recover through verified IPT", {
  d <- make_required_ipt_fixture()
  bad_primary <- list(
    list(converged = TRUE, argument = c(0.4, -0.1)),
    list(converged = TRUE, argument = c(NA, 0)),
    function(...) stop("forced primary exception")
  )
  for (primary in bad_primary) {
    for (code in c(1, 3, 4)) {
      routes <- make_required_ipt_routes(primary, list(estimate = d$gamma, code = code))
      result <- suppressWarnings(routes$fit(d$D, d$X, d$w, length(d$D)))
      expect_equal(result$flag, 1)
      expect_equal(result$pscore, rep(0.5, length(d$D)))
      expect_equal(routes$calls$secondary, 1)
    }
  }
})

test_that("unsuccessful candidates never silently select ordinary logit", {
  d <- make_required_ipt_fixture()
  primary_body <- body(DRDID:::pscore.cal)
  original_trust <- trust::trust
  invalid_secondary <- list(
    list(estimate = c(0.4, -0.1), code = 1),
    list(estimate = c(0.4, -0.1), code = 5),
    list(estimate = c(Inf, 0), code = 1),
    list(estimate = 0, code = 1),
    function(...) stop("forced secondary exception")
  )
  for (secondary in invalid_secondary) {
    routes <- make_required_ipt_routes(
      list(converged = TRUE, argument = c(0.4, -0.1)), secondary)
    error <- tryCatch(routes$fit(d$D, d$X, d$w, length(d$D)), error = identity)
    expect_s3_class(error, "error")
    if (inherits(error, "error")) expect_identical(conditionMessage(error), required_ipt_error)
    expect_equal(routes$calls$secondary, 1)
  }
  routes <- make_required_ipt_routes(
    list(converged = TRUE, argument = c(0.4, -0.1)),
    list(estimate = c(0.4, -0.1), code = 1))
  error <- tryCatch(routes$fit(d$D, d$X, d$w, length(d$D), bootstrap.index = 7),
                     error = identity)
  expect_s3_class(error, "error")
  if (inherits(error, "error")) {
    expect_identical(conditionMessage(error),
                     paste0(required_ipt_error, "\nWeighted-bootstrap replicate 7 failed."))
  }
  expect_identical(body(DRDID:::pscore.cal), primary_body)
  expect_identical(trust::trust, original_trust)
})

test_that("a continued-loss optimum is rejected unless it solves original IPT", {
  # One control of weight 1 and seven treated observations of weight 3 give
  # original root log(21), but continued-loss optimum log(7)+2.
  D <- c(0, rep(1, 7))
  X <- matrix(1, 8, 1)
  w <- c(1, rep(3, 7))
  tail_root <- log(7) + 2
  expect_gt(abs(sum(w * (D - (1 - D) * exp(tail_root))) / sum(w)), 1)
  routes <- make_required_ipt_routes(
    list(converged = FALSE, argument = tail_root),
    list(estimate = tail_root, code = 1))
  error <- tryCatch(routes$fit(D, X, w, 8), error = identity)
  expect_s3_class(error, "error")
  if (inherits(error, "error")) expect_identical(conditionMessage(error), required_ipt_error)

  routes <- make_required_ipt_routes(
    list(converged = FALSE, argument = tail_root),
    list(estimate = log(21), code = 3))
  result <- suppressWarnings(routes$fit(D, X, w, 8))
  expect_equal(result$flag, 1)
  expect_equal(result$pscore, rep(21 / 22, 8), tolerance = 1e-12)
})

test_that("failed logit initialization still permits a verified IPT fit", {
  d <- make_required_ipt_fixture()
  initializers <- list(function(...) stop("forced initial-logit exception"),
                       function(...) list(coefficients = c(NA, Inf)))
  for (initial in initializers) {
    routes <- make_required_ipt_routes(
      list(converged = FALSE, argument = c(0.4, -0.1)),
      list(estimate = d$gamma, code = 1), initial = initial)
    result <- suppressWarnings(routes$fit(d$D, d$X, d$w, length(d$D)))
    expect_equal(routes$calls$primary_start, d$gamma)
    expect_equal(routes$calls$secondary_start, d$gamma)
    expect_equal(result$flag, 1)
    expect_equal(result$pscore, rep(0.5, length(d$D)))
  }
})

test_that("real IPT solvers preserve their solution under common weight scaling", {
  set.seed(17)
  n <- 120
  X <- cbind(1, stats::rnorm(n), stats::rnorm(n))
  D <- stats::rbinom(n, 1, stats::plogis(-0.2 + 0.4 * X[, 2] - 0.2 * X[, 3]))
  w <- exp(0.4 * X[, 2] - 0.2 * X[, 3] + 0.2 * stats::rnorm(n))
  w <- w / mean(w)
  reference <- DRDID:::pscore.cal(D, X, w, n)
  # With unnormalized optimizer inputs, the small-weight fit can satisfy the
  # solver's stopping rule before meeting the original-IPT validation criterion.
  for (scale in c(1e-3, 1e3)) {
    result <- DRDID:::pscore.cal(D, X, scale * w, n)
    expect_equal(result$pscore, reference$pscore, tolerance = 1e-6)
    expect_equal(result$flag, reference$flag)
  }
  # Execute the actual nlm fallback: its stopping accuracy must be sufficient
  # for the original-IPT gate, rather than rejecting a regular solvable target.
  routes <- make_required_ipt_routes(
    list(converged = FALSE, argument = rep(NA_real_, ncol(X))), stats::nlm)
  fallback <- expect_warning(routes$fit(D, X, w, n),
                              "trust algorithm did not converge")
  expect_equal(fallback$flag, 1)
  expect_equal(fallback$pscore, reference$pscore, tolerance = 1e-6)
  expect_equal(routes$calls$secondary, 1)
})

test_that("all improved estimators stop when no finite original-IPT fit exists", {
  cells <- expand.grid(D = 0:1, post = 0:1, level = 0:1, replicate = 1:2)
  x <- ifelse(cells$D == 1, 2 + cells$level, -1 + 2 * cells$level)
  X <- cbind(1, x)
  y0 <- 0.5 + 0.3 * x
  y1 <- y0 + cells$D + 0.4
  y <- ifelse(cells$post == 1, y1, y0)
  # The treated covariate mean lies outside the control convex hull. There is
  # variation within every D/period cell, but no finite IPT calibration root.
  fits <- list(
    function() drdid_imp_panel(y1, y0, cells$D, X),
    function() drdid_imp_rc(y, cells$post, cells$D, X),
    function() drdid_imp_rc1(y, cells$post, cells$D, X)
  )
  for (fit in fits) {
    error <- tryCatch(suppressWarnings(fit()), error = identity)
    expect_s3_class(error, "error")
    if (inherits(error, "error")) expect_identical(conditionMessage(error), required_ipt_error)
  }
})

test_that("all improved weighted helpers stop on a failed draw without redrawing", {
  n <- 4
  D <- c(0, 0, 1, 1)
  x <- c(-1, 1, 0, 1.5)
  X <- cbind(1, x)
  w <- rep(1, n)
  post <- c(0, 1, 0, 1)
  # The original weighted target has a finite root; seed 2's bootstrap treated
  # mean instead exceeds the largest control x, so that draw cannot be fitted.
  root <- c(0.5 * log(1 - 0.75^2), atanh(0.75))
  expect_true(DRDID:::ipt.valid(root, D, X, w))
  set.seed(2)
  v <- stats::rexp(n)
  expected_rng <- .Random.seed
  expect_gt(sum(v * D * x) / sum(v * D), max(x[D == 0]))
  fits <- list(
    function() DRDID:::wboot.dr.imp.panel(7, n, x, D, X, w),
    function() DRDID:::wboot_drdid_imp_rc(7, n, x, post, D, X, w),
    function() DRDID:::wboot_drdid_imp_rc1(7, n, x, post, D, X, w)
  )
  for (fit in fits) {
    set.seed(2)
    error <- tryCatch(suppressWarnings(fit()), error = identity)
    expect_s3_class(error, "error")
    if (inherits(error, "error")) {
      expect_identical(conditionMessage(error),
                       paste0(required_ipt_error, "\nWeighted-bootstrap replicate 7 failed."))
    }
    expect_identical(.Random.seed, expected_rng)
  }
})
