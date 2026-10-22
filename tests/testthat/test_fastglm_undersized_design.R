make_undersized_glm_probe <- function(native) {
  # Replace only the qualified native dependency in a private clone. Rejected
  # inputs must never enter it; no package namespace is modified.
  probe <- DRDID:::fastglm_fit
  replacements <- 0L
  replace_native <- function(node) {
    if (identical(node, quote(fastglm::fastglmPure))) {
      replacements <<- replacements + 1L
      return(as.name(".undersized_native_fit"))
    }
    if (is.call(node)) node <- as.call(lapply(node, replace_native))
    node
  }
  body(probe) <- replace_native(body(probe))
  stopifnot(replacements == 1L)
  private <- new.env(parent = environment(probe))
  private$.undersized_native_fit <- native
  environment(probe) <- private
  probe
}

test_that("undersized GLM designs are rejected before the native dependency", {
  original_body <- body(DRDID:::fastglm_fit)
  original_native <- fastglm::fastglmPure
  entered <- 0L
  probe <- make_undersized_glm_probe(function(...) {
    entered <<- entered + 1L
    list(native_entered = TRUE)
  })
  for (n in 0:3) {
    X <- matrix(1, n, 4)
    for (method in c(0L, 3L)) {
      expect_error(probe(X, rep(0, n), stats::gaussian(), rep(1, n), method),
                   "singular: fewer observations than covariates")
    }
  }
  expect_identical(entered, 0L)
  expect_identical(body(DRDID:::fastglm_fit), original_body)
  expect_identical(fastglm::fastglmPure, original_native)
})

test_that("square and taller supported GLMs preserve direct native fits", {
  set.seed(20261022)
  square <- diag(4)
  regular <- cbind(1, matrix(stats::rnorm(120), 40, 3))
  cases <- list(
    list(X = square, y = 1:4, family = stats::gaussian(), method = 0L),
    list(X = regular, y = stats::rnorm(40), family = stats::gaussian(), method = 0L),
    list(X = regular, y = rep(c(0, 1), 20), family = stats::binomial(), method = 3L)
  )
  for (case in cases) {
    n <- nrow(case$X)
    varying <- if (case$method == 3L) 1 + seq_len(n) %% 3 else
      seq(0.5, 2, length.out = n)
    weights <- list(rep(1, n), varying)
    if (n > ncol(case$X)) {
      weights[[3]] <- weights[[2]]
      weights[[3]][c(1, 3)] <- 0
    }
    for (w in weights) {
      expected <- fastglm::fastglmPure(case$X, case$y, case$family, w,
                                       method = case$method, tol = 1e-8, maxit = 100L)
      actual <- DRDID:::fastglm_fit(case$X, case$y, case$family, w, case$method)
      expect_identical(actual, expected)
    }
  }
})

test_that("each undersized repeated-cross-section outcome cell gives an R error", {
  support <- as.matrix(expand.grid(a = c(-1, 1), b = c(-1, 1), c = c(-1, 1)))
  X <- cbind(1, support[rep(1:8, 4), ])
  D <- rep(c(0, 0, 1, 1), each = 8)
  post <- rep(c(0, 1, 0, 1), each = 8)
  y <- 0.3 * X[, 2] + 0.2 * X[, 3] + D * post + seq_len(32) / 100
  for (cell in 1:4) {
    in_cell <- ((cell - 1) * 8 + 1):(cell * 8)
    keep <- setdiff(seq_len(32), in_cell[-(1:3)])
    # The overall propensity design and other three outcome designs remain
    # full rank. Only this Gaussian outcome fit has three rows/four columns.
    for (inffunc in c(FALSE, TRUE)) {
      expect_error(drdid_rc(y[keep], post[keep], D[keep], X[keep, ],
                            rep(1, length(keep)), inffunc = inffunc),
                   "singular: fewer observations than covariates")
    }
  }
})
