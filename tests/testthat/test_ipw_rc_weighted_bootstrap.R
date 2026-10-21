context("Weighted bootstrap refits the repeated-cross-section IPW estimator")

make_ipw_rc_bootstrap_atoms <- function() {
  # A saturated binary propensity model with population probabilities .2 and .8.
  # Repeated atoms keep both retained and trimmed controls away from threshold .5.
  cells <- expand.grid(x = 0:1, D = 0:1, post = 0:1)
  prob_D <- ifelse(cells$x == 1, 0.8, 0.2)
  counts <- 100 * 0.5 * ifelse(cells$D == 1, prob_D, 1 - prob_D) *
    ifelse(cells$post == 1, 0.4, 0.6)
  d <- cells[rep(seq_len(nrow(cells)), round(counts)), ]
  rownames(d) <- NULL
  d$id <- seq_len(nrow(d))
  d$y <- 1 + 0.7 * d$x + 0.4 * d$D +
    d$post * (0.5 + 0.2 * d$x + 1.1 * d$D)
  d$w <- exp(0.2 * d$x + 0.1 * d$D - 0.15 * d$post)
  d
}

ipw_rc_bootstrap_refits <- function(d, X, w, nboot, trim.level) {
  # Reference each draw by the public point estimator, including its weight
  # normalization, nuisance refit, and trimming. This does not duplicate the
  # weighted-bootstrap helper's algebra.
  if (is.null(w)) w <- rep(1, nrow(d))
  vapply(seq_len(nboot), function(i) {
    ipw_did_rc(d$y, d$post, d$D, X, i.weights = w * stats::rexp(nrow(d)),
               trim.level = trim.level)$ATT
  }, numeric(1))
}

test_that("a unit-weight bootstrap draw agrees with an independent raw-moment target", {
  d <- expand.grid(D = 0:1, post = 0:1, replicate = 1:2)
  d$y <- d$D * d$post
  n <- nrow(d)
  set.seed(20261024)
  v <- stats::rexp(n)
  expected_rng <- .Random.seed
  # Controls and pre-period outcomes are zero, eliminating propensity nuisance
  # terms. The raw-weight functional is S * M_post / (A * L_post).
  expected <- sum(v) * sum(v * d$D * d$post) /
    (sum(v * d$D) * sum(v * d$post))
  set.seed(20261024)
  result <- DRDID:::wboot_ipw_rc(1, n, d$y, d$post, d$D,
                                 matrix(1, n, 1), rep(1, n))
  expect_equal(result, expected, tolerance = 1e-12)
  expect_identical(.Random.seed, expected_rng)
  expect_equal(ipw_did_rc(d$y, d$post, d$D, NULL, i.weights = v)$ATT,
               expected, tolerance = 1e-12)

  # A common scale change to the original weights must not change the refit.
  set.seed(20261024)
  scaled <- DRDID:::wboot_ipw_rc(1, n, d$y, d$post, d$D,
                                 matrix(1, n, 1), rep(7, n))
  expect_equal(scaled, expected, tolerance = 1e-12)
})

test_that("weighted bootstrap helper draws equal direct point-estimator refits", {
  d <- make_ipw_rc_bootstrap_atoms()
  zero_w <- d$w
  zero_w[c(1, 10, 30, 70)] <- 0
  weight_cases <- list(NULL, d$w, zero_w)
  n <- nrow(d)
  B <- 7
  for (X in list(NULL, cbind(1, d$x))) {
    # Active, separated control trimming is checked for the saturated model.
    trim_cases <- if (is.null(X)) 0.995 else c(0.995, 0.5)
    for (trim in trim_cases) {
      for (w in weight_cases) {
        set.seed(908)
        expected <- ipw_rc_bootstrap_refits(d, X, w, B, trim)
        expected_rng <- .Random.seed
        q <- if (is.null(w)) rep(1, n) else w / mean(w)
        int_cov <- if (is.null(X)) matrix(1, n, 1) else X
        set.seed(908)
        draws <- vapply(seq_len(B), DRDID:::wboot_ipw_rc, numeric(1),
                         n = n, y = d$y, post = d$post, D = d$D,
                         int.cov = int_cov, i.weights = q, trim.level = trim)
        expect_equal(draws, expected, tolerance = 1e-10)
        expect_identical(.Random.seed, expected_rng)
      }
    }
  }
})

test_that("the trimming fixture contains separated retained and excluded controls", {
  d <- make_ipw_rc_bootstrap_atoms()
  d$w[c(1, 10, 30, 70)] <- 0
  set.seed(908)
  probabilities <- vapply(seq_len(19), function(i) {
    b <- d$w * stats::rexp(nrow(d))
    # For a saturated binary logit, each fitted probability is exactly the
    # weighted treated share within x, providing an independent membership check.
    vapply(0:1, function(x) {
      rows <- d$x == x
      sum(b[rows] * d$D[rows]) / sum(b[rows])
    }, numeric(1))
  }, numeric(2))
  expect_lt(max(probabilities[1, ]), 0.4)
  expect_gt(min(probabilities[2, ]), 0.6)
})

test_that("public weighted-bootstrap draws, uncertainty and RNG match direct refits", {
  d <- make_ipw_rc_bootstrap_atoms()
  n <- nrow(d)
  B <- 19
  for (with_covariates in c(FALSE, TRUE)) {
    X <- if (with_covariates) cbind(1, d$x) else NULL
    trim <- if (with_covariates) 0.5 else 0.995
    w <- if (with_covariates) d$w else NULL
    if (with_covariates) w[c(1, 10, 30, 70)] <- 0
    set.seed(908)
    expected <- ipw_rc_bootstrap_refits(d, X, w, B, trim)
    expected_rng <- .Random.seed
    point <- ipw_did_rc(d$y, d$post, d$D, X, w, trim.level = trim)
    expected_se <- stats::IQR(expected - point$ATT) /
      (stats::qnorm(0.75) - stats::qnorm(0.25))
    expected_cv <- stats::quantile(abs((expected - point$ATT) / expected_se), 0.95)
    expected_ci <- point$ATT + c(-1, 1) * expected_cv * expected_se

    set.seed(908)
    direct <- ipw_did_rc(d$y, d$post, d$D, X, w, boot = TRUE,
                         boot.type = "weighted", nboot = B, trim.level = trim)
    direct_rng <- .Random.seed
    wrapped_data <- d
    wrapped_data$w <- if (is.null(w)) rep(1, n) else w
    set.seed(908)
    wrapped <- ipwdid(yname = "y", tname = "post", idname = "id", dname = "D",
                      xformla = if (with_covariates) ~x else ~1,
                      data = wrapped_data, panel = FALSE, normalized = FALSE,
                      weightsname = "w", boot = TRUE, boot.type = "weighted",
                      nboot = B, trim.level = trim)
    expect_identical(direct_rng, expected_rng)
    expect_identical(.Random.seed, expected_rng)
    for (result in list(direct, wrapped)) {
      expect_equal(result$ATT, point$ATT, tolerance = 1e-12)
      expect_equal(result$boots, expected, tolerance = 1e-10)
      expect_equal(result$se, expected_se, tolerance = 1e-10)
      expect_equal(as.numeric(c(result$lci, result$uci)),
                   as.numeric(expected_ci), tolerance = 1e-10)
    }
  }
})
