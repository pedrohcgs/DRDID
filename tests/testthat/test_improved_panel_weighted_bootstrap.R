context("Improved panel weighted bootstrap refits the public statistic")

make_trimmed_panel_atoms <- function(no_intercept = FALSE) {
  cells <- expand.grid(x = c(-1, 0, 1), D = 0:1)
  counts <- ifelse(cells$D == 1, c(4, 10, 16), c(16, 10, 4))
  if (no_intercept) {
    counts[cells$x == 0] <- ifelse(cells$D[cells$x == 0] == 1, 12, 8)
  }
  d <- cells[rep(seq_len(nrow(cells)), counts), ]
  rownames(d) <- NULL
  d$id <- seq_len(nrow(d))
  d$dy <- d$x^2 + d$D
  d$w <- exp(0.3 * d$x + 0.1 * d$D)
  d
}

fixed_panel_bootstrap_weights <- function(helper, weights) {
  # A private function copy fixes the multiplier vector for the exact target.
  # Actual propensity and outcome fitting remain unchanged; no namespace is mocked.
  replacements <- 0L
  replace_rexp <- function(node) {
    if (identical(node, quote(stats::rexp))) {
      replacements <<- replacements + 1L
      return(as.name(".fixed_panel_multipliers"))
    }
    if (is.call(node)) node <- as.call(lapply(node, replace_rexp))
    node
  }
  body(helper) <- replace_rexp(body(helper))
  stopifnot(replacements == 1L)
  private <- new.env(parent = environment(helper))
  private$.fixed_panel_multipliers <- function(n) {
    stopifnot(length(weights) == n)
    weights
  }
  environment(helper) <- private
  helper
}

panel_bootstrap_refits <- function(d, X, weights, multipliers, trim.level) {
  vapply(seq_len(ncol(multipliers)), function(j) {
    drdid_imp_panel(d$dy, rep(0, nrow(d)), d$D, X,
                    i.weights = weights * multipliers[, j],
                    trim.level = trim.level)$ATT
  }, numeric(1))
}

test_that("trimmed improved panel draws have an independently soluble target", {
  d <- make_trimmed_panel_atoms()
  n <- nrow(d)
  X <- cbind(1, d$x)
  # At unit weights the IPT root is (0, log(4)), giving p=(.2,.5,.8).
  # Odds-WLS coefficients for x^2 are (32/57, 5/19). The treated residual
  # sum is 30; after trimming x=1 controls, the control residual sum is -160/57.
  # Both residual contributions therefore use the treated weight sum 30.
  expected <- 187 / 171
  point <- drdid_imp_panel(d$dy, rep(0, n), d$D, X, trim.level = 0.6)
  helper <- fixed_panel_bootstrap_weights(DRDID:::wboot.dr.imp.panel, rep(1, n))
  draw <- helper(1, n, d$dy, d$D, X, rep(1, n), trim.level = 0.6)
  expect_equal(point$ATT, expected, tolerance = 1e-10)
  expect_equal(draw, expected, tolerance = 1e-10)
  ps <- DRDID:::pscore.cal(d$D, X, rep(1, n), n)$pscore
  expect_gt(min(abs(ps[d$D == 0] - 0.6)), 0.09)
  expect_true(any(ps[d$D == 0] < 0.6) && any(ps[d$D == 0] > 0.6))
})

test_that("a supplied no-intercept design refits its public panel statistic", {
  d <- make_trimmed_panel_atoms(no_intercept = TRUE)
  n <- nrow(d)
  X <- matrix(d$x, ncol = 1)
  # Without an intercept, IPT need not equate the two weight sums, even
  # without trimming. Here they are 32 treated versus 28 control-odds weights.
  # Direct unnormalized residual sums are 44.8 and 12.8, giving ATT=1.
  point <- drdid_imp_panel(d$dy, rep(0, n), d$D, X, trim.level = 0.995)
  helper <- fixed_panel_bootstrap_weights(DRDID:::wboot.dr.imp.panel, rep(1, n))
  draw <- helper(1, n, d$dy, d$D, X, rep(1, n), trim.level = 0.995)
  expect_equal(point$ATT, 1, tolerance = 1e-10)
  expect_equal(draw, 1, tolerance = 1e-10)
})

test_that("matched multipliers refit panel ATT with unit, positive and zero weights", {
  d <- make_trimmed_panel_atoms()
  n <- nrow(d)
  zero_weights <- d$w
  zero_weights[unlist(lapply(split(seq_len(n), interaction(d$x, d$D)),
                              utils::head, 1))] <- 0
  designs <- list(cbind(1, d$x), matrix(d$x, ncol = 1), NULL)
  B <- 7
  for (X in designs) {
    int_cov <- if (is.null(X)) matrix(1, n, 1) else X
    for (weights in list(rep(1, n), d$w, zero_weights)) {
      weights <- weights / mean(weights)
      for (trim in if (is.null(X)) c(0.995, 0.4) else c(0.995, 0.6)) {
        set.seed(121)
        multipliers <- matrix(stats::rexp(n * B), nrow = n)
        expected_rng <- .Random.seed
        expected <- panel_bootstrap_refits(d, X, weights, multipliers, trim)
        # Keep the fixture away from threshold rounding ambiguity. Trimming
        # membership may change between draws, exactly as in direct refitting.
        margins <- vapply(seq_len(B), function(j) {
          ps <- DRDID:::pscore.cal(d$D, int_cov, weights * multipliers[, j], n)$pscore
          min(abs(ps[d$D == 0] - trim))
        }, numeric(1))
        expect_gt(min(margins), 1e-3)
        set.seed(121)
        draws <- vapply(seq_len(B), DRDID:::wboot.dr.imp.panel, numeric(1),
                         n = n, deltaY = d$dy, D = d$D, int.cov = int_cov,
                         i.weights = weights, trim.level = trim)
        expect_equal(draws, expected, tolerance = 1e-10)
        expect_identical(.Random.seed, expected_rng)
      }
    }
  }
})

test_that("public improved panel bootstrap draws and uncertainty match direct refits", {
  d <- make_trimmed_panel_atoms()
  n <- nrow(d)
  X <- cbind(1, d$x)
  zero_weights <- d$w
  zero_weights[unlist(lapply(split(seq_len(n), interaction(d$x, d$D)),
                              utils::head, 1))] <- 0
  B <- 19
  for (weights in list(rep(1, n), d$w, zero_weights)) {
    q <- weights / mean(weights)
    for (trim in c(0.995, 0.6)) {
      point <- drdid_imp_panel(d$dy, rep(0, n), d$D, X, weights, trim.level = trim)
      set.seed(908)
      multipliers <- matrix(stats::rexp(n * B), nrow = n)
      expected_rng <- .Random.seed
      expected <- panel_bootstrap_refits(d, X, q, multipliers, trim)
      expected_se <- stats::IQR(expected - point$ATT) /
        (stats::qnorm(0.75) - stats::qnorm(0.25))
      expected_cv <- stats::quantile(abs((expected - point$ATT) / expected_se), 0.95)
      expected_ci <- point$ATT + c(-1, 1) * expected_cv * expected_se

      set.seed(908)
      direct <- drdid_imp_panel(d$dy, rep(0, n), d$D, X, weights,
                                boot = TRUE, boot.type = "weighted", nboot = B,
                                trim.level = trim)
      direct_rng <- .Random.seed
      long <- data.frame(id = rep(d$id, 2), post = rep(0:1, each = n),
                          D = rep(d$D, 2), x = rep(d$x, 2),
                          y = c(rep(0, n), d$dy), w = rep(weights, 2))
      set.seed(908)
      wrapped <- drdid(yname = "y", tname = "post", idname = "id", dname = "D",
                       xformla = ~x, data = long, panel = TRUE, estMethod = "imp",
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
  }
})
