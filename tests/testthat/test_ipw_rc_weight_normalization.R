context("Observation-weight normalization in repeated-cross-section IPW inference")

make_ipw_rc_weight_atoms <- function() {
  # These counts represent an iid row law with independent w, D and post:
  # P(w = 0.5) = P(w = 1.5) = 1/2, P(D = 1) = 4/5, P(post = 1) = 7/10.
  cells <- expand.grid(w = c(0.5, 1.5), D = 0:1, post = 0:1)
  counts <- 100 * 0.5 * ifelse(cells$D == 1, 0.8, 0.2) *
    ifelse(cells$post == 1, 0.7, 0.3)
  d <- cells[rep(seq_len(nrow(cells)), round(counts)), ]
  rownames(d) <- NULL
  d$y <- d$D * d$post * (2 - 1 / d$w)
  d
}

fit_ipw_rc_weight_atoms <- function(d, ...) {
  ipw_did_rc(d$y, d$post, d$D, covariates = NULL,
             i.weights = d$w, inffunc = TRUE, ...)
}

ipw_rc_weight_atom_influence <- function(d) {
  # Comparison outcomes and pre-period outcomes are zero, so propensity-score
  # estimation has no effect. Differentiate S/A * M_post/L_post directly, where
  # S = E[w], A = E[w D], M_post = E[w D post Y], L_post = E[w post].
  d$w * (d$D * d$post * (2 - 1 / d$w) / (0.8 * 0.7) -
           d$D / 0.8 - d$post / 0.7 + 1)
}

test_that("unequal observation weights have the influence of the normalized functional", {
  d <- make_ipw_rc_weight_atoms()
  expected <- ipw_rc_weight_atom_influence(d)
  result <- fit_ipw_rc_weight_atoms(d)

  expect_equal(result$ATT, 1, tolerance = 1e-12)
  expect_equal(as.vector(result$att.inf.func), expected, tolerance = 1e-12)
  # This is the IF second moment under the iid row law represented above,
  # rather than a sampling variance conditional on fixed cell counts.
  expect_equal(mean(expected^2), 71 / 112, tolerance = 1e-12)
  expect_equal(result$se, sqrt((71 / 112) / nrow(d)), tolerance = 1e-12)

  d$id <- seq_len(nrow(d))
  wrapped <- ipwdid(yname = "y", tname = "post", idname = "id", dname = "D",
                    xformla = ~1, data = d, panel = FALSE, normalized = FALSE,
                    weightsname = "w", inffunc = TRUE)
  expect_equal(wrapped$ATT, result$ATT, tolerance = 1e-12)
  expect_equal(wrapped$se, result$se, tolerance = 1e-12)
})

test_that("zero-weight rows and common weight scaling preserve IPW inference", {
  d <- make_ipw_rc_weight_atoms()
  result <- fit_ipw_rc_weight_atoms(d)
  scaled <- d
  scaled$w <- 17.3 * scaled$w
  scaled_result <- fit_ipw_rc_weight_atoms(scaled)
  expect_equal(scaled_result$ATT, result$ATT, tolerance = 1e-12)
  expect_equal(scaled_result$se, result$se, tolerance = 1e-12)
  expect_equal(scaled_result$att.inf.func, result$att.inf.func, tolerance = 1e-12)

  extra <- rbind(d, data.frame(w = 0, D = 1, post = 1, y = 999))
  extra_result <- fit_ipw_rc_weight_atoms(extra)
  expect_equal(extra_result$ATT, result$ATT, tolerance = 1e-12)
  expect_equal(extra_result$se, result$se, tolerance = 1e-12)
  expect_equal(as.vector(extra_result$att.inf.func[nrow(extra)]), 0,
               tolerance = 1e-12)
  # IF/n is each observation's contribution; adding a zero-weight row changes n.
  expect_equal(as.vector(extra_result$att.inf.func[seq_len(nrow(d))]) / nrow(extra),
               as.vector(result$att.inf.func) / nrow(d), tolerance = 1e-12)
})

test_that("constant weights and zero ATT retain their analytic influences", {
  d <- make_ipw_rc_weight_atoms()
  constant <- d
  constant$w <- 1
  result <- fit_ipw_rc_weight_atoms(constant)
  theta <- 2 / 3
  expected <- d$D * d$post * (2 - 1 / d$w) / (0.8 * 0.7) -
    theta * (d$D / 0.8 + d$post / 0.7 - 1)
  expect_equal(result$ATT, theta, tolerance = 1e-12)
  expect_equal(as.vector(result$att.inf.func), expected, tolerance = 1e-12)
  constant$w <- rep(7, nrow(constant))
  expect_equal(fit_ipw_rc_weight_atoms(constant)$att.inf.func,
               result$att.inf.func, tolerance = 1e-12)

  # Equal treated outcome levels in both periods give ATT = 0 but a nonzero IF.
  d$y <- d$D
  zero <- fit_ipw_rc_weight_atoms(d)
  expected_zero <- d$w * (d$D / 0.8 - 1) *
    (d$post / 0.7 - (1 - d$post) / 0.3)
  expect_equal(zero$ATT, 0, tolerance = 1e-12)
  expect_gt(mean(expected_zero^2), 0)
  expect_equal(as.vector(zero$att.inf.func), expected_zero, tolerance = 1e-12)
})

test_that("weighted IPW influence agrees with fully refitted empirical derivatives", {
  set.seed(53)
  n <- 72
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  # The fitted logit in (1, x, z) is misspecified, and comparison outcomes are
  # nonzero, so this checks the propensity estimation effect as well.
  D <- stats::rbinom(n, 1, stats::plogis(-0.25 + 0.45 * x - 0.25 * z +
                                         0.9 * (x^2 - 1) + 0.5 * sin(2 * z)))
  post <- stats::rbinom(n, 1, 0.5)
  y <- 0.7 + 0.6 * x - 0.4 * z + 0.6 * D +
    post * (0.5 + 0.4 * x + (1.1 + 0.3 * x) * D) + stats::rnorm(n)
  w <- exp(0.5 * x - 0.2 * z)
  X <- cbind(1, x, z)
  result <- ipw_did_rc(y, post, D, X, w, inffunc = TRUE)
  # A multiplicative perturbation of row i is empirical contamination with
  # eps(h) = h/(n+h). Each call renormalizes weights and refits the propensity.
  h <- 3e-5
  eps <- function(hh) hh / (n + hh)
  numerical <- vapply(seq_len(n), function(i) {
    plus <- minus <- w
    plus[i] <- w[i] * (1 + h)
    minus[i] <- w[i] * (1 - h)
    (ipw_did_rc(y, post, D, X, plus)$ATT -
       ipw_did_rc(y, post, D, X, minus)$ATT) / (eps(h) - eps(-h))
  }, numeric(1))
  expect_equal(as.vector(result$att.inf.func), numerical, tolerance = 1e-6)
})

test_that("the multiplier bootstrap uses the normalized-functional influence", {
  d <- make_ipw_rc_weight_atoms()
  expected <- ipw_rc_weight_atom_influence(d)
  set.seed(908)
  expected_draws <- DRDID:::mboot.did(expected, nboot = 49)
  set.seed(908)
  result <- fit_ipw_rc_weight_atoms(d, boot = TRUE, boot.type = "multiplier",
                                    nboot = 49)
  expect_equal(result$boots, expected_draws, tolerance = 1e-12)
  expected_se <- stats::IQR(expected_draws) /
    (stats::qnorm(0.75) - stats::qnorm(0.25))
  expect_equal(result$se, expected_se, tolerance = 1e-12)
})
