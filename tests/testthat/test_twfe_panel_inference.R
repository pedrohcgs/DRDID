context("TWFE panel inference uses paired individual observations")

make_twfe_inference_data <- function(n = 80, rho = 0.65) {
  set.seed(1439)
  X <- cbind(x = rnorm(n), z = rnorm(n))
  D <- rbinom(n, 1, plogis(0.2 + 0.3 * X[, 1] - 0.2 * X[, 2]))
  e0 <- rnorm(n)
  e1 <- rho * e0 + sqrt(1 - rho^2) * rnorm(n)
  y0 <- 1 + 0.7 * X[, 1] - 0.2 * X[, 2] + 0.4 * D + 2 * e0
  y1 <- 0.8 + 0.7 * X[, 1] - 0.2 * X[, 2] + 1.65 * D + 2 * e1
  list(y0 = y0, y1 = y1, D = D, X = X,
       w = exp(0.2 * X[, 1] + 0.1 * X[, 2]), n = n)
}

# With two balanced periods, common unit weights, and time-invariant additive
# covariates, the interaction coefficient is the weighted DID in outcome changes.
twfe_change_reference <- function(d, w = d$w) {
  dy <- d$y1 - d$y0
  a1 <- w * d$D
  a0 <- w * (1 - d$D)
  mu1 <- sum(a1 * dy) / sum(a1)
  mu0 <- sum(a0 * dy) / sum(a0)
  phi <- a1 * (dy - mu1) / mean(a1) - a0 * (dy - mu0) / mean(a0)
  list(ATT = mu1 - mu0, IF = phi,
       se = sqrt(sum((phi - mean(phi))^2)) / length(phi))
}

# Separate pooled-WLS construction retains the two row contributions, permitting
# checks of the old paired multiplier algorithm and within-unit covariance signs.
twfe_row_reference <- function(d, X = d$X) {
  n <- d$n
  dd <- rep(d$D, 2)
  post <- rep(0:1, each = n)
  Z <- cbind(1, dd, post, rbind(X, X), interaction = dd * post)
  w <- rep(d$w / mean(d$w), 2)
  fit <- stats::lm.wfit(Z, c(d$y0, d$y1), w)
  H <- crossprod(Z, w * Z) / (2 * n)
  row_IF <- as.vector((w * fit$residuals * Z) %*% solve(H)[, ncol(Z)])
  list(ATT = unname(tail(fit$coefficients, 1)), row_IF = row_IF,
       IF = (row_IF[seq_len(n)] + row_IF[n + seq_len(n)]) / 2)
}

test_that("constant within-group changes have zero panel TWFE sampling variance", {
  set.seed(901)
  n <- 120
  X <- cbind(x = rnorm(n), z = rnorm(n))
  D <- rbinom(n, 1, 0.45)
  y0 <- 3 * rnorm(n) + X[, 1] - 0.5 * D
  y1 <- y0 + 1.25 * D
  w <- exp(0.3 * X[, 1])
  fit <- twfe_did_panel(y1, y0, D, X, w, inffunc = TRUE)

  # Treating the two rows as independent previously gave an SE around 0.842.
  expect_equal(unname(fit$ATT), 1.25, tolerance = 1e-12)
  expect_length(fit$att.inf.func, n)
  expect_equal(fit$att.inf.func, rep(0, n), tolerance = 1e-10)
  expect_equal(fit$se, 0, tolerance = 1e-10)
  expect_equal(unname(fit$lci), 1.25, tolerance = 1e-10)
  expect_equal(unname(fit$uci), 1.25, tolerance = 1e-10)
})

test_that("panel TWFE influence values match weighted outcome-change DID", {
  for (rho in c(-0.8, 0.8)) {
    d <- make_twfe_inference_data(rho = rho)
    # These shapes avoid the separate intercept-plus-one-covariate input issue.
    designs <- list(matrix(d$X[, 1], ncol = 1), d$X)
    for (X in designs) {
      for (w in list(rep(1, d$n), d$w)) {
        fit <- twfe_did_panel(d$y1, d$y0, d$D, X, w, inffunc = TRUE)
        expected <- twfe_change_reference(d, w)
        expect_equal(unname(fit$ATT), expected$ATT, tolerance = 1e-11)
        expect_length(fit$att.inf.func, d$n)
        expect_equal(fit$att.inf.func, expected$IF, tolerance = 1e-10)
        expect_equal(mean(fit$att.inf.func), 0, tolerance = 1e-12)
        expect_equal(fit$se, expected$se, tolerance = 1e-11)
        expect_equal(fit$se, sd(fit$att.inf.func) * sqrt(d$n - 1) / d$n,
                     tolerance = 1e-12)
        expect_equal(as.numeric(c(fit$lci, fit$uci)),
                     expected$ATT + c(-1, 1) * 1.96 * expected$se, tolerance = 1e-11)
      }
    }
  }
})

test_that("panel TWFE influence values differentiate individual weights", {
  d <- make_twfe_inference_data(n = 48, rho = -0.6)
  h <- 1e-4
  for (X in list(matrix(d$X[, 1], ncol = 1), d$X)) {
    for (w in list(rep(1, d$n), d$w)) {
      fit <- twfe_did_panel(d$y1, d$y0, d$D, X, w, inffunc = TRUE)
      numerical_IF <- vapply(seq_len(d$n), function(i) {
        plus <- minus <- w
        plus[i] <- w[i] * (1 + h)
        minus[i] <- w[i] * (1 - h)
        a <- twfe_did_panel(d$y1, d$y0, d$D, X, plus)$ATT
        b <- twfe_did_panel(d$y1, d$y0, d$D, X, minus)$ATT
        as.numeric((a - b) / (h / (d$n + h) + h / (d$n - h)))
      }, numeric(1))
      expect_equal(fit$att.inf.func, numerical_IF, tolerance = 1e-6)
    }
  }
})

test_that("panel covariance accounts for both signs of serial outcome correlation", {
  d <- make_twfe_inference_data()
  d$w <- rep(1, d$n)
  e <- stats::lm.fit(cbind(1, d$D, d$X), sin(seq_len(d$n)))$residuals
  baseline <- 1 + 0.7 * d$X[, 1] - 0.2 * d$X[, 2] + 0.4 * d$D
  d$y0 <- baseline + e
  positive <- negative <- d
  positive$y1 <- baseline + 1.25 * d$D + e
  negative$y1 <- baseline + 1.25 * d$D - e
  fits <- lapply(list(positive, negative), function(a)
    twfe_did_panel(a$y1, a$y0, a$D, a$X, a$w, inffunc = TRUE))
  references <- lapply(list(positive, negative), twfe_row_reference)

  expect_equal(fits[[1]]$se, 0, tolerance = 1e-10)
  expect_gt(fits[[2]]$se, 0.05)
  expect_equal(fits[[2]]$se, twfe_change_reference(negative)$se, tolerance = 1e-11)
  row_covariance <- vapply(references, function(r)
    sum(r$row_IF[seq_len(d$n)] * r$row_IF[d$n + seq_len(d$n)]), numeric(1))
  expect_lt(row_covariance[1], 0)
  expect_gt(row_covariance[2], 0)
  # Marginal row variances alone cannot distinguish these cases.
  expect_equal(sum(references[[1]]$row_IF^2), sum(references[[2]]$row_IF^2),
               tolerance = 1e-10)
})

test_that("paired multiplier bootstrap draws retain the same fixed-seed values", {
  d <- make_twfe_inference_data()
  expected <- twfe_row_reference(d)
  B <- 99
  set.seed(902)
  paired_draws <- DRDID:::mboot.twfep.did(d$n, expected$row_IF, B)
  set.seed(902)
  fit <- twfe_did_panel(d$y1, d$y0, d$D, d$X, d$w, boot = TRUE,
                         boot.type = "multiplier", nboot = B, inffunc = TRUE)
  expected_se <- IQR(paired_draws) / (qnorm(0.75) - qnorm(0.25))

  expect_equal(unname(fit$ATT), expected$ATT, tolerance = 1e-11)
  expect_length(fit$att.inf.func, d$n)
  expect_equal(fit$att.inf.func, expected$IF, tolerance = 1e-10)
  expect_equal(as.vector(fit$boots), paired_draws, tolerance = 1e-11)
  expect_equal(fit$se, expected_se, tolerance = 1e-11)
})

test_that("weighted bootstrap still refits the same panel TWFE point estimator", {
  d <- make_twfe_inference_data()
  B <- 59
  set.seed(1291)
  expected_draws <- vapply(seq_len(B), function(b)
    twfe_change_reference(d, d$w * rexp(d$n))$ATT, numeric(1))
  for (X in list(matrix(d$X[, 1], ncol = 1), d$X)) {
    set.seed(1291)
    fit <- twfe_did_panel(d$y1, d$y0, d$D, X, d$w, boot = TRUE,
                           boot.type = "weighted", nboot = B, inffunc = TRUE)
    expect_equal(unname(fit$ATT), twfe_change_reference(d)$ATT, tolerance = 1e-11)
    expect_equal(as.vector(fit$boots), expected_draws, tolerance = 1e-11)
    expect_equal(fit$se, IQR(expected_draws) / (qnorm(0.75) - qnorm(0.25)),
                 tolerance = 1e-11)
    expect_length(fit$att.inf.func, d$n)
  }
})

test_that("the no-covariate panel fallback retains its existing inference", {
  d <- make_twfe_inference_data()
  expected <- reg_did_panel(d$y1, d$y0, d$D, NULL, d$w, inffunc = TRUE)
  for (X in list(NULL, rep(1, d$n), matrix(1, d$n, 1))) {
    fit <- twfe_did_panel(d$y1, d$y0, d$D, X, d$w, inffunc = TRUE)
    expect_equal(fit$ATT, expected$ATT, tolerance = 1e-11)
    expect_equal(fit$se, expected$se, tolerance = 1e-11)
    expect_equal(as.vector(fit$att.inf.func), as.vector(expected$att.inf.func),
                 tolerance = 1e-11)
  }
  for (method in c("weighted", "multiplier")) {
    set.seed(1813)
    expected <- reg_did_panel(d$y1, d$y0, d$D, NULL, d$w, inffunc = TRUE,
                               boot = TRUE, boot.type = method, nboot = 49)
    set.seed(1813)
    fit <- twfe_did_panel(d$y1, d$y0, d$D, NULL, d$w, inffunc = TRUE,
                           boot = TRUE, boot.type = method, nboot = 49)
    expect_equal(fit$ATT, expected$ATT, tolerance = 1e-11)
    expect_equal(fit$se, expected$se, tolerance = 1e-11)
    expect_equal(fit$boots, expected$boots, tolerance = 1e-11)
    expect_equal(as.vector(fit$att.inf.func), as.vector(expected$att.inf.func),
                 tolerance = 1e-11)
  }
  expect_null(twfe_did_panel(d$y1, d$y0, d$D, d$X, d$w)$att.inf.func)
})

test_that("direct panel TWFE calls reject incomplete individuals before fitting", {
  d <- make_twfe_inference_data()
  original <- list(y1 = d$y1, y0 = d$y0, D = d$D, covariates = d$X,
                   i.weights = d$w, inffunc = TRUE, nboot = 7)
  for (component in c("y0", "y1", "D", "i.weights", "covariates")) {
    args <- original
    # An outcome NA removes just one stacked row in lm; an X NA removes both
    # rows of that individual. Neither may silently change the pairing indices.
    args[[component]][3] <- NA_real_
    for (method in c("analytic", "weighted", "multiplier")) {
      args$boot <- method != "analytic"
      args$boot.type <- if (method == "analytic") "weighted" else method
      expect_error(do.call(twfe_did_panel, args), "Panel inputs must be complete")
    }
  }
})

test_that("zero-weight individuals cannot have missing panel outcomes or covariates", {
  d <- make_twfe_inference_data()
  original <- list(y1 = d$y1, y0 = d$y0, D = d$D, covariates = d$X,
                   i.weights = d$w, inffunc = TRUE)
  original$i.weights[3] <- 0
  for (component in c("y0", "y1", "covariates")) {
    args <- original
    args[[component]][3] <- NA_real_
    expect_error(do.call(twfe_did_panel, args), "Panel inputs must be complete")
  }
})

test_that("direct panel TWFE inputs must have one entry or row per individual", {
  d <- make_twfe_inference_data()
  original <- list(y1 = d$y1, y0 = d$y0, D = d$D, covariates = d$X,
                   i.weights = d$w, inffunc = TRUE)
  for (component in c("y0", "y1", "D")) {
    for (size in c(d$n - 1, d$n + 1)) {
      args <- original
      args[[component]] <- rep(original[[component]], length.out = size)
      expect_error(do.call(twfe_did_panel, args),
                    "y0 and y1 must each have one outcome per individual in D")
    }
  }
  for (size in c(d$n - 1, d$n + 1)) {
    args <- original
    args$i.weights <- rep(d$w, length.out = size)
    expect_error(do.call(twfe_did_panel, args),
                  "i.weights must have one weight per individual in D")

    args <- original
    args$covariates <- matrix(rep(d$X, length.out = 2 * size), nrow = size)
    expect_error(do.call(twfe_did_panel, args),
                  "covariates must have one row per individual in D")
    args$covariates <- rep(d$X[, 1], length.out = size)
    expect_error(do.call(twfe_did_panel, args),
                  "covariates must have one row per individual in D")
  }
})
