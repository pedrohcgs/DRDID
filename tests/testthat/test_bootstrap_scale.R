bootstrap_scale_cases <- function(zero = FALSE) {
  d <- expand.grid(j = 1:12, post = 0:1, D = 0:1)
  d$y <- if (zero) 0 else 0.3 * d$post + d$D * d$post + sin(d$j)
  panel <- d[d$post == 1, ]
  c(setNames(lapply(c("drdid_panel", "drdid_imp_panel", "reg_did_panel",
                     "ipw_did_panel", "std_ipw_did_panel", "twfe_did_panel"),
                   function(m) list(method = m, args = list(y1 = panel$y,
                     y0 = rep(0, nrow(panel)), D = panel$D, covariates = NULL))),
             paste0("panel", 1:6)),
    setNames(lapply(c("drdid_rc", "drdid_rc1", "drdid_imp_rc", "drdid_imp_rc1",
                     "reg_did_rc", "ipw_did_rc", "std_ipw_did_rc", "twfe_did_rc"),
                   function(m) list(method = m, args = list(y = d$y, post = d$post,
                     D = d$D, covariates = NULL))), paste0("rc", 1:8)))
}

test_that("all estimators return point intervals for exactly zero bootstrap errors", {
  for (case in bootstrap_scale_cases(TRUE)) for (bt in c("weighted", "multiplier")) {
    set.seed(710)
    fit <- do.call(get(case$method), c(case$args,
      list(boot = TRUE, boot.type = bt, nboot = 19, inffunc = TRUE)))
    expect_equal(unname(fit$ATT), 0, info = case$method)
    expect_identical(fit$se, 0, info = case$method)
    expect_equal(as.numeric(c(fit$lci, fit$uci)), c(0, 0), info = case$method)
    expect_true(all(fit$boots == 0), info = case$method)
    expect_length(fit$boots, 19)
  }
})

test_that("nondegenerate sparse multiplier errors with zero IQR fail clearly", {
  D <- rep(0:1, each = 20)
  y <- rep(0, 40)
  y[21:22] <- c(1, -1)
  fit <- reg_did_panel(y, rep(0, 40), D, NULL, inffunc = TRUE)
  expect_equal(fit$se, sqrt(0.005))
  set.seed(710)
  draws <- mboot.did(fit$att.inf.func, 999)
  expect_identical(IQR(draws), 0)
  expect_gt(sd(draws), 0)
  for (m in c("drdid_panel", "drdid_imp_panel", "reg_did_panel",
              "ipw_did_panel", "std_ipw_did_panel", "twfe_did_panel")) {
    set.seed(710)
    expect_error(do.call(get(m), list(y1 = y, y0 = rep(0, 40), D = D,
      covariates = NULL, boot = TRUE, boot.type = "multiplier", nboot = 999)),
      "zero IQR but are not all zero", info = m)
  }
})

test_that("positive-scale summaries retain the existing CI arithmetic", {
  for (case in bootstrap_scale_cases()) for (bt in c("weighted", "multiplier")) {
    set.seed(812)
    fit <- do.call(get(case$method), c(case$args,
      list(boot = TRUE, boot.type = bt, nboot = 29)))
    errors <- if (bt == "weighted") fit$boots - fit$ATT else fit$boots
    expected_se <- IQR(errors) / (qnorm(0.75) - qnorm(0.25))
    expected_cv <- quantile(abs(errors / expected_se), 0.95)
    expect_gt(expected_se, 0)
    expect_identical(fit$se, expected_se)
    expect_identical(fit$lci, fit$ATT - expected_cv * expected_se)
    expect_identical(fit$uci, fit$ATT + expected_cv * expected_se)
  }
})

test_that("nonfinite errors and nonzero constant errors are never hidden", {
  for (errors in list(c(0, NA_real_), c(0, NaN), c(0, Inf), c(0, -Inf))) {
    expect_error(bootstrap_inference(errors, 1), "finite bootstrap errors")
  }
  expect_error(bootstrap_inference(c(0, 1), Inf), "finite estimate")
  expect_error(bootstrap_inference(rep(2, 10), 1), "zero IQR but are not all zero")
  expect_identical(bootstrap_inference(rep(0, 10), 2), list(se = 0, lci = 2, uci = 2))
})


test_that("overflow in otherwise finite bootstrap errors fails explicitly", {
  expect_error(bootstrap_inference(c(seq(-1e-300, 1e-300, length.out = 19), 1e300), 0),
               "non-finite critical value or interval")
  expect_error(bootstrap_inference(rep(c(-1e308, 1e308), each = 10), 0),
               "non-finite IQR scale")
})
