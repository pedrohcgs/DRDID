bootstrap_count_cases <- function(zero = FALSE) {
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


test_that("invalid bootstrap counts fail before fitting or consuming RNG", {
  invalid <- list(0, -1, 1, 1.5, 2.5, Inf, -Inf, NA_real_, NaN,
                  numeric(), c(2, 3), "5", TRUE, 2 + 0i)
  for (case in bootstrap_count_cases()) for (bt in c("weighted", "multiplier")) {
    for (count in invalid) {
      set.seed(812)
      before <- .Random.seed
      # Missing data arguments deliberately ensure count validation runs first.
      expect_error(do.call(get(case$method), list(boot = TRUE, boot.type = bt,
        nboot = count)), "nboot must be a finite scalar integer of at least 2",
        info = case$method)
      expect_identical(.Random.seed, before)
    }
  }
})

test_that("valid bootstrap counts produce the requested draw count", {
  for (case in bootstrap_count_cases()) for (bt in c("weighted", "multiplier")) {
    for (count in c(2L, 5)) {
      set.seed(812)
      fit <- do.call(get(case$method), c(case$args,
        list(boot = TRUE, boot.type = bt, nboot = count)))
      expect_length(fit$boots, count)
    }
  }
  expect_identical(validate_bootstrap_count(NULL), 999)
  set.seed(812)
  case <- bootstrap_count_cases()[[1]]
  fit <- do.call(get(case$method), c(case$args,
    list(boot = TRUE, boot.type = "multiplier", nboot = NULL)))
  expect_length(fit$boots, 999)
})

test_that("inactive bootstrap ignores nboot as documented", {
  for (case in bootstrap_count_cases()) {
    expected <- do.call(get(case$method), c(case$args, list(boot = FALSE)))
    for (count in list(0, NA_real_, c(2, 3))) {
      set.seed(812)
      before <- .Random.seed
      actual <- do.call(get(case$method), c(case$args,
        list(boot = FALSE, nboot = count)))
      expect_identical(actual$ATT, expected$ATT)
      expect_identical(actual$se, expected$se)
      expect_identical(actual$lci, expected$lci)
      expect_identical(actual$uci, expected$uci)
      expect_null(actual$boots)
      expect_identical(.Random.seed, before)
    }
  }
})
