context("Panel outcome missingness preserves unit pairs")

make_panel_missing_data <- function(n = 40, stochastic = FALSE) {
  id <- seq_len(n)
  D <- as.numeric(id > n / 2)
  x <- sin(id)
  w <- 1 + (id %% 5) / 5
  y0 <- 100 * id
  y1 <- y0 + 2 * D
  if (stochastic) {
    set.seed(20261023)
    x <- rnorm(n)
    w <- exp(0.3 * x + 0.2 * rnorm(n))
    y0 <- 3 + x + rnorm(n)
    y1 <- y0 + 0.6 + 0.3 * x + 2 * D + rnorm(n)
  }
  data.frame(id = rep(id, 2), t = rep(0:1, each = n), D = rep(D, 2),
             y = c(y0, y1), x = rep(x, 2), w = rep(w, 2), unused = NA_real_)
}

panel_missing_methods <- c("dr_default", "dr_trad", "or", "ipw", "std_ipw")

fit_panel_missing <- function(data, method, panel = TRUE, xformla = ~1,
                              weightsname = NULL) {
  specification <- switch(method,
    dr_default = list(fun = drdid),
    dr_trad = list(fun = drdid, estMethod = "trad"),
    or = list(fun = ordid),
    ipw = list(fun = ipwdid, normalized = FALSE),
    std_ipw = list(fun = ipwdid, normalized = TRUE))
  fun <- specification$fun
  specification$fun <- NULL
  do.call(fun, c(list(yname = "y", tname = "t", idname = "id", dname = "D",
                     data = data, panel = panel, xformla = xformla,
                     weightsname = weightsname, inffunc = TRUE), specification))
}

complete_panel_outcome_units <- function(data) {
  paired_ids <- intersect(data$id[data$t == 0], data$id[data$t == 1])
  retained_ids <- setdiff(paired_ids, data$id[is.na(data$y)])
  data[data$id %in% retained_ids, ]
}

expect_same_panel_missing_result <- function(actual, expected) {
  expect_equal(actual$ATT, expected$ATT, tolerance = 1e-9)
  expect_equal(actual$se, expected$se, tolerance = 1e-9)
  expect_equal(as.vector(actual$att.inf.func), as.vector(expected$att.inf.func),
               tolerance = 1e-8)
}

test_that("missing outcomes cannot pair outcomes belonging to different units", {
  d <- make_panel_missing_data()
  d$y[d$id == 1 & d$t == 0] <- NA_real_
  d$y[d$id == 21 & d$t == 1] <- NA_real_
  complete <- complete_panel_outcome_units(d)

  # Dropping just the two rows previously produced ATT = 102 instead of 2.
  for (method in panel_missing_methods) {
    actual <- expect_warning(fit_panel_missing(d, method), "Some outcomes are missing")
    expected <- fit_panel_missing(complete, method)
    expect_equal(actual$ATT, 2, tolerance = 1e-9)
    expect_same_panel_missing_result(actual, expected)
    expect_length(actual$att.inf.func, 38)
  }
})

test_that("panel filtering keeps outcomes, covariates, and weights aligned", {
  cases <- list(
    equal_counts = list(pre = 4, post = 37),
    unequal_counts = list(pre = c(4, 12), post = 37),
    one_period_only = list(pre = 4, post = integer(0)),
    both_periods = list(pre = 4, post = 4),
    initially_unbalanced = list(pre = 4, post = 37)
  )

  for (case in names(cases)) {
    d <- make_panel_missing_data(60)
    missing <- cases[[case]]
    d$y[d$id %in% missing$pre & d$t == 0] <- NA_real_
    d$y[d$id %in% missing$post & d$t == 1] <- NA_real_
    if (case == "initially_unbalanced") d <- d[!(d$id == 7 & d$t == 0), ]
    d$id <- 100 + 7 * d$id
    set.seed(735)
    d <- d[sample(nrow(d)), ]
    complete <- complete_panel_outcome_units(d)
    pre <- complete[complete$t == 0, ]
    pre <- pre[order(pre$id), ]
    post <- complete[complete$t == 1, ]
    post <- post[order(post$id), ]

    missing_units <- length(unique(d$id[is.na(d$y)]))
    missing_warning <- paste0("Some outcomes are missing.*for ", missing_units,
                              " affected panel unit")
    dp <- expect_warning(DRDID:::pre_process_drdid(
      yname = "y", tname = "t", idname = "id", dname = "D", xformla = ~x,
      weightsname = "w", data = d, panel = TRUE), missing_warning)
    expect_equal(dp$y0, pre$y)
    expect_equal(dp$y1, post$y)
    expect_equal(dp$D, post$D)
    expected_X <- model.matrix(~x, post)
    expect_equal(dim(dp$covariates), dim(expected_X))
    expect_equal(as.vector(dp$covariates), as.vector(expected_X))
    expect_equal(dp$i.weights, post$w / mean(post$w))
    expect_equal(dp$y1 - dp$y0, 2 * dp$D)
    expect_length(dp$D, nrow(post))
  }
})

test_that("weighted panel inference matches explicitly retained complete units", {
  d <- make_panel_missing_data(120, stochastic = TRUE)
  d$y[d$id == 2 & d$t == 0] <- NA_real_
  d$y[d$id %in% c(17, 83) & d$t == 1] <- NA_real_
  complete <- complete_panel_outcome_units(d)
  for (method in panel_missing_methods) {
    actual <- expect_warning(fit_panel_missing(d, method, xformla = ~x, weightsname = "w"),
                              "Some outcomes are missing")
    expected <- fit_panel_missing(complete, method, xformla = ~x, weightsname = "w")
    expect_gt(expected$se, 0.01)
    expect_same_panel_missing_result(actual, expected)
    expect_length(actual$att.inf.func, 117)
  }
})

test_that("panel group-size checks use units retained after missing outcomes", {
  # Cover a small remaining group, either group disappearing, and no units left.
  for (removed in list(1:16, 1:20, 21:40, 1:40)) {
    d <- make_panel_missing_data()
    d$y[d$id %in% removed & d$t == 0] <- NA_real_
    for (method in panel_missing_methods) {
      expect_error(suppressWarnings(fit_panel_missing(d, method)), "very small")
    }
  }
})

test_that("panel designs with no independent columns are rejected before fitting", {
  d <- make_panel_missing_data()
  d$zero <- 0
  for (xformla in list(~0, ~0 + zero)) {
    for (method in panel_missing_methods) {
      expect_error(fit_panel_missing(d, method, xformla = xformla),
                    "covariate matrix has no linearly independent columns")
    }
  }
})

test_that("complete panel data still match the corresponding low-level estimators", {
  d <- make_panel_missing_data(120, stochastic = TRUE)
  pre <- d[d$t == 0, ]
  post <- d[d$t == 1, ]
  low_level <- list(dr_default = drdid_imp_panel, dr_trad = drdid_panel,
                    or = reg_did_panel, ipw = ipw_did_panel,
                    std_ipw = std_ipw_did_panel)
  for (method in panel_missing_methods) {
    actual <- fit_panel_missing(d, method, xformla = ~x, weightsname = "w")
    expected <- low_level[[method]](y1 = post$y, y0 = pre$y, D = post$D,
      covariates = model.matrix(~x, post), i.weights = post$w, inffunc = TRUE)
    expect_same_panel_missing_result(actual, expected)
    expect_length(actual$att.inf.func, 120)
  }
})

test_that("repeated cross sections still drop missing outcome rows individually", {
  d <- make_panel_missing_data(120, stochastic = TRUE)
  d$y[d$id == 2 & d$t == 0] <- NA_real_
  d$y[d$id %in% c(17, 83) & d$t == 1] <- NA_real_
  complete_rows <- d[!is.na(d$y), ]
  for (method in panel_missing_methods) {
    actual <- expect_warning(fit_panel_missing(d, method, panel = FALSE,
      xformla = ~x, weightsname = "w"), "Some outcomes are missing")
    expected <- fit_panel_missing(complete_rows, method, panel = FALSE,
                                  xformla = ~x, weightsname = "w")
    expect_same_panel_missing_result(actual, expected)
    expect_length(actual$att.inf.func, 237)
  }
})
