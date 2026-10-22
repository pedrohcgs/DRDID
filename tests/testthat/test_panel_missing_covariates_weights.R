# The helper fixture in test_panel_missing_outcomes.R uses a deliberately unused
# all-NA column: only variables needed for estimation should exclude a unit.
make_incomplete_panel_fixture <- function() {
  set.seed(623)
  n <- 120L
  x <- rnorm(n); w <- exp(.2 * x)
  D <- rep(0:1, each = n / 2L)
  y0 <- x + rnorm(n); y1 <- y0 + .4 + D + .2 * x + rnorm(n)
  data.frame(id = rep(100 + 7 * seq_len(n), each = 2L),
             t = rep(0:1, n), D = rep(D, each = 2L),
             y = as.vector(rbind(y0, y1)), x = rep(x, each = 2L),
             w = rep(w, each = 2L), unused = NA_real_)
}

fit_incomplete_panel <- function(data, method, boot = FALSE) {
  options <- switch(method, imp = list(fun = drdid),
                    trad = list(fun = drdid, estMethod = "trad"),
                    or = list(fun = ordid),
                    ipw = list(fun = ipwdid, normalized = FALSE),
                    std_ipw = list(fun = ipwdid, normalized = TRUE))
  fun <- options$fun; options$fun <- NULL
  do.call(fun, c(list(yname = "y", tname = "t", idname = "id", dname = "D",
                     xformla = ~x + I(x^2), data = data, weightsname = "w",
                     inffunc = TRUE, boot = boot, boot.type = "weighted", nboot = 9L), options))
}

test_that("panel covariate and weight missingness excludes matched units first", {
  for (case in c("pre_x", "post_x", "both_x", "pre_w", "post_w", "both_w", "mixed")) {
    d <- make_incomplete_panel_fixture()
    rows <- switch(case, pre_x = 3L, post_x = 4L, both_x = 3:4,
                   pre_w = 3L, post_w = 4L, both_w = 3:4, mixed = c(3L, 144L))
    if (case == "mixed") {
      d$x[3] <- NA; d$w[144] <- NA; d$y[83] <- NA
    } else d[rows, if (grepl("_x$", case)) "x" else "w"] <- NA
    # Include imbalance and a permutation so index-based matching cannot pass.
    d <- d[-9L, ]
    set.seed(912); d <- d[sample(nrow(d)), ]
    paired <- names(which(table(d$id) == 2L))
    incomplete <- d$id[!complete.cases(d[c("y", "x", "w")])]
    reference <- d[d$id %in% paired & !d$id %in% incomplete, ]
    for (method in c("imp", "trad", "or", "ipw", "std_ipw")) {
      actual <- expect_warning(fit_incomplete_panel(d, method), "missing.*affected panel unit")
      expected <- fit_incomplete_panel(reference, method)
      for (field in c("ATT", "se", "lci", "uci", "att.inf.func")) {
        expect_equal(unname(actual[[field]]), unname(expected[[field]]), tolerance = 1e-10,
                     info = paste(case, method, field))
      }
      expect_length(actual$att.inf.func, nrow(reference) / 2L)
    }
  }
})

test_that("retained observed covariate and weight changes are still rejected", {
  for (variable in c("x", "w")) {
    d <- make_incomplete_panel_fixture()
    d$x[3] <- NA_real_
    d[24, variable] <- d[24, variable] + .5
    for (method in c("imp", "trad", "or", "ipw", "std_ipw")) {
      expect_error(suppressWarnings(fit_incomplete_panel(d, method)),
                   if (variable == "x") "covariates should be time invariant" else "weights must be time invariant")
    }
  }
})

test_that("missing panel inputs preserve bootstrap draws and RNG after exclusion", {
  d <- make_incomplete_panel_fixture(); d$x[3] <- NA; d$w[144] <- NA
  excluded <- unique(d$id[!complete.cases(d[c("x", "w")])])
  reference <- d[!d$id %in% excluded, ]
  for (method in c("imp", "trad", "or", "ipw", "std_ipw")) {
    set.seed(71); actual <- suppressWarnings(fit_incomplete_panel(d, method, boot = TRUE)); rng <- .Random.seed
    set.seed(71); expected <- fit_incomplete_panel(reference, method, boot = TRUE)
    expect_identical(.Random.seed, rng)
    for (field in c("ATT", "se", "boots", "att.inf.func")) expect_equal(unname(actual[[field]]), unname(expected[[field]]), tolerance = 1e-10)
  }
})

test_that("missing panel inputs still trigger retained-group size checks", {
  for (variable in c("x", "w")) for (removed in list(1:55, 1:60, 1:120)) {
    d <- make_incomplete_panel_fixture()
    d[2L * removed - 1L, variable] <- NA_real_
    for (method in c("imp", "trad", "or", "ipw", "std_ipw")) {
      expect_error(suppressWarnings(fit_incomplete_panel(d, method)), "very small")
    }
  }
})
