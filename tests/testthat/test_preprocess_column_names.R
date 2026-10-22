column_name_methods <- function() {
  list(list(fun = drdid), list(fun = drdid, estMethod = "trad"),
       list(fun = ipwdid), list(fun = ipwdid, normalized = FALSE), list(fun = ordid))
}

make_column_name_data <- function(panel = TRUE) {
  set.seed(418)
  n <- 80L
  xx <- rep(runif(n, -.8, .8), each = 2L)
  gg <- rep(rep(0:1, each = n / 2L), each = 2L)
  tt <- rep(0:1, n)
  data.frame(ii = if (panel) rep(seq_len(n), each = 2L) else seq_len(2L * n),
             tt = tt, gg = gg, xx = xx,
             yy = .4 * xx + tt * (1 + gg + .3 * xx^2) + rnorm(2L * n),
             ww = exp(.3 * xx))
}

column_name_args <- function(data, panel, weighted) {
  list(yname = "yy", tname = "tt", idname = "ii", dname = "gg",
       xformla = ~xx + I(xx^2), data = data, panel = panel,
       weightsname = if (weighted) "ww" else NULL, inffunc = TRUE)
}

test_that("public estimates preserve outcome and transformed covariate aliases", {
  for (panel in c(TRUE, FALSE)) for (weighted in c(FALSE, TRUE)) {
    args <- column_name_args(make_column_name_data(panel), panel, weighted)
    for (specification in column_name_methods()) {
      fun <- specification$fun; specification$fun <- NULL
      reference <- do.call(fun, c(args, specification))
      for (role in c("yy", "xx")) for (alias in c("y", "D", "post", "w")) {
        renamed <- args
        names(renamed$data)[names(renamed$data) == role] <- alias
        if (role == "yy") renamed$yname <- alias else {
          renamed$xformla <- stats::as.formula(paste0("~", alias, "+I(", alias, "^2)"))
        }
        actual <- do.call(fun, c(renamed, specification))
        for (field in c("ATT", "se", "lci", "uci", "att.inf.func")) {
          expect_equal(actual[[field]], reference[[field]], tolerance = 1e-12,
                       info = paste(panel, weighted, role, alias, field))
        }
      }
    }
  }
})

test_that("preprocessing keeps role names and formula inputs on original data", {
  args <- column_name_args(make_column_name_data(), TRUE, TRUE)
  reference <- do.call(pre_process_drdid, args)
  for (role in c("yy", "tt", "ii", "gg", "ww", "xx")) {
    for (alias in c("y", "D", "post", "w")) {
      renamed <- args
      names(renamed$data)[names(renamed$data) == role] <- alias
      argument <- c(yy = "yname", tt = "tname", ii = "idname",
                    gg = "dname", ww = "weightsname")
      if (role %in% names(argument)) renamed[[argument[[role]]]] <- alias
      if (role == "xx") renamed$xformla <- stats::as.formula(paste0("~", alias, "+I(", alias, "^2)"))
      actual <- do.call(pre_process_drdid, renamed)
      for (field in c("y1", "y0", "D", "i.weights")) expect_identical(actual[[field]], reference[[field]])
      expect_equal(unname(actual$covariates), unname(reference$covariates), tolerance = 0)
    }
  }
})

test_that("alias renaming preserves weighted-bootstrap draws and RNG", {
  for (panel in c(TRUE, FALSE)) {
    args <- column_name_args(make_column_name_data(panel), panel, TRUE)
    args$boot <- TRUE; args$boot.type <- "weighted"; args$nboot <- 9L
    renamed <- args
    names(renamed$data)[names(renamed$data) == "yy"] <- "w"
    names(renamed$data)[names(renamed$data) == "xx"] <- "y"
    renamed$yname <- "w"; renamed$xformla <- ~y + I(y^2)
    for (specification in column_name_methods()) {
      fun <- specification$fun; specification$fun <- NULL
      set.seed(731); reference <- do.call(fun, c(args, specification)); rng <- .Random.seed
      set.seed(731); actual <- do.call(fun, c(renamed, specification))
      expect_identical(.Random.seed, rng)
      for (field in c("ATT", "se", "boots", "att.inf.func")) expect_equal(actual[[field]], reference[[field]], tolerance = 1e-12)
    }
  }
})
