weight_validation_methods <- c("drdid_panel", "drdid_imp_panel", "reg_did_panel",
  "ipw_did_panel", "std_ipw_did_panel", "twfe_did_panel", "drdid_rc",
  "drdid_rc1", "drdid_imp_rc", "drdid_imp_rc1", "reg_did_rc", "ipw_did_rc",
  "std_ipw_did_rc", "twfe_did_rc")
weight_validation_fit <- function(method, w) {
  D <- rep(0:1, each = 40); post <- rep(rep(0:1, each = 20), 2)
  y <- sin(seq_along(D)) + D * post
  args <- list(D = D, covariates = NULL, i.weights = w, inffunc = TRUE)
  if (grepl("panel", method)) args <- c(args, list(y1 = y, y0 = rep(0, 80))) else
    args <- c(args, list(y = y, post = post))
  do.call(get(method, mode = "function"), args)
}
test_that("all public low-level methods reject invalid raw weights", {
  bad <- list(rep(-1, 80), rep(0, 80), c(Inf, rep(1,79)),
              c(NA_real_, rep(1,79)), rep(1,79))
  for (method in weight_validation_methods) for (w in bad)
    expect_error(weight_validation_fit(method, w), "i.weights must be finite|Panel inputs must be complete|i.weights must have one weight")
})
test_that("zero weights and common positive rescaling preserve estimates", {
  w <- rep(1,80); w[c(2,22,42,62)] <- 0
  for (method in weight_validation_methods) {
    a <- weight_validation_fit(method,w); b <- weight_validation_fit(method, 3*w)
    expect_true(is.finite(a$ATT))
    expect_equal(a$ATT,b$ATT,tolerance=1e-10)
    expect_equal(a$att.inf.func,b$att.inf.func,tolerance=1e-9)
  }
})
test_that("wrapper normalization cannot turn all-negative weights positive", {
  d <- expand.grid(id=1:40, time=0:1)
  d$group <- as.numeric(d$id > 20); d$outcome <- sin(d$id)+d$time*d$group
  d$weight <- -1
  for (panel in c(TRUE,FALSE)) for (method in c("drdid","ordid","ipwdid")) {
    args <- list(yname="outcome",tname="time",idname="id",dname="group",
                 data=d,panel=panel,weightsname="weight")
    expect_error(do.call(get(method,mode="function"),args),"i.weights must be finite")
  }
})
test_that("RC missing weights retain complete-row exclusion", {
  d <- expand.grid(id=1:40,time=0:1)
  d$group <- as.numeric(d$id>20); d$outcome <- sin(d$id)+d$time*d$group
  d$weight <- 1; d$weight[1] <- NA_real_
  args <- list(yname="outcome",tname="time",idname="id",dname="group",
               data=d,panel=FALSE,weightsname="weight")
  expect_warning(a <- do.call(ordid,args),"weights are missing")
  args$data <- d[-1,]; b <- do.call(ordid,args)
  expect_equal(a$ATT,b$ATT)
})
