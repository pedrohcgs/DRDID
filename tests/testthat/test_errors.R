context("Test Error messages")

test_that("Error and Warning messages are working", {

  # Let us generate some panel data
  #-----------------------------------------------------------------------------
  # DGP 1 used by Sant'Anna and Zhao (2020) (Panel data case)
  # Sample size
  n <- 500
  # pscore index (strength of common support)
  Xsi.ps <- .75
  # Researcher always observes Z
  #-----------------------------------------------------------------------------
  # Mean and Std deviation of Z's without truncation
  mean.z1 <- exp(0.25/2)
  sd.z1 <- sqrt((exp(0.25) - 1) * exp(0.25))
  mean.z2 <- 10
  sd.z2 <- 0.54164
  mean.z3 <- 0.21887
  sd.z3 <-   0.04453
  mean.z4 <- 402
  sd.z4 <-  56.63891
  #-----------------------------------------------------------------------------
  set.seed(1234)
  # Gen covariates
  x1 <- stats::rnorm(n, mean = 0, sd = 1)
  x2 <- stats::rnorm(n, mean = 0, sd = 1)
  x3 <- stats::rnorm(n, mean = 0, sd = 1)
  x4 <- stats::rnorm(n, mean = 0, sd = 1)

  z1 <- exp(x1/2)
  z2 <- x2/(1 + exp(x1)) + 10
  z3 <- (x1 * x3/25 + 0.6)^3
  z4 <- (x1 + x4 + 20)^2

  z1 <- (z1 - mean.z1)/sd.z1
  z2 <- (z2 - mean.z2)/sd.z2
  z3 <- (z3 - mean.z3)/sd.z3
  z4 <- (z4 - mean.z4)/sd.z4

  x <- cbind(x1, x2, x3, x4)
  z <- cbind(z1, z2, z3, z4)
  #-----------------------------------------------------------------------------
  # Gen treatment groups
  # Propensity score
  pi <- stats::plogis(Xsi.ps * (- z1 + 0.5 * z2 - 0.25 * z3 - 0.1 * z4))
  d  <- as.numeric(runif(n) <= pi)
  #-----------------------------------------------------------------------------
  # Generate aux indexes for the potential outcomes
  index.lin <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)
  index.unobs.het <- d * (index.lin)
  index.att <- 0

  #This is the key for consistency of outcome regression
  index.trend <- 210 + 27.4*z1 + 13.7*(z2 + z3 + z4)

  #v is the unobserved heterogeneity
  v <- stats::rnorm(n, mean = index.unobs.het, sd = 1)

  #Gen realized outcome at time 0
  y0 <- index.lin + v + stats::rnorm(n)

  # gen outcomes at time 1
  # First let's generate potential outcomes: y_1_potential
  y10 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend #this is for the trend based on X

  y11 <- index.lin + v + stats::rnorm(n, mean = 0, sd = 1) +#This is the baseline
    index.trend + #this is for the trend based on X
    index.att # This is the treatment effects

  # Gen realized outcome at time 1
  y1 <- d * y11 + (1 - d) * y10
  #-----------------------------------------------------------------------------
  #Gen id
  id <- 1:n
  #-----------------------------------------------------------------------------
  # Put in a "wide" data frame
  dta_wide <- as.data.frame(cbind(id = id, y1 = y1, y0 = y0, d = d,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  # Make "long" data
  dta_long <- as.data.frame(cbind(id = id, y = y1, d = d, post = T,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- data.frame(rbind(dta_long,cbind(id = id, y = y0, d = d, post = F,
                                              x1 = z1, x2= z2, x3 = z3, x4 = z4)))
  dta_long <- dta_long[order(dta_long$id),]

  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Test errors
  #-----------------------------------------------------------------------------
  # Warning for estMethod

  expect_warning(drdid(yname="y",
                       tname = "post",
                       idname = "id",
                       dname = "d",
                       estMethod = "other",
                       xformla= ~ x1 + x2 + x3 + x4,
                       data = dta_long,
                       panel=T,
                       boot = F))

  # Warning for boot.type
  expect_warning(drdid(yname="y",
                       tname = "post",
                       idname = "id",
                       dname = "d",
                       xformla= ~ x1 + x2 + x3 + x4,
                       data = dta_long,
                       panel=T,
                       boot = T,
                       nboot = 100,
                       boot.type = "other"))

  # Warning for yname
  expect_error(drdid(yname="other",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long,
                     panel=T,
                     boot = F))
  # Warning for tname
  expect_error(drdid(yname="y",
                     tname = "other",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long,
                     panel=T,
                     boot = F))
  # Warning for dname
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "other",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long,
                     panel=T,
                     boot = F))

  # Warning for idname
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "other",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long,
                     panel=T,
                     boot = F))

  # Warning for weightsname
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     weightsname = "other",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long,
                     panel=T,
                     boot = F))

  # Warning for covariates
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4 + x5,
                     data = dta_long,
                     panel=T,
                     boot = F))

  # Warning for normalized
  expect_warning(ipwdid(yname="y",
                        tname = "post",
                        idname = "id",
                        dname = "d",
                        normalized = "W",
                        xformla= ~ x1 + x2 + x3 + x4,
                        data = dta_long,
                        panel=T,
                        boot = F))

  # Warning for more than two groups
  dta_long2 <- dta_long
  dta_long2$d <- base::round(stats::runif(n, 0,2))
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))

  # Warning for more than two periods
  dta_long2$d <- dta_long$d
  dta_long2$post <- base::round(stats::runif(n, 0,2))
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))

  # Warning for more than two periods
  dta_long2$post <- ifelse(dta_long$post==0, "pre", "post")
  expect_warning(drdid(yname="y",
                       tname = "post",
                       idname = NULL,
                       dname = "d",
                       xformla= ~ x1 + x2 + x3 + x4,
                       data = dta_long2,
                       panel=F,
                       boot = F))

  # Error for time-varying covariates
  dta_long2$post <- dta_long$post
  dta_long2$x1 <- rnorm(n)
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))

  # Error for time-varying weights
  dta_long2$x1 <- dta_long$x1
  dta_long2$ww <- abs(rnorm(n))
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     weightsname  = "ww",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))
  # Error for time-varying  groups
  dta_long2$d <- round(runif(n))
  expect_error(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))
  dta_long2$d <- dta_long$d

  # Error for small groups
  dta_long2$d <- 1
  dta_long2$d[1:4]<- 0

  expect_warning(drdid(yname="y",
                     tname = "post",
                     idname = "id",
                     dname = "d",
                     xformla= ~ x1 + x2 + x3 + x4,
                     data = dta_long2,
                     panel=T,
                     boot = F))



  #-----------------------------------------------------------------------------

})
