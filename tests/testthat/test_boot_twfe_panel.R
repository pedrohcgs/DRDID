context("Compare bootstrapped and analytical std errors: TWFE Panel")

test_that("Analytical and bootstrapped std errors are similar: TWFE Panel", {

  # Let us generate some panel data
  #-----------------------------------------------------------------------------
  # DGP 1 used by Sant'Anna and Zhao (2020) (Panel data case)
  # Sample size
  n <- 500
  # pscore index (strength of common support)
  Xsi.ps <- .75
  #-----------------------------------------------------------------------------
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
  dta_long <- as.data.frame(cbind(id = id, y = y1, d = d, post = TRUE,
                                  x1 = z1, x2= z2, x3 = z3, x4 = z4))
  dta_long <- data.frame(rbind(dta_long,cbind(id = id, y = y0, d = d, post = FALSE,
                                              x1 = z1, x2= z2, x3 = z3, x4 = z4)))
  dta_long <- dta_long[order(dta_long$id),]
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Use the different estimators to compute ATT
  #-----------------------------------------------------------------------------
  # Analytical std errors
  twfe.did_panel <- twfe_did_panel(dta_wide$y1,dta_wide$y0, dta_wide$d,
                                   covariates = dta_wide[,5:8], boot = FALSE)

  twfe.did_panel_n <- twfe_did_panel(dta_wide$y1, dta_wide$y0, dta_wide$d,
                                     covariates = NULL, boot = FALSE)
  #-----------------------------------------------------------------------------
  # Now with bootstrap (weighted)
  twfe.did_panel2 <- twfe_did_panel(dta_wide$y1,dta_wide$y0, dta_wide$d,
                                    covariates = dta_wide[,5:8], boot = TRUE)

  twfe.did_panel_n2 <- twfe_did_panel(dta_wide$y1, dta_wide$y0, dta_wide$d,
                                      covariates = NULL, boot = TRUE,
                                      boot.type = "weighted")
  # #-----------------------------------------------------------------------------
  # Now with bootstrap (multiplier)
  twfe.did_panel3 <- twfe_did_panel(dta_wide$y1,dta_wide$y0, dta_wide$d,
                                    dta_wide[,5:8], boot = TRUE, boot.type = "multiplier")

  twfe.did_panel_n3 <- twfe_did_panel(dta_wide$y1, dta_wide$y0, dta_wide$d,
                                      covariates = rep(1, n), boot = TRUE,
                                      boot.type = "multiplier")


  #-----------------------------------------------------------------------------
  # Check if all point estimates are equal
  expect_equal(twfe.did_panel$ATT, twfe.did_panel2$ATT)
  expect_equal(twfe.did_panel3$ATT, twfe.did_panel2$ATT)

  expect_equal(twfe.did_panel_n3$ATT, twfe.did_panel_n$ATT)
  expect_equal(twfe.did_panel_n2$ATT, twfe.did_panel_n$ATT)


  # Check if all standard errors are equal
  expect_equal(twfe.did_panel2$se, twfe.did_panel$se, tol = 0.05)
  expect_equal(twfe.did_panel3$se, twfe.did_panel$se, tol = 0.05)

  expect_equal(twfe.did_panel_n3$se, twfe.did_panel_n$se, tol = 0.15)
  expect_equal(twfe.did_panel_n2$se, twfe.did_panel_n$se, tol = 0.15)


})
