
## Simulate data Following Sant'Anna and Zhao (2020)'s DGP 1 (RC Case)

# Sample size
n <- 1000
# pscore index (strength of common support)
Xsi.ps <- .75
# Proportion in each period
lambda <- 0.5
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

# Generate "T"
post <- as.numeric(stats::runif(n) <= lambda)
# observed outcome
y <- post * y1 + (1 - post) * y0
#-----------------------------------------------------------------------------
#Gen id
id <- as.numeric(1:n)
#-----------------------------------------------------------------------------
# Put in a long data frame
sim_rc <- as.data.frame(cbind(id = id, post = post, y = y, d = d,
                                x1 = z1, x2= z2, x3 = z3, x4 = z4))
sim_rc <- sim_rc[order(sim_rc$id),]

usethis::use_data(sim_rc, overwrite = TRUE)
