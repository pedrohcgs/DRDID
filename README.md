# Doubly Robust Difference-in-Differences <img src="man/figures/logo.png" align="right" alt="" width="95" />
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/DRDID)](https://CRAN.R-project.org/package=DRDID)
[![Codecov test coverage](https://codecov.io/gh/pedrohcgs/DRDID/branch/master/graph/badge.svg)](https://app.codecov.io/gh/pedrohcgs/DRDID?branch=master)
[![R build status](https://github.com/pedrohcgs/DRDID/workflows/R-CMD-check/badge.svg)](https://github.com/pedrohcgs/DRDID/actions)
<!-- badges: end -->


The `DRDID` R package implements different estimators for the Average Treatment Effect on the Treated (ATT) in Difference-in-Differences (DID) setups where the parallel trends assumption holds after conditioning on a vector of pre-treatment covariates.


The main estimators implemented here are the locally efficient, doubly-robust DID estimators proposed by [Sant'Anna and Zhao (2020), Doubly Robust Difference-in-Differences Estimators](https://doi.org/10.1016/j.jeconom.2020.06.003). The package covers both panel data and repeated cross-section data setups with two treatment groups (treated and comparison group) and two time periods (pre-treatment and post-treatment).


See the [package manual](https://psantanna.com/DRDID/reference/index.html) for documentation of all package functions (with examples).


If you end up using this package, please cite our paper:
* Sant'Anna, Pedro H. C., and Zhao, Jun (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://doi.org/10.1016/j.jeconom.2020.06.003), *Journal of Econometrics*, Vol. 219 (1), pp. 101-122.


## Installation
To install the most recent version of the `DRDID` package from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("pedrohcgs/DRDID")
```

## Short example
The following is a portion of the empirical illustration considered by Sant'Anna and Zhao (2020)
that uses the LaLonde sample from the NSW experiment and considers data from the Current Population Survey (CPS) to form a non-experimental comparison group.

Let's first get the data ready:

``` r
library(DRDID)
# Load data in long format that comes in the DRDID package
data(nsw_long)
# Form the Lalonde sample with CPS comparison group
eval_lalonde_cps <- subset(nsw_long, nsw_long$treated == 0 | nsw_long$sample == 2)
```

Now, to estimate the ATT using the Improved Locally Efficient Doubly Robust DID estimator, we can use the **drdid** function:
```{r}
# Implement improved locally efficient DR DID:
out <- drdid(yname = "re", tname = "year", idname = "id", dname = "experimental",
             xformla= ~ age + educ + black + married + nodegree + hisp + re74,
             data = eval_lalonde_cps, panel = TRUE)
summary(out)

```

For additional details on the usage of the **drdid** function, see [check the manual](https://psantanna.com/DRDID/reference/drdid.html).
        

To implement IPW and outcome regression DID estimators, check [here](https://psantanna.com/DRDID/reference/ipwdid.html) and [here](https://psantanna.com/DRDID/reference/ordid.html), respectively.

