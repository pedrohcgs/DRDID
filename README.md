# DRDID: Doubly Robust Difference-in-Differences.
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/DRDID)](https://CRAN.R-project.org/package=DRDID)
[![Travis build status](https://travis-ci.com/pedrohcgs/DRDID.svg?branch=master)](https://travis-ci.com/pedrohcgs/DRDID)
[![Codecov test coverage](https://codecov.io/gh/pedrohcgs/DRDID/branch/master/graph/badge.svg)](https://codecov.io/gh/pedrohcgs/DRDID?branch=master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
<!-- badges: end -->

## Overview 

This `R` package implements different estimators for the Average Treatment Effect on the Treated (ATT) in Difference-in-Differences (DID) setups where the parallel trends assumption holds
after you condition on a vector of pre-treatment covariates. The main estimators implemented here are the locally efficient, doubly-robust DID estimators proposed by Sant'Anna and Zhao (2020), [Doubly Robust Difference-in-Differences Estimators](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315). In addition to these, the `DRDID` package also implements regression-based DID estimator in the spirit of Heckman, Ichimura and Todd (1997), and the inverse probability weighted (IPW) DID estimator proposed by Abadie (2005). 


See Sant'Anna and Zhao (2020) for additional details.

## Installing DRDID
This github website hosts the source code, and it always has the most updated version of the package.

To install the most recent version of the `DRDID` package from GitHub (this is what we recommend):

        library(devtools)
        devtools::install_github("pedrohcgs/DRDID")
        
## Authors 

Pedro H. C. Sant'Anna, Vanderbilt University, Nashville, TN. E-mail: pedro.h.santanna [at] vanderbilt [dot] edu.

Jun B. Zhao, Vanderbilt University, Nashville, TN. E-mail: jun.zhao [at] vanderbilt [dot] edu.


## References

* Abadie, Alberto (2005), "Semiparametric Difference-in-Differences Estimators", Review of Economic Studies, vol. 72(1), p. 1-19, <doi:10.1111/0034-6527.00321>.

* Graham, Bryan, Pinto, Cristine, and Egel, Daniel (2012), "Inverse Probability Tilting for Moment Condition Models with Missing Data", Review of Economic Studies, vol. 79(3), p. 1053-1079, <doi:10.1093/restud/rdr047>.

* Heckman, James J., Ichimura, Hidehiko, and Todd, Petra E. (1997), "Matching as an Econometric Evaluation Estimator: Evidence from Evaluating a Job Training Programme", Review of Economic Studies, vol. 64(4), p. 605–654, <doi:10.2307/2971733>.

* Sant'Anna, Pedro H. and Zhao, Jun B. (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315), Journal of Econometrics, Forthcoming.

