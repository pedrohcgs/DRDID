# DRDID: Doubly Robust Difference-in-Differences.
<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/DRDID)](https://CRAN.R-project.org/package=DRDID)
[![Travis build status](https://travis-ci.com/pedrohcgs/DRDID.svg?branch=master)](https://travis-ci.com/pedrohcgs/DRDID)
[![Codecov test coverage](https://codecov.io/gh/pedrohcgs/DRDID/branch/master/graph/badge.svg)](https://codecov.io/gh/pedrohcgs/DRDID?branch=master)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

## Overview 


This `R` package implements different DID procedures to estimate the ATT under the conditional parallel trends assumption, including a regression-based DID estimator in the spirit of Heckman, Ichimura and Todd (1997), the inverse probability weighted (IPW) DID estimator proposed by Abadie (2005), and the doubly robust (DR) DID estimators proposed
by Sant'Anna and Zhao (2020), [Doubly Robust Difference-in-Differences Estimators](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315).


On one hand, the regression DID estimator requires the regression model for the evolution of the potential outcomes for the comparison group to be correctly specified. On the other hand, the IPW DID estimator requires that the conditional probability of being in the treated group (the propensity score) is correctly specified. By noticing that these DID estimators rely on different, non-nested assumptions, [Sant'Anna and Zhao (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315) propose to combine them in a particular way such that the resulting DID estimators remain consistent if either (but not necessarily both) the outcome regression or the propensity score model is correctly specified. In addition, the DR DID estimators proposed by [Sant'Anna and Zhao (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315) attain the semiparametric efficiency bound when both the outcome regression and propensity score models are correctly specified.


As discussed in [Sant'Anna and Zhao (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315), a researcher can construct different DR DID estimators depending on how he/she estimates the nuisance functions. The `DRDID` package implements two DR DID estimators when panel data are available, and four DR DID estimators when repeated cross-section data are available:


**PANEL DATA ESTIMATORS**
        
* `drdid_panel` - This implements the "traditional" DR DID estimator proposed in Sant'Anna and Zhao (2020) when panel data are available, where one uses a logistic propensity score model and a linear regression model for the outcome growth of the comparison units. The propensity score parameters are estimated via maximum likelihood, and the outcome regression parameters are estimated via OLS - hence the label "traditional". When both working models are correctly specified, the DR DID estimator is semiparametrically efficient.

* `drdid_imp_panel` - This implements the "improved" DR DID estimator proposed in Sant'Anna and Zhao (2020) when panel data are available. On top of being DR consistent, and locally semiparametrically efficient as the estimator implemented via `drdid_panel`, it is also DR for inference - hence the label "improved". The propensity score parameters are estimated via the inverse probability tilting procedure proposed by Graham, Pinto and Egel (2012), and the outcome regression parameters via weighted least squares.


**REPEATED CROSS-SECTION DATA ESTIMATORS**
        
* `drdid_rc` - This implements the "traditional" DR and locally efficient DID estimator proposed in Sant'Anna and Zhao (2020) when repeated cross-section data are available. This estimator makes use of outcome regressions for both treated and comparison units. Like `drdid_panel`, the propensity score parameters are estimated via maximum likelihood, and the outcome regression parameters are estimated via OLS. 

* `drdid_rc1` - This implements the "Traditional" DR DID estimator proposed in Sant'Anna and Zhao (2020) when repeated cross-section data are available. This estimator makes use of outcome regressions only for the comparison units, though it does not attain the semiparametric efficient bound. Like `drdid_rc`, the propensity score parameters are estimated via maximum likelihood, and the outcome regression parameters are estimated via OLS. Given that `drdid_rc1` and `drdid_rc` enjoy the same double-robustness properties but `drdid_rc` is, in general, more efficient than `drdid_rc1`, we recommend practitioners to favor `drdid_rc` in detriment of `drdid_rc1`.

* `drdid_imp_rc` - This implements the "further improved" DR and locally efficient DID estimator proposed in Sant'Anna and Zhao (2020) when repeated cross-section data are available. This estimator makes use of outcome regressions for both treated and comparison units. Like `drdid_imp_panel`, the propensity score parameters are estimated using the procedure of Graham, Pinto and Egel (2012), and the outcome regression parameters for the comparison units are estimated via weighted least squares. On the other hand, the outcome regression parameters for the treated units are estimated via OLS.

* `drdid_imp_rc1` - This implements the "further improved" DR DID estimator proposed in Sant'Anna and Zhao (2020) when repeated cross-section data are available. This estimator makes use of outcome regressions only for the comparison units, though it does not attain the semiparametric efficient bound. Like `drdid_imp_panel`, the propensity score parameters are estimated using the procedure of Graham, Pinto and Egel (2012), and the outcome regression parameters for the comparison units are estimated via weighted least squares. Given that `drdid_imp_rc1` and `drdid_imp_rc` enjoy the same double-robustness properties but `drdid_imp_rc` is, in general, more efficient than `drdid_imp_rc1`, we recommend practitioners to favor `drdid_imp_rc` in detriment of `drdid_imp_rc`.


For further details, please see the paper Sant'Anna and Zhao (2020), [Doubly Robust Difference-in-Differences Estimators](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315). In case you have any comments and/or questions, please contact Pedro Sant'Anna (see email below).

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

* Sant'Anna, Pedro H. and Zhao, Jun B. (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315), Working paper.

