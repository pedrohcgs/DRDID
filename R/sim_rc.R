#' @title Simulated repeated cross-section data
#'
#' @description \code{sim_rc} contains a simulated dataset following the DGP1 in Sant'Anna
#' and Zhao (2020).
#'
#'
#' @format A data frame in "long" format with 1000 observations on the following and 8 variables:
#' \describe{
#'   \item{id}{unique identifier for each cross-sectional unit.}
#'   \item{post}{an indicator variable for post-treatment period (1 if post,
#'    0 if pre treatment period).}
#'   \item{y}{outcome of interest}
#'   \item{d}{an indicator variable for treatment group. Equal to 1 if
#'   experience treatment in the post-treatment period; equal to 0 if never experience
#'   treatment.}
#'   \item{x1}{Covariate z1 in Sant'Anna and Zhao(2020)}
#'   \item{x2}{Covariate z2 in Sant'Anna and Zhao(2020)}
#'   \item{x3}{Covariate z3 in Sant'Anna and Zhao(2020)}
#'   \item{x4}{Covariate z4 in Sant'Anna and Zhao(2020)}
#
#' }
#' @source Sant'Anna and Zhao (2020)
#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#' }
"sim_rc"
