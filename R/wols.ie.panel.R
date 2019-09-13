#' @import stats
NULL

###################################################################################
# Compute Weighted OLS regression parameters for the (Intrinsically-Efficient) Doubly-Robust DID estimator with Panel Data
#' Weighted OLS for Intrinsically Efficient DR-DID
#'
#' Compute Weighted OLS regression parameters for the (Intrinsically Efficient) Doubly-Robust DID estimator
#'
#' This function implements the Weighted OLS estimator for the Not-Treated group that is used to construct the
#' (Intrinsically Efficient) Doubly Robust DID estimator, when Panel Data is available
#'
#' @param deltaY An \eqn{n} x \eqn{1} vector of the difference of outcomes between post and pre treatment periods, i.e., Y_post - Y_pre.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param int.cov An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation (with intercept)
#' @param pscore An \eqn{n} x \eqn{1} vector of propensity score estimates.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used.
#'
#' @return   A list containing the following components:
#' \item{out.reg}{The fitted values of the outcome regression}
#'
#' @references
#'
#' Sant'Anna, Pedro H.C., and Zhao, Jun B. (2019) "Doubly Robust difference-in-difference estimators", \emph{Working Paper}.


wols.ie.panel <- function(deltaY, D, int.cov, pscore, i.weights){
  #-----------------------------------------------------------------------------
  #Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  i.weights <- as.vector( i.weights * (pscore/(1 - pscore))^2)
  #Run weighted OLS
  beta.ie. <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
                                    subset = D==0,
                                    weights = i.weights))

  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.ie., int.cov))

  # return pscore and flag
  return(list(out.reg = out.delta))

}
