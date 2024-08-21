#' @import stats
NULL
###################################################################################
# 'Traditional' Doubly Robust DiD estimator with Repeated Cross Section Data
#' Doubly robust DiD estimator for the ATT, with repeated cross-section data
#'
#' @description \code{drdid_rc1} is used to compute the doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with stationary repeated cross-sectional data. The resulting estimator is
#'  not locally efficient; see Section 3.2 of Sant'Anna and Zhao (2020).
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DiD point estimate}
#' \item{se}{ The DR DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, estMethod = "trad2", boot, boot.type, nboot, type="dr")}
#' @details
#'
#' The \code{drdid_rc1} function implements the doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.3)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of (separate) linear regression models for the outcome among the comparison units in both pre and
#' post-treatment time periods. Importantly, this estimator is not locally efficient for the ATT.
#'
#'
#' The propensity score parameters are estimated using maximum
#' likelihood, and the outcome regression coefficients are estimated using ordinary least squares.
#'
#' The resulting estimator is not not locally efficient; see Sant'Anna and Zhao (2020) for details.
#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#' }
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(sim_rc[,5:8])
#' # Implement the 'traditional' DR DiD estimator (not locally efficient!)
#' drdid_rc1(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'           covariates= covX)
#'
#' @export

drdid_rc1 <-function(y, post, D, covariates, i.weights = NULL,
                     boot = FALSE, boot.type =  "weighted", nboot = NULL,
                     inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all(as.matrix(covariates)[,1]==rep(1,n))){
      int.cov <- as.matrix(covariates)
    } else {
      int.cov <- as.matrix(cbind(1, covariates))
    }
  }

  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  # Normalize weights
  i.weights <- i.weights/mean(i.weights)
  #-----------------------------------------------------------------------------
  #Compute the Pscore by MLE
  pscore.tr <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  if(pscore.tr$converged == FALSE){
    warning(" glm algorithm did not converge")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- as.vector(pscore.tr$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  if(anyNA(reg.coeff.pre)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  if(anyNA(reg.coeff.post)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  # Combine the ORs
  out.y <- post * out.y.post + (1 - post) * out.y.pre
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y) / mean(w.cont.post)

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions

  # Asymptotic linear representation of OLS parameters in pre-period
  weights.ols.pre <- i.weights * (1 - D) * (1 - post)
  wols.x.pre <- weights.ols.pre * int.cov
  wols.eX.pre <- weights.ols.pre * (y - out.y.pre) * int.cov
  XpX.pre <- crossprod(wols.x.pre, int.cov)/n
  # Check if XpX is invertible
  if ( base::rcond(XpX.pre) < .Machine$double.eps) {
    stop("The regression design matrix for pre-treatment is singular. Consider removing some covariates.")
  }
  XpX.inv.pre <- solve(XpX.pre)
  asy.lin.rep.ols.pre <-  wols.eX.pre %*% XpX.inv.pre

  # Asymptotic linear representation of OLS parameters in post-period
  weights.ols.post <- i.weights * (1 - D) * post
  wols.x.post <- weights.ols.post * int.cov
  wols.eX.post <- weights.ols.post * (y - out.y.post) * int.cov
  XpX.post <- crossprod(wols.x.post, int.cov)/n
  # Check if XpX is invertible
  if ( base::rcond(XpX.post) < .Machine$double.eps) {
    stop("The regression design matrix for post-treatment is singular. Consider removing some covariates.")
  }
  XpX.inv.post <- solve(XpX.post)
  asy.lin.rep.ols.post <-  wols.eX.post %*% XpX.inv.post

  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(pscore.tr) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)

  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M1.post <- - base::colMeans(w.treat.post * post * int.cov)/mean(w.treat.post)
  M1.pre <- - base::colMeans(w.treat.pre * (1 - post) * int.cov)/mean(w.treat.pre)

  # Now get the influence function related to the estimation effect related to beta's
  inf.treat.or.post <- asy.lin.rep.ols.post %*% M1.post
  inf.treat.or.pre <- asy.lin.rep.ols.pre %*% M1.pre
  inf.treat.or <- inf.treat.or.post + inf.treat.or.pre

  # Influence function for the treated component
  inf.treat <- inf.treat.post - inf.treat.pre + inf.treat.or
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)

  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2.pre <- base::colMeans(w.cont.pre *(y - out.y - att.cont.pre) * int.cov)/mean(w.cont.pre)
  M2.post <- base::colMeans(w.cont.post *(y - out.y - att.cont.post) * int.cov)/mean(w.cont.post)
  # Now the influence function related to estimation effect of pscores
  inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)

  # Estimation effect from beta hat from post and pre-periods
  # Derivative matrix (k x 1 vector)
  M3.post <- - base::colMeans(w.cont.post * post * int.cov) / mean(w.cont.post)
  M3.pre <- - base::colMeans(w.cont.pre * (1 - post) * int.cov) / mean(w.cont.pre)

  # Now get the influence function related to the estimation effect related to beta's
  inf.cont.or.post <- asy.lin.rep.ols.post %*% M3.post
  inf.cont.or.pre <- asy.lin.rep.ols.pre %*% M3.pre
  inf.cont.or <- inf.cont.or.post + inf.cont.or.pre

  # Influence function for the control component
  inf.cont <- inf.cont.post - inf.cont.pre + inf.cont.ps + inf.cont.or
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + 1.96 * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - 1.96 * se.dr.att
    #Create this null vector so we can export the bootstrap draws too.
    dr.boot <- NULL
  }

  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      dr.boot <- mboot.did(dr.att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.dr.att <- stats::IQR(dr.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower doundary of 95% CI
      lci <- dr.att - cv * se.dr.att
    } else {
      # do weighted bootstrap
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_rc1,
                               n = n, y = y, post = post,
                               D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.dr.att <- stats::IQR((dr.boot - dr.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((dr.boot - dr.att)/se.dr.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower doundary of 95% CI
      lci <- dr.att - cv * se.dr.att

    }
  }


  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    estMethod = "trad2",
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "dr"
  )
  ret <- (list(ATT = dr.att,
               se = se.dr.att,
               uci = uci,
               lci = lci,
               boots = dr.boot,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
