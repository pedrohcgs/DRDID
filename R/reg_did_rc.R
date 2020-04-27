#' @import stats
NULL
###################################################################################
#' Regression-based DID estimator with Repeated Cross Section Data
#' Regression-based Difference-in-Differences Estimator for the ATT, with Repeated Cross Section Data
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation.
#' If covariates = NULL, this leads to an unconditional DID estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if boot = FALSE). Options are "weighted" and "multiplier".
#' If boot==T, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = FALSE). Default is 999 if boot = TRUE.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The Reg DID point estimate}
#'  \item{se}{The Reg DID standard error}
#'  \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#'  \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = F, boot, boot.type, nboot, type="or")}

#' @references Heckman, James J., Ichimura, Hidehiko, and Todd, Petra E. (1997),"Matching as an Econometric Evaluation Estimator: Evidence from Evaluating a Job Training Programme",
#' Review of Economic Studies, vol. 64(4), p. 605–654, <doi:10.2307/2971733>.
#' @references Sant'Anna, Pedro H. C. and Zhao, Jun (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315).

#' @export

reg_did_rc <-function(y, post, D, covariates, i.weights = NULL,
                      boot = F, boot.type = "weighted", nboot = NULL,
                      inffunc = F){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # post as vector
  post <- as.vector(post)
  # Sample size
  n <- length(D)
  # outcome of interested
  y <- as.vector(y)
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


  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                         subset = ((D==0) & (post==0)),
                                         weights = i.weights))
  out.y.pre <-   as.vector(tcrossprod(reg.coeff.pre, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group at the pre-treatment period, using ols.
  reg.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                          subset = ((D==0) & (post==1)),
                                          weights = i.weights))
  out.y.post <-   as.vector(tcrossprod(reg.coeff.post, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the OR DID estimators
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont <- i.weights * D

  reg.att.treat.pre <- w.treat.pre * y
  reg.att.treat.post <- w.treat.post * y
  reg.att.cont <- w.cont * (out.y.post - out.y.pre)

  eta.treat.pre <- mean(reg.att.treat.pre) / mean(w.treat.pre)
  eta.treat.post <- mean(reg.att.treat.post) / mean(w.treat.post)
  eta.cont <- mean(reg.att.cont) / mean(w.cont)

  reg.att <- (eta.treat.post - eta.treat.pre) - eta.cont

  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters in pre-period
  weights.ols.pre <- i.weights * (1 - D) * (1 - post)
  wols.x.pre <- weights.ols.pre * int.cov
  wols.eX.pre <- weights.ols.pre * (y - out.y.pre) * int.cov
  XpX.inv.pre <- solve(crossprod(wols.x.pre, int.cov)/n)
  asy.lin.rep.ols.pre <-  wols.eX.pre %*% XpX.inv.pre

  # Asymptotic linear representation of OLS parameters in post-period
  weights.ols.post <- i.weights * (1 - D) * post
  wols.x.post <- weights.ols.post * int.cov
  wols.eX.post <- weights.ols.post * (y - out.y.post) * int.cov
  XpX.inv.post <- solve(crossprod(wols.x.post, int.cov)/n)
  asy.lin.rep.ols.post <-  wols.eX.post %*% XpX.inv.post
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function
  inf.treat.pre <- (reg.att.treat.pre - w.treat.pre * eta.treat.pre) / mean(w.treat.pre)
  inf.treat.post <- (reg.att.treat.post - w.treat.post * eta.treat.post) / mean(w.treat.post)
  inf.treat <- inf.treat.post - inf.treat.pre
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.1 <- (reg.att.cont - w.cont * eta.cont)
  # Estimation effect from beta hat (OLS using only controls)
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w.cont * int.cov)
  # Now get the influence function related to the estimation effect related to beta's in post-treatment
  inf.cont.2.post <- asy.lin.rep.ols.post %*% M1
  # Now get the influence function related to the estimation effect related to beta's in pre-treatment
  inf.cont.2.pre <- asy.lin.rep.ols.pre %*% M1
  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2.post - inf.cont.2.pre) / mean(w.cont)
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  reg.att.inf.func <- (inf.treat - inf.control)
  #-----------------------------------------------------------------------------
  if (boot == F) {
    # Estimate of standard error
    se.reg.att <- stats::sd(reg.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- reg.att + 1.96 * se.reg.att
    # Estimate of lower doundary of 95% CI
    lci <- reg.att - 1.96 * se.reg.att
    #Create this null vector so we can export the bootstrap draws too.
    reg.boot <- NULL
  }

  if (boot == T) {
    if (is.null(nboot) == T) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      reg.boot <- mboot.did(reg.att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.reg.att <- stats::IQR(reg.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(reg.boot/se.reg.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- reg.att + cv * se.reg.att
      # Estimate of lower doundary of 95% CI
      lci <- reg.att - cv * se.reg.att
    } else {
      # do weighted bootstrap
      reg.boot <- unlist(lapply(1:nboot, wboot_reg_rc,
                                n = n, y = y, post = post, D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.reg.att <- stats::IQR((reg.boot - reg.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((reg.boot - reg.att)/se.reg.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- reg.att + cv * se.reg.att
      # Estimate of lower doundary of 95% CI
      lci <- reg.att - cv * se.reg.att

    }
  }


  if(inffunc==F) reg.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot==T, T, F)
  argu <- list(
    panel = F,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "or"
  )
  ret <- (list(ATT = reg.att,
               se = se.reg.att,
               uci = uci,
               lci = lci,
               boots = reg.boot,
               att.inf.func = reg.att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
