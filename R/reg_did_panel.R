#' @import stats
NULL
###################################################################################
#' Regression-based DID estimator with Panel Data
#' Regression-based Difference-in-Differences Estimator for the ATT, with Panel Data
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
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
#'  \item{uci}{Estimate of the upper boudary of a 95\% CI for the ATT}
#'  \item{lci}{Estimate of the lower boudary of a 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#' @export

reg_did_panel <-function(y1, y0, D, covariates,
                         i.weights = NULL,
                         boot = F,
                         boot.type = "weighted",
                         nboot = NULL,
                         inffunc = F){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
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
  #Compute the Outcome regression for the control group using ols.
  reg.coeff <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
                                     subset = D==0,
                                     weights = i.weights))
  out.delta <-   as.vector(tcrossprod(reg.coeff, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the OR-DID estimator
  # First, the weights
  w.treat <- i.weights * D
  w.cont <- i.weights * D

  reg.att.treat <- w.treat * deltaY
  reg.att.cont <- w.cont * out.delta

  eta.treat <- mean(reg.att.treat) / mean(w.treat)
  eta.cont <- mean(reg.att.cont) / mean(w.cont)

  reg.att <-   eta.treat - eta.cont
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters
  weights.ols <- i.weights * (1 - D)
  wols.x <- weights.ols * int.cov
  wols.eX <- weights.ols * (deltaY - out.delta) * int.cov
  XpX.inv <- solve(crossprod(wols.x, int.cov)/n)
  asy.lin.rep.ols <-  wols.eX %*% XpX.inv
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function
  inf.treat <- (reg.att.treat - w.treat * eta.treat) / mean(w.treat)
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.1 <- (reg.att.cont - w.cont * eta.cont)
  # Estimation effect from beta hat (OLS using only controls)
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w.cont * int.cov)
  # Now get the influence function related to the estimation effect related to beta's
  inf.cont.2 <- asy.lin.rep.ols %*% M1
  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)
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
      reg.boot <- unlist(lapply(1:nboot, wboot.reg.panel,
                                n = n, deltaY = deltaY, D = D, int.cov = int.cov, i.weights = i.weights))
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

  if(inffunc==F) att.inf.func <- NULL

  return(list(ATT = reg.att,
              se = se.reg.att,
              uci = uci,
              lci = lci,
              boots = reg.boot,
              att.inf.func = att.inf.func))
}
