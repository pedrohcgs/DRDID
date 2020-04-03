#' @import stats
NULL
###################################################################################
#' Intrinsically Efficient Doubly Robust DID estimator with Panel Data
#'
#' Intrinsically Efficient Doubly Robust Difference-in-Differences Estimator for the ATT, with Panel Data
#'
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation.
#' If covariates = NULL, this leads to an unconditional DID estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if boot = FALSE). Options are "weighted" and "multiplier".
#' If boot==T, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = FALSE). Default is 999 if boot = TRUE.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DID point estimate}
#' \item{se}{The DR DID standard error}
#' \item{uci}{Estimate of the upper boudary of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower boudary of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#' @export

drdid_ie_panel <-function(y1, y0, D, covariates,
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
  #Compute the Pscore by MLE
  pscore.ie <- stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights)
  ps.fit <- as.vector(pscore.ie$fitted.values)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group using wols.ie
  outcome.reg <- wols.ie.panel(deltaY, D, int.cov, ps.fit, i.weights = i.weights)
  out.delta <-  outcome.reg$out.reg
  #-----------------------------------------------------------------------------
  #Compute Intrinsically Efficient Doubly Robust DID estimators
  # First, the weights
  w.treat <- i.weights * D
  w.cont <- i.weights * ps.fit * (1 - D)/(1 - ps.fit)

  dr.att.treat <- w.treat * (deltaY - out.delta)
  dr.att.cont <- w.cont * (deltaY - out.delta)

  eta.treat <- mean(dr.att.treat) / mean(w.treat)
  eta.cont <- mean(dr.att.cont) / mean(w.cont)

  dr.att <-   eta.treat - eta.cont
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.1 <- (dr.att.treat - w.treat * eta.treat)
  # Estimation effect from beta hat
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w.treat * int.cov)

  # Asymptotic linear representation of OLS parameters (assuming weights fixed) (not assumming pscore is corre)
  weights.ols <- as.vector(i.weights * (1 - D) * ((ps.fit/(1 - ps.fit))^2))
  wols.x <- weights.ols * int.cov
  wols.eX <- weights.ols * (deltaY - out.delta) * int.cov
  XpX.inv <- solve(crossprod(wols.x, int.cov)/n)
  asy.lin.rep.wols1 <-  wols.eX %*% XpX.inv

  # term for the Asymptotic linear representation of OLS parameters that accounts for random weights
  # Not assuming pscore is correct
  wols.ps <- 2 * (1 - D) * (ps.fit/(1 - ps.fit)) * i.weights * (deltaY - out.delta) * int.cov
  M.ols.ps <- base::crossprod(wols.ps, int.cov)/n
  M.ols.ps <- M.ols.ps %*% XpX.inv

  # Now I need the Asymptotic linear representation of the Logit parameters
  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(pscore.ie) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps

  # Finally get the Asy lin rep of the OLS parameter that accounts for random weights
  asy.lin.rep.wols2 <- asy.lin.rep.ps %*% M.ols.ps

  # Asy lin rep of Weighted OLS
  asy.lin.rep.wols <- asy.lin.rep.wols1 + asy.lin.rep.wols2

  # Now get the influence function related to the estimation effect related to beta's
  inf.treat.2 <- asy.lin.rep.wols %*% M1

  # Influence function for the treated component
  inf.treat <- (inf.treat.1 - inf.treat.2) / mean(w.treat)
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.1 <- (dr.att.cont - w.cont * eta.cont)
  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2 <- base::colMeans(w.cont *(deltaY - out.delta - eta.cont) * int.cov)
  # Now the influence function related to estimation effect of pscores
  inf.cont.2 <- asy.lin.rep.ps %*% M2
  # Estimation Effect from beta hat (weighted OLS)
  M3 <-  base::colMeans(w.cont * int.cov)
  # Now the influence function related to estimation effect of regressions
  inf.cont.3 <- asy.lin.rep.wols %*% M3

  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2 - inf.cont.3) / mean(w.cont)

  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func <- inf.treat - inf.control
  #-----------------------------------------------------------------------------
  if (boot == F) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- dr.att + 1.96 * se.dr.att
    # Estimate of lower doundary of 95% CI
    lci <- dr.att - 1.96 * se.dr.att
    #Create this null vector so we can export the bootstrap draws too.
    dr.boot <- NULL
  }

  if (boot == T) {
    if (is.null(nboot) == T) nboot = 999
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
      dr.boot <- unlist(lapply(1:nboot, wboot.dr.ie.panel,
                               n = n, deltaY = deltaY, D = D, int.cov = int.cov, i.weights = i.weights))
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


  if(inffunc==F) att.inf.func <- NULL
  return(list(ATT = dr.att,
              se = se.dr.att,
              uci = uci,
              lci = lci,
              boots = dr.boot,
              att.inf.func = att.inf.func))
}
