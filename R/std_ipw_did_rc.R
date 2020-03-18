#' @import stats
NULL
###################################################################################
# Standardized version of Abadie's IPW DID estimator

#' Standardized Inverse probability weighted Difference-in-Differences Estimator with Repeated Cross Section Data
#'
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-tretment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if boot = FALSE). Options are "weighted" and "multiplier".
#' If boot==T, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = FALSE). Default is 999 if boot = TRUE.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The IPW DID point estimate.}
#' \item{se}{ The IPW DID standard error}
#' \item{uci}{Estimate of the upper boudary of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower boudary of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#' @export

std_ipw_did_rc <-function(y, post, D, covariates,
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
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Add constant to covariate vector
  int.cov <- as.matrix(cbind(1, covariates))

  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Pscore estimation (logit) and also its fitted values
  PS <- suppressWarnings(stats::glm(D ~ -1 + int.cov, family = "binomial", weights = i.weights))
  ps.fit <- as.vector(PS$fitted.values)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * y / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * y / mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * y / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * y / mean(w.cont.post)

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  # ATT estimator
  ipw.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre)
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  Hessian.ps <- stats::vcov(PS) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  inf.treat <- inf.treat.post - inf.treat.pre
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)
  inf.cont <- inf.cont.post - inf.cont.pre

  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2.pre <- base::colMeans(w.cont.pre *(y - att.cont.pre) * int.cov)/mean(w.cont.pre)
  M2.post <- base::colMeans(w.cont.post *(y - att.cont.post) * int.cov)/mean(w.cont.post)

  # Now the influence function related to estimation effect of pscores
  inf.cont.ps <- asy.lin.rep.ps %*% (M2.post - M2.pre)

  # Influence function for the control component
  inf.cont <- inf.cont + inf.cont.ps

  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  if (boot == F) {
    # Estimate of standard error
    se.att <- stats::sd(att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- ipw.att + 1.96 * se.att
    # Estimate of lower doundary of 95% CI
    lci <- ipw.att - 1.96 * se.att

    #Create this null vector so we can export the bootstrap draws too.
    ipw.boot <- NULL
  }
  if (boot == T) {
    if (is.null(nboot) == T) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      ipw.boot <- mboot.did(att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR(ipw.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(ipw.boot/se.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower doundary of 95% CI
      lci <- ipw.att - cv * se.att
    } else {
      # do weighted bootstrap
      ipw.boot <- unlist(lapply(1:nboot, wboot_std_ipw_rc,
                                n = n, y = y, post = post, D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR(ipw.boot - ipw.att) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((ipw.boot - ipw.att)/se.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower doundary of 95% CI
      lci <- ipw.att - cv * se.att

    }

  }
  if(inffunc==F) att.inf.func <- NULL

  return(list(ATT = ipw.att,
              se = se.att,
              uci = uci,
              lci = lci,
              boots = ipw.boot,
              att.inf.func = att.inf.func))
}
