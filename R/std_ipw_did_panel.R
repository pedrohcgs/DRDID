#' @import stats
NULL
###################################################################################
# Standardized version of Abadie's IPW DID estimator

#' Standardized Inverse probability weighted Difference-in-Differences Estimator with Panel Data
#'
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Deafault is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if boot = FALSE). Options are "weighted" and "multiplier".
#' If boot==T, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = FALSE). Deafault is 999 if boot = TRUE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The IPW DID point estimate.}
#' \item{se}{ The IPW DID standard error}
#' \item{uci}{Estimate of the upper boudary of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower boudary of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'
#' @export

std_ipw_did_panel <-function(y1, y0, D, covariates,
                         i.weights = NULL,
                         boot = F,
                         boot.type = "weighted",
                         nboot = NULL){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
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
  w.treat <- i.weights * D
  w.cont <- i.weights * ps.fit * (1 - D)/(1 - ps.fit)

  att.treat <- w.treat * deltaY
  att.cont <- w.cont * deltaY

  eta.treat <- mean(att.treat) / mean(w.treat)
  eta.cont <- mean(att.cont) / mean(w.cont)

  ipw.att <- eta.treat - eta.cont
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
  inf.treat <- (att.treat - w.treat * eta.treat)/mean(w.treat)
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.1 <- (att.cont - w.cont * eta.cont)
  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2 <- base::colMeans(w.cont *(deltaY - eta.cont) * int.cov)
  # Now the influence function related to estimation effect of pscores
  inf.cont.2 <- asy.lin.rep.ps %*% M2

  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)

  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- inf.treat - inf.control
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
      ipw.boot <- unlist(lapply(1:nboot, wboot.std.ipw.panel,
                                n = n, deltaY = deltaY, D = D, int.cov = int.cov, i.weights = i.weights))
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

  return(list(ATT = ipw.att,
              se = se.att,
              uci = uci,
              lci = lci,
              boots = ipw.boot))
}
