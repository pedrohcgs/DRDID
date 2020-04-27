#' @import stats trust
NULL
###################################################################################
# Improved Doubly Robust DID estimator with panel Data

#' Improved Doubly Robust Difference-in-Differences Estimator for the ATT, with Panel Data
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
#' @param nboot Number of bootstrap repetitions (not relevant if boot = FALSE). Default is 999 if boot = TRUE
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The DID point estimate.}
#'  \item{se}{The DID standard error.}
#'  \item{uci}{The upper bound of the 95\% CI for the ATT.}
#'  \item{lci}{The lower bound of the 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is \code{NULL}.}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation: =0 if \code{trust} algorithm converged,
#'    =1 if IPW algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = T, estMethod = "imp", boot, boot.type, nboot, type="dr")}
#'
#' @references Sant'Anna, Pedro H. C. and Zhao, Jun (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315).
#'
#' @export

drdid_imp_panel <-function(y1, y0, D, covariates, i.weights = NULL, boot = F, boot.type = "weighted",
                           nboot = NULL, inffunc = F){
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
  #Compute the Pscore using the pscore.cal
  pscore.br <- pscore.cal(D, int.cov, i.weights = i.weights, n = n)
  ps.fit <- as.vector(pscore.br$pscore)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  outcome.reg <- wols.br.panel(deltaY, D, int.cov, ps.fit, i.weights)
  out.delta <-  as.vector(outcome.reg$out.reg)

  #Compute Bias-Reduced Doubly Robust DID estimators
  dr.att.summand.num <- as.vector((1 - (1 - D)/(1 - ps.fit)) * (deltaY - out.delta))
  dr.att <- mean(i.weights * dr.att.summand.num)/mean(D * i.weights)

  #get the influence function to compute standard error
  dr.att.inf.func <- as.vector(i.weights * (dr.att.summand.num - D * dr.att) / mean(D * i.weights))

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
      dr.boot <- unlist(lapply(1:nboot, wboot.dr.imp.panel,
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


  if(inffunc==F) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot==T, T, F)
  argu <- list(
    panel = T,
    estMethod = "imp",
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
              ps.flag = pscore.br$flag,
              att.inf.func = dr.att.inf.func,
              call.param = call.param,
              argu = argu))

  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)

}
