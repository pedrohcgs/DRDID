#' @import stats
NULL
###################################################################################
# Standardized version of Abadie's IPW DID estimator

#' Standardized inverse probability weighted DiD estimator, with repeated cross-section data
#' @description \code{std_ipw_did_rc} is used to compute inverse probability weighted (IPW) estimators for the ATT
#'  in DID setups with stationary repeated cross-sectional data. IPW weights are normalized to sum up to one, that is,
#'  the estimator is of the Hajek type.
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation. Please add a column of ones if you want to include an intercept.
#' If covariates = NULL, this leads to an unconditional DID estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#' @param trim.level The level of trimming for the propensity score. Default is 0.995.
#'
#' @return A list containing the following components:
#' \item{ATT}{The IPW DID point estimate.}
#' \item{se}{ The IPW DID standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, normalized = TRUE, boot, boot.type, nboot, type="ipw")}

#' @references
#' \cite{Abadie, Alberto (2005), "Semiparametric Difference-in-Differences Estimators",
#' Review of Economic Studies, vol. 72(1), p. 1-19, \doi{10.1111/0034-6527.00321}.
#' }
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#'
#'
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(cbind(1, sim_rc[,5:8]))
#' # Implement normalized IPW DID estimator
#' std_ipw_did_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'                covariates= covX)
#'
#' @export

std_ipw_did_rc <-function(y, post, D, covariates, i.weights = NULL,
                          boot = FALSE, boot.type = "weighted", nboot = NULL,
                          inffunc = FALSE,
                          trim.level = 0.995){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # y as vector
  y <- as.vector(y)
  # post as vector
  post <- as.vector(post)
  # Covariate vector
  if(is.null(covariates)){
    int.cov <- as.matrix(rep(1,n))
  } else{
    int.cov <- as.matrix(covariates)
  }


  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  # Normalize weights
  i.weights <- i.weights/mean(i.weights)
  #-----------------------------------------------------------------------------
  # Pscore estimation (logit) and also its fitted values
  PS <- suppressWarnings(fastglm::fastglm(
                                  x = int.cov,
                                  y = D,
                                  family = stats::binomial(),
                                  weights = i.weights,
                                  intercept = FALSE,
                                  method = 3
  ))
  ps.fit <- fitted(PS)
  # Do not divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  W <- ps.fit * (1 - ps.fit) * i.weights
  # Trim the propensity score
  trim.ps <- (ps.fit < 1.01)
  trim.ps[D==0] <- (ps.fit[D==0] < trim.level)
  #-----------------------------------------------------------------------------
  #Compute IPW estimator
  # First, the weights
  w.treat.pre <- trim.ps * i.weights * D * (1 - post)
  w.treat.post <-  trim.ps * i.weights * D * post
  w.cont.pre <- trim.ps *  i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- trim.ps *  i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

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
  #Hessian.ps <- stats::vcov(PS) * n
  Hessian.ps <- chol2inv(chol(t(int.cov) %*% (W * int.cov))) * n
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
  if (boot == FALSE) {
    # Estimate of standard error
    se.att <- stats::sd(att.inf.func)*sqrt(n-1)/(n)
    # Estimate of upper boundary of 95% CI
    uci <- ipw.att + 1.96 * se.att
    # Estimate of lower boundary of 95% CI
    lci <- ipw.att - 1.96 * se.att

    #Create this null vector so we can export the bootstrap draws too.
    ipw.boot <- NULL
  }
  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      ipw.boot <- mboot.did(att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR(ipw.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmetric critical values
      cv <- stats::quantile(abs(ipw.boot/se.att), probs = 0.95)
      # Estimate of upper boundary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower boundary of 95% CI
      lci <- ipw.att - cv * se.att
    } else {
      # do weighted bootstrap
      ipw.boot <- unlist(lapply(1:nboot, wboot_std_ipw_rc,
                                n = n, y = y, post = post, D = D, int.cov = int.cov,
                                i.weights = i.weights,
                                trim.level = trim.level))
      # get bootstrap std errors based on IQR
      se.att <- stats::IQR(ipw.boot - ipw.att) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmetric critical values
      cv <- stats::quantile(abs((ipw.boot - ipw.att)/se.att), probs = 0.95)
      # Estimate of upper boundary of 95% CI
      uci <- ipw.att + cv * se.att
      # Estimate of lower boundary of 95% CI
      lci <- ipw.att - cv * se.att

    }

  }
  if(inffunc == FALSE) att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = FALSE,
    normalized = TRUE,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "ipw",
    trim.level = trim.level
  )
  ret <- (list(ATT = ipw.att,
               se = se.att,
               uci = uci,
               lci = lci,
               boots = ipw.boot,
               att.inf.func = att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
