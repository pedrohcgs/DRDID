#' @import stats
NULL
###################################################################################
# 'Improved' Doubly Robust DID estimator with Repeated Cross Section Data
#' 'Improved' Doubly Robust Difference-in-Differences Estimator for the ATT, with Repeated Cross Section Data
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
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
#' \item{ATT}{The DR DID point estimate}
#' \item{se}{ The DR DID standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation: =0 if \code{trust} algorithm converged,
#'    =1 if IPW algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = F, estMethod = "imp2", boot, boot.type, nboot, type="dr")}
#'
#' @references Sant'Anna, Pedro H. C. and Zhao, Jun (2020), ["Doubly Robust Difference-in-Differences Estimators"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3293315).
#' @export

drdid_imp_rc1 <- function(y, post, D, covariates, i.weights = NULL,
                          boot = F, boot.type =  "weighted", nboot = NULL,
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
  pscore.ipt <- pscore.cal(D, int.cov, i.weights = i.weights, n = n)
  ps.fit <- as.vector(pscore.ipt$pscore)
  ps.fit <- pmin(ps.fit, 1 - 1e-16)
  #Compute the Outcome regression for the control group
  out.y.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = T, treat = F)
  out.y.pre <-  as.vector(out.y.pre$out.reg)
  out.y.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = F, treat = F)
  out.y.post <-  as.vector(out.y.post$out.reg)
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
  #-----------------------------------------------------------------------------
  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.pre <- eta.treat.pre - w.treat.pre * att.treat.pre/mean(w.treat.pre)
  inf.treat.post <- eta.treat.post - w.treat.post * att.treat.post/mean(w.treat.post)
  # Influence function for the treated component
  inf.treat <- inf.treat.post - inf.treat.pre
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect from nuisance parameters
  inf.cont.pre <- eta.cont.pre - w.cont.pre * att.cont.pre/mean(w.cont.pre)
  inf.cont.post <- eta.cont.post - w.cont.post * att.cont.post/mean(w.cont.post)

  # Influence function for the control component
  inf.cont <- inf.cont.post - inf.cont.pre
  #-----------------------------------------------------------------------------
  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func <- inf.treat - inf.cont
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
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_imp_rc1,
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


  if(inffunc==F) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot==T, T, F)
  argu <- list(
    panel = F,
    estMethod = "imp2",
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
               ps.flag = pscore.ipt$flag,
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  #warning("This estimator is not locally efficient. Consider using 'drdid_imp_rc' instead.")
  # return the list
  return(ret)


}
