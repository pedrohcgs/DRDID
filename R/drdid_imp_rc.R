#' @import stats
NULL
###################################################################################
# 'Improved' and locally efficient Doubly Robust DID estimator with Repeated Cross Section Data
#' 'Improved' and locally efficient Doubly Robust Difference-in-Differences Estimator for the ATT, with Repeated Cross Section Data
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-tretment period,
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
#' \item{uci}{Estimate of the upper boudary of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower boudary of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation: =0 if \code{trust} algorithm converged,
#'    =1 if IPW algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'
#'
#' @export

drdid_imp_rc <- function(y, post, D, covariates,
                         i.weights = NULL,
                         boot = F,
                         boot.type =  "weighted",
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
  int.cov <- as.matrix(rep(1,n))
  if (!is.null(covariates)){
    if(all.equal(as.matrix(covariates)[,1], rep(1,n))){
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
  out.y.cont.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = T, treat = F)
  out.y.cont.pre <-  as.vector(out.y.cont.pre$out.reg)
  out.y.cont.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = F, treat = F)
  out.y.cont.post <-  as.vector(out.y.cont.post$out.reg)
  # Combine the ORs
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                               subset = ((D==1) & (post==0)),
                                               weights = i.weights))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                                subset = ((D==1) & (post==1)),
                                                weights = i.weights))
  out.y.treat.post <-   as.vector(tcrossprod(reg.treat.coeff.post, int.cov))
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

  w.d <- i.weights * D
  w.dt1 <- i.weights * D * post
  w.dt0 <- i.weights * D * (1 - post)

  # Elements of the influence function (summands)
  eta.treat.pre <- w.treat.pre * (y - out.y.cont) / mean(w.treat.pre)
  eta.treat.post <- w.treat.post * (y - out.y.cont)/ mean(w.treat.post)
  eta.cont.pre <- w.cont.pre * (y - out.y.cont) / mean(w.cont.pre)
  eta.cont.post <- w.cont.post * (y - out.y.cont) / mean(w.cont.post)

  # extra elements for the locally efficient DRDID
  eta.d.post <- w.d * (out.y.treat.post - out.y.cont.post)/mean(w.d)
  eta.dt1.post <- w.dt1 * (out.y.treat.post - out.y.cont.post)/mean(w.dt1)
  eta.d.pre <- w.d * (out.y.treat.pre - out.y.cont.pre)/mean(w.d)
  eta.dt0.pre <- w.dt0 * (out.y.treat.pre - out.y.cont.pre)/mean(w.dt0)

  # Estimator of each component
  att.treat.pre <- mean(eta.treat.pre)
  att.treat.post <- mean(eta.treat.post)
  att.cont.pre <- mean(eta.cont.pre)
  att.cont.post <- mean(eta.cont.post)

  att.d.post <- mean(eta.d.post)
  att.dt1.post <- mean(eta.dt1.post)
  att.d.pre <- mean(eta.d.pre)
  att.dt0.pre <- mean(eta.dt0.pre)

  # ATT estimator
  dr.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)
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
  dr.att.inf.func1 <- inf.treat - inf.cont
  #-----------------------------------------------------------------------------
  # Now, we only need to get the influence function of the adjustment terms
  # First, the terms as if all OR parameters were known
  inf.eff1 <- eta.d.post - w.d * att.d.post/mean(w.d)
  inf.eff2 <- eta.dt1.post - w.dt1 * att.dt1.post/mean(w.dt1)
  inf.eff3 <- eta.d.pre - w.d * att.d.pre/mean(w.d)
  inf.eff4 <- eta.dt0.pre - w.dt0 * att.dt0.pre/mean(w.dt0)
  inf.eff <- (inf.eff1 - inf.eff2) - (inf.eff3 - inf.eff4)
  #-----------------------------------------------------------------------------
  #get the influence function of the locally efficient DR estimator (put all pieces together)
  dr.att.inf.func <- dr.att.inf.func1 + inf.eff
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
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_rc,
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


  if(inffunc==F) att.inf.func <- NULL
  return(list(ATT = dr.att,
              se = se.dr.att,
              uci = uci,
              lci = lci,
              boots = dr.boot,
              ps.flag = pscore.ipt$flag,
              att.inf.func = att.inf.func))
}
