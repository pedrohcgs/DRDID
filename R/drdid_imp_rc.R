#' @import stats
NULL
###################################################################################
# 'Improved' and locally efficient Doubly Robust DiD estimator with Repeated Cross Section Data
#' Improved locally efficient doubly robust DiD estimator for the ATT, with repeated cross-section data
#'
#' @description \code{drdid_imp_rc} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with stationary repeated cross-sectional data. The resulting estimator is
#'  also doubly robust for inference; see Section 3.2 of Sant'Anna and Zhao (2020).
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation. Please add a vector of constants if you want to include an intercept in the models.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#' @param trim.level The level of trimming for the propensity score. Default is 0.995.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DiD point estimate}
#' \item{se}{ The DR DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation: =0 if \code{trust} algorithm converged,
#'    =1 if IPW algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = FALSE, estMethod = "imp", boot, boot.type, nboot, type="dr")}
#'
#' @references
#' \cite{Graham, Bryan, Pinto, Cristine, and Egel, Daniel (2012),
#' "Inverse Probability Tilting for Moment Condition Models with Missing Data."
#'  Review of Economic Studies, vol. 79 (3), pp. 1053-1079, \doi{10.1093/restud/rdr047}}
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#'
#' @details
#'
#' The \code{drdid_imp_rc} function implements the locally efficient doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.4)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of (separate) linear regression models for the outcome of both
#' treated and comparison units, in both pre and post-treatment periods.
#'
#'
#' The nuisance parameters (propensity score and outcome regression parameters) are estimated using the methods
#' described in Sections 3.2 of Sant'Anna and Zhao (2020). In short, the propensity score parameters are estimated
#' using the inverse probability tilting estimator proposed by Graham, Pinto and Pinto (2012), and the outcome
#' regression coefficients are estimated using weighted least squares,where the weights depend on
#' the propensity score estimates; see Sant'Anna and Zhao (2020) for details.
#'
#'
#' The resulting estimator is not only locally efficient and doubly robust for the ATT,
#' but it is also doubly robust for inference; see Sant'Anna and Zhao (2020) for details.
#'
#' @examples
#' # use the simulated data
#' covX = as.matrix(cbind(1, sim_rc[,5:8]))
#' # Implement the improved, locally efficient DR DiD estimator
#' drdid_imp_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'              covariates= covX)
#'
#' @export

drdid_imp_rc <- function(y, post, D, covariates, i.weights = NULL, boot = FALSE,
                         boot.type =  "weighted",  nboot = NULL, inffunc = FALSE,
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
  #Compute the Pscore using the pscore.cal
  pscore.ipt <- pscore.cal(D, int.cov, i.weights = i.weights, n = n)
  ps.fit <- as.vector(pscore.ipt$pscore)
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  trim.ps <- (ps.fit < 1.01)
  trim.ps[D==0] <- (ps.fit[D==0] < trim.level)
  #Compute the Outcome regression for the control group
  out.y.cont.pre <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = TRUE, treat = FALSE)
  out.y.cont.pre <-  as.vector(out.y.cont.pre$out.reg)
  out.y.cont.post <- wols_rc(y, post, D, int.cov, ps.fit, i.weights, pre = FALSE, treat = FALSE)
  out.y.cont.post <-  as.vector(out.y.cont.post$out.reg)
  # Combine the ORs
  out.y.cont <- post * out.y.cont.post + (1 - post) * out.y.cont.pre
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the treated group at the pre-treatment period, using ols.
  # reg.treat.coeff.pre <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                              subset = ((D==1) & (post==0)),
  #                                              weights = i.weights))
  treat_pre_filter <- (D == 1) & (post == 0)
  reg.treat.coeff.pre <- stats::coef(fastglm::fastglm(
                                      x = int.cov[treat_pre_filter, , drop = FALSE],
                                      y = y[treat_pre_filter],
                                      weights = i.weights[treat_pre_filter],
                                      family = gaussian(link = "identity")
  ))
  out.y.treat.pre <-   as.vector(tcrossprod(reg.treat.coeff.pre, int.cov))
  #Compute the Outcome regression for the treated group at the post-treatment period, using ols.
  # reg.treat.coeff.post <- stats::coef(stats::lm(y ~ -1 + int.cov,
  #                                               subset = ((D==1) & (post==1)),
  #                                               weights = i.weights))
  treat_post_filter <- (D == 1) & (post == 1)
  reg.treat.coeff.post <- stats::coef(fastglm::fastglm(
                                      x = int.cov[treat_post_filter, , drop = FALSE],
                                      y = y[treat_post_filter],
                                      weights = i.weights[treat_post_filter],
                                      family = gaussian(link = "identity")
  ))
  out.y.treat.post <-   as.vector(tcrossprod(reg.treat.coeff.post, int.cov))
  #-----------------------------------------------------------------------------
  # First, the weights
  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <-  i.weights * D * post
  w.cont.pre <- trim.ps * i.weights * ps.fit * (1 - D) * (1 - post)/(1 - ps.fit)
  w.cont.post <- trim.ps * i.weights * ps.fit * (1 - D) * post/(1 - ps.fit)

  w.d <-  i.weights * D
  w.dt1 <-  i.weights * D * post
  w.dt0 <-  i.weights * D * (1 - post)

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
  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)*sqrt(n-1)/(n)
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
      # get symmetric critical values
      cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
      # Estimate of upper bound of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower bound of 95% CI
      lci <- dr.att - cv * se.dr.att
    } else {
      # do weighted bootstrap
      dr.boot <- unlist(lapply(1:nboot, wboot_drdid_imp_rc,
                               n = n, y = y, post = post,
                               D = D, int.cov = int.cov, i.weights = i.weights,
                               trim.level = trim.level))
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
    estMethod = "imp",
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    type = "dr",
    trim.level = trim.level
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
  # return the list
  return(ret)


}
