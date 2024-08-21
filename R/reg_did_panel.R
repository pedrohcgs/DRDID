#' @import stats
NULL
###################################################################################
#' Outcome regression DiD estimator for the ATT, with panel data
#'
#' @description \code{reg_did_panel} computes the outcome regressions estimators for the average treatment effect on the
#' treated in difference-in-differences (DiD) setups with panel data.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The OR DiD point estimate}
#'  \item{se}{The OR DiD standard error}
#'  \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#'  \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, boot, boot.type, nboot, type="or")}
#'
#' @details
#'
#' The \code{reg_did_panel} function implements
#' outcome regression difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) defined in equation (2.2) of Sant'Anna and Zhao (2020) when panel data are available.
#' The estimator follows the same spirit of the nonparametric estimators proposed by Heckman, Ichimura and Todd (1997),
#' though here the the outcome regression models are assumed to be linear in covariates (parametric),
#'
#' The nuisance parameters (outcome regression coefficients) are estimated via ordinary least squares.


#' @references
#' \cite{Heckman, James J., Ichimura, Hidehiko, and Todd, Petra E. (1997),"Matching as an Econometric Evaluation Estimator: Evidence from Evaluating a Job Training Programme",
#' Review of Economic Studies, vol. 64(4), p. 605–654, \doi{10.2307/2971733}.
#' }
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#'
#'
#' @examples
#' # Form the Lalonde sample with CPS comparison group
#' eval_lalonde_cps <- subset(nsw, nsw$treated == 0 | nsw$sample == 2)
#' # Further reduce sample to speed example
#' set.seed(123)
#' unit_random <- sample(1:nrow(eval_lalonde_cps), 5000)
#' eval_lalonde_cps <- eval_lalonde_cps[unit_random,]
#' # Select some covariates
#' covX = as.matrix(cbind(eval_lalonde_cps$age, eval_lalonde_cps$educ,
#'                        eval_lalonde_cps$black, eval_lalonde_cps$married,
#'                        eval_lalonde_cps$nodegree, eval_lalonde_cps$hisp,
#'                        eval_lalonde_cps$re74))
#' # Implement OR DiD with panel data
#' reg_did_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'                D = eval_lalonde_cps$experimental,
#'                covariates = covX)
#'


#' @export

reg_did_panel <-function(y1, y0, D, covariates, i.weights = NULL,
                         boot = FALSE, boot.type = "weighted", nboot = NULL,
                         inffunc = FALSE){
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
  # Normalize weights
  i.weights <- i.weights/mean(i.weights)
  #-----------------------------------------------------------------------------
  #Compute the Outcome regression for the control group using ols.
  # reg.coeff <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
  #                                    subset = D==0,
  #                                    weights = i.weights))
  control_filter <- (D == 0)
  reg.coeff <- stats::coef(fastglm::fastglm(
                            x = int.cov[control_filter, , drop = FALSE],
                            y = deltaY[control_filter],
                            weights = i.weights[control_filter],
                            family = gaussian(link = "identity")
  ))
  if(anyNA(reg.coeff)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is probably the reason for it.")
  }
  out.delta <-   as.vector(tcrossprod(reg.coeff, int.cov))
  #-----------------------------------------------------------------------------
  #Compute the OR-DiD estimator
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
  XpX <- crossprod(wols.x, int.cov)/n
  # Check if XpX is invertible
  if ( base::rcond(XpX) < .Machine$double.eps) {
    stop("The regression design matrix is singular. Consider removing some covariates.")
  }
  XpX.inv <- solve(XpX)
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
  if (boot == FALSE) {
    # Estimate of standard error
    se.reg.att <- stats::sd(reg.att.inf.func)/sqrt(n)
    # Estimate of upper boudary of 95% CI
    uci <- reg.att + 1.96 * se.reg.att
    # Estimate of lower doundary of 95% CI
    lci <- reg.att - 1.96 * se.reg.att
    #Create this null vector so we can export the bootstrap draws too.
    reg.boot <- NULL
  }

  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
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

  if(inffunc == FALSE) reg.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = TRUE,
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
