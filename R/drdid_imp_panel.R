#' @import stats trust
NULL
###################################################################################
# Improved Doubly Robust DiD estimator with panel data

#' Improved locally efficient doubly robust DiD estimator for the ATT, with panel data
#'
#' @description \code{drdid_imp_panel} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with panel data. The resulting estimator is also doubly robust
#'  for inference; see Section 3.1 of Sant'Anna and Zhao (2020).
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The DiD point estimate.}
#'  \item{se}{The DiD standard error.}
#'  \item{uci}{The upper bound of the 95\% CI for the ATT.}
#'  \item{lci}{The lower bound of the 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is \code{NULL}.}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation: =0 if \code{trust} algorithm converged,
#'    =1 if IPW algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, estMethod = "imp", boot, boot.type, nboot, type="dr")}
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
#'
#' @details
#'
#' The \code{drdid_imp_panel} function implements the locally efficient doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.1)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of a linear regression model for the outcome evolution among the comparison units.
#'
#'
#' The nuisance parameters (propensity score and outcome regression parameters) are estimated using the methods
#' described in Sections 3.1 of Sant'Anna and Zhao (2020). In short, the propensity score parameters are estimated
#' using the inverse probability tilting estimator proposed by Graham, Pinto and Pinto (2012), and the outcome
#' regression coefficients are estimated using weighted least squares,where the weights depend on
#' the propensity score estimates; see Sant'Anna and Zhao (2020) for details.
#'
#'
#' The resulting estimator is not only locally efficient and doubly robust for the ATT,
#' but it is also doubly robust for inference; see Sant'Anna and Zhao (2020) for details.
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
#'                              eval_lalonde_cps$black, eval_lalonde_cps$married,
#'                              eval_lalonde_cps$nodegree, eval_lalonde_cps$hisp,
#'                              eval_lalonde_cps$re74))
#'
#' # Implement improved DR locally efficient DiD with panel data
#' drdid_imp_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'                 D = eval_lalonde_cps$experimental,
#'                 covariates = covX)
#'

#' @export

drdid_imp_panel <-function(y1, y0, D, covariates, i.weights = NULL, boot = FALSE, boot.type = "weighted",
                           nboot = NULL, inffunc = FALSE){
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
  #Compute the Pscore using the pscore.cal
  pscore.ipt <- pscore.cal(D, int.cov, i.weights = i.weights, n = n)
  ps.fit <- as.vector(pscore.ipt$pscore)
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  #Compute the Outcome regression for the control group
  outcome.reg <- wols.br.panel(deltaY, D, int.cov, ps.fit, i.weights)
  out.delta <-  as.vector(outcome.reg$out.reg)

  #Compute Bias-Reduced Doubly Robust DiD estimators
  dr.att.summand.num <- as.vector((1 - (1 - D)/(1 - ps.fit)) * (deltaY - out.delta))
  dr.att <- mean(i.weights * dr.att.summand.num)/mean(D * i.weights)

  #get the influence function to compute standard error
  dr.att.inf.func <- as.vector(i.weights * (dr.att.summand.num - D * dr.att) / mean(D * i.weights))

  if (boot == FALSE) {
    # Estimate of standard error
    se.dr.att <- stats::sd(dr.att.inf.func)/sqrt(n)
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


  if(inffunc == FALSE) dr.att.inf.func <- NULL
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  boot.type <- ifelse(argu$boot.type=="multiplier", "multiplier", "weighted")
  boot <- ifelse(argu$boot == TRUE, TRUE, FALSE)
  argu <- list(
    panel = TRUE,
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
              ps.flag = pscore.ipt$flag,
              att.inf.func = dr.att.inf.func,
              call.param = call.param,
              argu = argu))

  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)

}
