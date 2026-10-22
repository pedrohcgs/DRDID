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
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation. Please add a vector of constants if you want to include an intercept in the models.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' Exactly zero centered bootstrap draws give a zero standard error and point interval.
#' Nonfinite draws or a zero IQR with nonzero draws stop with an error.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot A finite scalar integer of at least 2 giving the number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The DiD point estimate.}
#'  \item{se}{The DiD standard error.}
#'  \item{uci}{The upper bound of the 95\% CI for the ATT.}
#'  \item{lci}{The lower bound of the 95\% CI for the ATT}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is \code{NULL}.}
#'  \item{ps.flag}{Source of the verified IPT propensity-score fit: =0 for the \code{trust}
#'    candidate, =1 for the secondary IPT candidate. A candidate must satisfy the original
#'    exponential-IPT equations to numerical tolerance; optimizer status alone is insufficient.}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, estMethod = "imp", boot, boot.type, nboot, type="dr")}
#' @param trim.level The level of trimming for the propensity score. Default is 0.995.
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
#' Propensity-score fitting must satisfy the original exponential-IPT equations; numerically
#' ill-conditioned designs are rejected. If neither
#' fitting route provides a verified solution, estimation stops rather than substituting
#' a logit MLE. To explicitly use the traditional estimator, call \code{drdid_panel()}
#' (or \code{drdid(..., estMethod = "trad")}). This changes both propensity-score and
#' outcome-regression fitting. Each weighted-bootstrap replicate must also have a verified
#' IPT fit; a failed replicate stops the bootstrap without skipping or redrawing it.
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
#' Under the conditions in Sant'Anna and Zhao (2020), with trimming and
#' propensity-score clipping inactive, the resulting estimator is locally efficient,
#' doubly robust for the ATT, and doubly robust for inference.
#' The weighted bootstrap refits the implemented point statistic. After active
#' trimming or clipping, the IPT/WLS balance identities need not hold for the
#' modified statistic; matching bootstrap refits alone does not establish these
#' double-robustness or influence-function guarantees.
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
#' covX = as.matrix(cbind(1, eval_lalonde_cps$age, eval_lalonde_cps$educ,
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
                           nboot = NULL, inffunc = FALSE, trim.level = 0.995){
  if (boot) nboot <- validate_bootstrap_count(nboot)

  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  validate_treatment_indicator(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
  # Covariate Vector
  if(is.null(covariates)){
    int.cov <- as.matrix(rep(1,n))
  } else{
    int.cov <- as.matrix(covariates)
  }

  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  }
  validate_sampling_weights(i.weights, n)

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
  outcome.reg <- wols.br.panel(deltaY, D, int.cov, ps.fit, i.weights)
  out.delta <-  as.vector(outcome.reg$out.reg)

  #Compute Bias-Reduced Doubly Robust DiD estimators
  dr.att.summand.num <- as.vector(trim.ps *(1 - (1 - D)/(1 - ps.fit)) * (deltaY - out.delta))
  dr.att <- mean(i.weights * dr.att.summand.num)/mean(D * i.weights)

  #get the influence function to compute standard error
  dr.att.inf.func <- as.vector(i.weights * trim.ps * (dr.att.summand.num - D * dr.att) / mean(D * i.weights))

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
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      dr.boot <- mboot.did(dr.att.inf.func, nboot)
      # get bootstrap std errors based on IQR
      boot.summary <- bootstrap_inference(dr.boot, dr.att)
      se.dr.att <- boot.summary$se
      uci <- boot.summary$uci
      lci <- boot.summary$lci
    } else {
      # do weighted bootstrap
      dr.boot <- unlist(lapply(seq_len(nboot),
                               wboot.dr.imp.panel,
                               n = n,
                               deltaY = deltaY,
                               D = D,
                               int.cov = int.cov,
                               i.weights = i.weights,
                               trim.level = trim.level))
      # get bootstrap std errors based on IQR
      boot.summary <- bootstrap_inference((dr.boot - dr.att), dr.att)
      se.dr.att <- boot.summary$se
      uci <- boot.summary$uci
      lci <- boot.summary$lci

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
