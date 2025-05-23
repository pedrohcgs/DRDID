#' @import stats
NULL
###################################################################################
# Standardized version of Abadie's IPW DiD estimator

#' Standardized inverse probability weighted DiD estimator, with panel data
#' @description \code{std_ipw_did_panel} is used to compute inverse probability weighted (IPW) estimators for the ATT
#'  in difference-in-differences (DiD) setups with panel data. IPW weights are normalized to sum up to one, that is,
#'  the estimator is of the Hajek type.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score estimation. Please add a column of ones if you want to include an intercept in the model.
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
#' \item{ATT}{The IPW DiD point estimate.}
#' \item{se}{ The IPW DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, normalized = TRUE, boot, boot.type, nboot, type="ipw")}

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
#' # Form the Lalonde sample with CPS comparison group
#' eval_lalonde_cps <- subset(nsw, nsw$treated == 0 | nsw$sample == 2)
#' # Further reduce sample to speed example
#' set.seed(123)
#' unit_random <- sample(1:nrow(eval_lalonde_cps), 5000)
#' eval_lalonde_cps <- eval_lalonde_cps[unit_random,]
#' # Select some covariates
#' covX = as.matrix(cbind(1, eval_lalonde_cps$age, eval_lalonde_cps$educ,
#'                        eval_lalonde_cps$black, eval_lalonde_cps$married,
#'                        eval_lalonde_cps$nodegree, eval_lalonde_cps$hisp,
#'                        eval_lalonde_cps$re74))
#' # Implement normalized IPW DiD with panel data
#' std_ipw_did_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'                 D = eval_lalonde_cps$experimental,
#'                 covariates = covX)
#'
#' @export

std_ipw_did_panel <-function(y1, y0, D, covariates, i.weights = NULL,
                             boot = FALSE, boot.type = "weighted", nboot = NULL,
                             inffunc = FALSE,
                             trim.level = 0.995){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # generate deltaY
  deltaY <- as.vector(y1 - y0)
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
  if(PS$converged == FALSE){
    warning("Propernsity score estimation did not converge.")
  }
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
  w.treat <- trim.ps * i.weights * D
  w.cont <- trim.ps * i.weights * ps.fit * (1 - D)/(1 - ps.fit)

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
  #Hessian.ps <- stats::vcov(PS) * n
  Hessian.ps <- chol2inv(chol(t(int.cov) %*% (W * int.cov))) * n
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
  # Batch multiple matrix multiplications for inf functions
  #batch_results <- batch_matrix_operations(NULL, NULL, score.ps, Hessian.ps, NULL, M2, NULL)
  inf.cont.2 <- asy.lin.rep.ps %*% M2
  #inf.cont.2 <- batch_results$inf_cont_2

  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2) / mean(w.cont)

  #get the influence function of the DR estimator (put all pieces together)
  att.inf.func <- inf.treat - inf.control
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.att <- stats::sd(att.inf.func)/sqrt(n)
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
      ipw.boot <- unlist(lapply(1:nboot, wboot.std.ipw.panel,
                                n = n, deltaY = deltaY, D = D, int.cov = int.cov,
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
    panel = TRUE,
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
