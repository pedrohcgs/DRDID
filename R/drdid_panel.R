#' @import stats
NULL
###################################################################################
#  Locally Efficient Doubly Robust DiD estimator with panel Data
#' Locally efficient doubly robust DiD estimator for the ATT, with panel data
#'
#' @description \code{drdid_panel} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups with panel data.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the propensity score and regression estimation. Please add a vector of constants if you want to include an intercept in the models.
#' If covariates = NULL, this leads to an unconditional DiD estimator.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DiD point estimate.}
#' \item{se}{ The DR DiD standard error.}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT.}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT.}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL.}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL.}
#'  \item{call.param}{The matched call.}
#'  \item{argu}{Some arguments used (explicitly or not) in the call (panel = TRUE, estMethod = "trad", boot, boot.type, nboot, type="dr")}
#'
#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#' \doi{10.1016/j.jeconom.2020.06.003}}
#' }
#' @details
#'
#' The \code{drdid_panel} function implements the locally efficient doubly robust difference-in-differences (DiD)
#' estimator for the average treatment effect on the treated (ATT) defined in equation (3.1)
#' in Sant'Anna and Zhao (2020). This estimator makes use of a logistic propensity score model for the probability
#' of being in the treated group, and of a linear regression model for the outcome evolution among the comparison units.
#'
#'
#' The propensity score parameters are estimated using maximum
#' likelihood, and the outcome regression coefficients are estimated using ordinary least squares.
#'
#'
#' @examples
#' # Form the Lalonde sample with CPS comparison group (data in wide format)
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
#'
#' # Implement traditional DR locally efficient DiD with panel data
#' drdid_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'              D = eval_lalonde_cps$experimental,
#'              covariates = covX)
#'
#' @export

drdid_panel <-function(y1, y0, D, covariates, i.weights = NULL,
                       boot = FALSE, boot.type =  "weighted", nboot = NULL,
                       inffunc = FALSE){
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

  # int.cov <- as.matrix(rep(1,n))
  # if (!is.null(covariates)){
  #   if(all(as.matrix(covariates)[,1]==rep(1,n))){
  #     int.cov <- as.matrix(covariates)
  #   } else {
  #     int.cov <- as.matrix(cbind(1, covariates))
  #   }
  # }

  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")

  # Normalize weights
  i.weights <- i.weights/mean(i.weights)
  #-----------------------------------------------------------------------------
  #Compute the Pscore by MLE
  pscore.tr <- suppressWarnings(fastglm::fastglm(
                                x = int.cov,
                                y = D,
                                family = stats::binomial(),
                                weights = i.weights,
                                intercept = FALSE,
                                method = 3
  ))
  class(pscore.tr) <- "glm" #this allow us to use vcov
  if(pscore.tr$converged == FALSE){
    warning("Propernsity score estimation did not converge.")
  }
  if(anyNA(pscore.tr$coefficients)){
    stop("Propensity score model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  ps.fit <- fitted(pscore.tr) #as.vector(pscore.tr$fitted.values)
  # Avoid divide by zero
  ps.fit <- pmin(ps.fit, 1 - 1e-6)
  W <- ps.fit * (1 - ps.fit) * i.weights
  # Compute the Outcome regression for the control group using wols
  control_filter <- (D == 0)
  reg.coeff <- stats::coef(fastglm::fastglm(
    x = int.cov[control_filter, , drop = FALSE],
    y = deltaY[control_filter],
    weights = i.weights[control_filter],
    family =  stats::gaussian(link = "identity")
  ))
  if(anyNA(reg.coeff)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason.")
  }
  out.delta <-   as.vector(tcrossprod(reg.coeff, int.cov))
  #-----------------------------------------------------------------------------
  #Compute Traditional Doubly Robust DiD estimators
  # First, the weights
  w.treat <- i.weights * D
  w.cont <- i.weights * ps.fit * (1 - D)/(1 - ps.fit)
  dr.att.treat <- w.treat * (deltaY - out.delta)
  dr.att.cont <- w.cont * (deltaY - out.delta)

  eta.treat <- mean(dr.att.treat) / mean(w.treat)
  eta.cont <- mean(dr.att.cont) / mean(w.cont)

  dr.att <-   eta.treat - eta.cont
  #-----------------------------------------------------------------------------
  #get the influence function to compute standard error
  #-----------------------------------------------------------------------------
  # First, the influence function of the nuisance functions
  # Asymptotic linear representation of OLS parameters
  weights.ols <- i.weights * (1 - D)
  wols.x <- weights.ols * int.cov
  wols.eX <- weights.ols * (deltaY - out.delta) * int.cov
  #XpX <- opt_crossprod(wols.x, int.cov, n)
  XpX <- base::crossprod(wols.x, int.cov)/n
  # Check if XpX is invertible
  if ( base::rcond(XpX) < .Machine$double.eps) {
    stop("The regression design matrix is singular. Consider removing some covariates.")
  }
  XpX.inv <- solve(XpX)
  asy.lin.rep.wols <-  wols.eX %*% XpX.inv

  # Asymptotic linear representation of logit's beta's
  score.ps <- i.weights * (D - ps.fit) * int.cov
  #Hessian.ps <- solve(t(int.cov) %*% (W * int.cov)) * n
  Hessian.ps <- chol2inv(chol(t(int.cov) %*% (W * int.cov))) * n
  asy.lin.rep.ps <-  score.ps %*% Hessian.ps


  # Now, the influence function of the "treat" component
  # Leading term of the influence function: no estimation effect
  inf.treat.1 <- (dr.att.treat - w.treat * eta.treat)
  # Estimation effect from beta hat
  # Derivative matrix (k x 1 vector)
  M1 <- base::colMeans(w.treat * int.cov)

  # Now get the influence function related to the estimation effect related to beta's
  inf.treat.2 <- asy.lin.rep.wols %*% M1

  # Influence function for the treated component
  inf.treat <- (inf.treat.1 - inf.treat.2) / mean(w.treat)
  #-----------------------------------------------------------------------------
  # Now, get the influence function of control component
  # Leading term of the influence function: no estimation effect
  inf.cont.1 <- (dr.att.cont - w.cont * eta.cont)
  # Estimation effect from gamma hat (pscore)
  # Derivative matrix (k x 1 vector)
  M2 <- base::colMeans(w.cont *(deltaY - out.delta - eta.cont) * int.cov)
  # Now the influence function related to estimation effect of pscores
  inf.cont.2 <- asy.lin.rep.ps %*% M2
  # Estimation Effect from beta hat (weighted OLS)
  M3 <-  base::colMeans(w.cont * int.cov)
  inf.cont.3 <- asy.lin.rep.wols %*% M3

  # # Batch multiple matrix multiplications for inf functions
  # batch_results <- batch_matrix_operations(wols.eX, XpX.inv, score.ps, Hessian.ps, M1, M2, M3)
  # # Now get the influence function related to the estimation effect related to beta's
  # inf.treat.2 <- batch_results$inf_treat_2
  # # Now the influence function related to estimation effect of pscores
  # inf.cont.2 <- batch_results$inf_cont_2
  # # Now the influence function related to estimation effect of regressions
  # inf.cont.3 <- batch_results$inf_cont_3

  # Influence function for the treated component
  inf.treat <- (inf.treat.1 - inf.treat.2) / mean(w.treat)

  # Influence function for the control component
  inf.control <- (inf.cont.1 + inf.cont.2 - inf.cont.3) / mean(w.cont)

  #get the influence function of the DR estimator (put all pieces together)
  dr.att.inf.func <- inf.treat - inf.control
  #-----------------------------------------------------------------------------
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
      # get symmetric critical values
      cv <- stats::quantile(abs(dr.boot/se.dr.att), probs = 0.95)
      # Estimate of upper boundary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower boundary of 95% CI
      lci <- dr.att - cv * se.dr.att
    } else {
      # do weighted bootstrap
      dr.boot <- unlist(lapply(1:nboot, wboot.dr.tr.panel,
                               n = n, deltaY = deltaY, D = D, int.cov = int.cov, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.dr.att <- stats::IQR((dr.boot - dr.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmetric critical values
      cv <- stats::quantile(abs((dr.boot - dr.att)/se.dr.att), probs = 0.95)
      # Estimate of upper boundary of 95% CI
      uci <- dr.att + cv * se.dr.att
      # Estimate of lower boundary of 95% CI
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
    estMethod = "trad",
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
               att.inf.func = dr.att.inf.func,
               call.param = call.param,
               argu = argu))
  # Define a new class
  class(ret) <- "drdid"

  # return the list
  return(ret)
}
