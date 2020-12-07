#' @import stats
NULL
###################################################################################
#' Two-way fixed effects DiD estimator, with repeated cross-section data
#' @description \code{twfe_did_rc} is used to compute linear two-way fixed effects estimators for the ATT
#'  in difference-in-differences (DiD) setups with stationary repeated cross-sectional data. As illustrated
#'  by Sant'Anna and Zhao (2020),this estimator generally do not recover the ATT. We encourage empiricists
#'  to adopt alternative specifications.
#'
#' @param y An \eqn{n} x \eqn{1} vector of outcomes from the both pre and post-treatment periods.
#' @param post An \eqn{n} x \eqn{1} vector of Post-Treatment dummies (post = 1 if observation belongs to post-treatment period,
#'             and post = 0 if observation belongs to pre-treatment period.)
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment period, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The TWFE DID point estimate}
#'  \item{se}{The TWFE DID standard error}
#'  \item{uci}{Estimate of the upper bound of a 95\% CI for the TWFE parameter.}
#'  \item{lci}{Estimate of the lower bound of a 95\% CI for the TWFE parameter.}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'
#' @examples
#' # use the simulated data provided in the package
#' covX = as.matrix(sim_rc[,5:8])
#' # Implement TWFE DID estimator (you probably should consider something else....)
#' twfe_did_rc(y = sim_rc$y, post = sim_rc$post, D = sim_rc$d,
#'            covariates= covX)
#'
#' @export

twfe_did_rc <- function(y, post, D, covariates = NULL, i.weights = NULL,
                        boot = FALSE, boot.type = "weighted", nboot = NULL,
                        inffunc = FALSE){
  #-----------------------------------------------------------------------------
  # D as vector
  D <- as.vector(D)
  # Sample size
  n <- length(D)
  # Weights
  if(is.null(i.weights)) {
    i.weights <- as.vector(rep(1, n))
  } else if(min(i.weights) < 0) stop("i.weights must be non-negative")
  #-----------------------------------------------------------------------------
  #Create dataset for TWFE approach
  if (is.null(covariates)) {
    x = NULL
  } else {
    if(all(as.matrix(covariates)[,1] == rep(1,n))) {
      # Remove intercept if you include it
      covariates <- as.matrix(covariates)
      covariates <- covariates[,-1]
      if(dim(covariates)[2]==0) {
        covariates = NULL
        x = NULL
      }
    }
  }

  if(!is.null(covariates))  x <- as.matrix(covariates)
  dd <- D
  post <- post
  i.weights <- as.vector(i.weights)
  #---------------------------------------------------------------------------
  #Estimate TWFE regression
  if(!is.null(x)){
    reg <- stats::lm(y ~  dd:post + post + dd + x, x = TRUE, weights = i.weights)
  }
  if(is.null(x)){
    reg <- stats::lm(y ~  dd:post + post + dd, x = TRUE, weights = i.weights)
  }
  twfe.att <- reg$coefficients["dd:post"]
  #-----------------------------------------------------------------------------
  #Elemenets for influence functions
  inf.reg <- (i.weights * reg$x * reg$residuals) %*%
    base::solve(crossprod(i.weights * reg$x, reg$x) / dim(reg$x)[1])

  sel.theta <- matrix(c(rep(0, dim(inf.reg)[2])))

  index.theta <- which(dimnames(reg$x)[[2]]=="dd:post",
                       arr.ind = TRUE)

  sel.theta[index.theta, ] <- 1
  #-----------------------------------------------------------------------------
  #get the influence function of the TWFE regression
  twfe.inf.func <- as.vector(inf.reg %*% sel.theta)
  #-----------------------------------------------------------------------------
  if (boot == FALSE) {
    # Estimate of standard error
    se.twfe.att <- stats::sd(twfe.inf.func)/sqrt(length(twfe.inf.func))
    # Estimate of upper boudary of 95% CI
    uci <- twfe.att + 1.96 * se.twfe.att
    # Estimate of lower doundary of 95% CI
    lci <- twfe.att - 1.96 * se.twfe.att
    #Create this null vector so we can export the bootstrap draws too.
    twfe.boot <- NULL
  }

  if (boot == TRUE) {
    if (is.null(nboot) == TRUE) nboot = 999
    if(boot.type == "multiplier"){
      # do multiplier bootstrap
      twfe.boot <- mboot.did(twfe.inf.func, nboot)
      # get bootstrap std errors based on IQR
      se.twfe.att <- stats::IQR(twfe.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs(twfe.boot/se.twfe.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- twfe.att + cv * se.twfe.att
      # Estimate of lower doundary of 95% CI
      lci <- twfe.att - cv * se.twfe.att
    } else {
      # do weighted bootstrap
      twfe.boot <- unlist(lapply(1:nboot, wboot_twfe_rc,
                                 n = n, y = y, dd = dd, post = post, x = x, i.weights = i.weights))
      # get bootstrap std errors based on IQR
      se.twfe.att <- stats::IQR((twfe.boot - twfe.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
      # get symmtric critival values
      cv <- stats::quantile(abs((twfe.boot - twfe.att)/se.twfe.att), probs = 0.95)
      # Estimate of upper boudary of 95% CI
      uci <- twfe.att + cv * se.twfe.att
      # Estimate of lower doundary of 95% CI
      lci <- twfe.att - cv * se.twfe.att

    }
  }


  if(inffunc == FALSE) att.inf.func <- NULL
  return(list(ATT = twfe.att,
              se = se.twfe.att,
              uci = uci,
              lci = lci,
              boots = twfe.boot,
              att.inf.func = att.inf.func))
}
