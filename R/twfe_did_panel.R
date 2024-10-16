#' @import stats
NULL
###################################################################################
#' Two-way fixed effects DiD estimator, with panel data
#'
#' @description \code{twfe_did_panel} is used to compute linear two-way fixed effects estimators for the ATT
#'  in difference-in-differences (DiD) setups with panel data. As illustrated by Sant'Anna and Zhao (2020),
#'  this estimator generally do not recover the ATT. We encourage empiricists to adopt alternative specifications.
#'
#' @param y1 An \eqn{n} x \eqn{1} vector of outcomes from the post-treatment period.
#' @param y0 An \eqn{n} x \eqn{1} vector of outcomes from the pre-treatment period.
#' @param D An \eqn{n} x \eqn{1} vector of Group indicators (=1 if observation is treated in the post-treatment, =0 otherwise).
#' @param covariates An \eqn{n} x \eqn{k} matrix of covariates to be used in the regression estimation. We will always include an intercept.
#' @param i.weights An \eqn{n} x \eqn{1} vector of weights to be used. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is FALSE.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if \code{boot = FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is FALSE.
#'
#' @return A list containing the following components:
#'  \item{ATT}{The TWFE DiD point estimate}
#'  \item{se}{The TWFE DiD standard error}
#'  \item{uci}{Estimate of the upper bound of a 95\% CI for the TWFE parameter.}
#'  \item{lci}{Estimate of the lower bound of a 95\% CI for the TWFE parameter.}
#'  \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#'  \item{att.inf.func}{Estimate of the influence function. Default is NULL}
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
#' # Implement TWFE DiD with panel data
#' twfe_did_panel(y1 = eval_lalonde_cps$re78, y0 = eval_lalonde_cps$re75,
#'                D = eval_lalonde_cps$experimental,
#'                covariates = covX)
#'
#' @export

twfe_did_panel <-function(y1, y0, D, covariates, i.weights = NULL,
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
  # Normalize weights
  i.weights <- i.weights/mean(i.weights)
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

  if (!is.null(covariates)){
    if (ncol(as.matrix(covariates)) == 1) {
      x = as.matrix(c(covariates, covariates))
    } else {
      x <- as.matrix(rbind(covariates, covariates))
    }
  }

  # Post treatment indicator
  post <- as.vector(c(rep(0, length(y0)), rep(1,length(y1))))
  # treatment group
  dd <- as.vector((c(D, D)))
  # outcome
  y <- as.vector(c(y0, y1))
  # weights
  i.weights <- as.vector(c(i.weights, i.weights))

  # If there are covariates, proceed like this
  if(!is.null(x)){
    #---------------------------------------------------------------------------
    #Estimate TWFE regression

    reg <- stats::lm(y ~  dd:post + post + dd + x, x = TRUE, weights = i.weights)


    twfe.att <- reg$coefficients["dd:post"]
    #-----------------------------------------------------------------------------
    #Elemenets for influence functions
    XpX <- crossprod(i.weights * reg$x, reg$x) / dim(reg$x)[1]
    # Check if XpX is invertible
    if ( base::rcond(XpX) < .Machine$double.eps) {
      stop("The regression design matrix is singular. Consider removing some covariates.")
    }
    XpX.inv <- solve(XpX)

    inf.reg <- (i.weights * reg$x * reg$residuals) %*% XpX.inv

    sel.theta <- matrix(c(rep(0, dim(inf.reg)[2])))

    index.theta <- which(dimnames(reg$x)[[2]]=="dd:post",
                         arr.ind = TRUE)

    sel.theta[index.theta, ] <- 1
    #-----------------------------------------------------------------------------
    #get the influence function of the TWFE regression
    att.inf.func <- as.vector(inf.reg %*% sel.theta)
    #-----------------------------------------------------------------------------
    if (boot == FALSE) {
      # Estimate of standard error
      se.twfe.att <- stats::sd(att.inf.func)/sqrt(length(att.inf.func))
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
        twfe.boot <- mboot.twfep.did(n, att.inf.func, nboot)
        # get bootstrap std errors based on IQR
        se.twfe.att <- stats::IQR(twfe.boot) / (stats::qnorm(0.75) - stats::qnorm(0.25))
        # get symmetric critical values
        cv <- stats::quantile(abs(twfe.boot/se.twfe.att), probs = 0.95)
        # Estimate of the upper boundary of 95% CI
        uci <- twfe.att + cv * se.twfe.att
        # Estimate of lower boundary of 95% CI
        lci <- twfe.att - cv * se.twfe.att
      } else {
        # do weighted bootstrap
        twfe.boot <- unlist(lapply(1:nboot, wboot.twfe.panel,
                                   n = n, y = y, dd = dd, post = post, x = x, i.weights = i.weights))
        # get bootstrap std errors based on IQR
        se.twfe.att <- stats::IQR((twfe.boot - twfe.att)) / (stats::qnorm(0.75) - stats::qnorm(0.25))
        # get symmetric critical values
        cv <- stats::quantile(abs((twfe.boot - twfe.att)/se.twfe.att), probs = 0.95)
        # Estimate of the upper boundary of 95% CI
        uci <- twfe.att + cv * se.twfe.att
        # Estimate of lower boundary of 95% CI
        lci <- twfe.att - cv * se.twfe.att

      }
    }
  }

  #If no covariates, call ordid
  if(is.null(x)){

    # Create dta_long
    dta_long <- as.data.frame(cbind( y = y, post = post, d = dd,
                                     id = rep(1:n,2), w = i.weights))

    reg <- ordid(yname="y",
                 tname = "post",
                 idname = "id",
                 dname = "d",
                 weightsname  = "w",
                 xformla= NULL,
                 data = dta_long,
                 panel = TRUE,
                 boot = boot, boot.type = boot.type, nboot = nboot,
                 inffunc = inffunc)
    twfe.att <- reg$ATT
    se.twfe.att <- reg$se
    uci <- reg$uci
    lci <- reg$lci
    twfe.boot <- reg$boots
    att.inf.func <- reg$att.inf.func
  }


  if(inffunc == FALSE) att.inf.func <- NULL
  return(list(ATT = twfe.att,
              se = se.twfe.att,
              uci = uci,
              lci = lci,
              boots = twfe.boot,
              att.inf.func = att.inf.func))
}
