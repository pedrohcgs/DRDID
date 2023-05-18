NULL
###################################################################################
#' Outcome regression DiD estimators for the ATT
#'
#' @description \code{ordid} computes the outcome regressions estimators for the average treatment effect on the
#' treated in difference-in-differences (DiD) setups. It can be used with panel or repeated cross section data.
#' See Sant'Anna and Zhao (2020) for details.
#'
#' @param yname The name of the outcome variable.
#' @param tname The name of the column containing the time periods.
#' @param idname The name of the column containing the unit id name.
#' @param dname The name of the column containing the treatment group (=1 if observation is treated in the post-treatment, =0 otherwise)
#' @param xformla A formula for the covariates to include in the model. It should be of the form \code{~ X1 + X2}.
#' (intercept should not be listed as it is always automatically included). Default is NULL which is equivalent to \code{xformla=~1}.
#' @param data The name of the data.frame that contains the data.
#' @param panel Whether or not the data is a panel dataset. The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a particular point in time.  The default is TRUE.
#'  When \code{panel = FALSE}, the data is treated
#'  as stationary repeated cross sections.
#' @param weightsname The name of the column containing the sampling weights. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is \code{FALSE} and analytical
#'  standard errors are reported.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = \code{FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is \code{FALSE}.
#'
#' @return A list containing the following components:
#' \item{ATT}{The OR DiD point estimate}
#' \item{se}{ The OR DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#' \item{call.param}{The matched call.}
#' \item{argu}{Some arguments used in the call (panel, normalized, boot, boot.type, nboot, type=="or")}
#'
#' @details
#'
#' The \code{ordid} function implements
#' outcome regression difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) defined in equation (2.2) of Sant'Anna and Zhao (2020). The estimator follows the same spirit
#' of the nonparametric estimators proposed by Heckman, Ichimura and Todd (1997), though here the the outcome regression
#' models are assumed to be linear in covariates (parametric).
#'
#' The nuisance parameters (outcome regression coefficients) are estimated via ordinary least squares.
#'
#'
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
#'
#' @examples
#' # -----------------------------------------------
#' # Panel data case
#' # -----------------------------------------------
#' # Form the Lalonde sample with CPS comparison group
#' eval_lalonde_cps <- subset(nsw_long, nsw_long$treated == 0 | nsw_long$sample == 2)
#' # Further reduce sample to speed example
#' set.seed(123)
#' unit_random <- sample(unique(eval_lalonde_cps$id), 5000)
#' eval_lalonde_cps <- eval_lalonde_cps[eval_lalonde_cps$id %in% unit_random,]
#' # Implement OR DiD with panel data
#' ordid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE)
#'
#' # -----------------------------------------------
#' # Repeated cross section case
#' # -----------------------------------------------
#' # use the simulated data provided in the package
#' # Implement OR DiD with repeated cross-section data
#' # use Bootstrap to make inference with 199 bootstrap draws (just for illustration)
#' ordid(yname="y", tname = "post", idname = "id", dname = "d",
#'       xformla= ~ x1 + x2 + x3 + x4,
#'       data = sim_rc, panel = FALSE,
#'       boot = TRUE, nboot = 199)
#'
#'
#'
#' @export

ordid <- function(yname, tname, idname, dname, xformla = NULL, data,
                  panel = TRUE, weightsname = NULL, boot = FALSE,
                  boot.type =  c("weighted", "multiplier"), nboot = 999,
                  inffunc = FALSE) {
  #-----------------------------------------------------------------------------
  # Pre-process data
  dp <- pre_process_drdid(
    yname = yname,
    tname = tname,
    idname = idname,
    dname = dname,
    xformla = xformla,
    data = data,
    panel = panel,
    weightsname = weightsname,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    inffunc = inffunc
  )
  #-----------------------------------------------------------------------------
  # Implement the methods
  # First panel data
  if (dp$panel == TRUE) {

    att_est <- reg_did_panel(
      y1 = dp$y1,
      y0 = dp$y0,
      D = dp$D,
      covariates = dp$covariates,
      i.weights = dp$i.weights,
      boot = dp$boot,
      boot.type = dp$boot.type,
      nboot = dp$nboot,
      inffunc = dp$inffunc
    )
  }

  # Now repeated cross section
  if (dp$panel == FALSE) {
    att_est <- reg_did_rc(
      y = dp$y,
      post = dp$post,
      D = dp$D,
      covariates = dp$covariates,
      i.weights = dp$i.weights,
      boot = dp$boot,
      boot.type = dp$boot.type,
      nboot = dp$nboot,
      inffunc = dp$inffunc
    )

  }

  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  argu <- list(
    panel = argu$panel,
    normalized = argu$normalized,
    boot = argu$boot,
    boot.type = argu$boot.type,
    nboot = argu$nboot,
    type = "or"
  )

  # Return these variables
  ret <- list(
    ATT = att_est$ATT,
    se = att_est$se,
    lci = att_est$lci,
    uci = att_est$uci,
    boots = att_est$boots,
    att.inf.func = att_est$att.inf.func,
    call.param = call.param,
    argu = argu
  )

  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)


}
