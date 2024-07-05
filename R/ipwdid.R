NULL
###################################################################################
#' Inverse probability weighted DiD estimators for the ATT
#'
#' @description \code{ipwdid} computes the inverse probability weighted estimators for the average treatment effect
#'  on the treated in difference-in-differences (DiD) setups. It can be used with panel or stationary repeated
#'  cross-sectional data, with or without normalized (stabilized) weights. See Abadie (2005) and Sant'Anna and
#'  Zhao (2020) for details.
#'
#' @param yname The name of the outcome variable.
#' @param tname The name of the column containing the time periods.
#' @param idname The name of the column containing the unit id name.
#' @param dname The name of the column containing the treatment group (=1 if observation is treated in the post-treatment, =0 otherwise)
#' @param xformla A formula for the covariates to include in the model. It should be of the form \code{~ X1 + X2}
#' (intercept should not be listed as it is always automatically included). Default is NULL which is equivalent to \code{xformla=~1}.
#' @param data The name of the data.frame that contains the data.
#' @param panel Whether or not the data is a panel dataset. The panel dataset should be provided in long format -- that
#'  is, where each row corresponds to a unit observed at a particular point in time.  The default is TRUE.
#'  When \code{panel = FALSE}, the data is treated
#'  as stationary repeated cross sections.
#' @param normalized Logical argument to whether IPW weights should be normalized to sum up to one. Default is \code{TRUE}.
#' @param weightsname The name of the column containing the sampling weights. If NULL, then every observation has the same weights. The weights are normalized and therefore enforced to have mean 1 across all observations.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is \code{FALSE} and analytical
#'  standard errors are reported.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = \code{FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is \code{FALSE}.
#'
#' @return A list containing the following components:
#' \item{ATT}{The IPW DiD point estimate}
#' \item{se}{ The IPW DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#' \item{call.param}{The matched call.}
#' \item{argu}{Some arguments used in the call (panel, normalized, boot, boot.type, nboot, type=="ipw")}

#' @references
#' \cite{Abadie, Alberto (2005), "Semiparametric Difference-in-Differences Estimators",
#' Review of Economic Studies, vol. 72(1), p. 1-19, \doi{10.1111/0034-6527.00321}}
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
#'
#' # Implement IPW DiD with panel data (normalized weights)
#' ipwdid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE)
#'
#' # -----------------------------------------------
#' # Repeated cross section case
#' # -----------------------------------------------
#' # use the simulated data provided in the package
#' #Implement IPW DiD with repeated cross-section data (normalized weights)
#' # use Bootstrap to make inference with 199 bootstrap draws (just for illustration)
#' ipwdid(yname="y", tname = "post", idname = "id", dname = "d",
#'       xformla= ~ x1 + x2 + x3 + x4,
#'       data = sim_rc, panel = FALSE,
#'       boot = TRUE, nboot = 199)
#'
#' @details
#'
#' The \code{ipwdid} function implements the
#' inverse probability weighted (IPW) difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) proposed by Abadie (2005) (\code{normalized = FALSE}) or Hajek-type version
#' defined in equations (4.1) and (4.2) in Sant'Anna and Zhao (2020), when either panel data or
#' stationary repeated cross-sectional data are available. This estimator makes use of
#' a logistic propensity score model for the probability of being in the treated group, and the propensity score
#' parameters are estimated via maximum likelihood.
#'
#'
#' @export

ipwdid <- function(yname, tname, idname, dname, xformla = NULL,
                   data, panel = TRUE, normalized = TRUE,
                   weightsname = NULL, boot = FALSE,
                   boot.type =  c("weighted", "multiplier"),
                   nboot = 999,
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
    normalized = normalized,
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
    if (dp$normalized) {
      att_est <- std_ipw_did_panel(
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

    } else {
      att_est <- ipw_did_panel(
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

  }
  # Now repeated cross section
  if (dp$panel == FALSE) {
    if (dp$normalized) {
      att_est <- std_ipw_did_rc(
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

    } else {
      att_est <- ipw_did_rc(
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
    type = "ipw"
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
