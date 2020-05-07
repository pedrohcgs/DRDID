NULL
###################################################################################
#'  	Locally Efficient Doubly Robust Difference-in-Differences Estimators for the ATT
#'
#' @description \code{drdid} computes the doubly robust estimators for the average treatment effect on the treated
#'  in DID setups. It can be used with panel or repeated cross section data. The nuisance parameters can be estimated
#'  using "traditional" methods (OLS for the outcome regressions and maximum likelihood for the propensity score), or
#'  "improved" methods (weighted least squares for the outcome regressions and
#'  inverse probability tilting for the propensity score). Estimators based on "improved methods" are also doubly
#'  robust for inference. See Sant'Anna and Zhao (2020) for a detailed description.
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
#'  Whens \code{panel=TRUE}, the variable \code{idname} must be set.  When \code{panel=FALSE}, the data is treated
#'  as stationary repeated cross sections.
#' @param estMethod the method to estimate the nuisance parameters.
#' The default is "imp" which uses weighted least squares to estimate the outcome regressions and
#' inverse probability tilting to the estimate the the propensity score, leading to the improved locally efficient  DR DID estimator
#' proposed by Sant'Anna and Zhao (2020). The other alternative is "trad",
#' which then uses OLS to estimate outcome regressions and maximum likelihood to estimate propensity score. This leads
#' to the "traditional" locally efficient DR DID estimator proposed by Sant'Anna and Zhao (2020).
#' @param weightsname The name of the column containing the sampling weights. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is \code{FALSE} and analytical
#'  standard errors are reported.
#' @param boot.type Type of bootstrap to be performed (not relevant if boot = FALSE). Options are "weighted" and "multiplier".
#' If \code{boot==TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = \code{FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is \code{FALSE}.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DID point estimate}
#' \item{se}{ The DR DID standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation
#'  (only active if \code{estMethod = "imp"}.): =0 if \code{trust} algorithm converged,
#'    =1 if IPT (original) algorithm converged (in case it was used), =2 if GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#' \item{call.param}{The matched call.}
#' \item{argu}{Some arguments used in the call (panel, estMethod, boot, boot.type, nboot, type="dr")}

#' @references{
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Forthcoming,
#' \url{https://arxiv.org/abs/1812.01723}}
#' }
#' @examples
#' # -----------------------------------------------
#' # -----------------------------------------------
#' # Panel data case
#' # -----------------------------------------------
#' # Data preparation using Lalonde sample with CPS comparison group
#' # Create "selection" treatment dummy: 1 if in experimental sample, 0 if in non-experimental
#' nsw_long$treated2 <- ifelse(is.na(nsw_long$treated), 0 , 1)
#'
#' # Form the Lalonde sample with CPS comparison group
#' eval_lalonde_cps <- subset(nsw_long, nsw_long$treated == 0 | nsw_long$sample == 2)
#' # -----------------------------------------------
#' # Implement Further improved DR locally efficient DID with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "treated2",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE)
#'
#' #Implement "traditional" DR locally efficient DID with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "treated2",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE, estMethod = "trad")
#' # -----------------------------------------------
#' # -----------------------------------------------
#' # Repeated cross section case
#' # -----------------------------------------------
#' # use the simulated data
#' #Implement "improved" DR locally efficient DID with repeated cross-section data
#' drdid(yname="y", tname = "post", idname = "id", dname = "d",
#'       xformla= ~ x1 + x2 + x3 + x4,
#'       data = sim_rc, panel = FALSE, estMethod = "imp")
#'
#' #Implement "traditional" DR locally efficient DID with repeated cross-section data
#' drdid(yname="y", tname = "post", idname = "id", dname = "d",
#'       xformla= ~ x1 + x2 + x3 + x4,
#'       data = sim_rc, panel = FALSE, estMethod = "trad")
#'

#' @export

drdid <- function(yname, tname, idname, dname, xformla = NULL, data,
                  panel = TRUE, estMethod = c("imp", "trad"), weightsname = NULL,
                  boot = FALSE, boot.type =  c("weighted", "multiplier"),
                  nboot = 999, inffunc = FALSE) {
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
    estMethod = estMethod,
    weightsname = weightsname,
    boot = boot,
    boot.type = boot.type,
    nboot = nboot,
    inffunc = inffunc
  )
  #-----------------------------------------------------------------------------
  # Implement the methods
  # First panel data
  if (dp$panel == T) {
    if (dp$estMethod == "imp") {
      att_est <- drdid_imp_panel(
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
      ps.flag <- att_est$ps.flag

    } else {
      att_est <- drdid_panel(
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
      ps.flag <- NULL

    }

  }
  # Now repeated cross section
  if (dp$panel == F) {
    if (dp$estMethod == "imp") {
      att_est <- drdid_imp_rc(
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
      ps.flag <- att_est$ps.flag

    } else {
      att_est <- drdid_rc(
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
      ps.flag <- NULL

    }

  }
  #---------------------------------------------------------------------
  # record the call
  call.param <- match.call()
  # Record all arguments used in the function
  argu <- mget(names(formals()), sys.frame(sys.nframe()))
  argu <- list(
    panel = argu$panel,
    estMethod = argu$estMethod,
    boot = argu$boot,
    boot.type = argu$boot.type,
    nboot = argu$nboot,
    type = "dr"
  )

  # Return these variables
  ret <- list(
    ATT = att_est$ATT,
    se = att_est$se,
    lci = att_est$lci,
    uci = att_est$uci,
    boots = att_est$boots,
    att.inf.func = att_est$att.inf.func,
    ps.flag = ps.flag,
    call.param = call.param,
    argu = argu
  )

  # Define a new class
  class(ret) <- "drdid"
  # return the list
  return(ret)


}
