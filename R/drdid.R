NULL
###################################################################################
#'  Locally efficient doubly robust DiD estimators for the ATT
#'
#' @description \code{drdid} is used to compute the locally efficient doubly robust estimators for the ATT
#'  in difference-in-differences (DiD) setups. It can be used with panel or stationary repeated cross section data.
#'  Data should be store in "long" format.
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
#' @param estMethod the method to estimate the nuisance parameters.
#' The default is "imp" which uses weighted least squares to estimate the outcome regressions and
#' inverse probability tilting to the estimate the the propensity score, leading to the improved locally efficient  DR DiD estimator
#' proposed by Sant'Anna and Zhao (2020). The other alternative is "trad",
#' which then uses OLS to estimate outcome regressions and maximum likelihood to estimate propensity score. This leads
#' to the "traditional" locally efficient DR DiD estimator proposed by Sant'Anna and Zhao (2020).
#' @param weightsname The name of the column containing the sampling weights. If NULL, then every observation has the same weights.
#' @param boot Logical argument to whether bootstrap should be used for inference. Default is \code{FALSE} and analytical
#'  standard errors are reported.
#' @param boot.type Type of bootstrap to be performed (not relevant if \code{boot = FALSE}). Options are "weighted" and "multiplier".
#' If \code{boot = TRUE}, default is "weighted".
#' @param nboot Number of bootstrap repetitions (not relevant if boot = \code{FALSE}). Default is 999.
#' @param inffunc Logical argument to whether influence function should be returned. Default is \code{FALSE}.
#'
#' @return A list containing the following components:
#' \item{ATT}{The DR DiD point estimate}
#' \item{se}{ The DR DiD standard error}
#' \item{uci}{Estimate of the upper bound of a 95\% CI for the ATT}
#' \item{lci}{Estimate of the lower bound of a 95\% CI for the ATT}
#' \item{boots}{All Bootstrap draws of the ATT, in case bootstrap was used to conduct inference. Default is NULL}
#' \item{att.inf.func}{Estimate of the influence function. Default is NULL}
#'  \item{ps.flag}{Convergence Flag for the propensity score estimation
#'  (only active if \code{estMethod = "imp"}.): =0 if \code{trust} algorithm converged,
#'    =1 if IPT (original) algorithm converged (in case it was used), =2 if
#'    GLM logit estimator was used (i.e., if both \code{trust} and IPT
#'    did not converged).}
#' \item{call.param}{The matched call.}
#' \item{argu}{Some arguments used in the call (panel, estMethod, boot, boot.type, nboot, type="dr")}
#'
#' @references
#' \cite{Graham, Bryan, Pinto, Cristine, and Egel, Daniel (2012),
#' "Inverse Probability Tilting for Moment Condition Models with Missing Data."
#'  Review of Economic Studies, vol. 79 (3), pp. 1053-1079, \doi{10.1093/restud/rdr047}}
#'
#'
#' \cite{Sant'Anna, Pedro H. C. and Zhao, Jun. (2020),
#' "Doubly Robust Difference-in-Differences Estimators." Journal of Econometrics, Vol. 219 (1), pp. 101-122,
#'  \doi{10.1016/j.jeconom.2020.06.003}}
#'
#'
#'
#' @details
#'
#' When panel data are available (\code{panel = TRUE}), the \code{drdid} function implements the
#' locally efficient doubly robust difference-in-differences (DiD) estimator for the average treatment effect
#' on the treated (ATT) defined in equation (3.1) in Sant'Anna and Zhao (2020). This estimator makes use of
#' a logistic propensity score model for the probability of being in the treated group,
#' and of a linear regression model for the outcome evolution among the comparison units.
#'
#'
#' When only stationary repeated cross-section data are available (\code{panel = FALSE}), the \code{drdid} function
#' implements the locally efficient doubly robust difference-in-differences (DiD) estimator for the
#' average treatment effect on the treated (ATT) defined in equation (3.4) in Sant'Anna and Zhao (2020).
#' This estimator makes use of a logistic propensity score model for the probability of being in the
#' treated group, and of (separate) linear regression models for the outcome of both treated and comparison units,
#' in both pre and post-treatment periods.
#'
#'
#' When one sets \code{estMethod = "imp"} (the default), the nuisance parameters (propensity score and
#' outcome regression parameters) are estimated using the methods described in Sections 3.1 and 3.2 of
#' Sant'Anna and Zhao (2020). In short, the propensity score parameters are estimated using the inverse
#' probability tilting estimator proposed by Graham, Pinto and Pinto (2012), and the outcome
#' regression coefficients are estimated using weighted least squares,where the weights depend on
#' the propensity score estimates; see Sant'Anna and Zhao (2020) for details.
#'
#' When one sets \code{estMethod = "trad"}, the propensity score parameters are estimated using maximum
#' likelihood, and the outcome regression coefficients are estimated using ordinary least squares.
#'
#' The main advantage of using \code{estMethod = "imp"} is that the resulting estimator is not only
#' locally efficient and doubly robust for the ATT, but it is also doubly robust for inference;
#' see Sant'Anna and Zhao (2020) for details.
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
#' # -----------------------------------------------
#' # Implement improved DR locally efficient DiD with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE)
#'
#' #Implement "traditional" DR locally efficient DiD with panel data
#' drdid(yname="re", tname = "year", idname = "id", dname = "experimental",
#'       xformla= ~ age+ educ+ black+ married+ nodegree+ hisp+ re74,
#'       data = eval_lalonde_cps, panel = TRUE, estMethod = "trad")
#'
#' # -----------------------------------------------
#' # Repeated cross section case
#' # -----------------------------------------------
#' # use the simulated data provided in the package
#' #Implement "improved" DR locally efficient DiD with repeated cross-section data
#' drdid(yname="y", tname = "post", idname = "id", dname = "d",
#'       xformla= ~ x1 + x2 + x3 + x4,
#'       data = sim_rc, panel = FALSE, estMethod = "imp")
#'
#' #Implement "traditional" DR locally efficient DiD with repeated cross-section data
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
  if (dp$panel == TRUE) {
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
  if (dp$panel == FALSE) {
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
