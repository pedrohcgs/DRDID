#' @title Print
#'
#' @description Prints a drdid Object
#'
#' @param x A drdid object
#' @param ... Other params (required as generic function, but not used)
#' @importFrom utils write.table
#' @export
#' @noRd
# Define new print function
print.drdid <- function(x, ...) {
  #-----------------------------------------------------------------------------
  # Preliminaries
  # Estimation methods
  if(x$argu$type == "dr"){
    if (x$argu$estMethod[1] == "imp") {
      estMeth1 <- "Further improved locally efficient DR DID estimator for the ATT:\n"
      estMeth2 <- "Outcome regression est. method: weighted least squares."
      estMeth3 <- "Propensity score est. method: inverse prob. tilting."
    }
    if (x$argu$estMethod[1] == "trad") {
      estMeth1 <- "Locally efficient DR DID estimator for the ATT:\n"
      estMeth2 <- "Outcome regression est. method: OLS."
      estMeth3 <- "Propensity score est. method: maximum likelihood."
    }
    if (x$argu$estMethod[1] == "imp2") {
      estMeth1 <- "Further improved DR (but not locally efficient) DID estimator for the ATT:\n"
      estMeth2 <- "Outcome regression est. method: weighted least squares."
      estMeth3 <- "Propensity score est. method: inverse prob. tilting."
    }
    if (x$argu$estMethod[1] == "trad2") {
      estMeth1 <- "DR (but not locally efficient) DID estimator for the ATT:\n"
      estMeth2 <- "Outcome regression est. method: OLS."
      estMeth3 <- "Propensity score est. method: maximum likelihood."
    }
  }
  if(x$argu$type == "ipw"){
    if (x$argu$normalized == T) {
      estMeth1 <- "IPW DID estimator for the ATT:\n"
      estMeth2 <- "Hajek-type IPW estimator (weights sum up to 1)."
      estMeth3 <- "Propensity score est. method: maximum likelihood."
    }
    if (x$argu$normalized == F) {
      estMeth1 <- "IPW DID estimator for the ATT:\n"
      estMeth2 <- "Horvitz-Thompson-type IPW estimator."
      estMeth3 <- "Propensity score est. method: maximum likelihood."
    }
  }

  if(x$argu$type == "or"){
    estMeth1 <- "Outcome-Regression DID estimator for the ATT:\n"
    estMeth2 <- ""
    estMeth3 <- "Outcome regression est. method: OLS."
  }

  # Panel vs. Repeated cross sections
  if (x$argu$panel == T)
    panel <- "Estimator based on panel data."
  if (x$argu$panel == F)
    panel <-
      "Estimator based on (stationary) repeated cross-sections data."

  #Creat parameters for the Table

  header <-
    c("ATT",
      "Std. Error",
      "t value",
      "Pr(>|t|)",
      "[95% Conf.",
      "Interval]")
  body <- cbind(
    round(x$ATT, digits = 4),
    round(x$se, digits = 4),
    round(x$ATT / x$se, digits = 4),
    round(2 * (stats::pnorm(-abs(
      x$ATT / x$se
    ))), digits = 4),
    round(x$lci, digits = 4),
    round(x$uci, digits = 4)
  )
  colnames(body) <- header
  #-----------------------------------------------------------------------------
  #Output
  cat(" Call:\n")
  print(x$call)
  cat("------------------------------------------------------------------")
  cat("\n", estMeth1, "\n")
  utils::write.table(format(rbind(header, body), justify= "centre", digits=2, nsmall=2),
                     row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  cat("------------------------------------------------------------------")
  # Panel data?
  cat("\n", panel)
  #Estimation Method
  if(x$argu$type != "or") cat("\n", estMeth2)
  cat("\n", estMeth3)
  # Analytical vs bootstrapped standard errors
  if (x$argu$boot == T) {
    boot1 <-
      cat(
        "\n Boostrapped standard error based on",
        x$argu$nboot,
        "bootstrap draws. \n Bootstrap method:",
        x$argu$boot.type[1],
        ". \n"
      )
  }
  if (x$argu$boot == F) {
    boot1 <- cat("\n Analytical standard error.\n")
  }
  cat("------------------------------------------------------------------")
  cat("\n See Sant'Anna and Zhao (2020) for details.")
}
