###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Panel Data


wols.br.panel <- function(deltaY, D, int.cov, pscore, i.weights){
  #-----------------------------------------------------------------------------

  i.weights <- as.vector(i.weights * pscore/(1 - pscore))

  #Run weighted OLS
  # beta.cal <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
  #                           subset = D==0,
  #                           weights = i.weights))
  control_filter <- (D == 0)
  beta.cal <- stats::coef(fastglm::fastglm(
                            x = int.cov[control_filter, , drop = FALSE],
                            y = deltaY[control_filter],
                            weights = i.weights[control_filter],
                            family = gaussian(link = "identity")
  ))

  if(anyNA(beta.cal)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }

  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.cal, int.cov))

  # return fitted values
  return(list(out.reg = out.delta))

}
