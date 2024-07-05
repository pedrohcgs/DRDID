###################################################################################
# Compute Weighted OLS regression parameters for the Improved Doubly-Robust DID estimator with Repeated Cross Section Data
###################################################################################
# Pre = T stands for pre-treatment period
# treat = F  stands for control group

wols_rc <- function(y, post, D, int.cov, pscore, i.weights, pre = NULL, treat = FALSE){
  #-----------------------------------------------------------------------------
  # Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  or.weights <- as.vector(i.weights * pscore/(1 - pscore))

  if((pre == T) & (treat==T)) {
    subs <- (D==1)*(post==0)
  }
  if((pre == F) & (treat==T)) {
    subs <- (D==1)*(post==1)
  }
  if((pre == T) & (treat==F)) {
    subs <- (D==0)*(post==0)
  }
  if((pre == F) & (treat==F)) {
    subs <- (D==0)*(post==1)
  }
  #Run weighted OLS
  beta.wls <- stats::coef(stats::lm(y ~ -1 + int.cov,
                                    subset = subs==1,
                                    weights = or.weights))
  if(anyNA(beta.wls)){
    stop("Outcome regression model coefficients have NA components. \n Multicollinearity (or lack of variation) of covariates is a likely reason")
  }

  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.wls, int.cov))

  # return fitted values
  return(list(out.reg = out.delta))

}
