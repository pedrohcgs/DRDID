
###################################################################################
# Weighted OLS for Intrinsically Efficient DR-DID


wols.ie.panel <- function(deltaY, D, int.cov, pscore, i.weights){
  #-----------------------------------------------------------------------------
  #Do not divide by zero
  #pscore <- pmin(pscore, 1 - 1e-16)
  i.weights <- as.vector( i.weights * (pscore/(1 - pscore))^2)
  #Run weighted OLS
  beta.ie. <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
                                    subset = D==0,
                                    weights = i.weights))

  #get fitted values
  out.delta <-  as.numeric(tcrossprod(beta.ie., int.cov))

  # return pscore and flag
  return(list(out.reg = out.delta))

}
