###################################################################################
# Compute Propensity Score using IPT
#
# Propensity Score estimator based on Inverse Probability of Tilting.
# flag is to understand convergence. =0 if trust algorithm converge, =1 if IPT algorithm converged, if needed,
#       = 2 if GLM logit estimator was used (both IPT and trust did not converge)}


pscore.cal <- function(D, int.cov, i.weights, n){
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Initial conditions for pscore
  init.gamma <- suppressWarnings(stats::coef(stats::glm(D ~ -1 + int.cov, family = "binomial"),
                                             weights = i.weights))

  #Compute IPT pscore
  pscore.cal <- suppressWarnings(trust::trust(loss.ps.cal, parinit = init.gamma, rinit=1,
                                              rmax=1000, iterlim=1000,
                                              D = D, int.cov = int.cov, iw = i.weights))

  flag <- ifelse(pscore.cal$converged, 0, 1)

  gamma.cal <- try(pscore.cal$argument)

  # If Algorithm doesn't converge, update like in Graham et al
  if(pscore.cal$converged==F) {

    pscore.IPT <- suppressWarnings(stats::nlm(loss.ps.IPT, init.gamma, iterlim = 10000, gradtol = 1e-06,
                                              check.analyticals = F,
                                              D = D, int.cov = int.cov, iw = i.weights,
                                              n = n))
    gamma.cal <- try(pscore.IPT$estimate)

    if(pscore.IPT$code>3) {
      gamma.cal <- init.gamma
      flag <- 2
    }
  }

  #Compute fitted pscore and weights for regression
  pscore.index <- tcrossprod(gamma.cal, int.cov)
  pscore <- as.numeric(stats::plogis(pscore.index))

  # return pscore and flag
  return(list(pscore = pscore,
              flag = flag))

}
