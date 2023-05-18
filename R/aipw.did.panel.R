###################################################################################
# DR DiD estimator for the ATT with panel Data


aipw.did.panel <- function(deltaY, D, ps, out.reg, i.weights){
  #-----------------------------------------------------------------------------
  # Compute the AIPW estimator
  w.treat <- i.weights * D
  w.cont <- i.weights * (1 - D) * ps / (1 - ps)

  aipw.1 <- mean(w.treat * (deltaY - out.reg)) / mean(w.treat)
  aipw.0 <- mean(w.cont * (deltaY - out.reg)) / mean(w.cont)

  aipw.att <- aipw.1 - aipw.0

  return(aipw.att)
}
