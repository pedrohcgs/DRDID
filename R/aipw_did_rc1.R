###################################################################################
# DR DiD estimator with Repeated Cross Section Data

aipw_did_rc1 <- function(y, post, D, ps, out.reg, i.weights){
  #-----------------------------------------------------------------------------
  # Compute the AIPW estimator
  # Compute  standardized IPW estimator
  i.weights <- i.weights/mean(i.weights)

  w.treat.pre <- i.weights * D * (1 - post)
  w.treat.post <- i.weights * D * post
  w.cont.pre <- i.weights * ps * (1 - D) * (1 - post)/ (1 - ps)
  w.cont.post  <- i.weights * ps * (1 - D) * post/ (1 - ps)

  aipw.1.pre <- mean(w.treat.pre * (y - out.reg)) / mean(w.treat.pre)
  aipw.1.post <- mean(w.treat.post * (y - out.reg)) / mean(w.treat.post)
  aipw.0.pre <- mean(w.cont.pre * (y - out.reg)) / mean(w.cont.pre)
  aipw.0.post <- mean(w.cont.post * (y - out.reg)) / mean(w.cont.post)

  aipw.att <- (aipw.1.post - aipw.1.pre) - (aipw.0.post - aipw.0.pre)

  return(aipw.att)
}
