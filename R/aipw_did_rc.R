###################################################################################
# Locally Efficient DR DiD estimator with Repeated Cross Section Data
###################################################################################

aipw_did_rc <- function(y, post, D, ps,
                        out.y.treat.post, out.y.treat.pre,
                        out.y.cont.post, out.y.cont.pre,
                        i.weights, trim.ps){
  #-----------------------------------------------------------------------------
  # First, the weights
  i.weights <- trim.ps * i.weights/mean(i.weights)
  w.treat.pre <- trim.ps * i.weights * D * (1 - post)
  w.treat.post <- trim.ps * i.weights * D * post
  w.cont.pre <- trim.ps * i.weights * ps * (1 - D) * (1 - post)/(1 - ps)
  w.cont.post <- trim.ps * i.weights * ps * (1 - D) * post/(1 - ps)

  #Extra weights for efficiency
  w.d <- trim.ps * i.weights * D
  w.dt1 <- trim.ps * i.weights * D * post
  w.dt0 <- trim.ps * i.weights * D * (1 - post)

  # Estimator of each component
  att.treat.pre <- mean(w.treat.pre * (y - out.y.cont.pre))/ mean(w.treat.pre)
  att.treat.post <- mean(w.treat.post * (y - out.y.cont.post))/ mean(w.treat.post)
  att.cont.pre <- mean(w.cont.pre * (y - out.y.cont.pre))/ mean(w.cont.pre)
  att.cont.post <- mean(w.cont.post * (y - out.y.cont.post))/ mean(w.cont.post)

  att.d.post <- mean(w.d * (out.y.treat.post - out.y.cont.post))/mean(w.d)
  att.dt1.post <- mean(w.dt1 * (out.y.treat.post - out.y.cont.post))/mean(w.dt1)
  att.d.pre <- mean( w.d * (out.y.treat.pre - out.y.cont.pre))/mean(w.d)
  att.dt0.pre <- mean(w.dt0 * (out.y.treat.pre - out.y.cont.pre))/mean(w.dt0)


  # ATT estimator
  aipw.att <- (att.treat.post - att.treat.pre) - (att.cont.post - att.cont.pre) +
    (att.d.post - att.dt1.post) - (att.d.pre - att.dt0.pre)

  return(aipw.att)
}
