#Loss function for estimation of the Bias reduced PS, based on Graham, Pinton and Egel (2012, 2016)

loss.ps.IPT <- function(gamma1, n, D, int.cov, iw){
  #Coefficients for quadratic extrapolation
  cn <- -(n - 1)
  bn <- -n + (n - 1) * log(n - 1)
  an <- -(n - 1) * (1 - log(n - 1) + 0.5 * (log(n - 1))^2)
  vstar <- log(n - 1)

  v <- gamma1 %*% t(int.cov)
  phi <- ifelse(v < vstar, - v - exp(v), an + bn * v + 0.5 * cn * (v^2))
  phi1 <- ifelse(v < vstar, - 1 - exp(v), bn + cn * v)
  phi2 <- ifelse(v < vstar, - exp(v), cn)

  #phi <- (v<vstar) * (- v - exp(v)) + (v>=vstar) * (an + bn * v + 0.5 * cn * (v^2))
  #phi1 <- (v<vstar) * (- 1 - exp(v)) + (v>=vstar) * (bn  + cn * v)
  #phi2 <- (v<vstar) * (- exp(v)) + (v>=vstar) * cn

  # Minus is because nlm minimizes functions, and we aim to maximize!
  res <- - sum(iw * (1 - D) * phi + v)

  attr(res, "gradient") <- - t(int.cov) %*% as.vector(iw * ((1-D) * phi1 + 1))
  attr(res, "hessian")  <-  - t(as.vector((1-D) * iw * phi2) * int.cov) %*% int.cov
  return(res)
}
