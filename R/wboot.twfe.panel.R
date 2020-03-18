# Bootrstapped TWFE Difference-in-Differences with panel data
# 2 periods and 2 groups

wboot.twfe.panel <- function(nn, n, y, dd, post, x, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  v <- as.matrix(c(v, v))
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the TWFE Regression
  reg.b <- stats::lm(y ~  dd:post + post + dd + x, weights = b.weights)
  twfe.att.b <- reg.b$coefficients["dd:post"]
  #-----------------------------------------------------------------------------
  return(twfe.att.b)
}
