# Bootrstapped TWFE Panel data
# 2 periods and 2 groups
mboot.twfep.did = function(n, linrep, nboot){
  # Use the Mammen(1993) binary V's
  k1 = 0.5 * (1 - 5^0.5)
  k2 = 0.5 * (1 + 5^0.5)
  pkappa = 0.5 * (1 + 5^0.5)/(5^0.5)


  bootapply <- function(nn, n=n, linrep=linrep) {
    v <- stats::rbinom(n, 1, pkappa)
    v <- ifelse(v == 1, k1, k2)
    v <- as.vector(c(v, v))
    b.did <- mean(linrep * v)
    return(b.did)
  }

  boot.did <- unlist(lapply(1:nboot, bootapply, n=n, linrep=linrep))


}
