levenetest <- function(dosage, p0, p1, p2, res) {
  if (length(p0) == 0) {
    return (list(f = NA,
                 numer = NA,
                 denom = NA,
                 df1 = NA,
                 df2 = NA))
  }
  n0 <- sum(p0)
  n1 <- sum(p1)
  n2 <- sum(p2)
  n <- length(res)
  k <- 3
  if (n1 < 5) {
    return (list(f = -1,
                 denom = 0,
                 numer = 0,
                 df1 = 0,
                 df2 = 0))
  }

  if (n2 < 5) {
    k <- 2
    muy0 <- sum(p0 * res) / n0
    muy1 <- sum((p1 + p2) * res) / (n1 + n2)
    z0 <- abs(res - muy0)
    z1 <- abs(res - muy1)
    muz0 <- sum(p0 * z0)
    muz1 <- sum((p1 + p2)*z1)
    muz <- (muz0 + muz1) / n
    muz0 <- muz0 / n0
    muz1 <- muz1 / (n1 + n2)
    denom <- (sum(p0*(z0 - muz0)^2) + sum((p1 + p2)*(z1 - muz1)^2)) / (n - k)
    numer <- (n0*sum((muz0 - muz)^2) + (n1 + n2)*sum((muz1 - muz)^2)) / (k - 1)
    return (list(f = numer / denom,
                 numer = numer,
                 denom = denom,
                 df1 = k - 1,
                 df2 = n - k))
  } else if (n0 < 5) {
    k <- 2
    muy0 <- sum((p0 + p1) * res) / (n0 + n1)
    muy1 <- sum(p2 * res) / n2
    z0 <- abs(res - muy0)
    z1 <- abs(res - muy1)
    muz0 <- sum((p0 + p1)*z0)
    muz1 <- sum(p2 * z1)
    muz <- (muz0 + muz1) / n
    muz0 <- muz0 / (n0 + n1)
    muz1 <- muz1 / n2
    denom <- (sum((p0 + p1)*(z0 - muz0)^2) + sum(p2*(z1 - muz1)^2)) / (n - k)
    numer <- ((n0 + n1)*sum((muz0 - muz)^2) + n2*sum((muz1 - muz)^2)) / (k - 1)
    return (list(f = numer / denom,
                 numer = numer,
                 denom = denom,
                 df1 = k - 1,
                 df2 = n - k))
  }
  
  k <- 3
  muy0 <- sum(p0 * res) / n0
  muy1 <- sum(p1 * res) / n1
  muy2 <- sum(p2 * res) / n2
  
  z0 <- abs(res - muy0)
  z1 <- abs(res - muy1)
  z2 <- abs(res - muy2)
  muz0 <- sum(p0*z0)
  muz1 <- sum(p1*z1)
  muz2 <- sum(p2*z2)
  muz <- (muz0 + muz1 + muz2) / n
  muz0 <- muz0 / n0
  muz1 <- muz1 / n1
  muz2 <- muz2 / n2
  denom <- (sum(p0*(z0 - muz0)^2) + sum(p1*(z1 - muz1)^2) + sum(p2*(z2 - muz2)^2)) / (n - k)
  numer <- (n0*(muz0 - muz)^2 + n1*(muz1 - muz)^2  + n2*(muz2 - muz)^2) / (k - 1)
  return (list(f = numer / denom,
               numer = numer,
               denom = denom,
               df1 = k - 1,
               df2 = n - k))
}