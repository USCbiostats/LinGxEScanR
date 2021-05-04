levinetest <- function(res, p0, p1, p2) {
  nk <- numeric(3)
  tau <- numeric(3)
  lz <- matrix(0, length(res), 3)
  mlz <- numeric(3)
  d <- numeric(3)
  nk[1] <- sum(p0)
  nk[2] <- sum(p1)
  nk[3] <- sum(p2)
  tau[1] <- sum(p0*res) / nk[1]
  tau[2] <- sum(p1*res) / nk[2]
  tau[3] <- sum(p2*res) / nk[3]
  lz[,1] <- abs(res - tau[1])
  lz[,2] <- abs(res - tau[2])
  lz[,3] <- abs(res - tau[3])
  mlz[1] <- sum(p0*lz[,1]) / nk[1]
  mlz[2] <- sum(p1*lz[,2]) / nk[2]
  mlz[3] <- sum(p2*lz[,3]) / nk[3]
  d[1] <- sum(p0*(lz[,1] - mlz[1])^2)
  d[2] <- sum(p1*(lz[,2] - mlz[2])^2)
  d[3] <- sum(p2*(lz[,3] - mlz[3])^2)
  denom <- sum(d) / (length(res) - 2)
  mz <- (sum(p0*lz[,1]) + sum(p1*lz[,2]) + sum(p2*lz[,3]))/(length(res) - 2)
  numer <- ((mlz[1] - mz)^2 + (mlz[2] - mz)^2 + (mlz[3] - mz)^2) / 2
  f <- numer/denom
  return(list(nk = nk,
              tau = tau,
              lz = lz,
              mlz = mlz,
              d = d,
              denom = denom,
              mz = mz,
              numer = numer,
              f = f))
}