# Determine which statistics will be output for a model
# Also verifies if inputs are valid
# Returns a logical vector

assignstats <- function(statlist, modelname) {
  if (is.character(statlist) == FALSE)
    stop(paste(modelname, "must be a character vector"))
  if (length(statlist) == 0)
    return (rep(FALSE, 4))
  if (length(statlist) == 1) {
    if(statlist == "all")
      return (rep(TRUE, 4))
    if (statlist == "" | statlist == "none")
      return (rep(FALSE, 4))
  }
  statsidx <- match(statlist, c("lrt", "score", "Wald", "WaldHW"))
  if (all(is.na(statsidx) == FALSE) == FALSE)
    stop(paste("Unknown values in", modelname))
  teststats <- rep(FALSE, 4)
  teststats[statsidx] <- TRUE
  return(teststats)
}

assigncolumnnames <- function(outfile, binarye, teststats, levene) {
  snpcols <- c("SNP", "CHR", "LOC", "REF", "ALT")
  nfcols <- c("n", "aaf")
  if (binarye == TRUE)
    nfcols <- c(nfcols, "n_e0", "aaf_e0", "n_e1", "aaf_e1")
  statcols <- character(0)
  modelname <- c("go", "ge", "gxe")
  statnameg <- vector("list", 4)
  statnameg[[1]] <- "lrt_chi_g"
  statnameg[[2]] <- "score_z_g"
  statnameg[[3]] <- c("wald_se_g", "wald_t_g", "wald_df_g")
  statnameg[[4]] <- c("waldhw_se_g", "waldhw_z_g")
  statnamegxe <- vector("list", 4)
  statnamegxe[[1]] <- "lrt_chi_gxe"
  statnamegxe[[2]] <- "score_z_gxe"
  statnamegxe[[3]] <- c("wald_se_gxe", "wald_t_gxe", "wald_df_gxe")
  statnamegxe[[4]] <- c("waldhw_se_gxe", "waldhw_z_gxe")
  statname2df <- vector("list", 4)
  statname2df[[1]] <- "lrt_chi2df"
  statname2df[[2]] <- "score_chi2df"
  statname2df[[3]] <- c("wald_chi2df", "wald_var_g", "wald_var_gxe", "wald_cov_ggxe")
  statname2df[[4]] <- c("waldhw_chi2df", "waldhw_var_g", "waldhw_var_gxe", "waldhw_cov_ggxe")
  for (i in 1:length(teststats)) {
    if (all(teststats[[i]] == FALSE) == TRUE)
      next
    statcols <- c(statcols, paste("bg", modelname[i], sep = '_'))
    for (j in 1:length(teststats[[i]])) {
      if (teststats[[i]][j] == TRUE)
        statcols <- c(statcols, paste(statnameg[[j]], modelname[i], sep = '_'))
    }
  }
  if (all(teststats[[3]] == FALSE) == FALSE) {
    statcols <- c(statcols, "bgxe")
    for (j in 1:length(teststats[[3]])) {
      if (teststats[[3]][j] == TRUE)
        statcols <- c(statcols, statnamegxe[[j]])
    }
    for (j in 1:length(teststats[[3]])) {
      if (teststats[[3]][j] == TRUE)
        statcols <- c(statcols, statname2df[[j]])
    }
  }
  levenecols <- character(0)
  if (levene[1] == TRUE)
    levenecols <- c(levenecols, "flevenee", "fchisqnume", "fchisqdene",
                    "dfnume", "dfdene")
  if (levene[2] == TRUE)
    levenecols <- c(levenecols, "flevenene", "fchisqnumne", "fchisqdenne",
                    "dfnumne", "dfdenne")
  out.columnnames <- paste0(c(snpcols, nfcols, statcols, levenecols), collapse = '\t')
  if (outfile != '') {
    filecon <- file(outfile, open = "w")
    writeLines(out.columnnames, filecon)
    close(filecon)
  }
  return(0)
}

gstats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2) {
  outstats <- numeric(0)

  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  bg <- lslinregout$bb[1]
  outstats <- bg
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)
  
  if (teststats[1] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
  
  if (teststats[2] == TRUE) {
    cols <- rep(TRUE,p+q)
    cols[p+1] <- FALSE
    i2.1 <- (lslinregout$xtx[p+1, p+1] - lslinregout$xtx[p+1, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+1]) / s2
    l2 <- sum((resids0*lslinregout$xr[,1])) / s2
    scoreg <- l2 / sqrt(i2.1)
    outstats <- c(outstats, scoreg)
  }
  
  if (teststats[3] == TRUE) {
    seg <- lslinregout$std_err[p + 1]
    tg <- bg / seg
    dfg <- n - p - q
    outstats <- c(outstats, seg, tg, dfg)
  }
  
  if (teststats[4] == TRUE) {
    seg <- sqrt(lslinregout$hws2[p + 1, p + 1])
    zg <- bg / seg
    outstats <- c(outstats, seg, zg)
  }
  return(outstats)
}

gxestats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2) {
  outstats <- numeric(0)
  
  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  bgxe <- lslinregout$bb[2]
  outstats <- bgxe
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)

  if (teststats[1] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
  
  if (teststats[2] == TRUE) {
    cols <- 1:(p+q-1)
    i2.1 <- (lslinregout$xtx[p+2, p+2] - lslinregout$xtx[p+2, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+2]) / s2
    l2 <- sum((resids0*lslinregout$xr[,2])) / s2
    scoreg <- l2 / sqrt(i2.1)
    outstats <- c(outstats, scoreg)
  }
  
  if (teststats[3] == TRUE) {
    segxe <- lslinregout$std_err[p + 2]
    tgxe <- bgxe / segxe
    dfgxe <- n - p - q
    outstats <- c(outstats, segxe, tgxe, dfgxe)
  }
  
  if (teststats[4] == TRUE) {
    segxe <- sqrt(lslinregout$hws2[p + 2, p + 2])
    zgxe <- bgxe / segxe
    outstats <- c(outstats, segxe, zgxe)
  }
  return(outstats)
}

ggxestats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2) {
  outstats <- numeric(0)
  
  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)
  cols <- (p+1):(p+q)
  
  if (teststats[1] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))

  if (teststats[2] == TRUE) {
    i2.1 <- (lslinregout$xtx[cols,cols] - lslinregout$xtx[cols, 1:p] %*% xtxinv0 %*% lslinregout$xtx[1:p, cols]) / s2
    l2 <- c(sum(resids0*lslinregout$xr[,1]), sum(resids0*lslinregout$xr[,2])) / s2
    score2df <- t(l2) %*% solve(i2.1) %*% l2
    outstats <- c(outstats, score2df)
  }
  
  if (teststats[3] == TRUE) {
    outstats <- c(outstats, lslinregout$chi2, lslinregout$xrs2[1,1],
                  lslinregout$xrs2[2,2], lslinregout$xrs2[1,2])
  }

  if (teststats[4] == TRUE) {
    hw2df <- lslinregout$bb %*% solve(lslinregout$hws2[cols,cols]) %*% t(lslinregout$bb)
    outstats <- c(outstats, hw2df, lslinregout$hws2[p+1,p+1],
                  lslinregout$hws2[p+q,p+q], lslinregout$hws2[p+1,p+q])
  }
  return(outstats)
}
