# Determine which statistics will be output for a model
# Also verifies if inputs are valid
# Returns a logical vector

assignstats <- function(statlist, modelname) {
  teststats <- rep(FALSE, 6)
  names(teststats) <- c("fit", "lrt", "score", "scoreHW", "Wald", "WaldHW")
  if (is.character(statlist) == FALSE)
    stop(paste(modelname, "must be a character vector"))
  if (length(statlist) == 0)
    return (teststats)
  if (length(statlist) == 1) {
    if(statlist == "all") {
      teststats[1:6] <- TRUE
      return (teststats)
    }
    if (statlist == "" | statlist == "none")
      return (teststats)
  }
  statsidx <- match(statlist, names(teststats))
  if (all(is.na(statsidx) == FALSE) == FALSE)
    stop(paste("Unknown values in", modelname))
  teststats[statsidx] <- TRUE
  teststats[1] <- teststats[1] || teststats[2] || teststats[5] || teststats[6]
  return(teststats)
}

assigncolumnnames <- function(outfile, binarye, teststats, levene) {
  snpcols <- c("SNP", "CHR", "LOC", "REF", "ALT")
  nfcols <- c("n", "aaf")
  if (binarye == TRUE)
    nfcols <- c(nfcols, "n_e0", "aaf_e0", "n_e1", "aaf_e1")
  statcols <- character(0)
  modelname <- c("go", "ge", "gxe")
  statnameg <- vector("list", 6)
#  names(statnameg) <- c("beta", "lrt", "score", "scoreHW", "Wald", "WaldHW")
  statnameg[[1]] <- "bg"
  statnameg[[2]] <- "lrt_chi_g"
  statnameg[[3]] <- c("score_z_g", "score_bg", "info_bg")
  statnameg[[4]] <- c("scorehw_z_g", "scorehw_bg", "infohw_bg")
  statnameg[[5]] <- c("wald_se_g", "wald_t_g", "wald_df_g")
  statnameg[[6]] <- c("waldhw_se_g", "waldhw_z_g")
  statnamegxe <- vector("list", 6)
  statnamegxe[[1]] <- "bgxe"
  statnamegxe[[2]] <- "lrt_chi_gxe"
  statnamegxe[[3]] <- c("score_z_gxe", "score_bgxe", "info_bgxe")
  statnamegxe[[4]] <- c("scorehw_z_gxe", "scorehw_bgxe", "infohw_bgxe")
  statnamegxe[[5]] <- c("wald_se_gxe", "wald_t_gxe", "wald_df_gxe")
  statnamegxe[[6]] <- c("waldhw_se_gxe", "waldhw_z_gxe")
  statname2df <- vector("list", 5)
  statname2df[[1]] <- NA
  statname2df[[2]] <- "lrt_chi2df"
  statname2df[[3]] <- c("score_chi2df", "score_bg_joint", "score_bgxe_joint", "info_bg_joint", "info_bgxe_joint", "info_bg_bgxe_joint")
  statname2df[[4]] <- c("scorehw_chi2df", "scorehw_bg_joint", "scorehw_bgxe_joint", "infohw_bg_joint", "infohw_bgxe_joint", "infohw_bg_bgxe_joint")
  statname2df[[5]] <- c("wald_chi2df", "wald_var_g", "wald_var_gxe", "wald_cov_ggxe")
  statname2df[[6]] <- c("waldhw_chi2df", "waldhw_var_g", "waldhw_var_gxe", "waldhw_cov_ggxe")
  for (i in 1:3) {
    if (all(teststats[[i]] == FALSE) == TRUE)
      next
    for (j in 1:length(teststats[[i]])) {
      if (teststats[[i]][j] == TRUE)
        statcols <- c(statcols, paste(statnameg[[j]], modelname[i], sep = '_'))
    }
  }
  if (all(teststats[[4]] == FALSE) == FALSE) {
    for (j in 1:length(teststats[[4]])) {
      if (teststats[[4]][j] == TRUE)
        statcols <- c(statcols, statnamegxe[[j]])
    }
  }
  if (all(teststats[[5]] == FALSE) == FALSE) {
    for (j in 1:length(teststats[[5]])) {
      if (teststats[[5]][j] == TRUE)
        statcols <- c(statcols, statname2df[[j]])
    }
  }
  
  levenecols <- character(0)
  if (levene[1] == TRUE)
    levenecols <- c(levenecols, "f_levene_e", "chisq_num_e", "chisq_den_e",
                    "df_num_e", "df_den_e")
  if (levene[2] == TRUE)
    levenecols <- c(levenecols, "f_levene_ne", "chisq_num_ne", "chisq_den_ne",
                    "df_num_ne", "df_den_ne")
  if (outfile != '') {
    out.columnnames <- paste0(c(snpcols, nfcols, statcols, levenecols), collapse = '\t')
    filecon <- file(outfile, open = "w")
    writeLines(out.columnnames, filecon)
    close(filecon)
  } else {
    out.columnnames <- c(snpcols, nfcols, statcols, levenecols)
  }
  return(out.columnnames)
}

gstats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2) {
  outstats <- numeric(0)

  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)
  
  if (teststats[1] == TRUE){
    bg <- lslinregout$bb[1]
    outstats <- bg
  }

  if (teststats[2] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
  
  if (teststats[3] == TRUE) {
    cols <- rep(TRUE,p+q)
    cols[p+1] <- FALSE
    i2.1 <- (lslinregout$xtx[p+1, p+1] - lslinregout$xtx[p+1, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+1]) * s2
    l2 <- sum(resids0*lslinregout$xr[,1])
    scoreg <- l2 / sqrt(i2.1)
    outstats <- c(outstats, scoreg, l2, i2.1)
  }
  
  if (teststats[4] == TRUE) {
    cols <- rep(TRUE,p+q)
    cols[p+1] <- FALSE
    lslinreguut(xl = lslinregout$xl,
                xr = lslinregout$xr,
                resids = resids0,
                s2a = lslinregout$s2a,
                uut = lslinregout$uut)
    cmat <- matrix(0, 1, p+q)
    cmat[1,cols] <- -lslinregout$xtx[p+1, cols] %*% xtxinv0
    cmat[1,p+1] <- 1
    testscore <- sum(resids0*lslinregout$xr[,1])
    testinfo <- cmat %*% lslinregout$uut %*% t(cmat)
    scorehwg <- testscore / sqrt(testinfo)
    outstats <- c(outstats, scorehwg, testscore, testinfo)
  }
  
  if (teststats[5] == TRUE) {
    seg <- lslinregout$std_err[p + 1]
    tg <- bg / seg
    dfg <- n - p - q
    outstats <- c(outstats, seg, tg, dfg)
  }
  
  if (teststats[6] == TRUE) {
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

  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)

  if (teststats[1] == TRUE){
    bgxe <- lslinregout$bb[2]
    outstats <- bgxe
  }
    
  if (teststats[2] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
  
  if (teststats[3] == TRUE) {
    cols <- 1:(p+q-1)
    i2.1 <- (lslinregout$xtx[p+2, p+2] - lslinregout$xtx[p+2, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+2]) * s2
    l2 <- sum((resids0*lslinregout$xr[,2]))
    scoreg <- l2 / sqrt(i2.1)
    outstats <- c(outstats, scoreg, l2, i2.1)
  }
  
  if (teststats[4] == TRUE) {
    lslinreguut(xl = lslinregout$xl,
                xr = lslinregout$xr,
                resids = resids0,
                s2a = lslinregout$s2a,
                uut = lslinregout$uut)
    cmat <- matrix(0, 1, p+q)
    cmat[1,1:(p+1)] <- -lslinregout$xtx[p+2, 1:(p+1)] %*% xtxinv0
    cmat[1,(p+q)] <- 1
    testscore <- sum(resids0*lslinregout$xr[,2])
    testinfo <- cmat %*% lslinregout$uut %*% t(cmat)
    scorehwgxe <- testscore / sqrt(testinfo)
    outstats <- c(outstats, scorehwgxe, testscore, testinfo)
  }
  
  if (teststats[5] == TRUE) {
    segxe <- lslinregout$std_err[p + 2]
    tgxe <- bgxe / segxe
    dfgxe <- n - p - q
    outstats <- c(outstats, segxe, tgxe, dfgxe)
  }
  
  if (teststats[6] == TRUE) {
    segxe <- sqrt(lslinregout$hws2[p + 2, p + 2])
    zgxe <- bgxe / segxe
    outstats <- c(outstats, segxe, zgxe)
  }
  
  return(outstats)
}

twodfstats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2) {
  outstats <- numeric(0)
  
  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)
  cols <- (p+1):(p+q)
  
  if (teststats[2] == TRUE)
    outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))

  if (teststats[3] == TRUE) {
    i2.1 <- (lslinregout$xtx[cols,cols] - lslinregout$xtx[cols, 1:p] %*% xtxinv0 %*% lslinregout$xtx[1:p, cols]) * s2
    l2 <- c(sum(resids0*lslinregout$xr[,1]), sum(resids0*lslinregout$xr[,2]))
    score2df <- t(l2) %*% solve(i2.1) %*% l2
    outstats <- c(outstats, score2df, l2[1], l2[2], i2.1[1,1], i2.1[2,2], i2.1[2,1])
  }
  
  if (teststats[4] == TRUE) {
    lslinreguut(xl = lslinregout$xl,
                xr = lslinregout$xr,
                resids = resids0,
                s2a = lslinregout$s2a,
                uut = lslinregout$uut)
    cmat <- matrix(0, 2, p+q)
    cmat[1:2,1:p] <- -lslinregout$xtx[(p+1):(p+2), 1:p] %*% xtxinv0
    cmat[1:2,(p+1):(p+2)] <- diag(1,2)
    uhwg <- matrix(0,2,1)
    uhwg[1,1] <- sum(resids0*lslinregout$xr[,1])
    uhwg[2,1] <- sum(resids0*lslinregout$xr[,2])
    testinfo <- cmat %*% lslinregout$uut %*% t(cmat)
    score2df <- t(uhwg) %*% solve(testinfo) %*% uhwg
    outstats <- c(outstats, score2df, uhwg[1,1], uhwg[2,1], testinfo[1,1], testinfo[2,2], testinfo[2,1])
  }
  
  if (teststats[5] == TRUE) {
    outstats <- c(outstats, lslinregout$chi2, lslinregout$xrs2[1,1],
                  lslinregout$xrs2[2,2], lslinregout$xrs2[1,2])
  }

  if (teststats[6] == TRUE) {
    hw2df <- t(lslinregout$bb) %*% solve(lslinregout$hws2[cols,cols]) %*% lslinregout$bb
    outstats <- c(outstats, hw2df, lslinregout$hws2[p+1,p+1],
                  lslinregout$hws2[p+q,p+q], lslinregout$hws2[p+1,p+q])
  }
  return(outstats)
}

assigntests <- function(teststats) {
  testmask <- vector("list", 8)
  testmask[[1]] <- c( TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  testmask[[2]] <- c( TRUE,  TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE)
  testmask[[3]] <- c(FALSE, FALSE, FALSE,  TRUE, FALSE, FALSE, FALSE)
  testmask[[4]] <- c(FALSE, FALSE, FALSE,  TRUE, FALSE, FALSE, FALSE)
  testmask[[5]] <- c( TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE, FALSE)
  testmask[[6]] <- c( TRUE,  TRUE,  TRUE,  TRUE,  TRUE, FALSE,  TRUE)
  testmask[[7]] <- c( TRUE,  TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE)
  testmask[[8]] <- c( TRUE,  TRUE,  TRUE,  TRUE,  TRUE, FALSE, FALSE)
  testtomask <- c(1L, 2L, 3L, 3L, 3L)
  
  codemask <- vector("list", 4)
  names(codemask) = c("gonly", "ge", "gxe", "g0gxe")
  for(i in 1L:4L) {
    codemask[[i]] <- rep(FALSE, 7)
    names(codemask[[i]]) <- c("fit", "resids", "sigma2",
                              "xtx", "xtxinv", "Wald", "WaldHW")
  }
  for (i in 1L:5L) {
    for (j in 1L:6L) {
      if (teststats[[i]][j] == TRUE) {
        codemask[[testtomask[i]]] <- codemask[[testtomask[i]]] | testmask[[j]]
      }
    }
  }
  
  if (teststats$ggxe["lrt"] == TRUE)
    codemask$g0gxe <- codemask$g0gxe | testmask[[7]]
  if (teststats$ggxe["score"] == TRUE)
    codemask$g0gxe <- codemask$g0gxe | testmask[[8]]
  if (teststats$ggxe["scoreHW"] == TRUE)
    codemask$g0gxe <- codemask$g0gxe | testmask[[8]]
  
  # 4 is the two df test
  if (teststats$gxe["lrt"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[7]]
  if (teststats$gxe["score"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[8]]
  if (teststats$gxe["scoreHW"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[8]]
  
  
  return(codemask)
}
