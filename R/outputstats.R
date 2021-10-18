# Determine which statistics will be output for a model
# Also verifies if inputs are valid
# Returns a logical vector

assignstats <- function(statlist, modelname) {
  teststats <- rep(FALSE, 6)
  names(teststats) <- c("fit", "lrt", "score", "robustscore", "Wald", "robustWald")
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

assigncolumnnames <- function(outfile, binarye, teststats, pout, statout, meta, levene) {
  snpcols <- c("SNP", "CHR", "LOC", "REF", "ALT")
  nfcols <- c("n", "aaf")
  if (binarye == TRUE)
    nfcols <- c(nfcols, "n_e0", "aaf_e0", "n_e1", "aaf_e1")
  statcols <- character(0)
  modelname <- c("go", "ge", "gxe")

  statnameg <- vector("list", 6)
  statnameg[[1]] <- "bg"

  statnameg[[2]] <- character(0)
  if (statout == TRUE)
    statnameg[[2]] <- c(statnameg[[2]], "lrt_chisq_bg")
  if (pout == TRUE)
    statnameg[[2]] <- c(statnameg[[2]], "lrt_p_bg")

  if (meta == TRUE)
    statnameg[[3]] <- c("score_bg", "info_bg")
  else
    statnameg[[3]] <- character(0)
  if (pout == TRUE)
    statnameg[[3]] <- c("score_p_bg", statnameg[[3]])
  if (statout == TRUE)
    statnameg[[3]] <- c("score_z_bg", statnameg[[3]])
  
  if (meta == TRUE)
    statnameg[[4]] <- c("robustscore_bg", "robustinfo_bg")
  else
    statnameg[[4]] <- character(0)
  if (pout == TRUE)
    statnameg[[4]] <- c("robustscore_p_bg", statnameg[[4]])
  if (statout == TRUE)
    statnameg[[4]] <- c("robustscore_z_bg", statnameg[[4]])
  
  statnameg[[5]] <- c("wald_se_bg")
  if (statout == TRUE)
    statnameg[[5]] <- c(statnameg[[5]], "wald_t_bg", "wald_df_bg")
  if (pout == TRUE)
    statnameg[[5]] <- c(statnameg[[5]], "wald_p_bg")
  
  statnameg[[6]] <- c("robustwald_se_bg")
  if (statout == TRUE)
    statnameg[[6]] <- c(statnameg[[6]], "robustwald_z_bg")
  if (pout == TRUE)
    statnameg[[6]] <- c(statnameg[[6]], "robustwald_p_bg")
  
  statnamegxe <- vector("list", 6)
  statnamegxe[[1]] <- "bgxe"
  
  statnamegxe[[2]] <- character(0)
  if (statout == TRUE)
    statnamegxe[[2]] <- c(statnamegxe[[2]], "lrt_chisq_bgxe")
  if (pout == TRUE)
    statnamegxe[[2]] <- c(statnamegxe[[2]], "lrt_p_bgxe")
  
  if (meta == TRUE)
    statnamegxe[[3]] <- c("score_bgxe", "info_bgxe")
  else
    statnamegxe[[3]] <- character(0)
  if (pout == TRUE)
    statnamegxe[[3]] <- c("score_p_bgxe", statnamegxe[[3]])
  if (statout == TRUE)
    statnamegxe[[3]] <- c("score_z_bgxe", statnamegxe[[3]])
  
  if (meta == TRUE)
    statnamegxe[[4]] <- c("robustscore_bgxe", "robustinfo_bgxe")
  else
    statnamegxe[[4]] <- character(0)
  if (pout == TRUE)
    statnamegxe[[4]] <- c("robustscore_p_bgxe", statnamegxe[[4]])
  if (statout == TRUE)
    statnamegxe[[4]] <- c("robustscore_z_bgxe", statnamegxe[[4]])
  
  statnamegxe[[5]] <- c("wald_se_bgxe")
  if (statout == TRUE)
    statnamegxe[[5]] <- c(statnamegxe[[5]], "wald_t_bgxe", "wald_df_bgxe")
  if (pout == TRUE)
    statnamegxe[[5]] <- c(statnamegxe[[5]], "wald_p_bgxe")
  
  statnamegxe[[6]] <- c("robustwald_se_bgxe")
  if (statout == TRUE)
    statnamegxe[[6]] <- c(statnamegxe[[6]], "robustwald_z_bgxe")
  if (pout == TRUE)
    statnamegxe[[6]] <- c(statnamegxe[[6]], "robustwald_p_bgxe")

  statname2df <- vector("list", 6)
  statname2df[[1]] <- NA
  
  statname2df[[2]] <- character(0)
  if (statout == TRUE)
    statname2df[[2]] <- c(statname2df[[2]], "lrt_chi2df_joint")
  if (pout == TRUE)
    statname2df[[2]] <- c(statname2df[[2]], "lrt_p_joint")

  if (meta == TRUE)
    statname2df[[3]] <- c("score_bg_joint", "score_bgxe_joint", "info_bg_joint", "info_bgxe_joint", "info_bg_bgxe_joint")
  else
    statname2df[[3]] <- character(0)
  if (pout == TRUE)
    statname2df[[3]] <- c("score_p_joint", statname2df[[3]])
  if (statout == TRUE)
    statname2df[[3]] <- c("score_chi2df_joint", statname2df[[3]])

  if (meta == TRUE)  
    statname2df[[4]] <- c("robustscore_bg_joint", "robustscore_bgxe_joint", "robustinfo_bg_joint", "robustinfo_bgxe_joint", "robustinfo_bg_bgxe_joint")
  else
    statname2df[[4]] <- character(0)
  if (pout == TRUE)
    statname2df[[4]] <- c("robustscore_p_joint", statname2df[[4]])
  if (statout == TRUE)
    statname2df[[4]] <- c("robustscore_chi2df_joint", statname2df[[4]])

  if (meta == TRUE)
    statname2df[[5]] <- c("wald_var_bg_joint", "wald_var_bgxe_joint", "wald_cov_bg_bgxe_joint")
  else
    statname2df[[5]] <- character(0)
  if (statout == TRUE)
    statname2df[[5]] <- c("wald_chi2df_joint", statname2df[[5]])
  if (pout == TRUE)
    statname2df[[5]] <- c("wald_p_joint", statname2df[[5]])
  
  if (meta == TRUE)
    statname2df[[6]] <- c("robustwald_var_bg_joint", "robustwald_var_bgxe_joint", "robustwald_cov_bg_bgxe_joint")
  else
    statname2df[[6]] <- character(0)
  if (statout == TRUE)
    statname2df[[6]] <- c("robustwald_chi2df_joint", statname2df[[6]])
  if (pout == TRUE)
    statname2df[[6]] <- c("robustwald_p_joint", statname2df[[6]])
  
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

gstats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2,
                   pout, statout, meta) {
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

  if (teststats[2] == TRUE) {
    if (statout == TRUE)
      outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(2*(lslinregout$loglike[2] - loglike0), 1, lower.tail = FALSE))
  }
  
  if (teststats[3] == TRUE) {
    cols <- rep(TRUE,p+q)
    cols[p+1] <- FALSE
    i2.1 <- (lslinregout$xtx[p+1, p+1] - lslinregout$xtx[p+1, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+1]) * s2
    l2 <- sum(resids0*lslinregout$xr[,1])
    scoreg <- l2 / sqrt(i2.1)
    if (statout == TRUE)
      outstats <- c(outstats, scoreg)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(scoreg), lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, l2, i2.1)
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
    if (statout == TRUE)
      outstats <- c(outstats, scorehwg)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(scorehwg), lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, testscore, testinfo)
  }
  
  if (teststats[5] == TRUE) {
    seg <- lslinregout$std_err[p + 1]
    tg <- bg / seg
    dfg <- n - p - q
    outstats <- c(outstats, seg)
    if (statout == TRUE)
      outstats <- c(outstats, tg, dfg)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pt(abs(tg), dfg, lower.tail = FALSE))
  }
  
  if (teststats[6] == TRUE) {
    seg <- sqrt(lslinregout$hws2[p + 1, p + 1])
    zg <- bg / seg
    outstats <- c(outstats, seg)
    if (statout == TRUE)
      outstats <- c(outstats, zg)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(zg), lower.tail = FALSE))
  }
  return(outstats)
}

gxestats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2,
                     pout, statout, meta) {
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
    
  if (teststats[2] == TRUE) {
    if (statout == TRUE)
      outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(2*(lslinregout$loglike[2] - loglike0), 1, lower.tail = FALSE))
  }

  if (teststats[3] == TRUE) {
    cols <- 1:(p+q-1)
    i2.1 <- (lslinregout$xtx[p+2, p+2] - lslinregout$xtx[p+2, cols] %*% xtxinv0 %*% lslinregout$xtx[cols, p+2]) * s2
    l2 <- sum((resids0*lslinregout$xr[,2]))
    scoreg <- l2 / sqrt(i2.1)
    if (statout == TRUE)
      outstats <- c(outstats, scoreg)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(scoreg), lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, l2, i2.1)
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
    if (statout == TRUE)
      outstats <- c(outstats, scorehwgxe)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(scorehwgxe), lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, testscore, testinfo)
  }
  
  if (teststats[5] == TRUE) {
    segxe <- lslinregout$std_err[p + 2]
    tgxe <- bgxe / segxe
    dfgxe <- n - p - q
    outstats <- c(outstats, segxe)
    if (statout == TRUE)
      outstats <- c(outstats, tgxe, dfgxe)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pt(abs(tgxe), dfgxe, lower.tail = FALSE))
  }
  
  if (teststats[6] == TRUE) {
    segxe <- sqrt(lslinregout$hws2[p + 2, p + 2])
    zgxe <- bgxe / segxe
    outstats <- c(outstats, segxe)
    if (statout == TRUE)
      outstats <- c(outstats, zgxe)
    if (pout == TRUE)
      outstats <- c(outstats, 2*pnorm(abs(zgxe), lower.tail = FALSE))
  }
  
  return(outstats)
}

twodfstats <- function(teststats, lslinregout, loglike0, resids0, xtxinv0, s2,
                       pout, statout, meta) {
  outstats <- numeric(0)
  
  if (all(teststats == FALSE) == TRUE)
    return(outstats)
  
  n <- nrow(lslinregout$xl)
  p <- ncol(lslinregout$xl)
  q <- ncol(lslinregout$xr)
  cols <- (p+1):(p+q)
  
  if (teststats[2] == TRUE) {
    if (statout == TRUE)
      outstats <- c(outstats, 2*(lslinregout$loglike[2] - loglike0))
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(2*(lslinregout$loglike[2] - loglike0), 2, lower.tail = FALSE))
  }

  if (teststats[3] == TRUE) {
    i2.1 <- (lslinregout$xtx[cols,cols] - lslinregout$xtx[cols, 1:p] %*% xtxinv0 %*% lslinregout$xtx[1:p, cols]) * s2
    l2 <- c(sum(resids0*lslinregout$xr[,1]), sum(resids0*lslinregout$xr[,2]))
    score2df <- t(l2) %*% solve(i2.1) %*% l2
    if (statout == TRUE)
      outstats <- c(outstats, score2df)
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(score2df, 2, lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, l2[1], l2[2], i2.1[1,1], i2.1[2,2], i2.1[2,1])
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
    if (statout == TRUE)
      outstats <- c(outstats, score2df)
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(score2df, 2, lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, uhwg[1,1], uhwg[2,1], testinfo[1,1], testinfo[2,2], testinfo[2,1])
  }
  
  if (teststats[5] == TRUE) {
    if (statout == TRUE)
      outstats <- c(outstats, lslinregout$chi2)
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(lslinregout$chi2,2,lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, lslinregout$xrs2[1,1], lslinregout$xrs2[2,2], lslinregout$xrs2[1,2])
  }

  if (teststats[6] == TRUE) {
    hw2df <- t(lslinregout$bb) %*% solve(lslinregout$hws2[cols,cols]) %*% lslinregout$bb
    if (statout == TRUE)
      outstats <- c(outstats, hw2df)
    if (pout == TRUE)
      outstats <- c(outstats, pchisq(hw2df,2,lower.tail = FALSE))
    if (meta == TRUE)
      outstats <- c(outstats, lslinregout$hws2[p+1,p+1], lslinregout$hws2[p+q,p+q], lslinregout$hws2[p+1,p+q])
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
                              "xtx", "xtxinv", "Wald", "robustWald")
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
  if (teststats$ggxe["robustscore"] == TRUE)
    codemask$g0gxe <- codemask$g0gxe | testmask[[8]]
  
  # 4 is the two df test
  if (teststats$gxe["lrt"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[7]]
  if (teststats$gxe["score"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[8]]
  if (teststats$gxe["robustscore"] == TRUE)
    codemask$ge <- codemask$ge | testmask[[8]]
  
  
  return(codemask)
}
