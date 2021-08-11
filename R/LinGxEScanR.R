#' @useDynLib LinGxEScanR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

getdosages <- function(dosage, p0, p1, p2) {
  return (dosage)
}

allocatelinregmem <- function(data, n, p, q) {
  y <- data[,1]
  xl <- as.matrix(data[,1:p])
  xl[,1] <- 1.
  xr <- matrix(0., n, q)
  xtx <- matrix(0., p + q, p + q)
  bt <- matrix(0., p, 1)
  bb <- matrix(0., q, 1)
  ql <- matrix(0., n, p)
  qr <- matrix(0., n, q)
  rtl <- matrix(0., p, p)
  rtr <- matrix(0., p, q)
  rbr <- matrix(0., q, q)
  h <- matrix(0., p, q)
  k <- matrix(0., p, 1)
  t <- matrix(0., n, q)
  zt <- matrix(0, p, 1)
  zb <- matrix(0., q, 1)
  resids0 <- numeric(n)
  resids <- numeric(n)
  sigma2 <- numeric(2)
  s2 <- numeric(2)
  hws2 <- matrix(0., p + q, p + q)
  xtxinv0 <- matrix(0., p, p)
  xtxinv <- matrix(0., p + q, p + q)
  std_err <- numeric(p + q)
  xrs2 <- matrix(0., q, q)
  chi2 <- numeric(1)
  loglike <- numeric(2)
  s2a <- matrix(0., n, n)
  
  return(list(y = y,
              xl = xl,
              xr = xr,
              xtx = xtx,
              bt = bt,
              bb = bb,
              ql = ql,
              qr = qr,
              rtl = rtl,
              rtr = rtr,
              rbr = rbr,
              h = h,
              k = k,
              t = t,
              zt = zt,
              zb = zb,
              resids0 = resids0,
              resids = resids,
              sigma2 = sigma2,
              s2 = s2,
              s2a = s2a,
              hws2 = hws2,
              xtxinv = xtxinv,
              xtxinv0 = xtxinv0,
              std_err = std_err,
              xrs2 = xrs2,
              chi2 = chi2,
              loglike = loglike))
}

initializelslinreg <- function(linregmem) {
  initlslinreg(y = linregmem$y,
               xl = linregmem$xl,
               xtx = linregmem$xtx,
               xtxinv0 = linregmem$xtxinv0,
               ql = linregmem$ql,
               rtl = linregmem$rtl,
               k = linregmem$k,
               zt = linregmem$zt,
               resids = linregmem$resids0,
               sigma2 = linregmem$sigma2,
               s2 = linregmem$s2,
               loglike = linregmem$loglike)
}

runlslinreg <- function(dosage, linregmem, codemask) {
  lslinreg(dosage = dosage,
           y = linregmem$y,
           xl = linregmem$xl,
           xr = linregmem$xr,
           xtx = linregmem$xtx,
           bt = linregmem$bt,
           bb = linregmem$bb,
           ql = linregmem$ql,
           qr = linregmem$qr,
           rtl = linregmem$rtl,
           rtr = linregmem$rtr,
           rbr = linregmem$rbr,
           h = linregmem$h,
           k = linregmem$k,
           t = linregmem$t,
           zb = linregmem$zb,
           resids = linregmem$resids,
           sigma2 = linregmem$sigma2,
           s2 = linregmem$s2,
           hws2 = linregmem$hws2,
           xtxinv = linregmem$xtxinv,
           std_err = linregmem$std_err,
           xrs2 = linregmem$xrs2,
           chi2 = linregmem$chi2,
           loglike = linregmem$loglike)
}

runlslinreg2 <- function(dosage, linregmem, codemask) {
  if (codemask[1] == TRUE) {
    lslinregfit(dosage = dosage,
                y = linregmem$y,
                xl = linregmem$xl,
                xr = linregmem$xr,
                bt = linregmem$bt,
                bb = linregmem$bb,
                ql = linregmem$ql,
                qr = linregmem$qr,
                rtl = linregmem$rtl,
                rtr = linregmem$rtr,
                rbr = linregmem$rbr,
                h = linregmem$h,
                k = linregmem$k,
                t = linregmem$t,
                zb = linregmem$zb)
  }
  if (codemask[2] == TRUE) {
    lslinregresiduals(y = linregmem$y,
                      xl = linregmem$xl,
                      xr = linregmem$xr,
                      bt = linregmem$bt,
                      bb = linregmem$bb,
                      resids = linregmem$resids)
  }
  if (codemask[3] == TRUE) {
    lslinregsigma2(xl = linregmem$xl,
                   xr = linregmem$xr,
                   resids = linregmem$resids,
                   sigma2 = linregmem$sigma2,
                   s2 = linregmem$s2,
                   loglike = linregmem$loglike)
  }
  if (codemask[4] == TRUE) {
    lslinregxtx(xl = linregmem$xl,
                xr = linregmem$xr,
                xtx = linregmem$xtx)
  }
  if (codemask[5] == TRUE) {
    lslinregxtxinv(xtx = linregmem$xtx,
                   xtxinv = linregmem$xtxinv)
  }
  if (codemask[6] == TRUE) {
    lslinregwaldtest(xl = linregmem$xl,
                     xr = linregmem$xr,
                     bb = linregmem$bb,
                     s2 = linregmem$s2,
                     xtxinv = linregmem$xtxinv,
                     std_err = linregmem$std_err,
                     xrs2 = linregmem$xrs2,
                     chi2 = linregmem$chi2)
  }
  if (codemask[7] == TRUE) {
    lslinreghwtest(xl = linregmem$xl,
                   xr = linregmem$xr,
                   resids = linregmem$resids,
                   xtxinv = linregmem$xtxinv,
                   s2a = linregmem$s2a,
                   hws2 = linregmem$hws2)
  }
}

runalllslinreg <- function(dosage, p0, p1, p2,
                           subindex, binarye, eindex0, eindex1,
                           minmac, maxmac,
                           lslinregg, lslinregge, lslinreggxe, lslinregg0gxe,
                           teststats, codemask, levene, snpinfo, snplist, snpnum, outfile) {
  subdose <- dosage[subindex]
  mac <- sum(subdose)
  if (mac < minmac || mac > maxmac) {
    increment(snpnum)
    return (NA)
  }
  if (binarye == TRUE) {
    mac <- sum(dosage[eindex0])
    if (mac < 2 || mac > 2*length(subdose) - 2) {
      increment(snpnum)
      return (NA)
    }
    mac <- sum(dosage[eindex1])
    if (mac < 2 || mac > 2*length(subdose) - 2) {
      increment(snpnum)
      return (NA)
    }
  }
  subp0 <- p0[subindex]
  subp1 <- p1[subindex]
  subp2 <- p2[subindex]
  snpid <- snpinfo$snpid[snplist[snpnum]]
  chromosome <- snpinfo$chromosome[snplist[snpnum]]
  location <- snpinfo$location[snplist[snpnum]]
  reference <- snpinfo$reference[snplist[snpnum]]
  alternate <- snpinfo$alternate[snplist[snpnum]]
  increment(snpnum)
  lslinregg$xr[,1] <- subdose
  runlslinreg2(subdose, lslinregg, codemask[[1]])
  lslinregge$xr[,1] <- subdose
  runlslinreg2(subdose, lslinregge, codemask[[2]])
  lslinreggxe$xr[,1] <- subdose
  lslinreggxe$xr[,2] <- lslinreggxe$xl[,ncol(lslinreggxe$xl)] * subdose
  runlslinreg2(subdose, lslinreggxe, codemask[[3]])
  runlslinreg2(lslinreggxe$xr[,2], lslinregg0gxe, codemask[[4]])
  n <- length(subindex)
  aaf <- mean(subdose) / 2
  if (binarye) {
    n0 <- length(eindex0)
    n1 <- length(eindex1)
    aaf0 <- mean(dosage[eindex0]) / 2
    aaf1 <- mean(dosage[eindex1]) / 2
  }
  gostats <- gstats(teststats[[1]],
                    lslinregg,
                    lslinregg$loglike[1],
                    lslinregg$resids0,
                    lslinregg$xtxinv0,
                    lslinregg$s2[1])
  gestats <- gstats(teststats[[2]],
                    lslinregge,
                    lslinregge$loglike[1],
                    lslinregge$resids0,
                    lslinregge$xtxinv0,
                    lslinregge$s2[1])
  ggxestatsout <- gstats(teststats[[3]],
                      lslinreggxe,
                      lslinregg0gxe$loglike[2],
                      lslinregg0gxe$resids,
                      lslinregg0gxe$xtxinv,
                      lslinregg0gxe$s2[2])
  gxestatsout <- gxestats(teststats[[4]],
                          lslinreggxe,
                          lslinregge$loglike[2],
                          lslinregge$resids,
                          lslinregge$xtxinv,
                          lslinregge$s2[2])
  twodfstats <- twodfstats(teststats[[5]],
                           lslinreggxe,
                           lslinreggxe$loglike[1],
                           lslinreggxe$resids0,
                           lslinreggxe$xtxinv0,
                           lslinreggxe$s2[1])
  levenestats <- numeric(0)
  if (levene[1] == TRUE) {
    w <- levenetest(dosage = subdose,
                    p0 = subp0,
                    p1 = subp1,
                    p2 = subp2,
                    res = lslinregge$resids0)
    levenestats <- c(levenestats, w$f, w$numer, w$denom, w$df1, w$df2)
  }
  if (levene[2] == TRUE) {
    w <- levenetest(dosage = subdose,
                    p0 = subp0,
                    p1 = subp1,
                    p2 = subp2,
                    res = lslinregg$resids0)
    levenestats <- c(levenestats, w$f, w$numer, w$denom, w$df1, w$df2)
  }
  if (outfile == '') {
    if (binarye) {
      return (list(snpid, chromosome, location, reference, alternate,
                   n, aaf, n0, aaf0, n1, aaf1, g0stats, gestats, ggxestatsout,
                   gxestatsout, twodfstats, levenestats))
    }
    return (list(snpid, chromosome, location, reference, alternate,
                 n, aaf, gostats, gestats, ggxestatsout,
                 gxestatsout, twodfstats, levenestats))
  }
  snpout <- paste(snpid, chromosome, location, reference, alternate, sep = '\t')
  if (binarye)
    nfout <- paste(n, aaf, n0, aaf0, n1, aaf1,  sep = '\t')
  else
    nfout <- paste(n, aaf, sep = '\t')

  outline <- paste0(c(snpout, nfout, gostats, gestats, ggxestatsout,
                   gxestatsout, twodfstats, levenestats),collapse = '\t')
  filecon <- file(outfile, open = "a")
  writeLines(outline, filecon)
  close(filecon)
    
  return (0)
}

subsetsnps <- function(snps, snplist) {
  if (length(snps) == 0) {
    stop("No SNPs selected")
  }
  if (is.character(snps) == TRUE) {
    if (length(snps) == 1 & snps[1] == "all")
      return(rep(TRUE, length(snplist)))
    snps2 <- match(snps, snplist)
    snps2 <- snps2[!is.na(snps2)]
    if (length(snps2) == 0)
      stop("No matching SNPs found")
    if (length(snps) != length(snps2))
      print(paste(length(snps) - length(snps2),
                  " of ",
                  length(snps),
                  " SNPs were not found"))
    snps <- snps2
  }
  if (is.numeric(snps) == FALSE)
    stop("snps must be a character or integer array")
  if (is.integer(snps) == FALSE) {
    if (all(floor(snps) == snps) == FALSE)
      stop("snps must be a character or integer array")
  }
  if (min(snps) < 1)
    stop("snp indices must be positive")
  if (max(snps) > length(snplist))
    stop("at least one snp index is greater than the number of SNPs available")
  snpstouse <- rep(FALSE, length(snplist))
  snpstouse[snps] <- TRUE
  
  return (snpstouse)
}

# Subsets the covariate data to complete cases
# Also gets subject indices in the genetic data
subsetdata <- function(subdata, ginfo, mincov) {
  # There must be at least two columns in the subject data
  if(ncol(subdata) < 2)
    stop("There must me at least two columns in the subject data")
  
  # Check if the first column is a character value
  if (is.character(subdata[,1]) == FALSE)
    stop("First column of subject data must be a character value")
  
  # Remove subjects without complete data
  covdata <- subdata[complete.cases(subdata),]
  if (nrow(covdata) == 0)
    stop("No subjects have complete phenotype/covariate data")
  
  colnames(covdata) <- colnames(subdata)
  # Determine if family IDs are used
  # and get the indices of the subjects in the genetic data that are
  # in the covariate data
  if (ginfo$usesfid == FALSE) {
    covdata$genindex <- match(covdata[,1], ginfo$samples$sid)
    phenocol <- 2
  } else {
    # When family ID is used, there must be at least 3 columns
    # in the subject data
    if (ncol(covdata) < 3)
      stop("When using family ID, subject data must have at least 3 columns")
    # When family ID is used, the second columns must be a character value
    if (is.character(covdata[,2]) == FALSE)
      stop("When using family ID, the first two columns must be character values")
    covdata$genindex <- row.match(covdata[,1:2], ginfo$samples)
    phenocol <- 3
  }
  if (ncol(subdata) < phenocol + mincov)
    stop("Subject data has no covariates")
  for (i in phenocol:ncol(subdata)) {
    if (is.numeric(covdata[,i]) == FALSE)
      stop("Phenotype and covariate values must be numeric")
  }
  
  # Drop subjects that don't have genetic values
  covdata <- covdata[complete.cases(covdata),]
  # Check if any subjects have complete data
  if (nrow(covdata) == 0)
    stop("No subjects have complete data")
  
  phenocov <- as.matrix(covdata[,phenocol:(ncol(covdata)-1)])
  dimnames(phenocov) <- list(covdata[,ncol(covdata)],
                             colnames(covdata)[phenocol:(ncol(covdata)-1)])
  eunique <- sort(unique(phenocov[,ncol(phenocov)]))
  if (length(eunique) == 2) {
    binarye <- all(eunique == c(0,1))
  } else {
    binarye <- FALSE
  }
    
  return (list(subdata = phenocov,
               genindex = covdata$genindex,
               binarye = binarye))
}


#####################################################
###          Check for valid input values
#####################################################
validateinput <- function(data, ginfo, outfile, skipfile,
                          minmaf, blksize) {
  # Check if input values are of correct type
  if (is.data.frame(data) == FALSE)
    stop("data must be a data frame")
  if (class(ginfo) != "genetic-info")
    stop("ginfo not a genetic-info class")
  if (class(ginfo$additionalinfo) != "bdose-info" &
      class(ginfo$additionalinfo) != "vcf-info")
    stop("ginfo does not have information about a binary dosage or a vcf file")
  if (is.character(outfile) == FALSE)
    stop("outfile must be a character value")
  if (length(outfile) != 1)
    stop("outfile must be a character vector of length 1")
  if (is.character(skipfile) == FALSE)
    stop("skipfile must be a character value")
  if (length(skipfile) != 1)
    stop("skipfile must be a character vector of length 1")
  if (is.numeric(minmaf) == FALSE)
    stop("minmaf must be a numeric value")
  if (length(minmaf) != 1)
    stop("minmaf must be a numeric vector of length 1")
  if (minmaf < 0 | minmaf > 0.25)
    stop("minmaf must be a value from 0 to 0.25, inclusive")
  if (is.numeric(blksize) == FALSE)
    stop("blksize must be an integer")
  if (length(blksize) != 1)
    stop("blksize must be an integer vector of length 1")
  if (blksize != floor(blksize))
    stop("blksize must be an integer")
  blksize <- as.integer(blksize)
  if (blksize < 0)
    stop("blksize must be greater than or equal to 0")
}

#####################################################
###                  LINGWEIS
#####################################################

#' lingweis
#'
#' Run a linear gweis using genetic data from a binary dosage or vcf file
#' @param data Data frame containing the subject ID, phenotype
#' and covariates
#' @param ginfo Information about the binary dosage or vcf file returned
#' from the BinaryDosage::getbdinfo or BinaryDosage::getvcfinfo routine
#' @param snps The SNPs to be used in the scan. This may be an integer
#' vector indicate which SNPs to use in the binary dosage file or a 
#' character vector of the SNP IDs to use. The value may also be "all",
#' indicating to use all SNPs. The default value is "all".
#' @param outfile The file name for the results Can be blank.
#' If the value is "", the results are returned as a data frame. Default
#' value is ""
#' @param skipfile The name of the file to write the SNPs that were not
#' used and the reason they weren't used. If the value is blank, there is
#' no output of the unused SNPs. Default value is "".
#' @param minmaf Minimum minor allele frequency of SNPs to include
#' in analysis. SNPS that have less than 20 minor alleles observed
#' will be excluded from the analysis regardless of the value of
#' minmaf. A value of 0 indicates to use all the SNPs that have 20
#' minor alleles observed. Default value is 0.
#' @param blksize Size of blocks of SNPs to read in at one time.
#' Larger blocks can improve overall speed but require larger
#' amounts of computer memory. A value of 0 indicates to use the
#' recommended block size. Default value is 0.
#' @return
#' 0
#' @export
#'
#' @examples
#' bdinfo <- readRDS(system.file("extdata/pdata_4_1.bdinfo", package = "GxEScanR"))
#' covdata <- readRDS(system.file("extdata/covdata.rds", package = "GxEScanR"))
#'
#' results <- gweis(data = covdata, ginfo = bdinfo)
lingweis <- function(data, ginfo, snps,
                     gonly, ge, ggxe, gxe, twodf, levene,
                     outfile, skipfile,
                     minmaf, blksize) {
  ## Set values to default for missing input values
  if (missing(data) == TRUE)
    stop("No subject data")
  if (missing(ginfo) == TRUE)
    stop("No genetic data, ginfo")
  if (missing(snps) == TRUE)
    snps <- "all"
  if (missing(gonly) == TRUE)
    gonly <- "Wald"
  if (missing(ge) == TRUE)
    ge <- "WaldHW"
  if (missing(ggxe) == TRUE)
    ggxe <- "WaldHW"
  if (missing(gxe) == TRUE)
    gxe <- "WaldHW"
  if (missing(twodf) == TRUE)
    twodf <- "WaldHW"
  if (missing(levene) == TRUE)
    levene <- c(TRUE, FALSE)
  if (missing(outfile) == TRUE)
    outfile <- ""
  if (missing(skipfile) == TRUE)
    skipfile <- ""
  if (missing(minmaf) == TRUE)
    minmaf <- 0.
  if (missing(blksize))
    blksize <- 0L

  ## Make sure input values are valid  
  validateinput(data, ginfo, outfile, skipfile, minmaf, blksize)
  
  ## Subset the SNPs for analysis
  snps <- subsetsnps(snps = snps,
                     snplist = ginfo$snps$snpid)
  snps <- (1:length(snps))[snps]
  
  ## Subset the subjects for analysis - complete covariate and SNP data
  ## Also match subjects with values in genetic file
  subsetinfo <- subsetdata(subdata = data,
                            ginfo = ginfo,
                            mincov = 1L)
  data <- subsetinfo$subdata
  subindex <- subsetinfo$genindex
  if (subsetinfo$binarye) {
    eindex0 <- subindex[data[,ncol(data)] == 0]
    eindex1 <- subindex[data[,ncol(data)] == 1]
  } else {
    eindex0 <- integer(0)
    eindex1 <- integer(0)
  }
  nsub <- nrow(data)
  ncov <- ncol(data)

  ## Determine test statistics to calculate
  teststats <- vector("list", 3)
  teststats[[1]] <- assignstats(gonly, "gonly")
  teststats[[2]] <- assignstats(ge, "ge")
  teststats[[3]] <- assignstats(ggxe, "ggxe")
  teststats[[4]] <- assignstats(gxe, "gxe")
  teststats[[5]] <- assignstats(twodf, "2df")
  teststats[[5]][1] <- FALSE
  if (is.logical(levene) == FALSE)
    stop("levene must be a logical value")
  if (length(levene) < 1 | length(levene) > 2)
    stop("leven must be a logical vector of length 1 or 2")
  if (length(levene) == 1)
    levene <- rep(levene[1], 2)

  #####################################################
  ###       Calculate the minimum number
  ###       of observed genes
  #####################################################
  minsum <- 2 * nrow(data) * minmaf
  if (minsum < 10)
    minsum <- 10
  maxsum <- 2*nrow(data) - minsum

  #####################################################
  ###       Do initial regressions
  #####################################################
  lslinregg <- allocatelinregmem(data = data,
                                 n = nrow(data),
                                 p = ncol(data) - 1L,
                                 q = 1L)
  lslinregge <- allocatelinregmem(data = data,
                                  n = nrow(data),
                                  p = ncol(data),
                                  q = 1L)
  lslinreggxe <- allocatelinregmem(data = data,
                                   n = nrow(data),
                                   p = ncol(data),
                                   q = 2L)
  lslinregg0gxe <- allocatelinregmem(data = data,
                                   n = nrow(data),
                                   p = ncol(data),
                                   q = 1L)
  initializelslinreg(lslinregg)
  initializelslinreg(lslinregge)
  initializelslinreg(lslinreggxe)
  initializelslinreg(lslinregg0gxe)

  columnnames <- assigncolumnnames(outfile, subsetinfo$binarye, teststats, levene)
  codemask <- assigntests(teststats)

  snpnumber <- integer(1)
  snpnumber[1] <- 1L
  if (class(ginfo$additionalinfo) == "vcf-info") {
    result <- vcfapply2(ginfo,
                        func = runalllslinreg,
                        snps = snps,
                        minmac = minsum,
                        maxmac = maxsum,
                        subindex = subindex,
                        binarye = subsetinfo$binarye,
                        eindex0 = eindex0,
                        eindex1 = eindex1,
                        lslinregg = lslinregg,
                        lslinregge = lslinregge,
                        lslinreggxe = lslinreggxe,
                        lslinregg0gxe = lslinregg0gxe,
                        teststats = teststats,
                        codemask = codemask,
                        levene = levene,
                        snpinfo = ginfo$snps,
                        snplist = snps,
                        snpnum = snpnumber,
                        outfile = outfile)
  } else {
    result <- bdapply2(ginfo,
                       func = runalllslinreg,
                       snps = snps,
                       minmac = minsum,
                       maxmac = maxsum,
                       subindex = subindex,
                       binarye = subsetinfo$binarye,
                       eindex0 = eindex0,
                       eindex1 = eindex1,
                       lslinregg = lslinregg,
                       lslinregge = lslinregge,
                       lslinreggxe = lslinreggxe,
                       lslinregg0gxe = lslinregg0gxe,
                       teststats = teststats,
                       codemask = codemask,
                       levene = levene,
                       snpinfo = ginfo$snps,
                       snplist = snps,
                       snpnum = snpnumber,
                       outfile = outfile)
  }
  return (result)
}