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
  bt <- matrix(0., 1, p)
  bb <- matrix(0., 1, q)
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
  resids <- numeric(n)
  s2 <- numeric(1)
  xtxinv <- matrix(0., p + q, p + q)
  std_err <- numeric(p + q)
  xrs2 <- matrix(0., q, q)
  chi2 <- numeric(1)
  
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
              resids = resids,
              s2 = s2,
              xtxinv = xtxinv,
              std_err = std_err,
              xrs2 = xrs2,
              chi2 = chi2))
}

initializelslinreg <- function(linregmem) {
  initlslinreg(y = linregmem$y,
               xl = linregmem$xl,
               xtx = linregmem$xtx,
               ql = linregmem$ql,
               rtl = linregmem$rtl,
               k = linregmem$k,
               zt = linregmem$zt,
               resids = linregmem$resids,
               s2 = linregmem$s2)
}

runlslinreg <- function(dosage, linregmem) {
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
           s2 = linregmem$s2,
           xtxinv = linregmem$xtxinv,
           std_err = linregmem$std_err,
           xrs2 = linregmem$xrs2,
           chi2 = linregmem$chi2)
}

runalllslinreg <- function(dosage, p0, p1, p2, subindex, minmac, maxmac,
                        lslinregg, lslinregge, lslinreggxe, leveneresids,
                        snpinfo, snplist, snpnum, outfile) {
  subdose <- dosage[subindex]
  if (sum(subdose) < minmac || sum(subdose) > maxmac)
    return (NA)
  snpid <- snpinfo$snpid[snplist[snpnum]]
  chromosome <- snpinfo$chromosome[snplist[snpnum]]
  location <- snpinfo$location[snplist[snpnum]]
  reference <- snpinfo$reference[snplist[snpnum]]
  alternate <- snpinfo$alternate[snplist[snpnum]]
  increment(snpnum)
  runlslinreg(subdose, lslinregg)
  runlslinreg(subdose, lslinregge)
  runlslinreg(subdose, lslinreggxe)
  n <- length(subindex)
  p_1 <- ncol(lslinregg$xl)
  q_1 <- ncol(lslinregg$xr)
  p_2 <- ncol(lslinregge$xl)
  q_2 <- ncol(lslinregge$xr)
  p_3 <- ncol(lslinreggxe$xl)
  q_3 <- ncol(lslinreggxe$xr)
  bg <- lslinregg$bb[1]
  seg <- lslinregg$std_err[p_1 + 1]
  tg <- lslinregg$bb[1] / lslinregg$std_err[p_1 + 1]
  dfg <- n - p_1 - q_1
  bge <- lslinregge$bb[1]
  seg <- lslinregge$std_err[p_2 + 1]
  tge <- lslinregge$bb[1] / lslinregge$std_err[p_2 + 1]
  dfge <- n - p_2 - q_2
  bggxe <- lslinreggxe$bb[1]
  seggxe <- lslinreggxe$std_err[p_3 + 1]
  tggxe <- lslinreggxe$bb[1] / lslinreggxe$std_err[p_3 + 1]
  dfggxe <- n - p_3 - q_3
  bgxe <- lslinreggxe$bb[2]
  segxe <- lslinreggxe$std_err[p_3 + 2]
  tgxe <- lslinreggxe$bb[2] / lslinreggxe$std_err[p_3 + 2]
  dfgxe <- n - p_3 - q_3
  varg <- lslinreggxe$xrs2[1,1]
  vargxe <- lslinreggxe$xrs2[2,2]
  covggxe <- lslinreggxe$xrs2[1,2]
  chisqggxe <- lslinreggxe$chi2
  dfchisqggxe <- 2L
  w <- levenetest(dosage = dosage,
                  p0 = p0,
                  p1 = p1,
                  p2 = p2,
                  res = leveneresids)
  if (outfile == '') {
    return (list(snpid, chromosome, location, reference, alternate,
                 n, bg, seg, tg, dfg, bge, seg, tge, dfge,
                 bggxe, seggxe, tggxe, dfggxe, bgxe, segxe, tgxe, dfgxe,
                 chisqggxe, dfchisqggxe, varg, vargxe, covggxe,
                 w$f, w$numer, w$denom, w$df1, w$df2))
  }
  outline <- paste(snpid, chromosome, location, reference, alternate,
                   n, bg, seg, tg, dfg, bge, seg, tge, dfge,
                   bggxe, seggxe, tggxe, dfggxe, bgxe, segxe, tgxe, dfgxe,
                   chisqggxe, dfchisqggxe, varg, vargxe, covggxe,
                   w$f, w$numer, w$denom, w$df1, w$df2, sep = '\t')
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
  return (phenocov)
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
lingweis <- function(data, ginfo, snps, outfile, skipfile,
                     minmaf, blksize) {
  if (missing(data) == TRUE)
    stop("No subject data")
  if (missing(ginfo) == TRUE)
    stop("No genetic data, ginfo")
  if (missing(snps) == TRUE)
    snps = "all"
  if (missing(outfile) == TRUE)
    outfile <- ""
  if (missing(skipfile) == TRUE)
    skipfile <- ""
  if (missing(minmaf) == TRUE)
    minmaf <- 0.
  if (missing(blksize))
    blksize <- 0L
  validateinput(data, ginfo, outfile, skipfile, minmaf, blksize)
  snps <- subsetsnps(snps = snps,
                     snplist = ginfo$snps$snpid)
  snps <- (1:length(snps))[snps]
  data <- subsetdata(subdata = data,
                     ginfo = ginfo,
                     mincov = 1L)
  nsub <- nrow(data)
  ncov <- ncol(data)
  
  subindex <- match(subdata[,1], ginfo$samples$sid)
  
  #####################################################
  ###       Calcualte the minimum number
  ###       of observed genes
  #####################################################
  if (minmaf == 0)
    minsum <- 2 * nrow(data) * 0.05
  else
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
  initializelslinreg(lslinregg)
  initializelslinreg(lslinregge)
  initializelslinreg(lslinreggxe)
  leveneresids <- numeric(nrow(data))
  leveneresids[1:nrow(data)] <- lslinregge$resids

  if (outfile != '') {
    columnnames <- paste("SNP", "CHR", "LOC", "REF", "ALT",
                         "n", "bg_g", "seg_g", "wtg_g", "dfg_g",
                         "bg_e", "seg_e", "wtg_e", "dfg_e",
                         "bg_gxe", "seg_gxe", "wtg_gxe", "dfg_gxe",
                         "bgxe", "segxe", "wtgxe", "dfgxe",
                         "wchisqggxe", "dfwchisqggxe", "varg", "vargxe", "covggxe",
                         "W", "wnumer", "wdenom", "wdf1", "wdf2",
                         sep = '\t')
    filecon <- file(outfile, open = "w")
    writeLines(columnnames, filecon)
    close(filecon)
  }
  snpnumber <- integer(1)
  snpnumber[1] <- 1L
  if (class(ginfo$additionalinfo) == "vcf-info") {
    result <- vcfapply2(ginfo,
                        func = runalllslinreg,
                        snps = snps,
                        minmac = minsum,
                        maxmac = maxsum,
                        subindex = subindex,
                        lslinregg = lslinregg,
                        lslinregge = lslinregge,
                        lslinreggxe = lslinreggxe,
                        leveneresids = leveneresids,
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
                       lslinregg = lslinregg,
                       lslinregge = lslinregge,
                       lslinreggxe = lslinreggxe,
                       leveneresids = leveneresids,
                       snpinfo = ginfo$snps,
                       snplist = snps,
                       snpnum = snpnumber,
                       outfile = outfile)
  }
  return (result)
}