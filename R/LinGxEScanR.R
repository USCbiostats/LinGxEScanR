#' @useDynLib LinGxEScanR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

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
  
  phenocov = as.matrix(covdata[,phenocol:(ncol(covdata)-1)])
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
  data <- subsetdata(subdata = data,
                     ginfo = ginfo,
                     mincov = 1)
  nsub <- nrow(data)
  ncov <- ncol(data)
  
  return (data)
}