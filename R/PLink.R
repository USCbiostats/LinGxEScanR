###########################################################
#                Binary Dosage                            #
###########################################################

#' Apply a function to each SNP in a
#' plink binary file
#'
#' A routine that reads in the SNP data serially
#' from a binary dosage file and applies a user
#' specified function to the data.
#'
#' @param plinkinfo List with information about the
#' binary dosage file returned from getplinkinfo
#' @param func A user supplied function to apply
#' to the data for each snp. The function must be
#' provide with the following parameters, dosage,
#' p0, p1, and p2, where dosage is the dosage values
#' for each subject and p0, p1, and p2 are the
#' probabilities that a subject has zero, one,
#' and two copies of the alternate allele,
#' respectively.
#' @param snps List of SNPs to perform the function
#' on. Null value indicates to perform function all
#' all SNPs. Default NULL.
#' @param ... Additional parameters needed by the
#' user supplied function
#'
#' @return A list with length equal to the number
#' of SNPs in the binary dosage file. Each element
#' of the list is the value returned by the user
#' supplied function
#' @export
#' @family Iterating functions
#' @examples
#' # Get information about a plink binary file
#'
#' famfile <- system.file("extdata", "test.fam", package = "LinGxEScanR")
#' bimfile <- system.file("extdata", "test.bim", package = "LinGxEScanR")
#' bedfile <- system.file("extdata", "test.bed", package = "LinGxEScanR")
#' plinkinfo <- getplinkinfo(famfile, bimfile, bedfile)
#'
#' # Apply the getaaf, get alternate allele frequency, function
#' # to all the SNPs in the binary dosage file
#'
#' aaf <- plinkapply2(plinkinfo = plinkinfo,
#'                func = BinaryDosage:::getaaf)
plinkapply2 <- function(plinkinfo, func, snps, ...) {
  if (missing(plinkinfo) == TRUE)
    stop("No plink binary file information specified")
  if (is.na(match("genetic-info", class(plinkinfo))) == TRUE)
    stop("plinkinfo does not contain information about a plink binary file")
  if (is.na(match("plink-info", class(plinkinfo$additionalinfo))) == TRUE)
    stop("plinkinfo does not contain information about a plink binary file")
  
  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")
  
  if (missing(snps) == TRUE | length(snps) == 0)
    snps <- 1:nrow(plinkinfo$snps)
  if (is.numeric(snps) == FALSE)
    stop("snps must be a numeric vector")
  if (all(as.integer(snps) == snps) == FALSE)
    stop ("snps must be an integer vector")
  snps <- as.integer(snps)
  if (min(snps) < 0 | max(snps) > nrow(plinkinfo$snps))
    stop ("snps contains value out of range")
  
  retval <- vector("list", length(snps))
  dosage <- numeric(nrow(plinkinfo$samples))
  p0 <- numeric(nrow(plinkinfo$samples))
  p1 <- numeric(nrow(plinkinfo$samples))
  p2 <- numeric(nrow(plinkinfo$samples))
  p0[1:nrow(plinkinfo$samples)] <- NA
  p1[1:nrow(plinkinfo$samples)] <- NA
  p2[1:nrow(plinkinfo$samples)] <- NA
  for (i in snps) {
    dosage[1:nrow(plinkinfo$samples)] <- as.numeric(getplinksnp(plinkinfo, i))
    retval[[i]] <- func(dosage, p0, p1, p2, ...)
  }
  return (retval)
}

#' Title
#'
#' @param plinkinfo Information about a binary plink data set returned from
#' getplinkinof
#' @param snp SNP name or number
#'
#' @return Integer vector with SNP values
#' @export
#'
#' @examples
#' famfile <- system.file("extdata", "test.fam", package = "LinGxEScanR")
#' bimfile <- system.file("extdata", "test.bim", package = "LinGxEScanR")
#' bedfile <- system.file("extdata", "test.bed", package = "LinGxEScanR")
#' plinkinfo <- getplinkinfo(famfile, bimfile, bedfile)
#'
#' snpvalue <- getplinksnp(plinkinfo, "snp2")
getplinksnp <- function(plinkinfo, snp) {
  if (class(plinkinfo) != "genetic-info")
    stop("plinkinfo not a genetic-info object")
  if (class(plinkinfo$additionalinfo) != "plink-info")
    stop("plinkinfo not a plink-info object")
  if (is.character(snp) == TRUE) {
    snp <- match(snp, plinkinfo$snps$snpid)
    if (is.na(snp) == TRUE)
      stop("Cannot find SNP")
  }
  infile <- file(plinkinfo$filename, "rb")
  seek(infile, 3 + (snp - 1)*plinkinfo$datasize)
  snpbytes <- readBin(infile, "raw", plinkinfo$datasize)
  close(infile)
  
  snpvalues <- integer(nrow(plinkinfo$samples))
  snprawtoint(snpbytes, snpvalues)
  
  return(snpvalues)
}

#' Title
#'
#' @param famfile name of plink family file
#' @param bimfile name of plink binary map file
#' @param bedfile anme of plink binary SNP file
#'
#' @return List with information about binary plink data set
#' @export
#'
#' @examples
#' famfile <- system.file("extdata", "test.fam", package = "PlinkBin")
#' bimfile <- system.file("extdata", "test.bim", package = "PlinkBin")
#' bedfile <- system.file("extdata", "test.bed", package = "PlinkBin")
#' plinkinfo <- getplinkinfo(famfile, bimfile, bedfile)
getplinkinfo <- function(famfile, bimfile, bedfile) {
  pheader <- as.raw(c(0x6c, 0x1b, 0x01))
  famcolclasses <- c(rep("character", 4), "integer", "numeric")
  famdf <- read.table(famfile, colClasses = famcolclasses)
  colnames(famdf) <- c("fid", "sid", "pid", "mid", "sex", "status")
  mapcolclasses <- c(rep("character", 2), "numeric", "integer", rep("character", 2))
  mapdf <- read.table(bimfile, colClasses = mapcolclasses)
  colnames(mapdf) <- c("chromosome", "snpid", "cm", "location", "reference", "alternate")
  plinkfilename <- normalizePath(bedfile, winslash = "/")
  plinkheader <- readplinkheader(plinkfilename)
  if (all(plinkheader == pheader) == FALSE)
    stop("bedfile is not a plink binary file")
  plinkinfo <- list(filename = plinkfilename,
                    usesfid = TRUE,
                    samples = famdf,
                    onechr = all(mapdf$chromosome == mapdf$chromosome[1]),
                    snpidformat = 0,
                    snps = mapdf[,c(1,4,2,3,5,6)],
                    snpinfo = list(),
                    additionalinfo = list(famfile = normalizePath(famfile, winslash = "/"),
                                          mapfile = normalizePath(bimfile, winslash = "/"),
                                          headersize = 3),
                    datasize = as.integer(floor((nrow(famdf) + 3)/4)),
                    indices = numeric(0))
  class(plinkinfo) <-  "genetic-info"
  class(plinkinfo$additionalinfo) <- "plink-info"
  return(plinkinfo)
}

readplinkheader <- function(bedfile) {
  infile <- file(bedfile, "rb")
  plinkheader <- readBin(infile, raw(), 3)
  close(infile)
  return(plinkheader)
}
