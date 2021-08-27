###########################################################
#                Binary Dosage                            #
###########################################################

#' Apply a function to each SNP in a
#' binary dosage file
#'
#' A routine that reads in the SNP data serially
#' from a binary dosage file and applies a user
#' specified function to the data.
#'
#' @param bdinfo List with information about the
#' binary dosage file returned from getbdinfo
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
#' # Get information about a binary dosage file
#'
#' vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
#'
#' # Apply the getaaf, get alternate allele frequency, function
#' # to all the SNPs in the binary dosage file
#'
#' aaf <- bdapply2(bdinfo = bdinfo,
#'                func = BinaryDosage:::getaaf)
bdapply2 <- function(bdinfo, func, snps, ...) {
  if (missing(bdinfo) == TRUE)
    stop("No binary dosage file information specified")
  if (is.na(match("genetic-info", class(bdinfo))) == TRUE)
    stop("bdinfo does not contain information about a binary dosage file")
  if (is.na(match("bdose-info", class(bdinfo$additionalinfo))) == TRUE)
    stop("bdinfo does not contain information about a binary dosage file")
  
  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")

  if (missing(snps) == TRUE | length(snps) == 0)
    snps <- 1:nrow(bdinfo$snps)
  if (is.numeric(snps) == FALSE)
    stop("snps must be a numeric vector")
  if (all(as.integer(snps) == snps) == FALSE)
    stop ("snps must be an integer vector")
  snps <- as.integer(snps)
  if (min(snps) < 0 | max(snps) > nrow(bdinfo$snps))
    stop ("snps contains value out of range")
  
  retval <- vector("list", length(snps))
  dosage <- numeric(nrow(bdinfo$samples))
  p0 <- numeric(nrow(bdinfo$samples))
  p1 <- numeric(nrow(bdinfo$samples))
  p2 <- numeric(nrow(bdinfo$samples))
  us <- integer(2 * nrow(bdinfo$samples))
  for (i in snps) {
    dosage[1:nrow(bdinfo$samples)] <- NA
    p0[1:nrow(bdinfo$samples)] <- NA
    p1[1:nrow(bdinfo$samples)] <- NA
    p2[1:nrow(bdinfo$samples)] <- NA
    BinaryDosage:::ReadBinaryDosageData(bdinfo, i, dosage, p0, p1, p2, us)
    retval[[i]] <- func(dosage, p0, p1, p2, ...)
  }
  return (retval)
}

###########################################################
#                        VCF                              #
###########################################################

#' Apply a function to each SNP in a vcf file
#'
#' A routine that reads in the SNP data serially
#' from a vcf file and applies a user specified
#' function to the data for the selected SNPs.
#'
#' @param vcfinfo List with information about the
#' vcf file returned from getvcfinfo
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
#' @param blksz Number of lines to read in on
#' each iteration. A larger number uses more memory
#' but can reduce run time. Maximum value 1000.
#' Default 100.
#' @param ... Additional parameters needed by the
#' user supplied function
#'
#' @return A list with length equal to the number
#' of SNPs in the vcf file. Each element
#' of the list is the value returned by the user
#' supplied function
#' @export
#' @family Iterating functions
#' @examples
#' # Get information about a vcf file
#'
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' vcfinfo <- getvcfinfo(vcffiles = vcf1afile)
#'
#' # Apply the getaaf, get alternate allele frequency, function
#' # to all the SNPs in the vcf file
#'
#' aaf <- vcfapply2(vcfinfo = vcfinfo,
#'                  func = BinaryDosage:::getaaf)
vcfapply2 <- function(vcfinfo, func, snps, blksz, ...) {
  if (missing(vcfinfo) == TRUE)
    stop("No vcf file information specified")
  if (is.na(match("genetic-info", class(vcfinfo))) == TRUE)
    stop("vcfinfo does not appear to contain information about a vcf file")
  if (is.na(match("vcf-info", class(vcfinfo$additionalinfo))) == TRUE)
    stop("vcfinfo does not appear to contain information about a vcf file")
  
  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")
  
  if (missing(snps) == TRUE | length(snps) == 0)
    snps <- 1:nrow(vcfinfo$snps)
  if (is.numeric(snps) == FALSE)
    stop("snps must be a numeric vector")
  if (all(as.integer(snps) == snps) == FALSE)
    stop ("snps must be an integer vector")
  snps <- as.integer(snps)
  if (min(snps) < 0 | max(snps) > nrow(vcfinfo$snps))
    stop ("snps contains value out of range")
  snpstouse <- 1:nrow(vcfinfo$snps) %in% snps
  maxsnp <- max(snps)
  
  if (missing(blksz) == TRUE)
    blksz <- 100L
  if (is.numeric(blksz) == FALSE)
    stop("blksz must be a numeric value")
  if (all(as.integer(blksz)) == FALSE)
    stop("blksz must be an integer value")
  blksz <- as.integer(blksz)
  if (length(blksz) != 1)
    stop ("blksz must be a single integer value")
  if (blksz > 1000L)
    blksz <- 1000L
  
  retval <- vector("list", length(snps))
  if (vcfinfo$additionalinfo$gzipped == FALSE)
    con <- file(vcfinfo$filename, "r")
  else
    con <- gzfile(vcfinfo$filename, "r")
  line <- readLines(con, n = vcfinfo$additionalinfo$headerlines)
  
  dosage <- numeric(nrow(vcfinfo$samples))
  p0 <- numeric(nrow(vcfinfo$samples))
  p1 <- numeric(nrow(vcfinfo$samples))
  p2 <- numeric(nrow(vcfinfo$samples))
  
  j <- 0
  firstline <- seq(1, maxsnp, blksz)
  lastline <- firstline + (blksz - 1)
  if (lastline[length(lastline)] > maxsnp)
    lastline[length(lastline)] <- maxsnp
  for (i in 1:length(firstline)) {
    if (i == length(firstline))
      blksz <- lastline[i] - firstline[i] + 1
    if (blksz < 1)
      break;
    line <- readLines(con, n = blksz)
    linenumber <- 0
    for (k in firstline[i]:lastline[i]){
      linenumber <- linenumber + 1
      if (snpstouse[k] == FALSE)
        next
      j <- j + 1
      x <- unlist(strsplit(line[linenumber], "\t"))
      y <- unlist(strsplit(x[10:length(x)], ":"))
      if (length(vcfinfo$additionalinfo$datacolumns$dosage) == 1) {
        dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage
        gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob
        numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns
      } else {
        dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage[i]
        gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob[i]
        numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns[i]
      }
      if(is.na(dosagecol) == FALSE) {
        dosage[1:nrow(vcfinfo$samples)] <- as.numeric(y[seq(dosagecol, length(y) - numcolumns + dosagecol, numcolumns)])
      }
      if(is.na(gpcol) == FALSE) {
        gpstring <- y[seq(gpcol, length(y) - numcolumns + gpcol, numcolumns)]
        z <- unlist(strsplit(gpstring, ","))
        p0[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(1, length(z) - 2, 3)])
        p1[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(2, length(z) - 1, 3)])
        p2[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(3, length(z), 3)])
      } else {
        p0[1:nrow(vcfinfo$samples)] <- NA
        p1[1:nrow(vcfinfo$samples)] <- NA
        p2[1:nrow(vcfinfo$samples)] <- NA
      }
      if(is.na(dosagecol) == FALSE) {
        dosage[1:nrow(vcfinfo$samples)] <- as.numeric(y[seq(dosagecol, length(y) - numcolumns + dosagecol, numcolumns)])
      } else {
        dosage[1:nrow(vcfinfo$samples)] <- p1 + p2 + p2
      }
      retval[[j]] <- func(dosage = dosage,
                          p0 = p0,
                          p1 = p1,
                          p2 = p2,
                          ...)
    }
  }
  close(con)
  return (retval)
}
