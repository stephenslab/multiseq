#' @useDynLib multiseq
#' @title A list containing counts, library size and covariates for 4 samples.
#' @description We generated this data as follows. We downloaded Encode ChipSeq
#' data and ran multiseq on 4 samples belonging to 2 different experiments to get
#' an estimated baseline and effect. Then we generated counts using a negative
#' binomial model using the null (the estimated baseline) and the alternative
#' (the estimated baseline + the estimated effect), using 2 samples for each case.
#' @docType data
#' @keywords datasets
#' @format A list with elements: \code{x} a 4 by 2^14 matrix of counts (4 is the
#' number of samples and 2^14 is the number of adjacent bases where reads have been
#' counted),\code{read.depth} a 4-dimensional vector specifiying the total number
#' of counts per sample, and \code{g} a 4-dimensional vector specifying a
#' binary covariate for each sample.
#' 
#' @name example1
NULL


#' @title A list containing counts, location of the counts for 4 samples, as well
#' as specifications about the samples.
#' @description We generated the data as follows. We downloaded Encode RNASeq
#' data in \code{bam} format. Then we simulated data using a betabinomial model
#' and thinned a region (specified in file \code{/extdata/sim/peaks.bb})
#' by a factor of 2 to get a differential expression pattern.
#' @format A list containing \code{x}, a matrix of counts for 4 samples,
#' \code{region} specifying the genomic location of the counts, and \code{samples}
#' with details about the data. 
#' @docType data
#' @keywords datasets
#' #' 
#' @name example2
NULL


#' Multiseq: analysis of sequence data from multiple samples
#' @docType package
#' @name multiseqr
#NULL
