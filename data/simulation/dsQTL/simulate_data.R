source("my.utils.R")
##' `simulate.data' simulate data sets by resampling reads from real data (real.read.counts) for given signals (sig0, sig1).
##' 
##' @param seed seed number to set up
##' @param numGroup0 number of samples for Group0
##' @param numGroup1 number of samples for Group1
##' @param real.read.counts a vector of size T (e.g., 1024); t-th element contains number of reads at t-th position; from which we will resample reads for simulation;
##' @param sig0 a vector of size T (e.g., 1024); t-th element contains probability of sampling read at t-th position in real.read.counts for Group0
##' @param sig1 a vector of size T (e.g., 1024); t-th element contains probability of sampling read at t-th position in real.read.counts for Group1
##' @param real.library.read.depth a vector of size M (>1); from whch we will sample library read depth
##' @param over.dispersion parameter used in sample.from.Binomial.with.Overdispersion; see `my.utils.R' for details. 
##' @return a list of data, group, library.read.depth; data contains simulated data; matrix of (numGroup0+numGroup1) by T; group contains group indicator for simulated data; a vector of size (numGroup0+numGroup1); library.read.depth contains simulated library read depth; a vector of size (numGroup0+numGroup1)
simulate.data <- function(seed = 1, numGroup0, numGroup1, sig0, sig1, real.read.counts, real.library.read.depth = NULL, over.dispersion=NULL){
  
  genoD = c(rep(0, numGroup0), rep(1, numGroup1))
  numSam = length(genoD)
  numBPs = length(sig0)
  
  ## phenotype data
  phenoD = matrix(data=NA, nr= length(genoD), nc = numBPs)
  
  ## let's sample!!!
  set.seed(seed)
  
  ## upper and lower bound!
  trunc.fun = function(x){
    x = max(0, x)
    return(min(1,x))
  }
  p.sig0 = sapply(sig0, trunc.fun)
  p.sig1 = sapply(sig1, trunc.fun)
  
  ## geno = 0
  wh0 = which(genoD == 0)
  if(length(wh0) > 0){
    phenoD[wh0,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh0), total.count = real.read.counts, mu.sig = p.sig0, over.dispersion = over.dispersion)
  }  
  ## geno = 1
  wh1 = which(genoD == 1)
  if(length(wh1) > 0){
    phenoD[wh1,] = sample.from.Binomial.with.Overdispersion(num.sam = length(wh1), total.count = real.read.counts, mu.sig = p.sig1, over.dispersion = over.dispersion)
  }
  if(is.null(real.library.read.depth)){
    library.read.depth=NULL
  }else{
    library.read.depth = sample(real.library.read.depth, numSam, replace = TRUE)
  }
  return(list(data = phenoD, group = genoD, library.read.depth = library.read.depth))
}
