## This file contains useful functions for analysis.
## Sometime I copied original functions and modified them. Then, I specified original sources. 
##
## Copyright (C) 2014 Heejung Shim
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.



##' 'compute.proportions.given.effect.size' takes effect sizes over a retion (as a ratio)
##' and sum of proportions (in binomial sampling) for two groups as input. And it computes
##' proportion parameters in binomial sampling for two groups for the given effect sizes
##' (as a ratio). Specifically, let $p_1 = p$ and $p_2 = p \times effect.ratio$
##' where $\log(effect.ratio)$ is effec.size in log space and $effect.ratio = 1$ means
##' there is no effect. Use $p$ which satisfies $p_1 + p_2 = sum.prop$.
##' So $p_1 = p = sum.prop\frac{1}{1+effect.ratio}$ and $p_2 = sum.propr\frac{effect.ratio}{1+effect.ratio}$.
##'
##'
##' for example,
##' sum.prop = rep(2/70, 1024)
##' effect.ratio = rep(1, 1024)
##' effect.ratio[100:200] = 1/2
##' res = compute.proportions.given.effect.size(sum.prop = sum.prop, effect.ratio = effect.ratio)
##' res$prop0
##' res$prop1
##' 
##' @param sum.prop a vector of sum of proportions (in binomial sampling) for two groups over a region
##' @param effect.ratio a vector of effect sizes over a region as a ratio; default value = NULL; if effect.ratio == NULL, we assign the same proportion to two groups (no effect).
##' @return a list of prop0 and prop1; prop0 (prop1) is a vector of proportions (in binomial sampling) for group0 (group1) over a region 
compute.proportions.given.effect.size <- function(sum.prop, effect.ratio=NULL) {
  
  if(is.null(effect.ratio)){
    effect.ratio = rep(1, length(sum.prop))
  }
  prop0 = sum.prop/(effect.ratio + 1)
  prop1 = sum.prop/(effect.ratio + 1)*effect.ratio
  
  return(params = list(prop0 = prop0, prop1 = prop1))
}



##' 'estBetaParams' takes mean and variance for beta distribution and returns
##' alpha and beta parameters in beta distribution.
##'
##'
##' I copied the original version from http://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
##' Then I modified to handle cases when somme elements of mu are 0 or 1 (then it returns NA if mu <= 0 or >=1) or when computed alpha and beta <= 0 (then it returns NA)
##'
##' @param mu.orig a vector of mean for beta distribution
##' @param var a vector (or scalar) of variance for beta distribution
##' @return a list of alpha and beta 
estBetaParams <- function(mu.orig, var) {
  
  del.ix = ((mu.orig <= 0) | (mu.orig >= 1))
  if(sum(del.ix) > 0){
    mu = mu.orig[!del.ix]
  }else{
    mu = mu.orig
  }
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  
  if(sum(del.ix) > 0){
    alpha.orig = beta.orig = rep(NA, length(mu.orig))
    alpha.orig[!del.ix] = alpha
    beta.orig[!del.ix] = beta
  }else{
    alpha.orig = alpha
    beta.orig = beta
  }
  
  ## handle non-positive alpha or beta 
  invalid.para = which((alpha.orig <= 0) | (beta.orig <= 0))
  if(length(invalid.para) > 0){
    alpha.orig[invalid.para] = rep(NA, length(invalid.para))
    beta.orig[invalid.para] = rep(NA, length(invalid.para))
  }
  
  return(params = list(alpha = alpha.orig, beta = beta.orig))
}





##' 'sample.from.Binomial.with.Overdispersion' simulates binomial samples with/without
##' over dispersion. 
##'
##' For a given overdispersion parameter, computed parameters for beta distribution can be invalid (e.g., mu.sig are 0 or 1). Then we sample read without overdispersion for those positions.
##' 
##' @param num.sam number of samples to be sampled
##' @param total.count a vector of non-negative counts;
##' @param mu.sig a vector of probabilities (we allow 0 or 1 as probablity)
##' @param over.dispersion if over.dispersion == NULL, simulate data from binomial. If over.dispersion is provided, simulate data binomial with over.dispersion.
##' @return a matrix of num.sam by L (length of total.count) containing simulated data. 
sample.from.Binomial.with.Overdispersion <- function(num.sam, total.count, mu.sig, over.dispersion=NULL){
  
  invalid.entry = ((mu.sig < 0) | (mu.sig > 1))
  if(sum(invalid.entry) > 0){ stop("ERROR, mu.sig have some values outside of valid range [0, 1]")}
  
  if(is.null(over.dispersion)){
    return(matrix(data=rbinom(length(mu.sig)*num.sam, total.count, mu.sig), nr = num.sam, byrow = TRUE))
  }else{
    
    
    final.dat = matrix(data=NA, nr = num.sam, nc = length(mu.sig))
    
    # get alpha and beta
    resBeta = estBetaParams(mu.sig, over.dispersion)
    alpha = resBeta$alpha
    beta = resBeta$beta
    
    # for valid alpha and beta, sample data 
    del.ix = is.na(alpha)
    p.sig = rbeta(sum(!del.ix)*num.sam, alpha[!del.ix], beta[!del.ix])
    dat.V = rbinom(sum(!del.ix)*num.sam, total.count[!del.ix], p.sig) 
    final.dat[,!del.ix] = matrix(data=dat.V, nr = num.sam, byrow = TRUE)
    
    # for invalid alpha and beta, sample without over dispersion
    if(sum(del.ix) > 0){
      dat.IV = matrix(data=rbinom(sum(del.ix)*num.sam, total.count[del.ix], mu.sig[del.ix]), nr = num.sam, byrow = TRUE)
      final.dat[,del.ix] = matrix(data=dat.IV, nr = num.sam, byrow = TRUE)
    }
    
    return(final.dat)
    
  }
  
}






##' 'get.pcr.artifacts.posi' checks if any sample has pcr artifacts in a given position and returns index of samples (a vector) that are identified to have pcr artifacts in the given position. 
##'
##' Here, we consider a window of left and right `win.half.size' bps around a given position. If read count at the given positon is more than `prop.thresh' of reads in the window, we consider read count at the given position is due to pcr artifacts.
##' 
##' 
##' @param posi a candidate position for pcr artifacts in data 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param win.half.size default=50; a half of window size
##' @param prop.thresh default=0.9;  
##' @return pcr.sample a vector of samples that are identified to have pcr artifacts in a given position. 
get.pcr.artifacts.posi <- function(posi, data, win.half.size = 50, prop.thresh = 0.9){
  
  st.win = max(1, posi - win.half.size)        
  en.win = st.win + win.half.size*2
  if(en.win > dim(data)[2]){
    en.win = dim(data)[2]
    st.win = en.win - win.half.size*2
  }
  prop = data[,posi]/apply(data[,st.win:en.win], 1, sum)
  pcr.sample = which(prop > prop.thresh)
  return(pcr.sample)
}



##' 'remove.pcr.artifacts' identifies positions in data where at least one sample has pcr artifacts (using the function `get.pcr.artifacts.posi') and replace read counts of all samples at that position with max(1, average of read counts of sample without pcr artifacts).
##'
##'
##' 
##' 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param win.half.size default=50; argument to get.pcr.artifacts.posi
##' @param prop.thresh default=0.9; argument to get.pcr.artifacts.posi
##' @return a list of data and posi.with.pcr.artifacts; data contains pcr artifacts-removed data; posi.with.pcr.artifacts contains a list of positions where at least one sample has pcr.artifacts. 
remove.pcr.artifacts <- function(data, win.half.size = 50, prop.thresh = 0.9){
  
  # only consider positions with at least 2 reads as a candidate position 
  max.val = apply(data, 2, max)
  candidate.posi = which(max.val > 1)
  
  num.sam = dim(data)[1]
  len = length(candidate.posi)
  if(len > 0){
    # for each candidate position, get which samples have pcr artifacts
    pcr.posi.list = lapply(candidate.posi, get.pcr.artifacts.posi, data = data, win.half.size = win.half.size, prop.thresh = prop.thresh)
    pcr.posi = which(sapply(pcr.posi.list, length) > 0)
    len.pcr.posi = length(pcr.posi)
    if(len.pcr.posi > 0){
      # for positions where at least one sample has pcr artifacts, replace read counts of all samples with max(1, average of read counts of sample without pcr artifacts)
      for(p in 1:len.pcr.posi){
        pcr.sam = pcr.posi.list[[pcr.posi[p]]]
        ix = candidate.posi[pcr.posi[p]]
        if(length(pcr.sam) == num.sam){
          data[,ix] = 1
        }else{
          data[, ix] = max(1, ceiling(mean(data[-pcr.sam,ix])))
        }
      }
    }
  }else{
    len.pcr.posi=0
  }
  
  if(len.pcr.posi > 0){
    return(list(data=data, posi.with.pcr.artifacts = candidate.posi[pcr.posi]))
  }else{
    return(list(data=data, posi.with.pcr.artifacts = NULL))
  }
}


##' 'remove.pcr.artifacts.in.known.posi' uses pre-identified positions where at least one sample has pcr artifacts and replace read counts of all samples at that position with max(1, average of read counts of sample without pcr artifacts).
##'
##'
##' 
##' 
##' @param data a matrix of num of samples by site size; original count data of all samples in a given site
##' @param known.pcr.posi a vector of positions 
##' @param win.half.size default=50; argument to get.pcr.artifacts.posi
##' @param prop.thresh default=0.9; argument to get.pcr.artifacts.posi
##' @return a list of data and posi.with.pcr.artifacts; data contains pcr artifacts-removed data; posi.with.pcr.artifacts contains a list of positions where at least one sample has pcr.artifacts. 
remove.pcr.artifacts.in.known.posi <- function(data, known.pcr.posi, win.half.size = 50, prop.thresh = 0.9){
  
  if(length(known.pcr.posi) > 0){
    num.sam = dim(data)[2]
    
    pcr.posi.list = lapply(known.pcr.posi, get.pcr.artifacts.posi, data = data, win.half.size = win.half.size, prop.thresh = prop.thresh)
    
    len.pcr.posi = length(pcr.posi.list)
    for(p in 1:len.pcr.posi){
      pcr.sam = pcr.posi.list[[p]]
      ix = known.pcr.posi[p]
      if(length(pcr.sam) == num.sam){
        data[,ix] = 1
      }else{
        data[, ix] = max(1, ceiling(mean(data[-pcr.sam,ix])))
      }
    }
  }
  return(list(data=data))
}




##' 'get.pval.from.empirical.null.dist.discrete' takes empirical distribution of statistic
##' under the null (statistic.null) and a series of observed test statistics (statistic.alt),
##' and returns a series of p-values corresponding to the observed test statistics.
##' When test statistic is not completely continuous, there are problems for the qvalue
##' package which uses the insignificant p values to estimate the proportion of nulls.
##' Here, this function uses randomization to produce continuous p-values that are uniformly
##' distributed under the null.
##'
##' Let t be observed test statistic. p-value is defined by P(T >= t | H0), but this is not
##' uniformly distributed unless test statistic is continuous. We randomize them using
##' P(T > t | H0) + U*P(T = t | H0), where U ~ uniform(0,1). The proposed p-value is
##' between P(T > t | H0) and P(T > t | H0) + P(T = t | H0). This p-value can be empirically
##' computed by #[T > t] + U*(#[T = t] + 1) / (total number of statistic under the null + 1).
##' Here "+1" in the denominator is due to observed t. 
##'
##' @param statistic.null a vector of test statistics under the null
##' @param statistic.alt a vector of observed test statistics
##' @param big.sig bool indicating whether bigger statistic is more significant.
##' @return a vector of p-values corresponding to the observed test statistics.
get.pval.from.empirical.null.dist.discrete <- function(statistic.null, statistic.alt, big.sig = TRUE){
  
  numNulltests = length(statistic.null)
  if(big.sig){
    numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null > x))}, statistic.null = statistic.null)
  }else{
    numSig = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null < x))}, statistic.null = statistic.null)
  }
  numEqual = sapply(statistic.alt, function(x, statistic.null){ return(sum(statistic.null == x))}, statistic.null = statistic.null)
  Uval = runif(length(statistic.alt))
  
  pval.list = (numSig + Uval*(numEqual + 1))/(numNulltests + 1)
  return(pval.list)
}
