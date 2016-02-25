setwd("~/multiseq/data/simulation/dsQTL/")
source("simulate_data.R")

seed = 1
numGroup0 = 10
numGroup1 = 10
tolerance = 1e-10
sig0 = scan("alt.sig0", what=double())
sig1 = scan("alt.sig1", what=double())
real.DNase.dat = read.table("pheno.dat", as.is = TRUE)
real.read.counts = ceiling(as.numeric(apply(real.DNase.dat, 2, sum)))
real.library.read.depth = scan("library.read.depth.dat", what=double())

res = list()

## data with sample size 20 with library read depth 
res[[1]] = simulate.data(seed = seed, numGroup0 = numGroup0, numGroup1 = numGroup1, sig0 = sig0, sig1= sig1, real.read.counts = real.read.counts, real.library.read.depth = real.library.read.depth, over.dispersion=NULL)

## data with sample size 20 with library read depth, but all samples in group 2 have zero read count.
res[[2]] = simulate.data(seed = seed, numGroup0 = numGroup0, numGroup1 = numGroup1, sig0 = sig0, sig1= rep(0, length(sig1)), real.read.counts = real.read.counts, real.library.read.depth = real.library.read.depth, over.dispersion=NULL)

## data with sample size 20 without library read depth
res[[3]] = simulate.data(seed = seed, numGroup0 = numGroup0, numGroup1 = numGroup1, sig0 = sig0, sig1= sig1, real.read.counts = real.read.counts, over.dispersion=NULL)

## data with sample size 20 without library read depth, but all samples in group 2 have zero read count.
res[[4]] = simulate.data(seed = seed, numGroup0 = numGroup0, numGroup1 = numGroup1, sig0 = sig0, sig1= rep(0, length(sig1)), real.read.counts = real.read.counts, over.dispersion=NULL)

# set.seed(seed)
# multiseq_out_base = list()
# for(i in 1:length(res)){
#   multiseq_out_base[[i]] = multiseq(res[[i]]$data, g = res[[i]]$group, read.depth = res[[i]]$library.read.depth, verbose = TRUE)
# }
# 
# save(multiseq_out_base, file = 'multiseq_out_base.Robj')

set.seed(seed)
multiseq_out = list()
for(i in 1:length(res)){
  multiseq_out[[i]] = multiseq(res[[i]]$data, g = res[[i]]$group, read.depth = res[[i]]$library.read.depth, verbose = TRUE)
}

load("multiseq_out_base.Robj")

all.equal(multiseq_out_base, multiseq_out, tolerance = tolerance)
