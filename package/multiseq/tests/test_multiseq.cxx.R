require(rbenchmark)
install.packages(file.path("~/src/ash/package/ashr_0.1.tar.gz"),repos=NULL,type="source")
install.packages(file.path("~/src/multiseq/package/multiseq.tar.gz"),repos=NULL,type="source")
library(multiseq)

#load data from gene OAS1 
data(OAS1,package="multiseq")

M <- OAS1$M
g <- OAS1$g
read.depth <- OAS1$read.depth

#run multiseq with c++ code (cxx=TRUE) or without c++ code (cxx=FALSE)
#using profiler
Rprof("multiseq.cxx.out")
res.cxx <- multiseq(M, g=g, minobs=3, read.depth=read.depth, cxx=TRUE, forcebin=FALSE, prior="nullbiased", lm.approx=TRUE)
Rprof(NULL)


Rprof("multiseq.out")
res <- multiseq(M, g=g, minobs=3, read.depth=read.depth, cxx=FALSE, forcebin=FALSE, prior="nullbiased", lm.approx=TRUE)
Rprof(NULL)

all.equal(res, res.cxx)

#ben <- benchmark(multiseq(M, g=g, minobs=3, read.depth=read.depth, cxx=TRUE, forcebin=FALSE, prior="nullbiased", lm.approx=TRUE), multiseq(M, g=g, minobs=3, read.depth=read.depth, cxx=FALSE, forcebin=FALSE, prior="nullbiased", lm.approx=TRUE))

#to see output run this in the terminal
#R CMD Rprof multiseq.out > multiseq.txt
#R CMD Rprof multiseqcxx.out > multiseqcxx.txt
