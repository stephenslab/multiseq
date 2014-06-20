args <- commandArgs(TRUE)
samplesheet <- args[1]
regions.file <- args[2]
dir.name <- file.path(args[3]) #output directory

#samplesheet <- "~/src/stat45800/tests/simulations/Encode/data/high/chr1.1179649.1310720/NonNullSamplesheet.txt"
#regions.file <- "tmp2"
#dir.name="tmp_dir"


library(multiseq)
prior="nullbiased"
samples <- read.table(samplesheet, stringsAsFactors=F, header=T)
g <- factor(samples$Tissue)
g <- match(g, levels(g)) - 1
regions <- as.matrix(read.table(regions.file, stringsAsFactors=F, header=F))
J=log2(as.numeric(regions[1,3])-as.numeric(regions[1,2]))
#check if J is integer
if (!(J%%1==0))
        stop(paste("Error: length of region in file", regions.file, "is not a power of 2"))
#check if all regions have the same length
if (sum(as.numeric(regions[,3])-as.numeric(regions[,2])!=2^J)!=0)
        stop(paste("Error: regions in file ", regions.file, "have different length"))
regions <- paste0(regions[,1],":",as.numeric(regions[,2])+1,"-",as.numeric(regions[,3]))
#save output in foleder fitted_g
dir.create(dir.name)
dir.name <- file.path(dir.name, "fitted_g")
dir.create(dir.name)

nregions=length(regions)

start=1
y.list=list()
y.rate=NULL
sample.size=sapply(1:J,function(x){ceiling(2^(J-x)/nregions)})
intervals=c(0,cumsum(sapply(J:1,function(x){2^(x-1)})))

for (i in 1:nregions){
        #M <- get.counts(samples, regions[i])
        #res <- multiseq(M, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth, smoothing=FALSE, cyclespin=FALSE, get.fitted.g=FALSE, prior=prior)
    y.rate=cbind(y.rate,res$y.rate)
    for (j in 1:J){
        sampled=sample((intervals[j]+1):intervals[j+1],size=sample.size[j],replace=FALSE)
        ind=as.vector(rbind(2*sampled-1,2*sampled))
        if (i==1)
            y.list[[j]]=res$y[,ind]
        else
            y.list[[j]]=cbind(y.list[[j]],res$y[,ind])
    }
}
intervals=c(0,cumsum(sample.size*nregions))
y=NULL
for (j in 1:J)
    y=cbind(y,y.list[[j]])
ret <- multiseq(listy=list(y=y, y.rate=y.rate, intervals=intervals), g=g, cyclespin=FALSE, get.fitted.g=TRUE) #specify
#qsub to each locus
reg <- multiseq(M, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth, prior=prior, set.fitted.g=ret$fitted.g, set.fitted.g.intercept=ret$fitted.g.intercept)
