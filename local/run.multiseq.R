library(multiseq)

#********************************************************************
#
#     Get arguments
#
#********************************************************************

args            <- commandArgs(TRUE)
samplesheet     <- args[1]
region          <- args[2]
dir.name        <- file.path(args[3])
fitted.g.file   <- file.path(args[4])


#set fitted g
set.fitted.g           -> NULL
set.fitted.g.intercept -> NULL
if (fitted.g.file!="NA"){
    load(fitted.g.file)
    set.fitted.g=ret$fitted.g
    set.fitted.g.intercept=ret$fitted.g.intercept
}

prior           <- "nullbiased" 
lm.approx       <- FALSE #=TRUE #so far we run things with approx=FALSE
                             
samples         <- read.table(samplesheet, stringsAsFactors=F, header=T)   
g               <- factor(samples$Type)
g               <- match(g, levels(g)) - 1  
split_region <- unlist(strsplit(region, "\\:|\\-"))
chr          <- split_region[1]
locus.start  <- as.numeric(split_region[2])
locus.end    <- as.numeric(split_region[3])
locus.name   <- paste(chr, locus.start, locus.end, sep=".")

#create output directories
dir.create(dir.name)
dir.name     <- file.path(dir.name, locus.name)
dir.create(dir.name)

#timing
ptm <- proc.time()

#get counts
load(file=file.path(dir.name, "data.RData"))
if (sum(x)<10){
    stop("Total number of reads over all samples is <10. Stopping.")
}
file.remove(file=file.path(dir.name, "data.RData"))


#run multiseq
ashparam=list(prior=prior)
res        <- multiseq(x,
                       g=g,
                       minobs=1,
                       lm.approx=lm.approx,
                       read.depth=samples$ReadDepth,
                       set.fitted.g=set.fitted.g,
                       set.fitted.g.intercept=set.fitted.g.intercept,
                       ashparam=list(prior=prior))
#timing
t <- proc.time() - ptm
write(t[[1]], file.path(dir.name,"time.txt"))
warnings()

#write output
res$region <- region
intervals  <- get.intervals(res, z.threshold=2)
write.bed(intervals, file.path(dir.name, "multiseq.bed"))
write.gz(res, file.path(dir.name, "multiseq.gz"))
write(res$logLR$value, file=file.path(dir.name,"logLR.txt"))

