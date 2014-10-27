library(multiseq)
library(rhdf5)

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
hub.name <- NULL
onlyoneend=TRUE #if bam files are paired end only take one of the reads in the pair

prior           <- "nullbiased" 
lm.approx       <- FALSE #=TRUE #so far we run things with approx=FALSE
                             
samples         <- read.table(samplesheet, stringsAsFactors=F, header=T)   

#get g
g               <- factor(samples$Type)
g               <- match(g, levels(g)) - 1  

#get region
split_region <- unlist(strsplit(region, "\\:|\\-"))
chr          <- split_region[1]
locus.start  <- as.numeric(split_region[2])
locus.end    <- as.numeric(split_region[3])
locus.name   <- paste(chr, locus.start, locus.end, sep=".")

#create output directories
dir.create(dir.name)
dir.name     <- file.path(dir.name, locus.name)
dir.create(dir.name)

#get counts
M <- get.counts(samplesheet, region, onlyoneend=TRUE)
if (sum(M)<10){
    stop("Total number of reads over all samples is <10. Stopping.")
}

#set fitted g
set.fitted.g=NULL
set.fitted.g.intercept=NULL
if (fitted.g.file!="NA"){
    load(fitted.g.file)
    set.fitted.g=ret$fitted.g
    set.fitted.g.intercept=ret$fitted.g.intercept
}

#run multiseq
res        <- multiseq(M,
                       g=g,
                       minobs=1,
                       lm.approx=lm.approx,
                       read.depth=samples$ReadDepth,
                       prior=prior,
                       set.fitted.g=set.fitted.g,
                       set.fitted.g.intercept=set.fitted.g.intercept)
warnings()

#write output
res$region <- region
intervals  <- get.intervals(res, z.threshold=2)
write.bed(intervals, file.path(dir.name, "multiseq.bed"))
write.gz(res, file.path(dir.name, "multiseq.gz"))
