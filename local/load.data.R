library(multiseq)
library(rhdf5)

#********************************************************************
#
#     Get arguments
#
#********************************************************************

args            <- commandArgs(TRUE)
samplesheet     <- args[1]
list_loci       <- args[2]
dir.name        <- file.path(args[3])

onlyoneend      <- TRUE #if bam files are paired end only take one of the reads in the pair
                             
samples         <- read.table(samplesheet, stringsAsFactors=F, header=T)   
regions         <- read.table(list_loci, stringsAsFactors=F, header=F)

for (i in 1:nrow(regions)){
    chr         <- regions[i,1]
    locus.start <- regions[i,2]
    locus.end    <- regions[i,3]
    locus.name   <- paste(chr, locus.start, locus.end, sep=".")

    #create output directories
    dir.create(dir.name)
    dir.name     <- file.path(dir.name, locus.name)
    dir.create(dir.name)

    #timing
    ptm <- proc.time()

    #get counts
    x <- get.counts(samplesheet, region, onlyoneend=TRUE)
    save(x, file=file.path(dir.name, "data.RData"))

    #timing
    t <- proc.time() - ptm
    write(t[[1]], file.path(dir.name,"time_loading.txt"))
}
