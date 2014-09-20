#GIVEN A LIST OF SITES WHERE WE WANT TO RUN MULTISEQ (IN BED FORMAT) SELECT 1000 SITES (IF THE BED FILE CONTAINS LESS THAN 1000 SITES THEN TAKE ALL OF THEM). FOR EACH j IN 1:J select 2^(J-j)/1000  
66 33 17  9  5  3  2  1  1  1  1  1  1  1  1  1  1


#ql 20g
#list_loci="list_loci_multiples"
#tmp_list_loci="tmp"
#sh write_list_loci.sh 131072 ../data/chromosome.lengths.hg19.txt > $tmp_list_loci 
#cat $tmp | awk '{if ($3-$2==131072) print}' > $list_loci
#rm $tmp #!!!
#sample_sheet="~/src/stat45800/data/Encode/samplesheetEncode.txt"
#OUT_DIR="./results_Encode"
#Rscript pipeline.R $sample_sheet $OUT_DIR $list_loci 
#sh qsub_run_multiseq.sh $samplesheet $OUT_DIR < $list_loci 
                                        #Rscript ~/src/stat45800/tests/simulations/Encode/data/high/chr1.1179649.1310720/NonNullSamplesheet.txt tmp_dir tmp2
args           <- commandArgs(TRUE)
samplesheet    <- args[1]
dir.name       <- file.path(args[2]) #output directory 
regions.file   <- args[3]

#samplesheet  <- "~/src/stat45800/tests/simulations/Encode/data/high/chr1.1179649.1310720/NonNullSamplesheet.txt"
#regions.file <- "tmp2"
#dir.name     <- "tmp_dir"


library(multiseq)
library(rhdf5)
prior          <- "nullbiased"
lm.approx      <- FALSE
samples        <- read.table(samplesheet, stringsAsFactors=F, header=T)
g              <- factor(samples$Tissue)
g              <- match(g, levels(g)) - 1
nsig           <- length(g)
regions        <- as.matrix(read.table(regions.file, stringsAsFactors=F, header=F))
J              <- log2(as.numeric(regions[1,3])-as.numeric(regions[1,2]))
#check if J is integer
if (!(J%%1==0))
    stop(paste("Error: length of region in file", regions.file, "is not a power of 2"))
#check if all regions have the same length
if (sum(as.numeric(regions[,3])-as.numeric(regions[,2])!=2^J)!=0)
    stop(paste("Error: regions in file ", regions.file, "have different length"))
regions        <- paste0(regions[,1],":",as.numeric(regions[,2])+1,"-",as.numeric(regions[,3]))
lregions       <- length(regions)
nregions       <- min(1000,lregions)
#save output in folder fitted_g
dir.create(dir.name)

#for each scale choose how many observations to sample
sample.size    <- sapply(1:J,function(x){ceiling(2^(J-x)/nregions)})
intervals      <- c(0,cumsum(sapply(J:1,function(x){2^(x-1)})))
y.intervals    <- c(0,cumsum(sample.size*nregions))

y              <- matrix()
length(y)      <- nsig*2*tail(y.intervals,1)
dim(y)         <- c(nsig,2*tail(y.intervals,1))
y.rate         <- matrix()
length(y.rate) <- nsig*2*nregions
dim(y.rate)    <- c(nsig,2*nregions)
i              <- 0
for (counter in 1:lregions){
    M          <- get.counts(samples, regions[counter])
    if (sum(M)>nsig*J/100){
        i <- i+1
        res        <- multiseq(M, g=g, minobs=1, read.depth=samples$ReadDepth, smoothing=FALSE, cyclespin=FALSE, get.fitted.g=FALSE, prior=prior)
        y.rate[1:nsig,(2*i-1):(2*i)] <- res$y.rate
        for (j in 1:J){
                                        #for each scale choose which observations to sample
            sampled <- sample((intervals[j]+1):intervals[j+1],size=sample.size[j],replace=FALSE)
            y[1:nsig,2*y.intervals[j]+1+((i-1)*2*sample.size[j]):(i*2*sample.size[j]-1)] <- res$y[,as.vector(rbind(2*sampled-1,2*sampled))]
        }
        rm(res)
        if (i==nregions)
            break
    }
    rm(M)    
}
if (i<nregions)
    stop(paste("error: not enough data (less than", nregions, "loci); provide a different list of loci"))
ret   <- multiseq(listy=list(y=y, y.rate=y.rate, intervals=y.intervals), g=g, cyclespin=FALSE, get.fitted.g=TRUE, lm.approx=lm.approx) 
save.image("Encode.RData")
save(ret, file=file.path(dir.name,"fitted.g.RData"))
#qsub to each locus
#load(file=file.path(dir.name,"fitted.g.RData")) 
#reg   <- multiseq(M, g=g, minobs=1, lm.approx=lm.approx, read.depth=samples$ReadDepth, prior=prior, set.fitted.g=ret$fitted.g, set.fitted.g.intercept=ret$fitted.g.intercept)








