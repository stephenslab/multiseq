#usage
#Rscript run_multiseq.R $chr $locus.start $locus.end $samplesheet $multiseq.output.folder 

#library(rhdf5)
#library(ashr)
library(multiseq)

#region,samples$Tissue,samples$

#course.repodir    <- scan(".course.repodir.txt", what=character())
#ash.repodir       <- scan(".ash.repodir.txt", what=character())
#RNAsplice.repodir <- scan(".RNAsplice.repodir.txt", what= character())
#source(file.path(course.repodir, "code", "Rcode", "utils.R"))

#********************************************************************
#
#     Get arguments
#
#********************************************************************
args            <- commandArgs(TRUE)
chr             <- args[1]
locus.start     <- as.numeric(args[2])+1
locus.end       <- as.numeric(args[3])
samplesheet     <- args[4]
dir.name        <- file.path(args[5])
hub.name        <- args[6]
chrom.file      <- file.path(args[7])

samples         <- read.table(samplesheet, stringsAsFactors=F, header=T)
fra             <- 2  #how many sd to plot when plotting effect size
#    merged.bed         <-
do.plot         <- FALSE
do.smooth       <- FALSE

#' @param dir.name: specify the directory where multiseq results will be saved. File effect_mean_var.txt.gz
#' is written, which has two columns: effect mean and effect variance. This parameter defaults to NULL in which case
#' results are returned but not saved in any directory
run.multiseq(samples, chr, locus.start, locus.end, dir.name)
run.multiseq <- function(samples, chr, start, end, dir.name=NULL){
    g               <- factor(samples$Tissue)
    g               <- match(g, levels(g)) - 1

    locus.name <- paste(chr, locus.start, locus.end, sep=".")
    dir.create(dir.name)
    dir.name   <- file.path(dir.name, locus.name)
    dir.create(dir.name)
    
    M <- get.counts(samples, chr, start, end)

    if (sum(M)<10){
        write.table(t(c(chr, locus.start, locus.end, NA, NA)),
                    quote = FALSE,
                    col.names=FALSE,
                    row.names=FALSE,
                    file=file.path(dir.name,"summary.txt"))
        stop("Total number of reads over all samples is <10. Stopping.")
    }

    print("Compute effect")
    #ptm      <- proc.time()
    res <- multiseq(M, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
    #my.time  <- proc.time() - ptm

    res$chr=chr
    res$locus.start=locus.start
    res$locus.end=locus.end
    res <- get.effect.intervals.bed(res,fra)
    
    if (!is.null(dir.name)){
        #save results in a compressed file 
        write.effect.mean.variance.gz(res,dir.name)
        write.effect.intervals.bed(res,dir.name,fra=2)
    }
}

#plotting
if (do.plot==TRUE){
    print("Plot Effect")
    #source(file.path(RNAsplice.repodir, "code", "R", "utils.R"))
                                        # Extract gene model from genePred annotation file
    annotation  <- "/mnt/lustre/home/epantaleo/data/annotations/hg19.ensGene.gp.gz"
    Transcripts <- get.Transcripts(annotation, chr, locus.start, locus.end)
                                        # Plot

                                        #centipede <- read.table("/mnt/lustre/data/share/DNaseQTLs/CentipedeAndPwmVar/CentipedeOverlapMaxP99.UCSC
                                        #sites     <- centipede[centipede[,1] == chr & centipede[,2] < (locus.end+1) & centipede[,3] > locus.start, ]
                                        #offset    <- -0.25 #0.0025

    pdf(file=file.path(dir.name, "Effect.pdf"))
    layout(matrix(c(1,2), 2, 1, byrow = FALSE))
                                        #if (dim(sites)[1] > 0){for (k in 1:dim(sites)[1]){offset <-  -offset
                                        #    text(x=(sites[k,2] + sites[k,3])/2, y=(ymax/2 -abs(offset) - offset), strsplit(as.character(sites[k,4]), split="=")[[1]][2])
                                        #    rect(sites[k,2], 0, sites[k,3], ymax/2 + 1, col=rgb(0,0,1,0.3), border='NA')}}
    par(mar=c(0,3,2,1))
    plotEffect(res$effect.mean, res$effect.var, title="Mean Effect")
    par(mar=c(4,3,0,1))
    plotTranscripts(Transcripts, locus.start, locus.end)
    dev.off()
}

#************************************************************
#
#    Apply cyclespin to each genotype class
#
#*************************************************************
if (do.smooth==TRUE){
    print("Smooth by group")
    res0 <- multiseq(M[g==0,], minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth[g==0])
    res1 <- multiseq(M[g==1,], minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth[g==1])
    
    if (do.plot==TRUE){
        print("Plot smoothed signals")
        pdf(file=file.path(dir.name,"Smoothing_by_group.pdf"))
        ymax <- max(c(as.vector(res0$est + fra*sqrt(res0$var)), as.vector(res1$est + fra*sqrt(res1$var))))
        ymin <- min(c(as.vector(res0$est - fra*sqrt(res0$var)), as.vector(res1$est - fra*sqrt(res1$var))))
        ylim <- c(ymin, ymax)
        layout(matrix(c(1,2,3,4), 4, 1, byrow = FALSE))
        par(mar=c(0,3,2,1))
        plotMean(res0$baseline.mean, res0$baseline.var, title=paste0("Smoothed ", samples$Tissue[g==0][1]), ylim=ylim)
        par(mar=c(0,3,2,1))
        plotMean(res1$baseline.mean, res1$baseline.var, title=paste0(samples$Tissue[g==1][1]), ylim=ylim)
        par(mar=c(0,3,2,1))
        plotTranscripts(Transcripts, locus.start, locus.end)
        dev.off()
    }
}

if (!is.null(hub.name)){
    source("~/src/NGS_utils/R/TrackHub.R")
    multiseq.folder=file.path(args[5])
    multiseqToTrackHub(paste0(chr, ":", locus.start, "-", locus.end), hub.name, multiseq.folder, chrom.file)
}
