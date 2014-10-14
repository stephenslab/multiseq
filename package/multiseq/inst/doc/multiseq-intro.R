## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center")

## ----load_data-----------------------------------------------------------
    library(multiseq)
    #load example data - R object
    #type ?example1 to get more information
    data(example1, package="multiseq")   
    x          <- dat$x
    g          <- dat$g
    read.depth <- dat$read.depth

## ----smoothing, echo=c(1,3,4,5,6,8)--------------------------------------
    #smoothing
    par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    for (i in which(g==1)) plot(log(x[i,]), type="l", xlab="Position", ylab="log(x)")
    unique(sort(log(x[1,])))
    res0         <- multiseq(x=x[g==1,], minobs=1, lm.approx=FALSE, read.depth=read.depth[g==1])
    #plot baseline mean +/- 2 posterior standard deviations
    invisible(dev.off())
    plot(res0, fra=2, what="baseline")

## ----estimate_effect-----------------------------------------------------
    #estimating an effect
    res           <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=read.depth)
    #plot estimated effect mean +/- 2 posterior standard deviations
    plot(res, fra=2)

    #print intervals where `multiseq` found a strong effect (zero is outside of +/- fra posterior standard deviations
    res$intervals <- get.effect.intervals(res, fra=2)
    res$intervals

## ----load_seq_data-------------------------------------------------------
    setwd(file.path(path.package("multiseq"),"extdata","sim"));
    samplesheet <- file.path(path.package("multiseq"),"extdata","sim","samplesheet.sim.txt")
    samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
    g <- factor(samples$Type)
    g <- match(g, levels(g))-1
    if (noExecutable("wigToBigWig")){
       data(example2, package="multiseq")
    }else{
       region      <- "chr1:154206209-154214400"
       x           <- get.counts(samples, region)
    }

## ----testing_smoothing, echo=c(2,3,4,6,7)--------------------------------
    par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    for (i in which(g==0)) plot(log(x[i,]), type="l", xlab="Position", ylab="log(x)")
    res0 <- multiseq(x=x[g==0,], minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth[g==0])
    #plot estimated log baseline +- 2 s.d.
    invisible(dev.off())
    res0$region <- region
    plot(res0, fra=2, what="baseline")

## ----testing_diff--------------------------------------------------------
    res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
    #plot estimated effect and s.d. 
    par(mfrow=c(2,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,1,1) + 0.1)
    res$region <- region
    plot(res, fra=2, is.xaxis=FALSE) 
    transcripts   <- get.transcripts(file.path(path.package("multiseq"),"extdata","sim","hg19.ensGene.part.gp"), region)
    plot(transcripts, region) 

## ----track_hub, results='hide'-------------------------------------------
    setwd(file.path(path.package("multiseq"),"extdata","sim"))
    hub_name <- "testMultiseq/sim"
    samplesheetToTrackHub(samplesheet, hub_name, chr="chr1")

## ----multiseqToTrackHub, results='hide'----------------------------------
    res$region    <- region
    res$intervals <- get.effect.intervals(res, fra=2)
    multiseqToTrackHub(res, fra=2, hub_name="testMultiseq/multiseq_sim")

