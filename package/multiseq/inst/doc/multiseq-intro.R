## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center")

## ----smoothing-----------------------------------------------------------
    #First load the package
    library(multiseq)
 
    spikes <- function(x){
	toreturn <- 0.75*exp(-500*(x-0.23)^2) +
		 1.5*exp(-2000*(x-0.33)^2) +
		 3*exp(-8000*(x-0.47)^2) +
		 2.25*exp(-16000*(x-0.69)^2) +
		 0.5*exp(-32000*(x-0.83)^2)
	return(toreturn)
    }

    n    <- 1024   
    t    <- 1:n/n    
    mu   <- 8/3*(3/16+spikes(t)) 

    #use multiseq to smooth a signal
    x    <- rpois(n, mu)
    res  <- multiseq(x)
    
    #use multiseq to smooth 6 signals simultaneously
    x6          <- t(replicate(6, rpois(n,mu)))
    res6        <- multiseq(x6)
   

    z.threshold <- 2
    ylim        <- c(0, max(c(mu, 
    		   	      x, 
                              exp(res$baseline.mean+z.threshold*sqrt(res$baseline.var)), 
                              exp(res6$baseline.mean+z.threshold*sqrt(res6$baseline.var)))))
    
    #plot
    par(mfrow=c(4,1), oma = c(0,0,0,0) + 0.1, mar = c(2,2,1,2), mgp=c(0,1,0))
    
    #plot underlying signal    		      
    plot(mu, type="l", col="red", ylim=ylim, xaxt="n", xlab="", ylab="")
    legend("topright", 
           legend =c("underlying signal", 
	             "simulated data", 
                     "estimated signal 1 sample", 
                     "estimated signal 6 samples"), 
           lty=c(1, 0, 1, 1), 
           pch=c(".","o",".","."), 
           col=c("red", "black", "dark green", "blue"))  
    title("underlying signal", line=-1)    

    #plot underlying signal and simulated data 
    plot(x, ylim=ylim, xaxt="n", xlab="", ylab="")
    lines(mu, col="red")
    title("simulated data",line=-1)

    #plot estimated baseline using 1 sample
    plot(mu, type="l", col="red", ylim=ylim, xaxt="n", xlab="", ylab="")
    lines(exp(res$baseline.mean), col="dark green")  
    title("estimated signal (1 sample)",line=-1) 

    #plot estimated baseline using 6 samples
    plot(mu, type="l", col="red", ylim=ylim, xlab="Position", ylab="")
    lines(exp(res6$baseline.mean), col="blue")
    title("estimated signal (6 samples)",line=-1)

## ----OAS1----------------------------------------------------------------
    data(OAS1, package="multiseq")
    
    res0         <- multiseq(x=OAS1$x[OAS1$g==0,], minobs=1, read.depth=OAS1$read.depth[OAS1$g==0])
    res1         <- multiseq(x=OAS1$x[OAS1$g==1,], minobs=1, read.depth=OAS1$read.depth[OAS1$g==1])
    res2         <- multiseq(x=OAS1$x[OAS1$g==2,], minobs=1, read.depth=OAS1$read.depth[OAS1$g==2])
    res          <- multiseq(x=OAS1$x, g=OAS1$g, minobs=1, read.depth=OAS1$read.depth)
    
    par(mfrow=c(5,1), oma = c(0,0,0,0) + 0.1, mar = c(2,2,1,2), mgp=c(0,1,0))
    res$region   <- OAS1$region
    M <- max(res0$baseline.mean+z.threshold*sqrt(res0$baseline.var),
		res1$baseline.mean+z.threshold*sqrt(res1$baseline.var), 
		res2$baseline.mean+z.threshold*sqrt(res2$baseline.var))
    m <- min(res0$baseline.mean-z.threshold*sqrt(res0$baseline.var),
		res1$baseline.mean-z.threshold*sqrt(res1$baseline.var),
		res2$baseline.mean-z.threshold*sqrt(res2$baseline.var)) 

    ylim         <- exp(c(m,M))
    plot(res0, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="(genotype AA)", ylim=ylim)
    plot(res1, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="(genotype AG)", ylim=ylim)
    plot(res2, z.threshold=2, highlight=FALSE, is.xaxis=FALSE, what="baseline", main="(genotype GG)", ylim=ylim)
    plot(res, z.threshold=2, is.xaxis=FALSE)
    transcripts <- get.transcripts(file.path(path.package("multiseq"),"extdata","hg19.OAS1.refGene.part.gp"), 
                                   OAS1$region)
    plot(transcripts, OAS1$region)

## ----plot----------------------------------------------------------------
    res$intervals <- get.intervals(res, what="effect")
    res$intervals 

## ----load_seq_data-------------------------------------------------------
    setwd(file.path(path.package("multiseq"), "extdata"));
    samplesheet <- file.path(path.package("multiseq"),"extdata","samplesheetEncode.txt")
    samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
    g <- factor(samples$Type)
    g <- match(g, levels(g))-1
    if (noExecutable("wigToBigWig")){
       data(dat, package="multiseq")
       region      <- dat$region
       x <- dat$x
    }else{
       region      <- "chr1:11740409-11756792"
       x           <- get.counts(samplesheet, region) 
    }

## ----testing_on_chipseq--------------------------------------------------
      #smooth data in each genotype class
      res0        <- multiseq(x=x[which(g==0),], minobs=1, read.depth=samples$ReadDepth[which(g==0)])
      res1  	  <- multiseq(x=x[which(g==1),], minobs=1, read.depth=samples$ReadDepth[which(g==1)]) 
      #find an effect given a covariate
      res         <- multiseq(x=x, g=g, minobs=1, read.depth=samples$ReadDepth)

      res0$region <- region
      res1$region <- region
      res$region <- region

      get.intervals(res0, what="baseline")
      get.intervals(res1, what="baseline")   
      get.intervals(res) 

      #plot
      par(mfrow=c(3,1), oma = c(0,0,0,0) + 0.1, mar = c(2,2,1,2), mgp=c(0,1,0))      
      z.threshold=2
      M          <- max(res0$baseline.mean+z.threshold*sqrt(res0$baseline.var), 
                        res1$baseline.mean+z.threshold*sqrt(res1$baseline.var))
      m          <- min(res0$baseline.mean-z.threshold*sqrt(res0$baseline.var), 
                        res1$baseline.mean-z.threshold*sqrt(res1$baseline.var))
      ylim       <- exp(c(m, M))
      plot(res0, z.threshold, is.xaxis=FALSE, what="baseline", main=samples$Type[g==0][1], ylim=ylim)
      plot(res1, z.threshold, is.xaxis=FALSE, what="baseline", main=samples$Type[g==1][1], ylim=ylim)
      plot(res, z.threshold)

## ----track_hub, results='hide'-------------------------------------------
    setwd(file.path(path.package("multiseq"),"extdata"))
    hub_name <- "testMultiseq/dat"
    samplesheetToTrackHub(samplesheet, hub_name, region=region)

## ----multiseqToTrackHub, results='hide'----------------------------------
    res$region    <- region
    multiseqToTrackHub(res, z.threshold=2, hub_name="testMultiseq/multiseq_dat")

