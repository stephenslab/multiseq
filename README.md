# Introduction

This repo contains an R package, *multiseq*, to analyze sequence data from multiple samples and it is ongoing work in the [Stephens lab](http://stephenslab.uchicago.edu/) at the University of Chicago. It also contains a set of R and bash scripts for data import and plot, including R scripts to plot input and output data in the [UCSC Genome Browser](http://genome.ucsc.edu/). 

This document summarizes how to install the package and use it and how to use the R and bash scripts for data import and plot. If you have any questions or run into any difficulty, please don't hesitate to contact us!


# Setup

## Downloading the repository

Let's assume you want to clone this repository into a directory named ~/src:

     mkdir ~/src
     cd ~/src
     git clone https://github.com/stephenslab/multiseq

     cd ~/src/multiseq
     git fetch
     git merge origin/master

## Installing the package

The (in development) *multiseq* package is in package/multiseq.tar.gz. 
It depends on *ashr* a package in development that you can download and install following instructions at [https://github.com/stephens999/ash/README](https://github.com/stephens999/ash/blob/master/README)

After installing *ashr* from within R use:

> install.packages("tools")

> install.packages("rhdf5")

> install.packages("data.table")

> install.packages("~/src/multiseq/package/multiseq.tar.gz",repos=NULL,type="source")

Some functions for sequencing data extraction/manipulation require additional executables to be in the user's PATH. The required executables are: `samtools`, `wigToBigWig`, `bigWigInfo`, and `bedToBigBed`.

## Setting up a mountpoint to visualize results in the UCSC Genome Browser

To visualize data in the UCS Genome Browser some shell environment variables must be set. We recommend setting these variables in your ~/.bashrc or ~/.profile files as follows:

    # specify the mountpoint and the http address associated with the mountpoint
    export MOUNTPOINT_PATH="/some/path"
    export MOUNTPOINT_HTTP_ADDRESS="https:some/address"

where you have to replace "/some/path" with the path to the mountpoint (i.e., "/data/internal/solexa_mountpoint/$USER on the PPS cluster") and "https:some/address" with the http address of the mountpoint. If you have access to the PPS cluster replace "/some/path" with "/data/internal/solexa_mountpoint/$USER" where $USER is your username in the cluster; the http address associated to this mountpoint is password protected (ask Ester).

Make sure that you remember to set these variables after adding them to your .bashrc for the first time. 
You can login again, or do source ~/.bashrc

Remember that when you submit jobs to a compute cluster (e.g. using SGE's qsub), they run in "batch mode" 
and may not execute your ~/.bashrc. To ensure that your jobs have the correct environment variables set, 
you should be able to pass a flag to your cluster submission command (e.g. the -V flag to qsub).

### Testing multiseq

> library(multiseq)

> data(OAS1,package="multiseq")

> M <- OAS1$M

> g <- OAS1$g

> read.depth <- OAS1$read.depth

> res <- multiseq(M, g=g, minobs=1, lm.approx=FALSE, read.depth=read.depth)

> fra=2 #fraction of sd

> plotResults(res,fra)

To print intervals where multiseq found an effect at 2 sd:

> get.effect.intervals(res,fra)

Smooth by group

> res0=multiseq(M[g==0,], minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth[g==0])

> plotResults(res0, fra, type="baseline")


The samplesheet should have the following format (see file ~/src/multiseq/data/sim/samplesheet.sim.txt):

    SampleID Type Tissue Replicate Peaks ReadDepth bigWigPath
    055A RNASeq 24hrShamControl 8 - 1550735 ./data/sim/055A.bw
    055B RNASeq 24hr2uMsimvastatinLPDS 8 - 2350343 ./data/sim/055BNonNull.bw
    056A RNASeq 24hrShamControl 9 - 1320166 ./data/sim/056A.bw
    056B RNASeq 24hr2uMsimvastatinLPDS 9 - 1723647 ./data/sim/056BNonNull.bw

### Testing multiseq on sequencing data

Need a sample sheet (samplesheet) a sequence name (chr), sequence start (start) and end (end) position:
 
> samplesheet="~/src/multiseq/data/sim/samplesheet.sim.txt"

> chr="chr5"

> start=131989505

> end=132120576
Rrun multiseq on all samples in samplesheet or select a subset of samples

> samples=read.table(samplesheet, stringsAsFactors=F, header=T) 

> g=factor(samples$Tissue) 

> g=match(g,levels(g))-1

> M=get.counts(samples, chr, start, end) 

> res=multiseq(M, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)

> fra=2 #fraction of sd

> get.effect.intervals(res,fra)

> plotResults(res,fra)

#to save results in dir.name

> dir.name="~/src/multiseq/data/multiseq_sim/"

> write.effect.mean.variance.gz(res,dir.name)

> write.effect.intervals(res,dir.name,fra)

Smooth by group

> res0=multiseq(M[g==0,], minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth[g==0]) 

> plotResults(res0, fra, type="baseline")

### Visualizing input and output in the UCSC Genome Browser

> hub_name="testMultiseq/sim"

> samplesheetToTrackHub(samplesheet,hub_name)

will create a track hub in "/some/path/testMultiseq/sim/" and will print the following message:

    go to http://genome.ucsc.edu/cgi-bin/hgHubConnect and click on the My Hubs window    
    copy paste the following string in the URL field
    https:some/address/testNGS/sim/hub.txt
    center the genome browser on the region of interest
    if the track is hidden click on show and then refresh

If the read tracks or the bed files are large, make sure enough memory is available to run simulationToTrackHub (e.g., use ql 10g on the PPS cluster).

This is a screenshot of the data in the Genome Browser:
![Image](data/sim/sim.png?raw=true)


As of now we can run multiseq region by region. Output of multiseq on a specific region consists of 3 files:

- *effect_mean_var.txt.gz*:  a file with two columns, first column a mean effect a
nd second column squared standard error of the effect) 
- *multiseq.effect.2sd.bed* and *multiseq.effect.3sd.bed*: two bed files containin
g significant intervals (as computed by multiseq) at 2 and 3 sd, respectively. 

multiseqToTrackHub wil create a track hub with 
- the effect +- 2 standard errors 
- the significant intervals at 2 sd 

in the UCSC Genome Browser. If multiseq output is in folder ~/src/multiseq/data/multiseq_sim/ and ~/src/multiseq/data/chromosome.lengths.hg19.txt is a file with chromosome names and lengths, then:

    region="chr5:131989505-132120576"
    hub_name="testMultiseq/multiseq_sim"
    multiseq_folder="~/src/multiseq/data/multiseq_sim"
    chrom_file="~/src/multiseq/data/chromosome.lengths.hg19.txt"
    multiseqToTrackHub(region, hub_name, multiseq_folder, chrom_file)

will create a track hub named *multiseq* in the "https:some/address/test/" folder and will print the following message:
  
    go to http://genome.ucsc.edu/cgi-bin/hgHubConnect and click on the My Hubs window
    copy paste the following string in the URL field
    https:some/address/test/multiseq/hub.txt
    note: center your genome browser around chr5:131989505-132120576 and make track visible

This is a screenshot of the track hub from the Genome Browser:
![Image](data/sim/plots/multiseq.png?raw=true)

### Running multiseq on multiple loci on the PPS cluster

If you want to run multiseq in the list of loci in list_loci.bed (a bed file) then we can use the script in the folder "local" from spudhead:
   sh qsub_run_multiseq.sh < list_loci.bed
This script will submit to the cluster as many jobs as there are lines in list_loci.bed.

[Note: the bash script qsub_run_multiseq.R runs the R script run.multiseq.R on each locus in list_loci.bed. run.multiseq.R calls multiseq with the appropriate arguments.]
 
Use 
    window_size=131072
    chrom_file=$HOME/src/multiseq/data/chromosome.lengths.hg19.txt
    sh write_list_loci.sh $window_size $chrom_file | head -n 5000 > list_loci.bed
to create a bed file list_loci.bed with 5000 adjacent intervals of size 131072 from chr1
