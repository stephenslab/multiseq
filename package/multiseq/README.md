## Installing the package

After downloading the package from https://github.com/stephenslab/multiseq, follow the steps below in order to install the package.

**Multiseq** depends on **ashr**, a package in development that you can download and install following instructions at [https://github.com/stephens999/ash/README](https://github.com/stephens999/ash/blob/master/README). It also depends on **Rcpp**, and therefore it requires an **R** version higher or equal to **3.0.0**. After installing **ashr** type the following commands in R:

```
    install.packages("tools")
    biocLite("rhdf5")
    install.packages("data.table")
    install.packages("Rcpp")
    install.packages("multiseq.tar.gz",repos=NULL,type="source")
```

## Optional: link to executables needed to read data in `hdf5`, `bam`, or `bigWig` format
This package contains a useful function, `get.counts`, that converts `rhdf5`, `bam` or `bigWig` input files into R objects that can then be used as input to the main `multiseq` function.

In order to use this function make sure you have `samtools`, `bigWigToWig` (part of the USCS tools) and `awk` in your path. If `/data/tools/ucsctools/` is your path to the UCSC tools and `/usr/local/bin/` is your path to `samtools` then add the following lines to your `~/.bashrc` file:

    export PATH=$PATH:"/data/tools/ucsctools/":"/usr/local/bin/"

Also see Note below.

## Optional: set some shell environmental variables to visualize results in the UCSC Genome Browser
This package contains two functions, namely `samplesheetToTrackHub` and `multiseqToTrackHub`, that can be used to display `multiseq` input and output in the UCSC Genome Browser. To display results in the UCSC Genome Browser you need a a folder that is accessible via the Internet. If `/some/path` is the path to a folder that is accessible via the Internet - the "mountpoint" - and `https:some/address` is the http address of the mountpoint, then you should set the following shell environmental variables in your `~/.bashrc` file as follows:

    # the following lines specify the mountpoint and the http address associated with the mountpoint
    export MOUNTPOINT_PATH="/some/path"
    export MOUNTPOINT_HTTP_ADDRESS="https:some/address"

[Note for PPS cluster users in the stephenslab group: replace `/some/path` with `/data/internal/solexa_mountpoint/$USER` where `$USER` is your username in the cluster; the http address associated to this mountpoint is password protected (ask John or Ester).]

Also see Note below.

## Note
After adding those lines to your `.bashrc` for the first time, remember to either login again, or do 
    source ~/.bashrc 

Remember that when you submit jobs to a compute cluster (e.g. using SGE's qsub), they run in "batch mode" and may not execute your `~/.bashrc`. To ensure that your jobs have the correct environmental variables set, you should be able to pass a flag to your cluster submission command (e.g. the `-V` flag to qsub).
