#' @title Split a region string into sequence name, locus start, and locus end
#'
#' @param region: a string specifying a genomic region: reference sequence name, start position (locus.start), end position (locus.end)
#' @examples
#' split_region("chr1:2345-234567")
#' print(split_region$chr)
#' #' @export
#' @return a list with elements chr, locus.start, locus.end
split_region <- function(region){
    split <- unlist(strsplit(region, "\\:|\\-"))
    if (length(split) != 3)
        stop("invalid region: example of a valid region is chr1:2345-234567 ")
    chr = split[1]
    locus.start=as.numeric(split[2])
    locus.end=as.numeric(split[3])

    locus.length=locus.end-locus.start+1
    
    if (locus.start%%1 | locus.end%%1 | locus.length<1) #check that locus.start and locus.end are integers
        stop("Incorrect parameters locus.start and/or locus.end")

    return(list(chr=chr, locus.start=locus.start, locus.end=locus.end))
}

#' @title Prepare input for multiseq function extracting counts in a genomic region from bam, hdf5, or bigWig files
#'
#' @description This functions extracts read counts from *bam*, *hdf5*, or *bigWig* files in a genomic region, preparing input for multiseq. If samples$bamReads is specified then this function extracts reads from the bam files in samples$bamReads using `samtools` (which must be in the USER's path) (no filter applied). Else if samples$h5FullPath is specified this function extracts reads from the hdf5 files in samples$h5FullPath using the R package rhdf5. If samples$bigWigPath is specified this function extracts reads from bigWig files using the executable `bigWigToWig` (which must be in the USER's path).
#'
#' @param samplesheet: a string or a data frame of size N equal to the number of samples. If a string it should be the path to a samplesheet; the samplesheet should contain a column with header SampleID and a column with header either bigWigPath or h5FullPath or bigWigPath. If a data frame, samples$SampleID must be specified and either samples$bigWigPath or samples$h5FullPath or samples$bamReads must be specified.
#' @param region: a string specifying a genomic region: reference sequence name, start position (locus.start), end position (locus.end)
#' @param onlyonend: a bool, defaults to FALSE; use TRUE if input is in bam format and only first end of the paired end read should be used.
#' @export 
#' @return a matrix with N rows and locus.end-locus.start+1 columns containing the number of reads that start at each base in the specified region in each sample. Rownames correspond to samples$SampleID
#' @examples
#'\dontrun{
#' setwd(file.path(path.package("multiseq"),"extdata","sim"))
#' samplesheet="samplesheet.sim.txt"
#' region="chr1:87297710-87305901"
#' M=get.counts(samplesheet, region)
#' }
get.counts <- function(samplesheet, region, onlyoneend=FALSE){
    region       <- split_region(region)
    locus.length <- region$locus.end-region$locus.start+1
    #load counts
    if (is.character(samplesheet)){
        samples <- read.table(samplesheet, stringsAsFactors=F, header=T)
    }else{
        samples <- samplesheet
    }
    M <- NULL
    if (is.null(samples$bigWigPath)){
        if (is.null(samples$h5FullPath)){
            if (is.null(samples$bamReads)){
                stop("no input file provided: provide paths to input files (in hdf5, bigWig or bam format) in the sample sheet file.")
            }else{
                #read bam files into matrix
                command.line <- paste0(region$chr,
                                       ":",
                                       region$locus.start,
                                       "-",
                                       region$locus.end)
                if (onlyoneend==TRUE)
                    command.line <- paste0(command.line,"| awk '{if (!($7=="=" && $4>$8)) print}'")
                command.line <- paste0(command.line,                           
                                       " | awk -v s='",
                                       region$locus.start,
                                       "' 'BEGIN{start=0; count=0}",
                                       "{st=$4; if ($4>=s){if (start==st) count+=1; ",
                                       "else {if (start>0) print start, count; start=st; count=1} }}' ")
                
                for (bamfile in samples$bamReads){
                    command <- paste0("samtools view ",
                                      bamfile,
                                      " ",
                                      command.line)
                    print(command)
                    con <- pipe(command, open = "r")
                    v   <- rep(0, locus.length)
                    while(length(oneLine <- readLines(con, n = 1)) > 0){
                        oneLine <- as.numeric(unlist(strsplit(oneLine," ")))
                        v[oneLine[1]-region$locus.start+1] <- oneLine[2]
                    }
                    close(con)
                    M <- rbind(M, v)
                }
            }
        }else{
            #read hdf5 files into R matrix
            for (h5file in samples$h5FullPath){
                print(paste0("Loading ", h5file))
                print(paste0("h5read(", h5file, ",", chr, ", index=list(", region$locus.start, ":", region$locus.end, ")"))
                M <- rbind(M, h5read(h5file, region$chr, index=list(region$locus.start:region$locus.end)))
            }
        }
    }else{
            #read bigWig files into R matrix
        for (bigWigfile in samples$bigWigPath){
            print(paste0("Loading ", bigWigfile))
            wigfile <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".wig")
            command=paste0("bigWigToWig ", bigWigfile, " -chrom=", region$chr, " -start=", region$locus.start, " -end=", region$locus.end, " stdout | grep -v fixed > " , wigfile)
            print(command)
            system(command)
            M <- rbind(M, as.numeric(as.matrix(read.table(wigfile, stringsAsFactors=F, header=FALSE))))
            file.remove(wigfile)
        }
    }
    
    row.names(M) <- samples$SampleID
    return(M)
}


#' @title getTranscripts
#' @description This function extracts transcript annotation from a file in GenePred format ithat overlap a specified region
#' @param GenePredIn: a file in GenePred format containing gene annotation
#' @param region: a string specifying a genomic region
#' @export
#' @examples
#'\dontrun{
#' GenePredIn="./knownGene_hg18.txt"
#' region="chr5:131989505-132120576"
#' getTranscripts(GenePredIn, region)
#' }
getTranscripts <- function(GenePredIn, region){
    genePred = data.frame(lapply(read.table(GenePredIn,
        fill=1,
        comment.char="",
        header=FALSE),
        as.character),
        stringsAsFactors=FALSE)
    return(genePred[genePred[,2]==chr &((genePred[,4]<=locus.start & genePred[,5]>=locus.start)|genePred[,4]<=locus.end & genePred[,5]>=locus.start),])
}

get.exons.start.end <- function(transcript){	
    exst <- as.numeric(strsplit(as.character(transcript[9]), ",")[[1]]) + 1
    exen <- as.numeric(strsplit(as.character(transcript[10]), ",")[[1]])
    return(list(exst=exst, exen=exen))
}

#' @title plotTranscripts
#' @param Transcripts: output of getTranscripts
#' @param plotStart: minimum value for x axis
#' @param plotEnd: maximum value for x axis
#' @param is.xaxis: bool, if TRUE plot x axis otherwise don't plot x axis
#' @param main: string, title
#' @param cex: number indicating the amount by which plotting text and symbols should be scaled relative to the default. 1=default, 1.5 is 50% larger, 0.5 is 50% smaller, etc. 
#' @export
#' @examples
#'\dontrun{
#' GenePredIn  <- "./knownGene_hg18.txt"
#' region      <- "chr5:131989505-132120576"
#' Transcripts <- getTranscripts(GenePredIn, region)
#' plotTranscripts(Transcripts)
#' }
plotTranscripts <- function(Transcripts, plotStart=NULL, plotEnd=NULL, is.xaxis=1, main=NULL, cex=1){
    if (is.null(Transcripts)){
        plot(1, type="n", axes=F, xlab="", ylab="")
        return(NULL)
    }
    nr  <- nrow(Transcripts)
    chr <- Transcripts[1,2]
    if (is.null(plotStart)==1)
        plotStart <- min(as.numeric(Transcripts[,4]))
    if (is.null(plotEnd)==1)
        plotEnd   <- max(as.numeric(Transcripts[,5]))
    plot("NA",
         axes="F",
         yaxt="n",
         main=main,
         xlim=c(plotStart,plotEnd),
         ylim=c(0,1),
         xlab=paste0("Position (Mb) on ", chr),
         ylab="",
         font.main=1,
         cex.main=cex,
         cex.lab=cex)
    if (is.xaxis==1){
        tck=axTicks(1)
        tcklab=format(tck/1000000)
        axis(1,
             at=tck,
             lab=tcklab,
             cex.lab=cex,
             cex.axis=cex)
    }
    
    for (i in 1:nr){#set colors
        if (Transcripts[i,3] == "-") fill.col="red" else fill.col="blue"
        border.col=fill.col;
        
        if (nr==1){y=0.5} else y=(nr+1-i)*(1/(nr+1))
        
        rect(Transcripts[i, 4],
             y-0.001,
             Transcripts[i, 5],
             y+0.001,
             col="black",
             border="black")
        trans <- get.exons.start.end(Transcripts[i,])
        
        ll=length(trans$exst)
        for (j in 1:ll)
            rect(trans$exst[j],
                 y-0.02,
                 trans$exen[j],
                 y+0.02,
                 col=fill.col,
                 border=border.col )
    }    
}

#' @title plotResults
#' @description This function plots the output of multiseq (either the effect or the baseline) and *fra* * posterior standard deviation
#' @param res: multiseq output
#' @param fra: a multiplier of the standard deviation
#' @param type: if type is "baseline" plot multiseq baseline output; if type is "effect" plot multiseq efefct output; defaults to "effect"
#' @param title: the title of the plot; defaults to NULL
#' @param ylim: plot ylim; defaults to default ylim
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr5:131989505-132120576"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#'
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#' plotResults(res)
#' plotResults(res,"baseline")
#' }   

plotResults <- function(res, fra=2, type="effect", title=NULL, ylim=NULL){
    if ((is.null(res$baseline.mean) | is.null(res$baseline.var)) & type=="baseline")
        stop("Error: no baseline in multiseq output res")
    if ((is.null(res$effect.mean) | is.null(res$effect.var)) & type=="effect")
        stop("Error: no effect in multiseq output res")
           
    paste0(title, " ", type, " (", fra, " s.d.)")
    if (type=="effect"){
        ybottom     <- res$effect.mean - fra*sqrt(res$effect.var)
        ytop        <- res$effect.mean + fra*sqrt(res$effect.var)
        N=length(res$effect.mean)
    }else if (type=="baseline"){
        ybottom     <- res$baseline.mean - fra*sqrt(res$baseline.var)
        ytop        <- res$baseline.mean + fra*sqrt(res$baseline.var)
        N=length(res$baseline.mean)
    }
    ymax        <- max(ytop) + 0.0000000001
    ymin        <- min(ybottom) - 0.0000000001
    wh.bottom   <- which(ybottom > 0)
    wh.top      <- which(ytop < 0)
    high.wh     <- sort(unique(union(wh.bottom, wh.top)))
    xval        <- 1:N
    col.posi    <- xval[high.wh]
    if(is.null(ylim)) ylim=c(ymin,ymax)
    plot(ybottom, pch=".", ylim=ylim, main=title, col="green", xlab='', xaxt='n')
    points(ytop, pch=".", col="green")
    points(res$effect.mean, pch=".", col="dark green")
    abline(h=0, col="red")

    if (type=="effect"){
        N.polygons  <- length(col.posi)
        if (N.polygons > 0)
            for(j in 1:N.polygons)
                rect(col.posi[j]-0.5, ymin-2, col.posi[j]+0.5, ymax+2,
                          col=rgb(1, 0, 0,0.5), border=NA, lty=NULL)
    }
}

#' @title get.effect.intervals
#' @description This function outputs intervals where multiseq found strong effect (zero is outside of +/- *fra* * posterior standard deviation). Output interval is in bed format (start is 0-based, end is 1-based)
#' @param res: multiseq output
#' @param fra: a multiplier of the standard deviation; this function will output intervals where there is an effect at plus or minus fra * standard deviation.
#' @param region:  a string specifying a genomic region: reference sequence name, start position (locus.start), end position (locus.end); defaults to NULL if provided, the function will output the interval in genomic coordinates
#' @export
#' @return a list with elements chr, start, end, sign (of the effect), fra, type (type can be "local" or "sequence" and specifies whether start and end are relative to the genomic sequence
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr5:131989505-132120576"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' get.effect.intervals(res, fra, region))
#' }
get.effect.intervals <- function(res, fra, region=NULL){
    if (is.null(res$effect.mean) | is.null(res$effect.var))
        stop("Error: no effect or effect var in multiseq output res")
    
    effect.start <- NULL
    effect.end   <- NULL
    effect.sign  <- NULL
    toreturn     <- NULL
    
    x=(res$effect.mean + fra * sqrt(res$effect.var) < 0)
    y=(res$effect.mean - fra * sqrt(res$effect.var) > 0)

    xs=sum(x)
    if (xs>0){ #if res$effect.mean + fra * sqrt(res$effect.var) < 0 at at least one position
        rlex=rle(x)
        boundaries=c(0,cumsum(rlex$lengths))
        for (i in 1:(length(boundaries)-1)){
            if (rlex$values[i]){ #is TRUE
                effect.start=c(effect.start,boundaries[i]) #0-based
                effect.end=c(effect.end,boundaries[i+1]) #1-based
                effect.sign=c(effect.sign,"-")
            }
        }
    }

    ys=sum(y)
    if (ys>0){
        rley=rle(y)
        boundaries=c(0,cumsum(rley$lengths))
        for (i in 1:(length(boundaries)-1)){
            if (rley$values[i]){ #is TRUE
                effect.start=c(effect.start,boundaries[i]) #0-based
                effect.end=c(effect.end,boundaries[i+1]) #1-based
                effect.sign=c(effect.sign,"-")
            }
        }
    }
    
    type="local"
    if (!is.null(effect.start)){
        if (!is.null(region)){
            region       <- split_region(region)
            if (!is.null(region$locus.start)&!is.null(region$locus.end)&!is.null(region$chr)){
                effect.start       <- region$locus.start+effect.start-1
                effect.end         <- region$locus.start+effect.end
                type               <- "sequence"
                toreturn$chr       <- region$chr
            }else
                "WARNING: missing locus start or locus end; effect start and end are local and not relative to the sequence"
        }
    }
    toreturn$start       <- effect.start
    toreturn$end         <- effect.end
    toreturn$sign        <- effect.sign
    toreturn$fra         <- fra
    toreturn$type        <- type

    return(toreturn)
}

#' @title write.effect.intervals
#' @description Write intervals with strong effect (zero is outside of +/- *fra* * posterior standard deviation) to bed file
#' @param intervals: output of get.effect.intervals
#' @param dir.name: output directory
#' @param fra: a multiplier of the standard deviation; this function will write intervals where there is an effect at plus or minus fra * standard deviation.
#' @param region:  a string specifying a genomic region: reference sequence name, start position (locus.start), end position (locus.end); defaults to NULL
 
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr5:131989505-132120576"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' write.effect.intervals(res, "~/output_folder", fra=2, region))
#' }
write.effect.intervals <- function(res, dir.name, fra=2, region=NULL){
    dir.create(dir.name, showWarnings = FALSE, recursive=TRUE)
    if (is.null(res$intervals) & is.null(region))
    if (is.null(res$intervals)){
        res$intervals <- get.effect.intervals(res, fra, region)
    }else{
         if (res$intervals$fra != fra)
             res$intervals <- get.effect.intervals(res, fra, region)
    }
             
    bedfile   <- file.path(dir.name, paste0("multiseq.effect.", fra, "sd.bed"))
    if (res$intervals$type=="sequence"){
        if (!is.null(res$intervals$start)){
            write(paste(res$intervals$chr,
                        res$intervals$start,
                        res$intervals$end,
                        ".",
                        "1000",
                        res$intervals$sign,
                        sep="\t"),
                  file=bedfile)
        }
    }else
        stop(paste("ERROR: missing multiseq output intervals; cannot generate", bedfile))
}

#' @title get.effect.length
#' @description This function returns the total lenght of intervals where multiseq found an effect at fra * posterior standard deviation
#' @param res: multiseq output
#' @param fra:  a multiplier of the standard deviation; this function will return the length of intervals where there is an effect at plus or minus fra * standard deviation.
#' @export   
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr5:131989505-132120576"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' get.effect.length(res, fra=2))
#' }
get.effect.length <- function(res, fra){
    if (is.null(res$intervals))
        res$intervals <- get.effect.intervals(res, fra)
    
    return(sum(res$intervals$end-res$intervals$start))
}

#' @title write.effect.mean.variance.gz
#' @description This function saves multiseq results in file effect_mean_var.txt.gz, a file with two columns: first column is effect (baseline) mean and second column is effect (baseline) variance
#' @param res: multiseq output
#' @param fra: a multiplier of the standard deviation; this function will return the length of intervals where there is an effect at plus or minus fra * standard deviation.
#' @param type: if type is "baseline" plot multiseq baseline output; if type is "effect" plot multiseq efefct output; defaults to "effect"
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr5:131989505-132120576"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' write.effect.mean.variance.gz(res, "~/output_folder")
#' }
write.effect.mean.variance.gz <- function(res, dir.name, type="effect"){
    dir.create(dir.name, showWarnings = FALSE, recursive=TRUE)
    gz1      <- gzfile(file.path(dir.name, "effect_mean_var.txt.gz"), "w")
    if (type=="baseline")
        if (is.null(res$baseline.mean) | is.null(res$baseline.var) ){
            close(gz1)
            stop("no baseline in multiseq output res")
        }else{
            write.table(cbind(res$baseline.mean, res$baseline.var), col.names=FALSE, row.names=FALSE, file=gz1)
        }
    if (type=="effect")
        if (is.null(res$effect.mean) | is.null(res$effect.var)){
           close(gz1)
           stop("no effect in multiseq output res")
       }else{
           write.table(cbind(res$effect.mean, res$effect.var), col.names=FALSE, row.names=FALSE, file=gz1)
       }
    close(gz1)
}

#' @title write.summary (to be removed)
#' @description File summary has 6 columns: chr, locus start, locus end, number of bases with effect at 2sd, number of bases with effect at 3sd, running time of multiseq
#' @param res: multiseq output
#' @param dir.name: output directory
#' @param timing: an integer indicating the running time
#' @export
write.summary <- function(res, dir.name, timing, region){
    region   <- split_region(region)
    Neffect2 <- get.effect.length(res, 2)
    Neffect3 <- get.effect.length(res, 3)
    write.table(t(c(region$chr, region$locus.start, region$locus.end, Neffect2, Neffect3, timing)),
                quote = FALSE,
                col.names=FALSE,
                row.names=FALSE, file=file.path(dir.name,"summary"))
}
