#' @title Split a region string into sequence name, region start position, and region end position.
#'
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position
#' @examples
#' region=split_region("chr1:2345-234567")
#' print(region)
#' print(region$chr)
#' print(region$start)
#' @export
#' @keywords internal
#' @return a list with elements \code{chr}, \code{start}, \code{end}. 
split_region <- function(region){
    split <- unlist(strsplit(region, "\\:|\\-"))
    if (length(split) != 3)
        stop("invalid region: example of a valid region is chr1:2345-234567 ")
    chr = split[1]
    start=as.numeric(split[2])
    end=as.numeric(split[3])

    locus.length=end-start+1
    
    if (start%%1 | end%%1 | locus.length<1) #check that start and end are integers
        stop("Incorrect parameters start and/or end")

    return(list(chr=chr, start=start, end=end))
}

#' @title Prepare input for multiseq function extracting counts in a genomic region from \code{bam}, \code{hdf5}, or \code{bigWig} files.
#'
#' @description This functions extracts read counts from \code{bam}, \code{hdf5}, or \code{bigWig} files in a genomic region, preparing input for \code{\link{multiseq}}. If \code{samples$bamReads} is specified then this function extracts reads from the bam files in \code{samples$bamReads} using \code{samtools} (which must be in the USER's path) (no filter applied). Else if \code{samples$h5FullPath} is specified this function extracts reads from the \code{hdf5} files in \code{samples$h5FullPath} using the R package \code{rhdf5}. If \code{samples$bigWigPath} is specified this function extracts reads from \code{bigWig} files using the executable \code{bigWigToWig} (which must be in the USER's path).
#'
#' @param samplesheet: a string or a data frame of size N equal to the number of samples. If a string it should be the path to a samplesheet; the samplesheet should contain a column with header \code{SampleID} and a column with header either \code{bigWigPath} or \code{h5FullPath} or \code{bigWigPath}. If a data frame, \code{samples$SampleID} must be specified and either \code{samples$bigWigPath} or \code{samples$h5FullPath} or \code{samples$bamReads} must be specified.
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position
#' @param onlyonend: a bool, defaults to FALSE; use TRUE if input is in bam format and only first end of the paired end read should be used.
#' @export 
#' @return a matrix with \code{N} rows and \code{end-start+1} columns containing the number of reads that start at each base in the specified region in each sample. Rownames correspond to \code{samples$SampleID}
#' @examples
#'\dontrun{
#' setwd(file.path(path.package("multiseq"),"extdata","sim"))
#' samplesheet="samplesheet.sim.txt"
#' region="chr1:154206209-154214400"
#' x=get.counts(samplesheet, region)
#' }
get.counts <- function(samplesheet=NULL, region=NULL, onlyoneend=FALSE){
    if (is.null(region) | is.null(samplesheet))
        stop("Invalid argument")
    #load samplesheet
    region       <- split_region(region)
    locus.length <- region$end-region$start+1
    #load counts
    if (is.character(samplesheet)){
        samples <- read.table(samplesheet, stringsAsFactors=F, header=T)
    }else{
        samples <- samplesheet
    }
    #load data
    M <- NULL
    if ("bamReads" %in% colnames(samples) & ! noExecutable("samtools")){
        #read bam files into matrix
        cmd <- paste0(region$chr,
                      ":",
                      region$start,
                      "-",
                      region$end)
        if (onlyoneend==TRUE)
            cmd <- paste0(cmd,"| awk '{if (!($7==\"=\" && $4>$8)) print}'")
        cmd <- paste0(cmd,                           
                      " | awk -v s='",
                      region$start,
                      "' 'BEGIN{start=0; count=0}",
                      "{st=$4; if ($4>=s){if (start==st) count+=1; ",
                      "else {if (start>0) print start, count; start=st; count=1} }}' ")
        
        for (bamfile in samples$bamReads){
            command <- paste0("samtools view ",
                              bamfile,
                              " ",
                              cmd)
            print(command)
            con <- pipe(command, open = "r")
            v   <- rep(0, locus.length)
            while(length(oneLine <- readLines(con, n = 1)) > 0){
                oneLine <- as.numeric(unlist(strsplit(oneLine," ")))
                v[oneLine[1]-region$start+1] <- oneLine[2]
            }
            close(con)
            M <- rbind(M, v)
        }
    }else if ("h5FullPath" %in%  colnames(samples)) {
                                        #read hdf5 files into R matrix
        for (h5file in samples$h5FullPath){
            print(paste0("Loading ", h5file))
            print(paste0("h5read(", h5file, ",", region$chr, ", index=list(", region$start, ":", region$end, ")"))
            M <- rbind(M, h5read(h5file, region$chr, index=list(region$start:region$end)))
        }
    }else if ("bigWigPath" %in%  colnames(samples) & !noExecutable("wigToBigWig")){  
        #read bigWig files into R matrix
        for (bigWigfile in samples$bigWigPath){
            print(paste0("Loading ", bigWigfile))
            wigfile <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".wig")
            command=paste0("bigWigToWig ", bigWigfile,
                " -chrom=", region$chr,
                " -start=", region$start,
                " -end=", region$end,
                " stdout | grep -v fixed > " ,
                wigfile)
            print(command)
            system(command)
            M <- rbind(M, as.numeric(as.matrix(read.table(wigfile, stringsAsFactors=F, header=FALSE))))
            file.remove(wigfile)
        }
    }else
    stop(paste("no input file provided and/or no executables available to read input:",
               "provide paths to input files (in hdf5, bigWig or bam format) in the samplesheet",
               "file and/or add required executables to user's PATH."))
    
    row.names(M) <- samples$SampleID
    return(M)
}

#' @title Extract transcript annotation (in a specified genomic region) from a file in \code{GenePred}.
#' format (a standard gene annotation format).
#' @param GenePredIn: a file in \code{GenePred} format containing gene annotation
#' @param region: a string specifying a genomic region
#' @export
#' @examples
#' GenePredIn <- file.path(path.package("multiseq"),"extdata","sim","hg19.ensGene.part.gp")
#' region     <- "chr1:154206209-154214400"
#' get.transcripts(GenePredIn, region)
get.transcripts <- function(GenePredIn, region){
    genePred    <- data.frame(lapply(read.table(GenePredIn,
                                                fill=1,
                                                comment.char="",
                                                header=FALSE),
                                     as.character),
                              stringsAsFactors=FALSE)
    region      <- split_region(region)
    Transcripts <-  genePred[genePred[,2]==region$chr &((genePred[,4]<=region$start & genePred[,5]>=region$start)|genePred[,4]<=region$end & genePred[,5]>=region$start),]
    t           <- list(t=Transcripts)
    return(structure(t,class="transcripts"))
}

get.exons.start.end <- function(transcript){	
    exst <- as.numeric(strsplit(as.character(transcript[9]), ",")[[1]]) + 1
    exen <- as.numeric(strsplit(as.character(transcript[10]), ",")[[1]])
    return(list(exst=exst, exen=exen))
}


#' @title Plot transcripts; if \code{region} is specified then plot only the specified region.
#' @description The only required field for this function is Transcripts.
#' @param Transcripts:  output of \code{\link{get.transcripts}}
#' @param region: region to be plotted
#' @param is.xaxis: bool, if TRUE plot \code{x} axis otherwise don't plot \code{x} axis
#' @export
#' @examples
#' GenePredIn  <- file.path(path.package("multiseq"),"extdata","sim","hg19.ensGene.part.gp") 
#' Transcripts <- get.transcripts(GenePredIn)
#' plot.transcripts(Transcripts)
plot.transcripts <- function(Transcripts, region=NULL, is.xaxis=TRUE, fra=2, what="effect", axes=F, xlim=NULL, type=NULL, cex=NULL, font.main=1, ...){
    if (!is.null(region)){
        region <- split_region(region)
        chr    <- region$chr
        if (is.null(xlim))
            xlim <- c(region$start, region$end)
    }

    #draw x axis
    if (is.null(Transcripts$t)){
        if (is.null(region$chr))
            chr="chr ?"
        y=1
        type="n"
    }else{    
        nr   <- nrow(Transcripts$t)
        tchr <- Transcripts$t[1,2]
        if (!is.null(region$chr)){
            if (tchr!=chr)
                stop("Chromosomes are not consistent.")
        }else{
            chr=tchr
        }
        if (is.null(xlim))
            xlim=c(min(as.numeric(Transcripts$t[,4])), max(as.numeric(Transcripts$t[,5])))
        y="NA"
    }
    if (is.xaxis)
        xlab=paste("Position (Mb) on", chr)
    else
        xlab=""
    plot(y,
                 type=type,
                 xlim=xlim,
                 axes=axes,
                 ylim=c(0,1),
                 yaxt="n",
                 xlab=xlab,
                 ylab="",
                 font.main=1,
                 cex.main=cex,
                 cex.lab=cex,
                 ...)
   
    if (is.xaxis){
        tck=axTicks(1)
        tcklab=format(tck/1000000)
        axis(1,
             at=tck,
             labels=tcklab,
             cex.lab=cex,
             cex.axis=cex)
    }
    #draw transcripts
    if (!is.null(Transcripts$t)){
        for (i in 1:nr){#set colors
            if (Transcripts$t[i,3] == "-") fill.col="red" else fill.col="blue"
            border.col=fill.col;
            
            if (nr==1){y=0.5} else y=(nr+1-i)*(1/(nr+1))
            
            rect(Transcripts$t[i,4],
                 y-0.001,
                 Transcripts$t[i,5],
                 y+0.001,
                 col="black",
                 border="black")
            trans <- get.exons.start.end(Transcripts$t[i,])
            
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
}

#' @title Plot the output of \code{\link{multiseq}} (either the effect or the baseline) and \code{fra} * posterior standard deviation.
#' @param x: multiseq output.
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position; defaults to NULL; if provided (or if not provided but x$region is defined) the x axis label and tickers will contain genomic information.
#' @param is.xaxis: bool, if TRUE plot \code{x} axis otherwise don't plot \code{x} axis.
#' @param fra: a multiplier of the standard deviation.
#' @param what: a string, it can be either "baseline" or "effect"; if "baseline", this function plots multiseq baseline output; if "baseline" this function plots multiseq effect output; defaults to "effect".
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr1:154206209-154214400"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#'
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#' plot.multiseq(res)
#' plot.multiseq(res, what="baseline")
#' }
plot.multiseq <- function(x, region=NULL, is.xaxis=TRUE, fra=2, what="effect",axes=F, type="l", col="green", main=NULL, ylim=NULL, xlab=NULL, ylab=NULL, cex=NULL, ...){
    if ((is.null(x$baseline.mean) | is.null(x$baseline.var)) & what=="baseline")
        stop("Error: no baseline in multiseq output")
    if ((is.null(x$effect.mean) | is.null(x$effect.var)) & what=="effect")
        stop("Error: no effect in multiseq output x")
    
    main <- paste0(what, " (", fra, " s.d.) ", main)
    if (what=="effect"){
        ybottom     <- x$effect.mean - fra*sqrt(x$effect.var)
        ytop        <- x$effect.mean + fra*sqrt(x$effect.var)
        N           <- length(x$effect.mean)
        y           <- x$effect.mean
    }else if (what=="baseline"){
        ybottom     <- x$baseline.mean - fra*sqrt(x$baseline.var)
        ytop        <- x$baseline.mean + fra*sqrt(x$baseline.var)
        N           <- length(x$baseline.mean)
        main        <- paste("log", main)
        y           <- x$baseline.mean
    }
    ymax        <- max(ytop) + 0.0000000001
    ymin        <- min(ybottom) - 0.0000000001
    wh.bottom   <- which(ybottom > 0)
    wh.top      <- which(ytop < 0)
    high.wh     <- sort(unique(union(wh.bottom, wh.top)))
    xval        <- 1:N
    col.posi    <- xval[high.wh]
    if (is.null(ylim))
        ylim=c(ymin,ymax)
    if (is.xaxis){
        if (is.null(xlab))
            xlab="Position"
        if (!is.null(region) | !is.null(x$region)){
            if (!is.null(region)){
                region <- split_region(region)
            }else if (!is.null(x$region)){
                region <- split_region(x$region)
            }
            xlab=paste("Position (Mb) on", region$chr)
        }
    }else{
        xlab=""
    }
    if (is.null(ylab))
        ylab=""
    plot(ybottom,
                 type=type,
                 ylim=ylim,
                 main=main,
                 col=col,
                 xlab=xlab,
                 ylab=ylab,
                 axes=axes)
    #draw axes
    axis(2)
    if (is.xaxis){
        if (!is.null(region) | !is.null(x$region)){
            axis(1,
                 at=seq(1,N,ceiling(N/10)),
                 labels=format(seq(region$start,region$end,ceiling(N/10))/1000000),
                 cex.lab=cex,
                 cex.axis=cex)
        }else
            axis(1)
    }
    points(ytop, type=type, col=col)
    points(y, type=type, col=paste("dark",col))
    
    abline(h=0, col="red")

    #draw intervals with effect
    if (what=="effect"){
        N.polygons  <- length(col.posi)
        if (N.polygons > 0)
            for(j in 1:N.polygons)
                rect(col.posi[j]-0.5, ymin-2, col.posi[j]+0.5, ymax+2,
                          col=rgb(1, 0, 0,0.5), border=NA, lty=NULL)
    }
}

#' @title Print intervals where \code{\link{multiseq}} found strong effect (zero is outside of +/- \code{fra} * posterior standard deviation).
#'
#' @description Output interval is in \code{bed} format (\code{start} is 0-based, \code{end} is 1-based).
#' @param res: \code{\link{multiseq}} output.
#' @param fra: a multiplier of the standard deviation; this function will output intervals where multiseq found strong effect (zero is outside of +/- \code{fra} * posterior standard deviation).
#' @param region: a string specifying a genomic region: reference sequence name, start position, end position; defaults to NULL if provided, the function will output the interval in genomic coordinates.
#' @export
#' @return a list with elements \code{chr}, \code{start}, \code{end}, \code{sign} (of the effect), \code{fra}, \code{type} (type can be "local" or "sequence" and specifies whether start and end are relative to a genomic sequence).
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr1:154206209-154214400"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' fra         <- 2
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' get.effect.intervals(res, fra, region))
#' }
get.effect.intervals <- function(res, fra, region=NULL){
    if (is.null(res$effect.mean) | is.null(res$effect.var))
        stop("ERROR: no effect or effect var in multiseq output")
    
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
        if (!is.null(region)| !is.null(res$region)){
            if (!is.null(region)){
                region       <- split_region(region)
            }else{
                region       <- split_region(res$region)
            }
            if (!is.null(region$start)&!is.null(region$end)&!is.null(region$chr)){
                effect.start       <- region$start+effect.start-1
                effect.end         <- region$start+effect.end
                type               <- "sequence"
                toreturn$chr       <- region$chr
            }else
                warning("WARNING: missing region start or region end; effect start and end are local and not relative to the sequence")
        }
    }
    toreturn$start       <- effect.start
    toreturn$end         <- effect.end
    toreturn$sign        <- effect.sign
    toreturn$fra         <- fra
    toreturn$type        <- type
    
    return(toreturn)
}

#' @title Write effect intervals to a \code{bed} file.
#' @param intervals: output of \code{\link{get.effect.intervals}}
#' @param bedfile: output \code{bed} file
#' @examples
#' \dontrun{
#' fra=2
#' region="chr1:154206209-154214400"
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#' intervals <- get.effect.intervals(res, fra))
#' write.bed(intervals, "out.bed")
#' }
#' @export
write.bed <- function(intervals, bedfile){
    if (is.null(intervals)){
        stop("Invalid input intervals")
    }else if (intervals$type!="sequence"){
        stop("Invalid input intervals")
    }else if (is.null(intervals$start) | is.null(intervals$end) | is.null(intervals$chr)){
        stop("Invalid input intervals")
    }
    dir.create(dirname(bedfile), showWarnings=FALSE, recursive=TRUE)
    cat(paste0(intervals$chr,
               "\t",
               intervals$start,
               "\t",
               intervals$end,
               "\t",
               ".",
               "\t",
               "1000",
               "\t",
               intervals$sign,
               "\n",
               collapse=""),
        sep="",
        file=bedfile)
}
       
    
#' @title Return the total lenght of intervals where \code{\link{multiseq}} found strong effect (zero is outside of +/- \code{fra} * posterior standard deviation).
#' @param res: \code{\link{multiseq}} output
#' @param fra:  a multiplier of the standard deviation; this function will return the length of intervals where \code{\link{multiseq}} found strong effect (zero is outside of +/- \code{fra} * posterior standard deviation).
#' @export   
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr1:154206209-154214400"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Tissue)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' fra=2
#' get.effect.length(res, fra))
#' }
get.effect.length <- function(res, fra){
    if (is.null(res$intervals))
        res$intervals <- get.effect.intervals(res, fra)
    
    return(sum(res$intervals$end-res$intervals$start))
}

#' @title Write \code{\link{multiseq}} results to a compressed file.
#' @description The output file has two columns: first column the effect (or the baseline) mean and second column is the effect (or the baseline) posterior standard deviation ^2.
#' @param res: \code{\link{multiseq}} output.
#' @param file: path to the output file.
#' @param what: if \code{what} is "baseline" then plot \code{\link{multiseq}} baseline output; if \code{what} is "effect" then plot \code{\link{multiseq}} effect output; defaults to "effect".
#' @export
#' @examples
#'\dontrun{
#' #run multiseq on sample data
#' samplesheet <- file.path(path.package("multiseq"), "extdata", "sim", "samplesheet.sim.txt")
#' region      <- "chr1:154206209-154214400"
#' x           <- get.counts(samplesheet, region)
#' samples     <- read.table(samplesheet, stringsAsFactors=F, header=T)
#' g           <- factor(samples$Type)
#' g           <- match(g, levels(g))-1
#' 
#' res <- multiseq(x=x, g=g, minobs=1, lm.approx=FALSE, read.depth=samples$ReadDepth)
#'
#' write.gz(res, "./results.mean.sd2.gz")
#' }
write.gz <- function(res, file="results.mean.sd2.gz", what="effect"){
    dir.create(dirname(file), showWarnings = FALSE, recursive=TRUE)
    gz1      <- gzfile(file, "w")
    if (what=="baseline")
        if (is.null(res$baseline.mean) | is.null(res$baseline.var) ){
            close(gz1)
            stop("no baseline in multiseq output res")
        }else{
            write.table(cbind(res$baseline.mean, res$baseline.var), col.names=FALSE, row.names=FALSE, file=gz1)
        }
    if (what=="effect")
        if (is.null(res$effect.mean) | is.null(res$effect.var)){
           close(gz1)
           stop("no effect in multiseq output res")
       }else{
           write.table(cbind(res$effect.mean, res$effect.var), col.names=FALSE, row.names=FALSE, file=gz1)
       }
    close(gz1)
}

#' @title write.summary (to be removed)
#' @description File summary has 6 columns: chr, region start, region end, number of bases where \code{multiseq}
#' found strong effect (zero is outside of +/- 2 or 3 * posterior standard deviation), running time of \code{multiseq}.
#' @param res: multiseq output
#' @param dir.name: output directory
#' @param timing: an integer indicating the running time
#' @keywords internal
#' @export
write.summary <- function(res, dir.name, timing, region){
    region   <- split_region(region)
    Neffect2 <- get.effect.length(res, 2)
    Neffect3 <- get.effect.length(res, 3)
    write.table(t(c(region$chr, region$start, region$end, Neffect2, Neffect3, timing)),
                quote = FALSE,
                col.names=FALSE,
                row.names=FALSE,
                file=file.path(dir.name,"summary"))
}
