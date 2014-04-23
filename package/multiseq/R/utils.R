#' Get sample counts in a genomic region from bam, hdf5, or bigWig file
#'
#' If samples$bamReads is specified then this function extracts reads from the bam files in samples$bamReads using `samtools` (which must be in the USER's path). Else if samples$h5FullPath is specified this function extracts reads from the hdf5 files in samples$h5FullPath using the R package rhdf5. Else if samples$bigWigPath is specified this function extracts reads from bigWig files using the executable `bigWigToWig` (which must be in the USER's path).
#'
#' @param samples: a data frame extracted from the samplesheet with N rows
#' @param chr: a string representing a reference sequence name
#' @param locus.start: an integer representing a position on chr
#' @param locus.end: an integer representing a position on chr (must be locus.end>locus.start)
#'
#' @return a matrix with N rows and locus.end-locus.start+1 columns containing the number of reads that start at each base in the specified region in each sample 
get.counts <- function(samples, chr, locus.start, locus.end){
    if (locus.start%%1 | locus.end%%1 | locus.end-locus.start<1) #check that locus.start and locus.end are integers
        stop("Incorrect parameters locus.start and/or locus.end")
    M <-  NULL
    if (is.null(samples$bigWigPath)){
        if (is.null(samples$h5FullPath)){
            if (is.null(samples$bamReads)){
                stop("no input file provided: provide paths to input files (in hdf5, bigWig or bam format) in the sample sheet file.")
            }else{
                #read bam files into matrix
                locus.length <- locus.end - locus.start + 1
                command.line <- paste0(chr,
                                       ":",
                                       locus.start,
                                       "-",
                                       locus.end,
                                       " | awk -v s='",
                                       locus.start,
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
                        v[oneLine[1]-locus.start+1] <- oneLine[2]
                    }
                    close(con)
                    M <- rbind(M, v)
                }
            }
        }else{
                                        #read hdf5 files into R matrix
            for (h5file in samples$h5FullPath){
                print(paste0("Loading ", h5file))
                print(paste0("h5read(", h5file, ",", chr, ", index=list(", locus.start, ":", locus.end, ")"))
                M <- rbind(M, h5read(h5file, chr, index=list(locus.start:locus.end)))
            }
        }
    }else{
                                        #read bigWig files into R matrix
        for (bigWigfile in samples$bigWigPath){
            print(paste0("Loading ", bigWigfile))
            wigfile <- paste0(bigWigfile,".wig")
            system(paste("bigWigToWig", bigWigfile, "stdout | grep -v fixed >" , wigfile))
            M <- rbind(M, as.numeric(as.matrix(read.table(wigfile, stringsAsFactors=F, header=FALSE))))
            file.remove(wigfile)
        }
    }
    
    row.names(M) <- samples$SampleID
    return(M)
}



getTranscripts <- function(GenePredIn, chr, locus.start, locus.end){
    genePred = data.frame(lapply(read.table(GenePredIn,
                                        fill=1,
                                        comment.char="",
                                        header=FALSE),
                             as.character),
                      stringsAsFactors=FALSE)
    return(genePred[genePred[,2]==chr,((genePred[,4]<=locus.start & genePred[,5]>=locus.start)|genePred[,4]<=locus.end & genePred[,5]>=locus.start)])
}

get.exons.start.end <- function(transcript){	
    exst <- as.numeric(strsplit(as.character(transcript[9]), ",")[[1]]) + 1
    exen <- as.numeric(strsplit(as.character(transcript[10]), ",")[[1]])
    return(list(exst=exst, exen=exen))
}

plotTranscripts <- function(Transcripts, plotStart=NULL, plotEnd=NULL, is.xaxis=1, main=NULL, expressions=NULL, cex=1){
    if (is.null(Transcripts)){
        plot(1, type="n", axes=F, xlab="", ylab="")
        return(NULL)
    }
    nr  <- nrow(Transcripts)
    chr <- Transcripts[1,2]
    if (is.null(plotStart)==1)
        plotStart <- min(as.numeric(OutTranscripts[,4]))
    if (is.null(plotEnd)==1)
        plotEnd   <- max(as.numeric(OutTranscripts[,5]))
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
        
        if (!is.null(expressions))
            text(x=Transcripts[i,5],
                 y=y+0.07,
                 labels=signif(expressions[i],1),
                 cex=0.7,
                 pos=3,
                 offset=0.1)
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

plotResults <- function(res, fra=2, title=NULL, ylim=NULL, intervals=TRUE, type="effect"){
    if (is.null(title))
        paste0(title, " ", type, " (", fra, " standard deviations)")
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
                polygon(c(polygon(c(col.posi[j]-0.5, col.posi[j]-0.5, col.posi[j]+0.5, col.posi[j]+0.5),
                                  c(ymin-2, ymax+2, ymax+2, ymin-2),
                                  col=rgb(1, 0, 0,0.5),
                                  border = NA)))
        return(N.polygons)
    }
}


get.effect.intervals <- function(res,fra){
    if (is.null(res$effect.mean))
        stop("Error: no effect in multiseq output")
    
    effect.start=NULL
    effect.end=NULL
    effect.sign=NULL
    
    x=(res$effect.mean + fra * sqrt(res$effect.var) < 0)
    y=(res$effect.mean - fra * sqrt(res$effect.var) > 0)

    x=cumsum(rle(x)$lengths)
    xl=length(x)
    if (xl>0){
        effect.start=c(effect.start,x[seq(1,xl,2)]+1)
        effect.end=c(effect.end,x[seq(2,xl-1,2)])
        effect.sign=c(effect.sign,rep("-",xl))
    }

    y=cumsum(rle(y)$lengths)
    yl=length(y)
    if (yl>0){
        effect.start=c(effect.start,y[seq(1,yl,2)]+1)
        effect.end=c(effect.end,y[seq(2,yl-1,2)])
        effect.sign=c(effect.sign,rep("+",yl))
    }
    effect.coordinates="local"
    if (!is.null(res$chr)){
        if (!is.null(res$locus.start)&!is.null(res$locus.end)){
            effect.start=locus.start+effect.start-1
            effect.end=locus.start+effect.end
            effect.coordinates="sequence"
        }else
            "WARNING: missing locus start or locus end; effect start and end are local and not relative to the sequence"
    }
    return(list(effect.start=effect.start,effect.end=effect.end,effect.sign=effect.sign,effect.fra=fra,effect.coordinates="sequence"))
}

write.effect.intervals <- function(res,dir.name,fra=2){
    if (res$effect.coordinates=="sequence"){
        bedfile <- file.path(dir.name, paste0("multiseq.effect.",fra, "sd.bed"))
        write(paste(res$chr,res$effect.start,res$effect.end,".","1000",res$effect.sign,sep="\t"),
              file=file.path(dir.name, paste0("multiseq.effect.",fra, "sd.bed")))
    }else
        warning(paste("WARNING: missing effect.start or effect.end; cannot generate",bedfile))
}

get.effect.length <- function(res,fra){
    if (is.null(res$effect.start)|res$fra!=fra)
        res <- get.effect.intervals(res,fra)
    return(sum(res$effect.end-res$effect.start))
}

write.effect.mean.variance.gz <- function(res,dir.name){
    gz1      <- gzfile(file.path(dir.name, "effect_mean_var.txt.gz"), "w")
    write.table(cbind(res$effect.mean, res$effect.var), col.names=FALSE, row.names=FALSE, file=gz1)
    close(gz1)
}

#this file has 6 columns
#chr, locus start, locus end, number of bases with effect at 2sd, number of bases with effect at 3sd, running time of multiseq 
write.summary <- function(res,dir.name,timing){
    Neffect2=get.effect.length(res,2)
    Neffect3=get.effect.length(res,3)
    write.table(t(c(chr, locus.start, locus.end, Neffect2, Neffect3, timing)),
                quote = FALSE,
                col.names=FALSE,
                row.names=FALSE, file=file.path(dir.name,"summary"))
}
