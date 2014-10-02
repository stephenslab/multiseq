MOUNTPOINT_PATH=Sys.getenv("MOUNTPOINT_PATH")
MOUNTPOINT_HTTP_ADDRESS=Sys.getenv("MOUNTPOINT_HTTP_ADDRESS")


#' @title Write genomes.txt and hub.txt files
#' @export
writeTrackHubSkeleton <- function(hub_dir, assembly="hg19", hub_name_string, email="esterpantaleo@gmail.com"){
    assembly_dir = file.path(hub_dir, assembly)

    #write genomes file
    print("write genomes file")
    cat(paste0("genome ", assembly, "\n",
               'trackDb ', assembly, "/trackDbFile.txt\n"),
        file=file.path(hub_dir, 'genomes.txt'), append=FALSE)
    
    #write hub file
    print("write hub file")
    cat(paste0('hub ', hub_name_string, "\n",
               'shortLabel ', hub_name_string, "\n",
               'longLabel ', hub_name_string, "\n",
               "genomesFile genomes.txt\n",
               'email ', email, "\n"),
        file=file.path(hub_dir, 'hub.txt'), append=FALSE)
}


#' @title Write a superTrack with BigBed files in the Track Hub trackdb_file
appendBedSuperTrack <- function(track, shortLabel, longLabel, tracks, shortLabels, longLabels, color, priority=2, out_file){
    cat(paste0('track ', track, "\n",
               'shortLabel ', shortLabel, "\n",
               'longLabel ', longLabel, "\n",
               "superTrack on none\n",
               'priority ', priority, "\n",
               "dragAndDrop subtracks\n\n"), file=out_file, append=TRUE)

    for (i in 1:length(tracks))
            cat(paste0('track signal', shortLabels[i], "\n",
                       "type bigBed\n",
                       'shortLabel ', shortLabels[i], "\n",
                       'longLabel ', longLabels[i], "\n",
                       'parent ', track, "\n",
                       "visibility full\n",
                       'bigDataUrl ', tracks[i], "\n",
                       'color ', color, "\n"), file=out_file, append=TRUE)
}

#' @title Write a message containing instructions on how to visualize the Track Hub
printGoToMessage <- function(hub_name, hub_dir, http_address, region=NULL){
    toprint=paste0("go to http://genome.ucsc.edu/cgi-bin/hgHubConnect and click on the My Hubs window\ncopy paste the following string in the URL field\n",http_address, "/", hub_name, "/hub.txt\nsubmit and ")
    if (!is.null(region))
        toprint=paste0(toprint," center your genome browser around ", region, " and make track visible\n")
    else
        toprint=paste0(toprint," center the genome browser on the region of interest\n")
    toprint=paste0(toprint,"(track has been saved in folder ", hub_dir, ")\n")

    cat(toprint)
}

#' @title convert hdf5 file (with data from a locus) to bigwig file
#'
#' @description This function requires the executable "wigToBigWig" to be in the user's PATH 
#' 
#' @param h5_track: path to the hdf5 track
#' @param chrom_file: path to the file containing chromosome names and lengths
#' @param bigWig_track: name of new bigWig track
#' @param assembly: genome assembly that reads were mapped to; default="hg19"

#' @return no return; it prints a Track Hub folder
#' @examples
#' hdf5ToBigWig(in.h5, chromosome_file, out)
#' converts in.h5 into out.bw using chromosome names and lengths in chromosome_file
#' 
hdf5ToBigWig <- function(h5_track, chrom_file, bigWig_track, assembly="hg19"){
    chromosomes = read.table(chrom_file, stringsAsFactors=FALSE)
    #write wiggle file
    tmpfile_wig = paste0(args.bigWig_track, "_tmp")
    if (file.exists(tmpfile_wig))
        file.remove(tmpfile_wig)
    for (i in 1:nRow(chromosomes)){
        myvalues = h5read(h5file, chromosomes[i,1], index=list(1:chromosomes[i,2]))
        print(paste("processing chromosome ", chromosomes[i,1]))
        cat("fixedStep chrom=%s start=1 step=1", myvalues, sep="\n", file=tmpfile_wig, append=TRUE)
        break
    }
    system(paste("wigToBigWig", tmpfile_wig, chrom_file, bigWig_track))
    file.remove(tmpfile_wig)
}

#' @title Create a UCSC Genome Browser "Track Hub" from read tracks and bed tracks (of significant intervals) listed in a samplesheet.
#'
#' @description This function requires the executables "wigToBigWig", "bedToBigBed", and "bigWigInfo" to be in the user's PATH.
#' Read tracks can be in bam, hdf5, wig or bigwig format and significant intervals can only be in bed format.
#'
#' @param samplesheet: a string specifying the path to the samplesheet; the samplesheet must contain a column with header sampleID and: either a column with header h5FullPath, or a column with header bigWigPath containing the path to the hdf5 files or the bigWig files, respectively. If the samplesheet has a column with header Peaks it must also have a column with header Tissue. Depending on the size of the hdf5 files this code might require a lot of memory
#' @param hub_name: name of the Track Hub. This string can be set to any value; it could contain a path, in which case the path will be relative to the mountpoint (see below); default=paste0(basename(samplesheet),".TrackHub")
#' @param chrom_file: path to the file containing chromosome names and lengths; default=file.path(path.package("multiseq"),"data","chromosome.lengths.hg19.txt")
#' @param assembly: genome assembly that reads were mapped to; default="hg19" 
#' @param mountpoint: path to the directory where the track hub folder will be saved in. This directory should be associated with an http address or an ftp address; default=Sys.getenv("MOUNTPOINT_PATH")
#' @param http_address: http or ftp address associated with the mountpoint; default=Sys.getenv("MOUNTPOINT_PATH")
#' @export
#' @examples
#'\dontrun{
#' setwd(file.path(path.package("multiseq"),"extdata","sim"))
#' samplesheet="samplesheet.sim.txt"
#' samplesheetToTrackHub(samplesheet)
#' }  
samplesheetToTrackHub <- function(samplesheet, hub_name=NULL, chrom_file=file.path(path.package("multiseq"),"data","chromosome.lengths.hg19.txt"), assembly="hg19", mountpoint=MOUNTPOINT_PATH, http_address=MOUNTPOINT_HTTP_ADDRESS){
    if (is.null(hub_name)) hub_name=paste0(basename(samplesheet),".TrackHub")
    samples         <- read.table(samplesheet, stringsAsFactors=F, header=T)
    hub_dir         <- file.path(mountpoint, hub_name)
    hub_name_string <- gsub("/", ".", hub_name)
    dir.create(hub_dir, showWarnings = FALSE, recursive=TRUE)
    assembly_dir    <- file.path(hub_dir, assembly)
    dir.create(assembly_dir, showWarnings = FALSE, recursive=TRUE) 
    sampleids = samples$SampleID

    bigwig_tracks <- NULL
    if ("bigWigPath" %in%  colnames(samples)){
        for (bigwig_track in samples$bigWigPath){
            track_name    <- basename(bigwig_track)
            dir.create(file.path(assembly_dir, dirname(track_name)), showWarnings = FALSE, recursive=TRUE)
            file.copy(from=bigwig_track,
                      to=file.path(assembly_dir, track_name))
            bigwig_tracks <- c(bigwig_tracks, track_name)
        }
    }else if ("h5FullPath" %in%  colnames(samples)){
        track_name <- paste0(basename(h5_track),".bw")
        h52bigwig(h5_track, chrom_file, bigWig_track, assembly)
    }
    
    
    writeTrackHubSkeleton(hub_dir, assembly, hub_name_string)
    cat(paste0("track reads\n",
               "shortLabel reads\n",
               "longLabel reads\n",
               "superTrack on none\n",
               "priority 1\n"), 
        file=file.path(assembly_dir, "trackDbFile.txt"),
        append=FALSE)
        
    #if bigwig files cover a region smaller than 2^20
    #use viewLimits
    command      <- paste0("bigWigInfo ", file.path(assembly_dir,bigwig_tracks[1]), " | grep basesCovered | tr -d \",\"" )
    print(command)
    bigWigLength <- unlist(strsplit(system(command, intern=TRUE)[1], " "))[2]
    print(bigWigLength)
    if (bigWigLength<2^20){
        #find ymax over all bigwig files
        bigWigM=0
        for (bigwig_track in bigwig_tracks){
            command   <- paste("bigWigInfo -minMax", file.path(assembly_dir, bigwig_track))
            bigWigMax <- unlist(strsplit(system(command, intern=TRUE)[1], " "))[2] 
            bigWigM   <- max(bigWigM, bigWigMax)
        }
        cat(paste0("autoScale off\n",
                   "viewLimits 0:",
                   bigWigM, "\n"),
            file=file.path(assembly_dir, "trackDbFile.txt"),
            append=TRUE)
    }else{
        cat("autoScale on\n",
            file=file.path(assembly_dir, "trackDbFile.txt"),
            append=TRUE)
    }
    cat("dragAndDrop subtracks\n\n",
        file=file.path(assembly_dir, "trackDbFile.txt"),
        append=TRUE)

    
    counter=1
    for (bigwig_track in bigwig_tracks){
        cat(paste0("track ", sampleids[counter], "\n",
                   "parent reads\n",
                   "type bigWig\n",
                   "graphType points\n",
                   "visibility full\n",
                   "color 0,0,0\n",
                   "bigDataUrl ", bigwig_track, "\n",
                   "track ", sampleids[counter], "\n",
                   "parent reads\n",
                   "type bigWig\n",
                   "graphType points\n",
                   "visibility full\n",
                   "color 0,0,0\n",
                   "bigDataUrl ", bigwig_track, "\n",
                   "shortLabel ", sampleids[counter], "\n",
                   "longLabel ", hub_name_string, " ", sampleids[counter], "\n\n"),
            file=file.path(assembly_dir, "trackDbFile.txt"),
            append=TRUE)
        counter=counter+1
    }
    
    #if bed files are available in the samplesheet
    #convert bed files into bigBed
    if ("Peaks" %in% colnames(samples)){
        peaks_files = unique(samples$Peaks)
        if (peaks_files[1]!="-"){
                tissues       <- unique(samples$Tissue)
                bigbed_tracks <- NULL
                for (peaks_track in peaks_files){
                        track_name <- basename(peaks_track)
                        if (file_ext(track_name)=="bb"){
                            dir.create(file.path(assembly_dir, dirname(track_name)), showWarnings = FALSE, recursive=TRUE)
                            file.copy(from=peaks_track, to=file.path(assembly_dir, track_name))
                            bigbed_tracks <- c(bigbed_tracks, track_name)
                        }else{
                            #convert bed file into bigBed file
                            bigbed_track <- paste0(track_name, '.bb')
                            cmd          <- paste("bedToBigBed -type=bed6+4",
                                                  peaks_track,
                                                  chrom_file,
                                                  file.path(hub_dir, bigbed_track))
                            print(cmd)
                            system(cmd)
                            bigbed_tracks <- c(bigbed_tracks, bigbed_track)
                        }
                } 
                #write the track hub                    
                appendBedSuperTrack("peaks", "peaks", "peaks", bigbed_tracks, tissues, paste0(hub_name_string,tissues), "0,0,0", out_file=file.path(assembly_dir, "trackDbFile.txt"))
        }
    }
    printGoToMessage(hub_name, hub_dir, http_address)
}

#' @title Create a UCSC genome browser "Track Hub" from multiseq output.
#'
#' @description This function requires the executables "wigToBigWig" and "bedToBigBed" to be in the user's PATH
#'
#' @param region: a region (e.g. chr1:2345-234567)
#' @param hub_name: name of the Track Hub. This string can be set to any value; it could contain a path, in which case the path will be relative to the mountpoint; default="multiseq"
#' @param multiseq_folder: path to the folder containing results from multiseq; this script requires output from multiseq to be in the format effect_mean_var.txt.gz where first column is effect and second column variance
#' @param chrom_file: path to the file containing chromosome names and lengths; default=file.path(path.package("multiseq"),"data","chromosome.lengths.hg19.txt")
#' @param assembly: genome assembly that reads were mapped to; default="hg19"
#' @param mountpoint: path to the directory where the Track Hub folder will be saved in. This directory should be associated with an http address or an ftp address; default=Sys.getenv("MOUNTPOINT_PATH")
#' @param http_address: http or ftp address associated with the mountpoint; default=Sys.getenv("MOUNTPOINT_PATH")
#' @export
#' @examples
#' \dontrun{
#' region="chr1:87297710-87305901"
#' multiseq_folder=file.path(path.package("multiseq"), "extdata", "multiseq_sim")
#' multiseqToTrackHub(region=region, multiseq_folder=multiseq_folder)
#' } 
multiseqToTrackHub <- function(region, hub_name="multiseq", multiseq_folder="./results_run_multiseq/", chrom_file=file.path(path.package("multiseq"),"data","chromosome.lengths.hg19.txt"), assembly="hg19", mountpoint=MOUNTPOINT_PATH, http_address=MOUNTPOINT_HTTP_ADDRESS){
    split_region = unlist(strsplit(region, "\\:|\\-"))
    if (length(split_region) != 3)
           stop("invalid region: example of a valid region is chr1:2345-234567 ")
    chrom = split_region[1]
    locus_start = as.numeric(split_region[2])
    locus_end = as.numeric(split_region[3])

    if (locus_start%%1 | locus_end%%1 | locus_end-locus_start<1) #check that locus.start
         stop("Incorrect parameters locus_start and/or locus_end") 
    hub_dir = file.path(mountpoint, hub_name)
    hub_name_string = gsub("/", ".", hub_name)
    dir.create(hub_dir, showWarnings = FALSE, recursive=TRUE)
    assembly_dir = file.path(hub_dir, assembly)
    dir.create(assembly_dir, showWarnings = FALSE, recursive=TRUE)
 
    multiseq_folder = file.path(multiseq_folder, paste0(c(chrom, locus_start, locus_end), collapse="."))
    multiseq_file = file.path(multiseq_folder, "effect_mean_var.txt.gz")
    #create bigWig file with mean
    mean_track = file.path(assembly_dir, "mean_track.bw")

    cmd = paste0("(echo fixedStep chrom=",
                 chrom,
                 " start=",
                 locus_start,
                 " step=1 ; zcat ",
                 multiseq_file,
                 " | awk '{print $1}' ) | wigToBigWig stdin ",
                 chrom_file,
                 " ",
                 mean_track)
    print(cmd)
    system(cmd)
    
    #create bigWig file with mean+2sd
    mean_plus_2sd_track = file.path(assembly_dir, "mean_plus_2sd_track.bw")
    cmd = paste0("(echo fixedStep chrom=",
                 chrom,
                 " start=",
                 locus_start,
                 " step=1 ; zcat ",
                 multiseq_file,
                 " | awk '{print $1+2*sqrt($2)}' ) | wigToBigWig stdin ",
                 chrom_file,
                 " ",
                 mean_plus_2sd_track)
    print(cmd)
    system(cmd)
            
    #create bigWig file with mean-2sd
    mean_minus_2sd_track = file.path(assembly_dir, "mean_minus_2sd_track.bw")
    cmd = paste0("(echo fixedStep chrom=",
                 chrom,
                 " start=",
                 locus_start,
                   " step=1 ; zcat ",
                 multiseq_file,
                 " | awk '{print $1-2*sqrt($2)}' ) | wigToBigWig stdin ",
                 chrom_file,
                 " ",
                 mean_minus_2sd_track)
    print(cmd)
    system(cmd)
        
    #create bigBed file with significant regions
    multiseq_bed_file = file.path(multiseq_folder, "multiseq.effect.2sd.bed")
    #check if bed file is empty
    no_bed=TRUE
    if (!is.null(file.info(multiseq_bed_file)$size)){
        no_bed=FALSE
        multiseq_peaks_track = file.path(assembly_dir, "multiseq_bed_file.bb")
        cmd = paste("bedToBigBed",
                     multiseq_bed_file,
                     chrom_file,
                     multiseq_peaks_track)
        print(cmd)
        system(cmd)
    }   

    #make track hub
    print("writeTrackHubSkeleton")
    writeTrackHubSkeleton(hub_dir, assembly, hub_name_string)
    #write trackdb_file
    print("write trackdb_file")
    cat("track SuperTrack\n",
        "shortLabel multiseq\n",
        "longLabel Plot of multiseq effect 2 sd\n",
        "superTrack on none\n",
        "priority 1\n\n",
        "track CompositeTrack\n",
        "container multiWig\n",
        "configurable on\n",
        "shortLabel Effect\n",
        "longLabel multiseq\n",
        "visibility full\n",
        "type bigWig\n",
        "autoScale on\n",
        "aggregate transparentOverlay\n",
        "windowingFunction mean\n",
        "superTrack SuperTrack full\n",
        "showSubtrackColorOnUi on\n",
        "smoothingWindow off\n",
        "dragAndDrop subtracks\n\n",
        file=file.path(assembly_dir, 'trackDbFile.txt'))
    shortLabel=c("Mean", "MeanPlus2Sd", "MeanMinus2Sd")
    longLabel=c("multiseq mean effect", "multiseq effect mean - 2 * se", "multiseq effect mean + 2 * se")
    bigDataUrl=c("mean_track.bw", "mean_plus_2sd_track.bw", "mean_minus_2sd_track.bw")
    color=c("0,0,0", "0,255,0", "0,255,0")
    for (i in c(1,2,3))
        cat(paste0("track Subtrack", i, "\n",
                   "type bigWig\n",
                   "shortLabel ", shortLabel[i], "\n",
                   "longLabel ", longLabel[i], "\n",
                   "parent CompositeTrack\n",
                   "graphType points\n",
                   "visibility full\n",
                   "bigDataUrl ", bigDataUrl[i], "\n",
                   "color ", color[i], "\n\n"), file=file.path(assembly_dir, 'trackDbFile.txt'), append=TRUE)
                   
    #plot significant region
    if (no_bed==FALSE)
        appendBedSuperTrack("SuperBedTrack", "multiseq_signal", "multiseq signal", "multiseq_bed_file.bb", "multiseq_signal", "multiseq signal 2sd", "255,0,0", out_file=file.path(assembly_dir, "trackDbFile.txt"))

    printGoToMessage(hub_name, hub_dir, http_address, region)
}
