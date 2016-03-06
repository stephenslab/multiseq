## About this repo

This repository contains an R package, **multiseq**, and a folder `local` with scripts to run **multiseq** on multiple loci on an SGE cluster.

### R package **multiseq**

This repository contains an R package, **multiseq**, that is ongoing work in the [Stephens lab](http://stephenslab.uchicago.edu/) at the University of Chicago. 

To download the package click on [link](https://github.com/stephenslab/multiseq/releases).

After uncompressing the downloaded archive into folder `multiseq`, you can find installation instructions in `multiseq/README.Rmd`.


### Instructions for developers

If you are developing the package, you can find useful (but not yet documented) instructions in file `multiseq/package/build.R` (how to build, document, archive the package).


### Running multiseq on multiple loci (mostly for internal use in the Stephens Lab)

This repository contains a folder `local` with scripts to run **multiseq** on multiple loci on an SGE cluster. In what follows we assume that you downloaded the repository to `~/src/multiseq`.

To run **multiseq** on a list of loci specified in a [bed](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) file `bed_file` using data specified in `samplesheet.txt` (the samplesheet format is explained in the **multiseq** vignette) and write output to `OUT_DIR` use script `~/src/multiseq/local/qsub_run_multiseq.sh`:

    #specify output folder OUT_DIR
    OUT_DIR="~/output_run_multiseq/"

    #run script
    samplesheet="/path/to/samplesheet.txt"
    bed_file="/path/to/list_loci.bed"
    cd ~/src/multiseq/local/ 
    sh qsub_run_multiseq.sh ${samplesheet} ${OUT_DIR} ${bed_file}

This script will first load data into R objects using 5 parallel jobs. When loading is complete (see option -hold_jid in file qsub_run_multiseq.sh) this script will submit (in parallel) as many jobs as the number of lines in the bed file. 

Note: To generate a bed file `list_loci.bed` with 5000 adjacent intervals of size 2^12=4096 on chr1 using hg19 annotation, use script `~/src/multiseq/local/write_list_loci.sh` as follows:
 
    #specify path to output listo of loci file
    bed_file="/path/to/list_loci.bed"

    #run script
    window_size=4096
    number_intervals=5000
    #if using hg19
    chrom_file="~/src/multiseq/package/multiseq/inst/extdata/chromosome.lengths.hg19.txt"
    cd ~/src/multiseq/local/
    sh write_list_loci.sh ${window_size} ${chrom_file} | head -n ${number_intervals}  > ${bed_file}



Instructions to run parse_bw_bb.sh

    #parse=~/src/NGS_utils/scripts/parse_bw_bb.sh                                                                                                                                                 #file_chrom_full_len=$course_repodir"/data/chromosome.lengths.hg19.txt"                                                                                                                        
    #cd "/data/internal/solexa_mountpoint/epantaleo/simulations/"$DATA_NAME"/"                                                                                                                     
    #sh $parse $file_chrom_full_len   
