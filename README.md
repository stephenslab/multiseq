## Content

R package **multiseq**
----------------------
This repository contains an R package, **multiseq**, that is ongoing work in the [Stephens lab](http://stephenslab.uchicago.edu/) at the University of Chicago. 

To download the package click on [link](https://github.com/stephenslab/multiseq/blob/master/package/multiseq.tar.gz?raw=true).

After uncompressing the downloaded archive into folder `multiseq`, you can find installation instructions in `multiseq/README.md`.


Running multiseq on multiple loci (mostly for internal use in the Stephens Lab)
-------------------------------------------------------------------------------
This repository contains a folder `local` with scripts to run **multiseq** on multiple loci on an SGE cluster. In what follows we assume that you downloaded the repository to `~/src/multiseq`.

To run **multiseq** on a list of loci specified in a [bed](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) file (e.g.: `~/src/multiseq/data/list_loci.bed`) using data specified in `samplesheet.txt` (the samplesheet format is explained in **multiseq** manual) use script `multiseq/local/qsub_run_multiseq.sh`:

    #specify output folder OUT_DIR
    OUT_DIR="~/output_run_multiseq/"

    #run script
    samplesheet="/path/to/samplesheet.txt"
    bed_file="/path/to/list_loci.bed"
    cd ~/src/multiseq/local/ 
    sh qsub_run_multiseq.sh ${samplesheet} ${OUT_DIR} < ${bed}

This script will submit as many jobs as there are lines in the bed file. By default the script sets memory usage to 10g but you can modify that option in the script if you are running **multiseq** on a smal region.

Note: To generate a bed file `list_loci.bed` with 5000 adjacent intervals of size 2^12=4096 on chr1 using hg19 annotation, use script `~/src/multiseq/local/write_list_loci.sh` as follows:
 
    #specify path to output listo of loci file
    bed_file="/path/to/list_loci.bed"

    #run script
    window_size=4096
    number_intervals=5000
    chrom_file="~/src/multiseq/package/data/chromosome.lengths.hg19.txt"
    cd ~/src/multiseq/local/
    sh write_list_loci.sh ${window_size} ${chrom_file} | head -n ${number_intervals}  > ${bed_file}



