Download the R package **multiseq**
-----------------------------------
This repository contains a R package, **multiseq**, to ..........smooth and or analyze functional data from one or multiple samples. This is ongoing work in the [Stephens lab](http://stephenslab.uchicago.edu/) at the University of Chicago. 

To download the package go to [link](https://github.com/stephenslab/multiseq/tree/master/package/multiseq.tar.gz) and click on "View Raw" or, alternatively, click on [link](https://github.com/stephenslab/multiseq/blob/master/package/multiseq.tar.gz?raw=true).

After uncompressing the downloaded archive into folder `multiseq`, you can find installation instructions in `multiseq/README.md`.


Additional material in the repository
-------------------------------------
This repo also contains a folder `local` with scripts to run **multiseq** on multiple loci on a SGE cluster. If you want to run **multiseq** on a list of loci specified in a bed file (e.g.: `list_loci.bed`) using data contained in sample sheet `.........samplesheet` and submit jobs to the cluster, then you can use the script `qsub_run_multiseq.sh` in folder `local`:

   #specify output folder OUT_DIR
   OUT_DIR="."
   sh qsub_run_multiseq.sh ..........saplesheet OUT_DIR < list_loci.bed

This script will submit as many jobs as there are lines in `list_loci.bed`.

To generate a bed file `list_loci.bed` with 5000 adjacent intervals of size 131072 on chr1 you could do:
 
    window_size=131072
    chrom_file=..........$HOME/src/multiseq/data/chromosome.lengths.hg19.txt
    sh ........write_list_loci.sh $window_size $chrom_file | head -n 5000 > list_loci.bed



