#!/bin/bash                         
#create bed file with adjacent intervals of size window_size
#usage: sh write_list_loci.sh 131072
                                                                                 
window_size=$1 #8192=2^13
file_with_chromosome_lengths=$HOME/src/stat45800/data/chromosome.lengths.hg19.txt
if [ ! -e  $file_with_chromosome_lengths ];then
    echo "ERROR: check path of the file containing chromosome lengths" >&2
else
    bedtools makewindows -g $file_with_chromosome_lengths -w $window_size | head -n 5000
fi
