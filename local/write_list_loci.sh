#!/bin/bash                         
#this script requires bedtools to be in the user's path
#create bed file with adjacent intervals of size window_size
#usage: sh write_list_loci.sh 131072 $HOME/src/multiseq/data/chromosome.lengths.hg19.txt | head 5000
                                                                                 
window_size=$1 #8192=2^13
file_with_chromosome_lengths=$2
if [ ! -e  $file_with_chromosome_lengths ];then
    echo "ERROR: specify/check path of the file containing chromosome lengths" >&2
else
    bedtools makewindows -g $file_with_chromosome_lengths -w $window_size 
fi
