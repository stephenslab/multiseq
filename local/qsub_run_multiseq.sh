#!/bin/sh
#
#
#**************************************************************************************************#    
#            FILE:       qsub_run_multiseq.sh                            
#           USAGE:       qsub_run_multiseq.sh sample_sheet OUT_DIR list_loci
#                        sample_sheet is a text file containing sample specification
#                        OUT_DIR is the directory where you want to save results            
#                        (If g is known you can give the path to g as a third argument to the script) 
#     DESCRIPTION:       use this script to run run.multiseq.R on the list of loci in list_loci
#         EXAMPLE:       qsub_run_multiseq.sh samplesheet.txt ./results_run_multiseq/ list_loci
#**************************************************************************************************#  

sample_sheet=$1
OUT_DIR=$2"/"
list_loci=$3

fitted_g=$OUT_DIR"/fitted.g.RData"

#check arguments
if [ -z ${OUT_DIR} ] || [ -z ${sample_sheet} ]; then
    echo "ERROR: incorrect number of arguments. Aborting. "
    exit 1
fi

#create output directory (if it doesn't already exist) 
if [ ! -d "${OUT_DIR}" ]; then
    mkdir ${OUT_DIR}
fi

#create directory for logs (if it doesn't already exist)
LOG_DIR=${OUT_DIR}"logs/"
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi

#load data into R objects using 5 parallel jobs
nlines=`cat $list_loci | wc -l`
let nn=${nlines}/5 
split -l ${nn} --suffix-length=5 ${list_loci} list_loci_part_
for ll in `ls list_loci_part_*`; do
    #load data
    PROCESS_NAME=`read_${ll}`
    echo "Rscript load.data.R ${sample_sheet} ${ll} ${OUT_DIR}" | \
	qsub -l h_vmem=2g -v PATH -cwd -N ${PROCESS_NAME} \
	-o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err" 
    
    #wait for data loading (option -hold_jid to qsub)
    #run multiseq on all regions in file ${ll} where data was loaded 
    while read region; do
	locus="`echo ${region} | awk '{print $1':'$2+1'-'$3}'`"
	PROCESS_NAME=`echo ${region} | awk '{print $1"."$2"."$3}'` 
	echo "Rscript run.multiseq.R ${sample_sheet} ${locus} ${OUT_DIR} ${fitted_g}" | \
	    qsub -hold_jid read_${ll} -l h_vmem=2g -v PATH -cwd -N ${PROCESS_NAME} \
	    -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
    done < ${ll}
done
