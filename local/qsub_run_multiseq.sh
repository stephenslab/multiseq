#!/bin/sh
#
#
#**************************************************************************************************#    
#            FILE:       qsub_run_multiseq.sh                            
#           USAGE:       qsub_run_multiseq.sh sample_sheet OUT_DIR < list_loci
#                        note: make sure that the bed file list_loci ends with a new line
#                        sample_sheet is a text file containing sample specification
#                        OUT_DIR is the directory where you want to save results            
#                        (If g is known you can give the path to g as a third argument to the script) 
#     DESCRIPTION:       use this script to run run.multiseq.R on the list of loci in list_loci
#         EXAMPLE:       qsub_run_multiseq.sh samplesheet.txt ./results_run_multiseq/ < list_loci
#**************************************************************************************************#  

MEM=2g
sample_sheet=$1
OUT_DIR=$2"/"
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

#run multiseq on each locus in list_loci using qsub
while read region; do
    locus="`echo ${region} | awk '{print $1':'$2+1'-'$3}'`"
    PROCESS_NAME=`echo ${region} | awk '{print $1"."$2"."$3}'`
    echo "Rscript run.multiseq.R ${sample_sheet} ${locus} ${OUT_DIR} $fitted_g" | \
	qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
        -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
done 
