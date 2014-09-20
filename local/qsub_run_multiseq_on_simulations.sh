#!/bin/sh
#
#
#**************************************************************************************************#
#    This is a bash script containing code to run multiseq on the simulated reads for each locus
#
#    usage: sh qsub_run_multiseq_on_simulations.sh Encode
#    results will be saved in ~/src/stat45800/tests/simulations/"$DATA_NAME"/results/"$sim_type"/"$peak_type
#**************************************************************************************************#  
DATA_NAME=$1
#data="data"
data="data_2fold_betabinomial_updated_0"
logs="logs_2fold_betabinomial_updated_0"
results="results_2fold_betabinomial_updated_0"
list_loci="list_loci" #"list_loci_13"
#results="results_squarem"
#logs="logs_squarem"



if [ -z $DATA_NAME ]; then
    echo "Error: specify data name"
    exit 1
fi


#seq=hg19
course_repodir=$HOME"/src/stat45800/"
MEM=10g

#if [ $DATA_NAME == "RonHepg2sim" ]; then
#    chrom_file=$course_repodir"/data/chromosome_short.lengths.hg19.txt"
#else if [ $DATA_NAME == "Encode" ]; then
#    chrom_file=$course_repodir"/data/chromosome.lengths.hg19.txt"  
#fi
#fi

LOG_DIR=$course_repodir"/tests/simulations/"$DATA_NAME"/"$results"/"$logs"/"
if [ ! -d $course_repodir"/tests/simulations/"$DATA_NAME"/"$results"/" ]; then
    mkdir $course_repodir"/tests/simulations/"$DATA_NAME"/"$results"/"
fi
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi
for peak_type in "all"; do
#for peak_type in "high" "low" "average"; do
    while read region; do
	locus_dot=`echo ${region} | awk '{start=$2+1; print $1"."start"."$3}'`
	locus="`echo ${region} | awk -v d=':' -v m='-' '{print $1,d,$2+1,m,$3}' OFS='' `"
	for sim_type in "Null" "NonNull"; do
	    OUT_DIR="~/src/stat45800/tests/simulations/"$DATA_NAME"/"$results"/"$sim_type"/"$peak_type
	    #hub_name="simulations/"$DATA_NAME"/"$sim_type"/"$peak_type"/"$locus_dot"/multiseq_new"
	    PROCESS_NAME=$peak_type"."$sim_type"."`echo ${region} | awk '{print $1"."$2"."$3}'`
	    echo $PROCESS_NAME
	    sim_samplesheet=$course_repodir"/tests/simulations/"$DATA_NAME"/"$data"/"$peak_type"/"$locus_dot"/"$sim_type"Samplesheet.txt"
	    echo Rscript run.multiseq.R $sim_samplesheet $locus $OUT_DIR $hub_name $chromfile $seq from $course_repodir"/tests/simulations/"$DATA_NAME"/"$list_loci"/list_loci"
	    echo "mkdir -p $OUT_DIR; Rscript run.multiseq.R $sim_samplesheet $locus $OUT_DIR $hub_name $chrom_file $seq" | \
                qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
                -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err" 
	done
    done < $course_repodir"/tests/simulations/"$DATA_NAME"/"$list_loci"/list_loci"
done

#######################################################################################
#
# ADD CHROMOSOME STRING
# if chrom does not contain the string "chr" add it (needed by the UCSC Genome Browser)
#
#######################################################################################
#parse=~/src/NGS_utils/scripts/parse_bw_bb.sh

#file_chrom_full_len=$course_repodir"/data/chromosome.lengths.hg19.txt" 
#cd "/data/internal/solexa_mountpoint/epantaleo/simulations/"$DATA_NAME"/" 
#sh $parse $file_chrom_full_len
