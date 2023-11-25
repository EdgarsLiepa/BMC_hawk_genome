
# ---
# --- About: Create FASTA QC report for raw reads FASTA.
# ---       
# ---
# --- In: 
# ---       Path to directory with FASTQ read files - .fq.gz
# ---      
# --- Out: 
# ---       result summary HTML in ${FILEPATH}/fastqc folder
# ---       result summary .zip in ${FILEPATH}/fastqc folder
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 25.11.22
# --- Modified: 16.01.23
#

#!/bin/bash
#PBS -N FASTQC_report
#PBS -l procs=4
#PBS -l walltime=25:00:00
#PBS -A bmc_pl_bior_covid
#PBS -j oe

# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes


# Directory with trimmed reads files.
# FILEPATH=''


# check if sample name is empty
if [ -z "$FILEPATH" ]; then
    echo "---------------ERROR--------------------"
    echo "----------------------------------------"
    echo "Sample name parameter \${FILEPATH} is empty"
    echo "Please specify \${FILEPATH} in cmd with option -v ARG_NAME1=\"ARG_VALUE1\""
    echo "Exmple:"
    echo "  qsub -v FILEPATH="../path_to_fastq_files" <Path_to_script_directory>/create_fastqc_report.sh"
    echo "----------------------------------------"
    exit 1
fi

echo "File directory $FILEPATH"

for SAMPLE in $(ls ${FILEPATH}/*_paired.fq.gz)
do

date
echo "Analyze $SAMPLE"
/home_beegfs/edgars01/tools/FastQC/fastqc --noextract --nogroup -o ${FILEPATH}/fastqc_trimmedReads $SAMPLE

done;

date