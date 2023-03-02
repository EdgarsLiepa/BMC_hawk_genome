# ---
# --- About: Index reference file using BWA and samtools
# ---       
# ---
# --- In: 
# ---      Pass Reference file without file extension as a parameter with option -v REFERENCE="<REFERENCE_PATH>"
# ---      Reference sequence should be in fasta format
# --- Out: 
# ---       <REFERENCE>.fai
# ---       <REFERENCE>.dict
# ---       <REFERENCE>.fa.gc
# ---       <REFERENCE>.fai
# ---       <REFERENCE>.fa.amb
# ---       <REFERENCE>.fa.ann
# ---       <REFERENCE>.fa.bwt
# ---       <REFERENCE>.fa.pac
# ---       <REFERENCE>.fa.sa
# ---       <REFERENCE>fna.fai
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 2.02.23

#!/bin/bash
#PBS -N index_reference
#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:00:00
#PBS -A bmc_pl_bior_covid
#PBS -j oe

# --Calculate the number of processors and nodes allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`


# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes

# check if sample name is empty
if [ -z "$REFERENCE" ]; then
    echo "---------------ERROR--------------------"
    echo "----------------------------------------"
    echo "Reference name parameter \${REFERENCE} is empty"
    echo "Please specify path \${REFERENCE} without file extension in cmd with option -v REFERENCE=\"ARG_VALUE1\""
    echo "Exmple:"
    echo "  qsub -v REFERENCE="../path_to_reference_files" <Path_to_script_directory>/index_reference.sh"
    echo "  qsub -v REFERENCE="<REF_DIR>/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic" <Path_to_script>/index_reference.sh"
    echo "----------------------------------------"
    exit 1
fi

echo "Index reference file: $REFERENCE"
echo "BWA (Burrows-Wheeler Aligner) is needed to be installed in the system."

bwa index $REFERENCE.fna

samtools faidx $REFERENCE.fna --fai-idx $REFERENCE.fai 

### Output job resurces used
echo "---------------------"
echo "Job resources used:"
/usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources
