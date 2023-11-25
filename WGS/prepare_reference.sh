# ---
# --- About: Index reference file using BWA and samtools.
# ---        Create dictionary file using picard tools.
# ---
# --- In: 
# ---      Pass Reference file as a parameter with option -i REFERENCE_PATH
# ---      Reference sequence should be in fasta format with one of the falowing extensions:
# ---      .fasta, .fna, .ffn, .faa, .frn, .fa
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
# --- Date: 3.02.23

#!/bin/bash
#PBS -N prepare_reference
#PBS -l nodes=1:ppn=4
#PBS -l walltime=1:00:00
#PBS -A bmc_pl_bior_covid
#PBS -j oe

# picard tools path
PICARD=/home_beegfs/edgars01/tools/picard.jar

# if job is running on cluster should PBS_JOBID be non-empty
# Compatible with qsub from Torque and PBSPro
if [ ! -z "$PBS_JOBID" ]; then
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
fi

helpFunction()
{
    echo "---------------Description--------------------"
    echo "----------------------------------------"
    echo "-i input reference in fasta format."
    echo "Following extensions are acceptable for fasta file: .fasta, .fna, .ffn, .faa, .frn, .fa"
    echo "----------------------------------------"
    echo "Exmple:"
    echo "  ./prepare_reference.sh -i ../path_to_reference_files"
    echo "  ./prepare_reference.sh -i \<REF_DIR\>/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fasta"
    echo "Exmple with qsub in Torque:"
    echo "  qsub -F "-i ../path_to_reference_files/ref.fasta" <Path_to_script_directory>/prepare_reference.sh"
    echo "  qsub -F "-i \<REF_DIR\>/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fasta" <Path_to_script>/prepare_reference.sh"
    echo "----------------------------------------"
    exit 1
}

while getopts "i:" opt
do
   case "$opt" in 
      i ) REFERENCE="$OPTARG" ;;
      h ) helpFunction ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# check if sample name is empty
if [ -z "$REFERENCE" ]; then
    echo "---------------ERROR--------------------"
    echo "----------------------------------------"
    echo "Reference name parameter \${REFERENCE} is empty"
    echo "Please specify path \${REFERENCE} in cmd with option -i inputPath"
    helpFunction
fi

# check if REFERENCE file has on of the extensions .fasta, .fna, .ffn, .faa, .frn, .fa
if [[ $REFERENCE != *.fasta ]] && [[ $REFERENCE != *.fna ]] && [[ $REFERENCE != *.ffn ]] && [[ $REFERENCE != *.faa ]] && [[ $REFERENCE != *.frn ]] && [[ $REFERENCE != *.fa ]]; then
    echo "---------------ERROR--------------------"
    echo "----------------------------------------"
    echo "Reference file \${REFERENCE} has wrong extension"
    echo "Please specify path \${REFERENCE} with one of the following extensions: .fasta, .fna, .ffn, .faa, .frn, .fa"
    helpFunction
fi


echo "Index reference file: ${REFERENCE##*/}"
echo "Reference file directory: ${REFERENCE%/*}"
echo ""
echo "Start bwa index"

bwa index $REFERENCE

echo ""
echo "Start samtools faidx "
samtools faidx $REFERENCE --fai-idx ${REFERENCE%.*}.fai 

echo ""
echo "Create Sequence Dictionary with java -jar ${PICARD}"




java -jar ${PICARD} CreateSequenceDictionary \
-R $REFERENCE \
-O ${REFERENCE%.*}.dict






# if job is running on cluster should PBS_JOBID be non-empty
# Compatible with qsub from Torque and PBSPro
if [ ! -z "$PBS_JOBID" ]; then

    ### Output job resurces used
    echo "---------------------"
    echo "Job resources used:"
    /usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources

fi