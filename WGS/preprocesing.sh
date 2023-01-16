# ---
# --- About: Cut adapter sequences and trim.
# ---       
# ---
# --- In: 
# ---       Fastq file
# ---      
# --- Out: 
# ---       Trimmed fastq file
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 14.12.22

#!/bin/bash
#PBS -l nodes=1:ppn=16,mem=25G
#PBS -l walltime=240:59:59
#PBS -j oe
#PBS -N preProcess 
#PBS -A bmc_pl_bior_covid
#PBS -q long


# path to in file directory
FPATH="/home_beegfs/edgars01/Ineta/WGS/starpfaili/combined"
# `path to save results in
OUTPATH="/home_beegfs/edgars01/Ineta/WGS/starpfaili/trimmed_sequences"

ADAPTER_FWD="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
ADAPTER_REV="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"



cat /proc/meminfo | grep MemTotal
cat /proc/cpuinfo | grep processor | wc -l

# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes

echo "Start PreProcessing: $SAMPLE_ID"


    
    module load bio/cutadapt/3.5


for SAMPLE_ID in $(ls /home_beegfs/edgars01/Ineta/WGS/starpfaili/combined/*.fq.gz | sed -e 's/_1.fq.gz//' -e 's/_2.fq.gz//' | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/combined\///'  | sort -u)
do
    echo "Start CutAdapt ${SAMPLE_ID} $date"
    cutadapt --cores=16 -a ${ADAPTER_FWD} -A ${ADAPTER_REV} -o ${OUTPATH}/${SAMPLE_ID}_1_cutadapt.fq.gz -p ${OUTPATH}/${SAMPLE_ID}_2_cutadapt.fq.gz ${FPATH}/${SAMPLE_ID}_1.fq.gz ${FPATH}/${SAMPLE_ID}_2.fq.gz
    
    echo "Start trimmomatic ${SAMPLE_ID} $date"
    java -jar /home_beegfs/edgars01/tools/trimmomatic/trimmomatic-0.39.jar PE -threads 16 ${OUTPATH}/${SAMPLE_ID}_1_cutadapt.fq.gz ${OUTPATH}/${SAMPLE_ID}_2_cutadapt.fq.gz ${OUTPATH}/${SAMPLE_ID}_F_paired.fq.gz ${OUTPATH}/${SAMPLE_ID}_F_unpaired.fq.gz ${OUTPATH}/${SAMPLE_ID}_R_paired.fq.gz ${OUTPATH}/${SAMPLE_ID}_R_unpaired.fq.gz LEADING:30 TRAILING:30 MINLEN:36

    echo "Preprocesing ${SAMPLE_ID} done! $date"   
    

done;


date