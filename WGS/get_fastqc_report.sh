
# ---
# --- About: Gather statistics about BAM files after 
# ---       
# ---
# --- In: 
# ---       Raw reads - .fq.gz 
# ---      
# --- Out: 
# ---       result summary HTML in ${FILEPATH}/fastqc folder
# ---       result summary .zip in ${FILEPATH}/fastqc folder
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 25.11.22

#!/bin/bash
#PBS -N count_reads
#PBS -l procs=1
#PBS -l walltime=42:00:00
#PBS -A bmc_pl_bior_covid
#PBS -j oe

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combined'

for f in $(ls ${FILEPATH}/* | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\/combined\///'  | sort -u)
do

echo $f
/home_beegfs/edgars01/tools/FastQC/fastqc --noextract --nogroup -o ${FILEPATH}/fastqc ${FILEPATH}/$f

done