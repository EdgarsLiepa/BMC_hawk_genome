
# ---
# --- About: Calculate Read depth, coverage, proportion of mapped reads
# ---       
# ---
# --- In: 
# ---       .markdup.fixedRG.bam
# ---      
# --- Out: 
# ---       Read mapping statistics in standart Output
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 24.11.22

#!/bin/bash
#PBS -N read_mapping_stat
#PBS -l procs=1
#PBS -l walltime=94:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'

module load bio/samtools/1.10

for f in $(ls ${FILEPATH}/*.markdup.fixedRG.bam sort -u)
do

    echo "$f "

    echo "Mean Read Depth: "

    samtools depth -a $f| awk '{c++;s+=$3}END{print s/c}'

    echo "Breadth of Coverage: "

    samtools depth -a $f | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'

    echo "Proportion of the Reads that Mapped to the Reference: "

    samtools flagstat $f | awk -F "[(|%]" 'NR == 3 {print $2}'

done