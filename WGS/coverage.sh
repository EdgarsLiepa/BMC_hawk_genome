
# ---
# --- About: Calculate genome coverage
# ---       
# ---
# --- In: 
# ---       BAM file with marked duplicates - .markdup.fixedRG.bam
# ---      
# --- Out: 
# ---       
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 24.11.22

#!/bin/bash
#PBS -N coverage
#PBS -l procs=1
#PBS -l walltime=94:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'

module load bio/bedtools/2.29.2

for f in $(ls ${FILEPATH}/*.markdup.fixedRG.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\///'  | sed -e 's/.markdup.fixedRG.bam//' | sort -u) 
do

bedtools genomecov -ibam ${FILEPATH}/$f.markdup.fixedRG.bam -pc > ${FILEPATH}/${f}_bedtools.txt

done