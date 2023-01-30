
# ---
# --- About: BAM files statistics with samtools stats
# ---       
# ---
# --- In: 
# ---       bam file
# ---      
# --- Out: 
# ---       
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 25.11.22

#!/bin/bash
#PBS -N samtoolsStats
#PBS -l procs=4
#PBS -l walltime=94:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM'
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/Results/sam_stat'


module load bio/samtools/1.10


echo "Results will be saved at  ${OUTPATH}"

for f in $(ls ${FILEPATH}/*.fixmate.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/mappedBAMandSAM\///'  | sed -e 's/.fixmate.bam//' | sort -u)
do

    date 
    echo "Start: $f"
    echo "analyze: ${FILEPATH}/$f.fixmate.bam "

    echo "Samtools stats: "
    samtools stats -r ${HG19FASTAPATH} ${FILEPATH}/$f.fixmate.bam > ${OUTPATH}/${f}_samStats_withDuplicates.txt
    
    echo "Samtools stats with no duplicates: "
    samtools stats -r ${HG19FASTAPATH} ${FILEPATH}/$f.markdup.fixedRG.bam > ${OUTPATH}/${f}_samStats.txt

    echo "Finished: $f"
    
done

