
# ---
# --- About: Gather statistics about BAM files after 
# ---       
# ---
# --- In: 
# ---       .markdup.fixedRG.bam
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

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/sam_stat'


module load bio/samtools/1.10


echo "Results will be saved at  ${OUTPATH}"

for f in $(ls ${FILEPATH}/*.markdup.fixedRG.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\///'  | sed -e 's/.markdup.fixedRG.bam//' | sort -u)
do

    date 
    echo "Start mapping: $f"
    echo "analyze: ${FILEPATH}/$f.markdup.fixedRG.bam "

    echo "Samtools stats: "
    samtools stats -r ${HG19FASTAPATH} ${FILEPATH}/$f.markdup.fixedRG.bam > ${OUTPATH}/${f}_samStats.txt

    echo "Finished mapping: $f"
    date 
done

