# ---
# --- About: 
# ---       CollectInsertSizeMetrics and get read counts from BAM file
# ---
# --- In: 
# ---       BAM file with marked duplicats .markdup.fixedRG.bam
# ---      
# --- Out: 
# ---       
# ---       Results txt - <Sample name>_insert_size_metrics.txt 
# ---       PDF histogramm - <Sample name>_insert_size_histogram.pdf
# ---       
# ---       Sequence count Results are printed in standart out.
# ---       TODO: add print to file.
# ---
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 22.11.22

#!/bin/bash
#PBS -N collect_insert_size
#PBS -l procs=4
#PBS -l walltime=44:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

module load R

module load bio/samtools/1.10

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'

for f in $(ls ${FILEPATH}/*.markdup.fixedRG.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\///'  | sort -u)
do

    echo "Analyze $f "

    java -jar /home_beegfs/edgars01/tools/picard.jar CollectInsertSizeMetrics -I ${FILEPATH}/$f -O ${FILEPATH}/${f}_insert_size_metrics.txt -H ${FILEPATH}/${f}_insert_size_histogram.pdf -M 0.5

    echo "Reads: "
    samtools view -c ${FILEPATH}/$f

    echo "Mapped Reads: "
    samtools view -c -F 260 ${FILEPATH}/$f
done
