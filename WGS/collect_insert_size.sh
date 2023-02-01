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
#PBS -l procs=4,mem=50G
#PBS -l walltime=44:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

# -- Switching to the directory from which the "qsub" command was run, moving to actual working directory
# -- Make sure any symbolic links are resolved to absolute path 
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)  ##ag
cd $PBS_O_WORKDIR
echo Working directory is: $PBS_O_WORKDIR
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


module load R

module load bio/samtools/1.10

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/sam_stat'

for f in $(ls ${FILEPATH}/*.fixmate.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\///'  | sort -u)
do

    # echo "Analyze $f "

    # java -jar /home_beegfs/edgars01/tools/picard.jar CollectInsertSizeMetrics -I ${FILEPATH}/$f -O ${FILEPATH}/${f}_insert_size_metrics.txt -H ${FILEPATH}/${f}_insert_size_histogram.pdf -M 0.5
    echo ${f} &>> DuplicatReadsCounts.txt

    echo "Reads: " &>> DuplicatReadsCounts.txt
    samtools view -c ${FILEPATH}/$f &>> DuplicatReadsCounts.txt

    echo "Mapped Reads: " DuplicatReadsCounts.txt
    samtools view -c -F 260 ${FILEPATH}/$f DuplicatReadsCounts.txt
done
