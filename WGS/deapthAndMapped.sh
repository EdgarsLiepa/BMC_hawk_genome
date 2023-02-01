
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

# -- Load modules

module load bio/samtools/1.10

# -- Set File paths

FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'


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