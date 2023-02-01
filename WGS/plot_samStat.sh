
# ---
# --- About: Plot BAM file statistics "samtools stats" results
# ---       
# ---
# --- In: 
# ---       "samtools stats" result file 
# ---       reference stats file with expected GC content (created with plot-bamstats -s from goshawk reference)  
# ---      
# --- Out: 
# ---       Result folder with HTML, png plots and text files with generated statistics
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 28.11.22

#!/bin/bash
#PBS -N PlotsamtoolsStats
#PBS -l procs=4
#PBS -l walltime=5:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

# -- Displaying job information  -- #

# Calculate the number of processors and nodes allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`
echo Time is `date`
echo Running on host `hostname`
cat /proc/meminfo | grep MemTotal
echo cores used 'cat /proc/cpuinfo | grep processor | wc -l'
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes



FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combined'
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fa.gc'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combine_broken/BAM_stats'


module load bio/samtools/1.10


date 

for f in $(ls ${OUTPATH}/*_samStats.txt | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\/starpfaili\/combine_broken\/BAM_stats//'  | sed -e 's/_samStats.txt//' | sort -u)
do

    echo "analyze: ${FILEPATH}/${f}_samStats.txt "
    echo "Results will be saved at  ${OUTPATH}/${f}"

    echo "plot-bamstats: "
    plot-bamstats -r ${HG19FASTAPATH} -p ${OUTPATH}/${f}/ ${OUTPATH}/${f}_samStats.txt

done

date

### Output job resurces used
/usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources