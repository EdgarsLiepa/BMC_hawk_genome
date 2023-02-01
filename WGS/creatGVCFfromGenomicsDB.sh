# ---
# --- About: 
# ---       Export GVCF file fom genomicsDB
# ---
# --- In: 
# ---        GenomicsDB 
# ---           
# ---      
# --- Out: 
# ---       A combined multi-sample gVCF. 
# ---           
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 19.10.22


#!/bin/bash
#PBS -N Export_GVCFs 
#PBS -l nodes=1:ppn=4,pmem=10,mem=40g
#PBS -l walltime=335:59:59
#PBS -A bmc_pl_bior_covid
#PBS -q long
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

#  --  Add file paths   -- #

#Mape, kur notiks visas starpdarbiibas
FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combine_broken/mapped'
FILEPATH2='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'
FILEPATHPARALLEL='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combine_broken/mapped_parallel'
DBPATH='/home_beegfs/edgars01/Ineta/WGS/CombinedGVCFs_Vanagi'
INTERVALLIST='/home_beegfs/edgars01/Ineta/WGS'
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'



/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk SelectVariants \
    -R ${HG19FASTAPATH} \
    -V gendb://${INTERVALLIST}/gatk.database_vanagi_ALL \
    -O ${DBPATH}/database_vanagi_ALL_cohort.g.vcf.gz



date

### Output job resurces used
/usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources
