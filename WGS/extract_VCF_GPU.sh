# ---
# --- About: 
# ---       Extract VCF file for gosh hawk from genomics DB
# ---       This script uses GPU acceleration with Nvidia Clara Parabricks toolkit 
# ---
# --- In: 
# ---       GenomicsDB workspace created by GenomicsDBImport:
# ---           gatk.database_vanagi
# ---      
# --- Out: 
# ---       vcf.gz file: 
# ---           VCF_ALL/vanagi.genotype.vcf.gz
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 30.11.22
# --- Modified: 26.06.25


#!/bin/bash
#PBS -N Export_VCF_GPU
#PBS -l nodes=1:ppn=32:gpus=2,mem=512G
#PBS -l feature=l40s
#PBS -q long
#PBS -l walltime=320:59:59
#PBS -A bmc_flpp_2022_0299
#PBS -j oe

cd $PBS_O_WORKDIR


cat /proc/meminfo | grep MemTotal
cat /proc/cpuinfo | grep processor | wc -l

# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes

    
echo "Export VCF file from GenomicsDB with L40 GPU and 32 cores"
date

module load singularity
module load cuda/cuda-12.4 


nvidia-smi

DBPATH='/home_beegfs/edgars01/Ineta/WGS/Results/CombinedVCF'
# OUTPATH='/home_beegfs/edgars01/Ineta/WGS/Results/CombinedVCF'

HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'

# export VCF from DB 
singularity run \
    --nv \
    /home_beegfs/edgars01/tools/singularity_containers/clara-parabricks_4.4.0-1.sif \
    pbrun genotypegvcf \
    --num-threads 32 \
    --ref $HG19FASTAPATH \
    --in-gvcf ${DBPATH}/vanagi_cohort.g.vcf.gz \
    --out-vcf ${DBPATH}/vanagi_26_06_25.vcf


if test -f "${DBPATH}/vanagi_26_06_25.vcf.gz"; then
    echo "VCF file ${DBPATH}/vanagi_26_06_25.vcf.gz created"
fi

date
