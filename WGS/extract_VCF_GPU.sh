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


#!/bin/bash
#PBS -N Export_VanagiVCF_GPU
#PBS -l nodes=1:ppn=6:gpus=2,mem=100G
#PBS -l walltime=20:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

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

    
echo "Export VCF file from GenomicsDB with v100 GPU and 6 cores"
date

module load singularity
module load cuda


DBPATH='/home_beegfs/edgars01/Ineta/WGS/CombinedGVCFs_Vanagi'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/VCF_ALL_GPU'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'

# export VCF from DB 
singularity run \
    --nv \
    /mnt/beegfs2/home/groups/bmc/test_parabricks/clara-parabricks_4.0.0-1.sif \
    pbrun genotypegvcf \
    --ref $HG19FASTAPATH \
    --in-gvcf ${DBPATH}/vanagi_cohort.g.vcf.gz \
    --out-vcf ${OUTPATH}/vanagi.joint.genotype.full.output_8_12.vcf.gz


if test -f "${OUTPATH}/vanagi.joint.genotype.full.output.vcf.gz"; then
    echo "VCF file ${OUTPATH}/vanagi.joint.genotype.full.output.vcf.gz exported"
fi

date