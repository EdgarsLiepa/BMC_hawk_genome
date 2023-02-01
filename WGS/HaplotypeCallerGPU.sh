# ---
# --- About: 
# ---       Call haplotypes (SNP/Indels)
# ---       with Nvidia Clara Parabricks GPU acceleration toolkit 
# ---       
# ---
# --- In: 
# ---       SAMPLENAME - Sample name of bam file with exstention .markdup.fixedRG.bam
# ---      
# --- Out: 
# ---       VCF file with haplotypes - 
# ---          ${OUTPATH}/<sample name>.genome.raw.snps.indels.g.vcf
# ---       
# --- Example:
# ---       qsub -v SAMPLENAME="sample1" HaplotypeCallerGPU.sh
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 14.12.22


#!/bin/bash
#PBS -N Haplotypecaller_GPU
#PBS -l nodes=1:ppn=6:gpus=2,mem=100G
#PBS -l walltime=20:59:59
#PBS -l feature=V100
#PBS -A bmc_pl_bior_covid
#PBS -j oe


# Path to bam file
FPATH='/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp'


cat /proc/meminfo | grep MemTotal
cat /proc/cpuinfo | grep processor | wc -l

# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo Using ${NPROCS} processors across ${NNODES} nodes


# check if sample name is empty
if [ -z "$SAMPLENAME" ]; then
    echo "----------------------------------------"
    echo "Sample name parameter \${SAMPLENAME} is empty"
    echo "Please specify \${SAMPLENAME} in cmd with option -v ARG_NAME1=\"ARG_VALUE1\""
    echo "Exmple:"
    echo "  qsub -v SAMPLENAME="sample1" <Path to script>/HaplotypeCallerGPU.sh"
    echo "----------------------------------------"
    exit 1
fi

module load singularity
module load cuda

echo "Call variants for ${SAMPLENAME}.markdup.fixedRG.bam"

# Run haplotypecaller on singularity Nvidia clara-parabricks container
singularity run \
    --nv \
    /mnt/beegfs2/home/groups/bmc/test_parabricks/clara-parabricks_4.0.0-1.sif \
    pbrun haplotypecaller \
    --ref $HG19FASTAPATH \
    --in-bam ${FPATH}/${SAMPLENAME}.markdup.fixedRG.bam \
    --gvcf \
    --out-variants ${FPATH}/${SAMPLENAME}.genome.raw.snps.indels.gpu.g.vcf


if test -f "${FPATH}/vanagi.joint.genotype.full.output.vcf.gz"; then
    echo "VCF file ${FPATH}/vanagi.joint.genotype.full.output.vcf.gz exported"
else
    echo "Result file not found at ${FPATH}/"
fi

date