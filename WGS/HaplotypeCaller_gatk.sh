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
#PBS -N Haplotypecaller_gatk
#PBS -l nodes=1:ppn=12,mem=100G
#PBS -l walltime=120:59:59
#PBS -A bmc_pl_bior_covid
#PBS -q long
#PBS -j oe


# Path to bam file
FPATH='/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/starpfaili/sampleVCFs_with_default_caller'
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

FILEPATH=/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM

module load singularity
module load cuda

# sed filepath from sample name 
sed -i 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/mappedBAMandSAM\///g' ${FILEPATH}/V350018459*.markdup.fixedRG.bam

for SAMPLENAME in $(ls ${FILEPATH}/V350018459*.markdup.fixedRG.bam |sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/mappedBAMandSAM\///'| sed -e 's/.markdup.fixedRG.bam//' | sort -u)
do

echo "Call variants for ${SAMPLENAME}.markdup.fixedRG.bam"

# Run haplotypecaller on singularity Nvidia clara-parabricks container
#Calling haplotypes (SNP/Indels)
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home_beegfs/edgars01/Ineta/WGS/tmp
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk HaplotypeCaller \
--reference $HG19FASTAPATH \
--input ${FPATH}/${SAMPLENAME}.markdup.fixedRG.bam \
-ERC GVCF \
-O ${OUTPATH}/${SAMPLENAME}.genome.raw.snps.indels.g.vcf;



if test -f "${OUTPATH}/${SAMPLENAME}.genome.raw.snps.indels.g.vcf"; then
    echo "VCF file ${SAMPLENAME}.genome.raw.snps.indels.g.vcf exported at ${OUTPATH}"
else
    echo "Result file not found at ${OUTPATH}/"
fi

done

date