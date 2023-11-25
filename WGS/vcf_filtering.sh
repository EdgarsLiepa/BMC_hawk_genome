# ---
# --- About: 
# ---       
# ---
# --- In: 
# ---       
# ---      
# --- Out: 
# ---       
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 19.12.22

#!/bin/bash
#PBS -N filterVCF
#PBS -l procs=8,mem=10GB
#PBS -l walltime=120:59:59
#PBS -A bmc_pl_bior_covid
#PBS -q long
#PBS -j oe

VCF=/home_beegfs/edgars01/Ineta/WGS/starpfaili/sampleVCFs
OUT=/home_beegfs/edgars01/Ineta/WGS/starpfaili/filtered_Sample_VCFs
REZ=/home_beegfs/edgars01/Ineta/WGS/Results/bcfTools
REF=/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna

# if job is running on cluster should PBS_JOBID be non-empty
# Compatible with qsub from Torque and PBSPro
if [ ! -z "$PBS_JOBID" ]; then

    # -- Displaying job information
    echo Running on host `hostname`
    echo Time is `date`
    echo Current working directory is `pwd`
    echo "Node file: $PBS_NODEFILE :"
    echo "---------------------"
    cat $PBS_NODEFILE
    echo "---------------------"
    echo Using ${NPROCS} processors across ${NNODES} nodes
fi 

# load vcftools from HPC program library
module load bio/vcftools/0.1.17

for SAMPLENAME in $(ls ${VCF}/*.genome.raw.snps.indels.g.vcf | sed -e 's/.genome.raw.snps.indels.g.vcf//' | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/sampleVCFs\///'  | sort -u)
do

    echo "START ${SAMPLENAME}"
    date 

    # filter cromosome Z (OV839371.1) and W (OV839360.1) and mitochondia OV839400.1 
    if test ! -s "$OUT/${SAMPLENAME}_chr_filtered.recode.vcf"; then
        vcftools --gzvcf $VCF/${SAMPLENAME}.genome.raw.snps.indels.g.vcf --not-chr OV839371.1 --not-chr OV839360.1 --not-chr OV839400.1 --recode --out ${OUT}/${SAMPLENAME}_chr_filtered 
    fi
    # filter out samples with quality > 20
    if test ! -s "$OUT/${SAMPLENAME}_quality_filtered.recode.vcf"; then
        vcftools --vcf ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf --minQ 20 --recode --out ${OUT}/${SAMPLENAME}_quality_filtered 
    fi
    # filter out samples with deapth > 10
    if test ! -s "$OUT/${SAMPLENAME}_deapth_filtered.recode.vcf"; then
        vcftools --vcf ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf --minDP 10 --recode --out ${OUT}/${SAMPLENAME}_deapth_filtered 
    fi
    # filter out samples quality > 20 && deapth > 10
    
    if test ! -s "$OUT/${SAMPLENAME}_deapthAndQuality_filtered.recode.vcf"; then
        vcftools --vcf ${OUT}/${SAMPLENAME}_quality_filtered.recode.vcf --minDP 10 --recode --out ${OUT}/${SAMPLENAME}_deapthAndQuality_filtered 
    fi


    # bcftools stat on chr_filtered
    module load bio/bcftools/1.10.2

    
    if test -s "$REZ/${SAMPLENAME}_chr_filtered.stats"; then
        echo "VCF file statistics without intervals [OV839371.1 OV839360.1 OV839400.1] $REZ/${SAMPLENAME}_chr_filtered.stats already created"
    else
        bcftools stats -v -s- --threads 8 -F $REF ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf > $REZ/${SAMPLENAME}_chr_filtered.stats
    fi

    if test -s "$REZ/${SAMPLENAME}_quality_filtered"; then
        echo "VCF file statistics for minQ 20 filter  $REZ/${SAMPLENAME}_quality_filtered already created"
    else
        bcftools stats -v -s- --threads 8 -F $REF ${OUT}/${SAMPLENAME}_quality_filtered.recode.vcf > $REZ/${SAMPLENAME}_quality_filtered
    fi
    
    if test -s "$REZ/${SAMPLENAME}_deapth_filtered.stats"; then
        echo "VCF file statistics for minDP 10 filter $REZ/${SAMPLENAME}_deapth_filtered.stats already created"
    else
        bcftools stats -v -s- --threads 8 -F $REF ${OUT}/${SAMPLENAME}_deapth_filtered.recode.vcf > $REZ/${SAMPLENAME}_deapth_filtered.stats
    fi

    if test -s "$REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats"; then
        echo "VCF file statistics for filtered deapth and quality $REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats already created"
    else
        bcftools stats -v -s- --threads 8 -F $REF ${OUT}/${SAMPLENAME}_deapthAndQuality_filtered.recode.vcf > $REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats
    fi

    echo "${SAMPLENAME} done"


done

date

if [ ! -z "$PBS_JOBID" ]; then
    ### Output job resurces used
    echo "---------------------"
    echo "Job resources used:"
    /usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources
fi