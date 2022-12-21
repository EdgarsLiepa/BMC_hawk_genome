# ---
# --- About: Cut adapter sequences and trim.
# ---       
# ---
# --- In: 
# ---       Fastq file
# ---      
# --- Out: 
# ---       Trimmed fastq file
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 19.12.22

#!/bin/bash
#PBS -N filterVCF_V350018543
#PBS -l procs=8,mem=10GB
#PBS -l walltime=30:59:59
#PBS -A bmc_pl_bior_covid
#PBS -j oe

VCF=/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM
OUT=/home_beegfs/edgars01/Ineta/WGS/starpfaili/sample_vcf
REZ=/home_beegfs/edgars01/Ineta/WGS/Results/bcfTools

date

for SAMPLENAME in $(ls /home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM/V350018543*.genome.raw.snps.indels.g.vcf | sed -e 's/.genome.raw.snps.indels.g.vcf//' | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/mappedBAMandSAM\///'  | sort -u)
do

    echo "START ${SAMPLENAME}"

    # filter cromosome Z (OV839371.1) and W (OV839360.1) and mitochondia OV839400.1 
    /home_beegfs/edgars01/tools/vcftools_0.1.13/bin/vcftools --gzvcf $VCF/${SAMPLENAME}.genome.raw.snps.indels.g.vcf --not-chr OV839371.1 --not-chr OV839360.1 --not-chr OV839400.1 --recode --out ${OUT}/${SAMPLENAME}_chr_filtered 
    # filter out samples with quality > 20
    /home_beegfs/edgars01/tools/vcftools_0.1.13/bin/vcftools --vcf ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf --minQ 20 --recode --out ${OUT}/${SAMPLENAME}_quality_filtered 
    # filter out samples with deapth > 10
    /home_beegfs/edgars01/tools/vcftools_0.1.13/bin/vcftools --vcf ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf --minDP 10 --recode --out ${OUT}/${SAMPLENAME}_deapth_filtered 
    # filter out samples quality > 20 && deapth > 10
    /home_beegfs/edgars01/tools/vcftools_0.1.13/bin/vcftools --vcf ${OUT}/${SAMPLENAME}_quality_filtered.recode.vcf --minDP 10 --recode --out ${OUT}/${SAMPLENAME}_deapthAndQuality_filtered 



    # bcftools stat on chr_filtered
    module load bio/bcftools/1.10.2

    bcftools stats -v -s- --threads 8 -F /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna ${OUT}/${SAMPLENAME}_chr_filtered.recode.vcf > $REZ/${SAMPLENAME}_chr_filtered.stats
    if test -f "$REZ/${SAMPLENAME}_chr_filtered.stats"; then
        echo "VCF file statistics without intervals [OV839371.1 OV839360.1 OV839400.1] $REZ/${SAMPLENAME}_chr_filtered.stats created"
    fi
    bcftools stats -v -s- --threads 8 -F /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna ${OUT}/${SAMPLENAME}_quality_filtered.recode.vcf > $REZ/${SAMPLENAME}_quality_filtered
    if test -f "$REZ/${SAMPLENAME}_quality_filtered"; then
        echo "VCF file statistics for minQ 20 filter  $REZ/${SAMPLENAME}_quality_filtered created"
    fi
    bcftools stats -v -s- --threads 8 -F /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna ${OUT}/${SAMPLENAME}_deapth_filtered.recode.vcf > $REZ/${SAMPLENAME}_deapth_filtered.stats
    if test -f "$REZ/${SAMPLENAME}_deapth_filtered.stats"; then
        echo "VCF file statistics for minDP 10 filter $REZ/${SAMPLENAME}_deapth_filtered.stats"
    fi
    bcftools stats -v -s- --threads 8 -F /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna ${OUT}/${SAMPLENAME}_deapthAndQuality_filtered.recode.vcf > $REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats
    if test -f "$REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats"; then
        echo "VCF file statistics for filtered deapth and quality $REZ/${SAMPLENAME}_deapthAndQuality_filtered.stats created"
    fi

echo "${SAMPLENAME} done"

done

date


