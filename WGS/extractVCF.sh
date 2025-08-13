# ---
# --- About: 
# ---       Extract VCF file for gosh hawk from genomics DB
# ---       Requires large ammount of memory for large sample set
# ---       4 example - 100 GB of ram was used for 80 samples that was run for 8 days
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
# --- Date: 19.10.22
# --- Modified: 14-07-25


#!/bin/bash
#PBS -N Export_Vanagi
#PBS -l nodes=1:ppn=12,mem=512g
#PBS -l walltime=335:59:59
#PBS -A bmc_flpp_2024_0109
#PBS -W x=HOSTLIST:wn67,wn68,wn69,wn70
#PBS -q long
#PBS -j oe
 
echo "Export VCF file from GenomicsDB"
date




#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'


SCRATCH=/scratch/$PBS_JOBID
mkdir -m 700 $SCRATCH

# copy data to scratch
echo "Copy ref"
cp "/mnt/beegfs2/home/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna" $SCRATCH/GCA_929443795.1_bAccGen1.1_genomic.fna
echo "Copy g.vcf"
cp /home_beegfs/edgars01/Ineta/WGS/Results/CombinedVCF/vanagi_cohort.g.vcf.gz $SCRATCH/vanagi_cohort.g.vcf.gz


date

echo "run gatk"
# export VCF from DB 
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk --java-options "-Xms20G -Xmx512G -XX:ParallelGCThreads=12 -Djava.io.tmpdir=$TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
-R $HG19FASTAPATH \
-V /home_beegfs/edgars01/Ineta/WGS/Results/CombinedVCF/vanagi_cohort.g.vcf.gz \
-O $SCRATCH/vanagi_genotype.vcf \
--tmp-dir $SCRATCH 

cp -r $SCRATCH/vanagi_genotype.vcf $HOME/Ineta/WGS/Results/CombinedVCF/
rm -rf $SCRATCH


if test -f "$HOME/Ineta/WGS/Results/CombinedVCF/vanagi_genotype.vcf"; then
    echo "VCF file $HOME/Ineta/WGS/Results/CombinedVCF/vanagi_genotype.vcf exported"
fi

date
