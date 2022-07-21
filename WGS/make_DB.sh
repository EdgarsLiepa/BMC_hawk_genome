#!/bin/bash
#PBS -N PL_BIOR_COVID_GenomicsDBI_vanagi
#PBS -l procs=32
#PBS -l walltime=35:59:59
#PBS -j oe
 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
source activate /home_beegfs/ditagu/miniconda3/envs/java11
 
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/VCF_ALL'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'
 
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk --java-options "-Xms26G -XX:ParallelGCThreads=2" GenomicsDBImport \
-V ${OUTPATH}/LTA1.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA2.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA3.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA4.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Gusts.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Irbis.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Rausis.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Vucens.genome.raw.snps.indels.g.vcf \
--genomicsdb-workspace-path ${OUTPATH}/gatk.database_vanagi \
--tmp-dir ${TMP} \
-L ${OUTPATH}/intervals.list \
--reader-threads 1;
 
chmod u=rwx,g=rx,o=r ${OUTPATH}/gatk.database_vanagi;
