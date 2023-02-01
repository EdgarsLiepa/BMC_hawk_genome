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


#!/bin/bash
#PBS -N Export_VanagiDB
#PBS -l nodes=1:ppn=6,mem=100g
#PBS -l walltime=329:59:59
#PBS -A bmc_pl_bior_covid
#PBS -W x=HOSTLIST:wn62,wn63,wn64
#PBS -q long
#PBS -j oe
 
echo "Export VCF file from GenomicsDB"
date




DBPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combine_broken/Vanagi_GenomicsDB'
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili/combine_broken/GVCF_out'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'


# export VCF from DB 
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk --java-options "-Xms2G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=$TMP -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs \
-R $HG19FASTAPATH \
-V gendb://${DBPATH} \
-O ${OUTPATH}/vanagi.joint.genotype.full.output.vcf.gz




if test -f "${OUTPATH}/vanagi.joint.genotype.full.output.vcf.gz"; then
    echo "VCF file ${OUTPATH}/vanagi.joint.genotype.full.output.vcf.gz exported"
fi

date