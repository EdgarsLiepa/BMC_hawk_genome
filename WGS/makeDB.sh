# ---
# --- About: 
# ---       Create genomics DB from gosh hawk mapped read VCF files
# ---
# --- In: 
# ---        mapped reads - __.genome.raw.snps.indels.g.vcf
# ---           
# ---      
# --- Out: 
# ---       gatk genomics Data Base - gatk.database_vanagi_ALL
# ---           
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 19.10.22


#!/bin/bash
#PBS -N makeVanagiDB_ALL
#PBS -l nodes=1:ppn=32,mem=180g
#PBS -l walltime=335:59:59
#PBS -A bmc_pl_bior_covid
#PBS -q long
#PBS -j oe
 
date

#Mape, kur notiks visas starpdarbiibas
FILEPATH='/home_beegfs/edgars01/Ineta/WGS'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'
 
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home_beegfs/edgars01/Ineta/WGS/tmp
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk --java-options "-Xms64G -Xmx170g -XX:ParallelGCThreads=32" GenomicsDBImport \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L01_ACG1.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L01_ACG2.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L01_ACG3.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L01_ACG4.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L01_ACG5.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L02_AV13_1.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L02_AV13_2.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L02_AV13_3.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L02_F75_2.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V300082518_L02_F78.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L01_ACG22.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L01_ACG23.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L01_ACG24.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L01_ACG25.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L02_ACG10.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L02_ACG11.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L02_ACG12.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L02_ACG13.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L03_ACG18.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L03_ACG19.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L03_ACG20.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L03_ACG21.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L04_ACG14.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L04_ACG15.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018459_L04_ACG17.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L03_AV-5.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L03_AV-6.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L03_ZL-17.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L04_ACG6.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L04_ACG7.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L04_ACG8.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350018543_L04_ACG9.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023344_L01_ACG-WGS.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L01_ACG26.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L01_ACG27.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L01_ACG28.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L01_ACG29.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L02_ACG30.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L02_ACG31.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L02_ACG32.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L02_ACG33.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L03_ACG46.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L03_ACG47.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L03_ACG48.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L03_ACG49.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L04_ACG42.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L04_ACG43.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023369_L04_ACG45.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L01_ACG34.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L01_ACG35.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L01_ACG36.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L01_ACG37.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L02_ACG38.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L02_ACG39.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L02_ACG40.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L02_ACG41.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L03_ACG54.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L03_ACG55.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L03_ACG56.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L03_ACG57.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L04_ACG50.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L04_ACG51.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L04_ACG52.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023385_L04_ACG53.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L01_ACG71.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L01_ACG72.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L01_ACG73.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L01_ACG74.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L02_ACG59.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L02_ACG60.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L02_ACG61.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L02_ACG62.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L03_ACG67.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L03_ACG68.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L03_ACG69.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L03_ACG70.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L04_ACG63.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L04_ACG64.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L04_ACG65.genome.raw.snps.indels.g.vcf \
-V ${FILEPATH}/raw_reads/starpfaili/V350023490_L04_ACG66.genome.raw.snps.indels.g.vcf \
--genomicsdb-workspace-path ${FILEPATH}/gatk.database_vanagi_ALL \
--tmp-dir ${TMP} \
-L ${FILEPATH}/intervals.list \
--reader-threads 32;
 
chmod u=rwx,g=rx,o=r ${FILEPATH}/gatk.database_vanagi_ALL;

date


