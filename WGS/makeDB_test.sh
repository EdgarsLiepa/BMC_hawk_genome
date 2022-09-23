#!/bin/bash
#!/usr/bin/env python3
#PBS -N makeVanagiDB
#PBS -l nodes=1:ppn=16
#PBS -l walltime=59:59:59
#PBS -A bmc_pl_bior_covid
#PBS -q long
#PBS -j oe

date

#Mape, kur notiks visas starpdarbiibas
DIR='/home_beegfs/edgars01/Ineta/WGS'
#Temporary dir
TEMP='/home_beegfs/groups/bmc/tmp/ditagu'
 
# export _JAVA_OPTIONS=-Djava.io.tmpdir=/home_beegfs/edgars01/Ineta/WGS/tmp
cp -r ${DIR}/gatk.database_vanagi ${DIR}/gatk.database_vanagi_backUp2 

/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk --java-options "-Xms26G -XX:ParallelGCThreads=16 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
-V ${DIR}/raw_reads/starpfaili/V350018543_L04_ACG6.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018543_L04_ACG7.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018543_L04_ACG8.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018543_L03_AV-5.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018543_L03_AV-6.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018543_L03_ZE-6.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L03_ACG20.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L03_ACG21.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L01_ACG22.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L01_ACG23.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L01_ACG24.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L01_ACG25.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L02_ACG10.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L02_ACG11.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L02_ACG12.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L02_ACG13.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L04_ACG14.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L04_ACG15.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L04_ACG17.genome.raw.snps.indels.g.vcf \
-V ${DIR}/raw_reads/starpfaili/V350018459_L03_ACG19.genome.raw.snps.indels.g.vcf \
--genomicsdb-update-workspace-path ${DIR}/gatk.database_vanagi \
--tmp-dir ${TEMP} \
-L ${DIR}/intervals.list \
--reader-threads 16;
 
# -V ${DIR}/raw_reads/starpfaili/V350018459_L04_ACG16.genome.raw.snps.indels.g.vcf \


chmod u=rwx,g=rx,o=r ${DIR}/gatk.database_vanagi;

date