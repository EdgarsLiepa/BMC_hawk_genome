#!/bin/bash
# nodes=1:ppn=12,mem=64g
#PBS -l nodes=1:ppn=12
#PBS -l walltime=48:59:59
#PBS -q long
#PBS -j oe
#PBS -N PL_BIOR_COVID_map_reads
 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
# # vai man šo vajag
# # source activate /home_beegfs/ditagu/miniconda3/envs/java11
 
#Mape, kur atrodas izejas sekvenatora outputs
FILEPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads'
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/raw_reads/starpfaili'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Faila pamatnosaukums
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp'

for f in $(ls /home_beegfs/edgars01/Ineta/WGS/raw_reads/V300082518_L02_84* | sed -e 's/_1.fq.gz//' | sed -e 's/_2.fq.gz//'| sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/raw_reads\///'  | sort -u)
do
date +'%Y-%m-%d'
date +'%r'


#Mapping reads to reference genome
/home/groups/bmc/projects/Innas_Eksomi/nikita/bwa-0.7.17/bwa mem -M -t 16 $HG19FASTAPATH \
${FILEPATH}/${f}_1.fq.gz ${FILEPATH}/${f}_2.fq.gz\
> ${OUTPATH}/${f}.mapped.sam;
 
# Convert SAM to BAM
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools view -bS \
${OUTPATH}/${f}.mapped.sam > ${OUTPATH}/${f}.mapped.bam;
 
# fills in mate coordinates and insert size fields
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools fixmate -m \
${OUTPATH}/${f}.mapped.bam ${OUTPATH}/${f}.fixmate.bam --threads 16;
 
# Sort BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools sort \
-o ${OUTPATH}/${f}.fixmate.sorted.bam ${OUTPATH}/${f}.fixmate.bam;
 
# Index BAM file 
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.fixmate.sorted.bam;
 
rm ${OUTPATH}/${f}.mapped.sam
rm ${OUTPATH}/${f}.mapped.bam
 
# mark duplicate alignments in a coordinate sorted file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools markdup -rs \
${OUTPATH}/${f}.fixmate.sorted.bam ${OUTPATH}/${f}.markdup.bam --threads 16;
 
rm ${OUTPATH}/${f}.fixmate.sorted.bam*
 
# index BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.markdup.bam;
 
# calculate depth
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools depth -a -d 0 \
${OUTPATH}/${f}.markdup.bam > ${OUTPATH}/${f}_samtools_depth.txt;
 
#Replacing readgroupinfo (probably gets lost)/unnecessary step, but conflicts downstream may arise
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools addreplacerg \
-m overwrite_all \
-r ID:1 -r SM:SLTA_8340 -r PL:bgi \
-o ${OUTPATH}/${f}.markdup.fixedRG.bam \
${OUTPATH}/${f}.markdup.bam;
 
rm ${OUTPATH}/${f}.markdup.bam*
 
# index BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.markdup.fixedRG.bam;
 
#Calling haplotypes (SNP/Indels)
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home_beegfs/edgars01/Ineta/WGS/tmp
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk HaplotypeCaller \
--reference $HG19FASTAPATH \
--input ${OUTPATH}/${f}.markdup.fixedRG.bam \
-ERC GVCF \
-O ${OUTPATH}/${f}.genome.raw.snps.indels.g.vcf;

date +'%Y-%m-%d'
date +'%r'

done