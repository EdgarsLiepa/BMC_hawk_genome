index-bwa.sh - tiek indeksēts references genoms, lai ir saderīgs ar bwa tālāk
#!/bin/bash
#PBS -N PL_BIOR_COVID
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -j oe
 
/home/groups/bmc/projects/Innas_Eksomi/nikita/bwa-0.7.17/bwa index /home_beegfs/ditagu/projekti/LZP_WGS/sheep_genome_consortium_ref/ncbi-genomes-2022-05-25/GCA_000298735.2_Oar_v4.0_genomic.fna
 
#!/bin/bash
#PBS -N PL_BIOR_COVID
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:10:00
#PBS -j oe

# pieinstalēju samtools 

/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools faidx /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna --fai-idx /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fai
 
 
#!/bin/bash
#PBS -N dict_ref
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -j oe

# šito vajag lejuplādēt  
java -jar /home_beegfs/ditagu/programs/picard.jar CreateSequenceDictionary \
-R /home_beegfs/ditagu/projekti/LZP_WGS/sheep_genome_consortium_ref/GCA_000298735.2_Oar_v4.0_genomic.fna \
-O /home_beegfs/ditagu/projekti/LZP_WGS/sheep_genome_consortium_ref/GCA_000298735.2_Oar_v4.0_genomic.dict
 
 
Mapēšana un co
#!/bin/bash
#PBS -l nodes=1:ppn=16,mem=64g
#PBS -l walltime=335:59:59
#PBS -j oe
#PBS -N map_reads
 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
source activate /home_beegfs/ditagu/miniconda3/envs/java11
 
#Mape, kur atrodas izejas sekvenatora outputs
FILEPATH='/home_beegfs/ditagu/projekti/LZP_WGS/aitas/SLTA_8340'
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/ditagu/projekti/LZP_WGS/aitas/SLTA_8340/starpfaili'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/ditagu/projekti/LZP_WGS/sheep_genome_consortium_ref/GCA_000298735.2_Oar_v4.0_genomic.fna'
#Faila pamatnosaukums
FILE='SLTA_8340'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp'
 
#Mapping reads to reference genome
/home/groups/bmc/projects/Innas_Eksomi/nikita/bwa-0.7.17/bwa mem -M -t 16 $HG19FASTAPATH \
${FILEPATH}/${FILE}_F_paired.fq.gz ${FILEPATH}/${FILE}_R_paired.fq.gz\
> ${OUTPATH}/${FILE}.mapped.sam;
 
# Convert SAM to BAM
/home_beegfs/ditagu/programs/samtools-1.14/samtools view -bS \
${OUTPATH}/${FILE}.mapped.sam > ${OUTPATH}/${FILE}.mapped.bam;
 
# fills in mate coordinates and insert size fields
/home_beegfs/ditagu/programs/samtools-1.14/samtools fixmate -m \
${OUTPATH}/${FILE}.mapped.bam ${OUTPATH}/${FILE}.fixmate.bam --threads 16;
 
# Sort BAM file
/home_beegfs/ditagu/programs/samtools-1.14/samtools sort \
-o ${OUTPATH}/${FILE}.fixmate.sorted.bam ${OUTPATH}/${FILE}.fixmate.bam;
 
# Index BAM file 
/home_beegfs/ditagu/programs/samtools-1.14/samtools index \
${OUTPATH}/${FILE}.fixmate.sorted.bam;
 
rm ${OUTPATH}/${FILE}.mapped.sam
rm ${OUTPATH}/${FILE}.mapped.bam
 
# mark duplicate alignments in a coordinate sorted file
/home_beegfs/ditagu/programs/samtools-1.14/samtools markdup -rs \
${OUTPATH}/${FILE}.fixmate.sorted.bam ${OUTPATH}/${FILE}.markdup.bam --threads 16;
 
rm ${OUTPATH}/${FILE}.fixmate.sorted.bam*
 
# index BAM file
/home_beegfs/ditagu/programs/samtools-1.14/samtools index \
${OUTPATH}/${FILE}.markdup.bam;
 
# calculate depth
/home_beegfs/ditagu/programs/samtools-1.14/samtools depth -a -d 0 \
${OUTPATH}/${FILE}.markdup.bam > ${OUTPATH}/${FILE}_samtools_depth.txt;
 
#Replacing readgroupinfo (probably gets lost)/unnecessary step, but conflicts downstream may arise
/home_beegfs/ditagu/programs/samtools-1.14/samtools addreplacerg \
-m overwrite_all \
-r ID:1 -r SM:SLTA_8340 -r PL:bgi \
-o ${OUTPATH}/${FILE}.markdup.fixedRG.bam \
${OUTPATH}/${FILE}.markdup.bam;
 
rm ${OUTPATH}/${FILE}.markdup.bam*
 
# index BAM file
/home_beegfs/ditagu/programs/samtools-1.14/samtools index \
${OUTPATH}/${FILE}.markdup.fixedRG.bam;
 
#Calling haplotypes (SNP/Indels)
# vajag java laikam 11 versiju
gatk --java-options "-Djava.io.tmpdir=$TMP" HaplotypeCaller \
--reference $HG19FASTAPATH \
--input ${OUTPATH}/${FILE}.markdup.fixedRG.bam \
-ERC GVCF \
-O ${OUTPATH}/${FILE}.genome.raw.snps.indels.g.vcf;
 

# 1. izveido genomics datubāzi dziļajām aitiņām
#!/bin/bash
#PBS -N GenomicsDBI_aita
#PBS -l procs=32
#PBS -l walltime=335:59:59
#PBS -j oe
 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
source activate /home_beegfs/ditagu/miniconda3/envs/java11
 
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/ditagu/projekti/LZP_WGS/aitas/VCF_ALL'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/ditagu/projekti/LZP_WGS/sheep_genome_consortium_ref/GCA_000298735.2_Oar_v4.0_genomic.fna'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'
 
gatk --java-options "-Xms26G -XX:ParallelGCThreads=2" GenomicsDBImport \
-V ${OUTPATH}/LTA1.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA2.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA3.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTA4.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Gusts.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Irbis.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Rausis.genome.raw.snps.indels.g.vcf \
-V ${OUTPATH}/LTT_Vucens.genome.raw.snps.indels.g.vcf \
--genomicsdb-workspace-path ${OUTPATH}/gatk.database_aitas \
--tmp-dir ${TMP} \
-L ${OUTPATH}/intervals.list \
--reader-threads 1;
 
chmod u=rwx,g=rx,o=r ${OUTPATH}/gatk.database_aitas;


# izveido kopēju VCF genotipu failu brūnajiem un zilajiem (kopējais)
#!/bin/bash
#PBS -N DB_vcf_deepseq
#PBS -l procs=6
#PBS -l walltime=95:59:59
#PBS -j oe
 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
source activate /home_beegfs/ditagu/miniconda3/envs/java11
 
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/ditagu/projekti/LZP_WGS/gotinas/VCF_ALL'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/ditagu/projekti/LZP_WGS/1000_bull_genomes/ARS-UCD1.2_Btau5.0.1Y.fa'
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp/ditagu'
 
gatk --java-options "-Djava.io.tmpdir=$TMP" GenotypeGVCFs \
 -R $HG19FASTAPATH \
 -V gendb://${OUTPATH}/gatk.database_deepseqonly \
 -O ${OUTPATH}/deepseqonly.joint.genotype.full.output.vcf.gz

