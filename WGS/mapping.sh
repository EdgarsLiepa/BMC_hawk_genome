# ---       Create Haplotaype (SNP/Indels) Genomic Variant Call file from sequencing reads
# ---
# --- About: 
# ---       map reads to reference, mate coordinates, sort, index, mark duplicates, 
# ---       calculate deapth, index and call haplotypes. 
# ---       
# ---
# --- In: 
# ---        Forward and reverse read FastQ files - _1.fq.gz;  _2.fq.gz
# ---           
# ---      
# --- Out: 
# ---        mate coordinates and size fields - .fixmate.bam
# ---        Haplotaype (SNP/Indels) Genomic Variant Call file - .genome.raw.snps.indels.g.vcf
# ---        Haplotaype (SNP/Indels) Genomic Variant Call indexes - .genome.raw.snps.indels.g.vcf.idx
# ---        Duplicate alignments with updated read group info - .markdup.fixedRG.bam
# ---        index file (BAI) for duplicate alignment BAM file - .markdup.fixedRG.bam.bai
# ---        Depth at each position or region. _samtools_depth.txt
# ---           
# ---       
# --- Author: Edgars Liepa edgars.liepa@biomed.lu.lv
# --- Date: 31.10.22

#!/bin/bash
#PBS -l nodes=1:ppn=20,mem=120G
#PBS -l walltime=320:59:59
#PBS -q long
#PBS -j oe
#PBS -N map_reads_V300082518_L02_F78
#PBS -A bmc_pl_bior_covid

 
module load conda
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8

# -- Set file paths -- #

#Mape, kur atrodas izejas sekvenatora outputs
FILEPATH="/home_beegfs/edgars01/Ineta/WGS/starpfaili/trimmed_sequences"
#Mape, kur notiks visas starpdarbiibas
OUTPATH='/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM'
OUTPATH_VCF='/home_beegfs/edgars01/Ineta/WGS/starpfaili/sampleVCFs'
#hg19 references genoms
HG19FASTAPATH='/home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna'
#Faila pamatnosaukums
#Temporary dir
TMP='/home_beegfs/groups/bmc/tmp'


# -- Displaying job information

cat /proc/meminfo | grep MemTotal
cat /proc/cpuinfo | grep processor | wc -l

echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo Using ${NPROCS} processors across ${NNODES} nodes


# -- Start processing -- #

for f in $(ls ${FILEPATH}/V300082518_L02_F78*_paired.fq.gz | sed -e 's/_F_paired.fq.gz//' | sed -e 's/_R_paired.fq.gz//'| sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/trimmed_sequences\///'  | sort -u)
do

date
echo "Start mapping: $f"

#Mapping reads to reference genome
/home/groups/bmc/projects/Innas_Eksomi/nikita/bwa-0.7.17/bwa mem -M -t 16 $HG19FASTAPATH \
${FILEPATH}/${f}_F_paired.fq.gz ${FILEPATH}/${f}_R_paired.fq.gz\
> ${OUTPATH}/${f}.mapped.sam;

echo "Convert SAM to BAM $date"

# Convert SAM to BAM
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools view -bS \
${OUTPATH}/${f}.mapped.sam > ${OUTPATH}/${f}.mapped.bam;
 
echo "fills in mate coordinates and insert size $date"

# fills in mate coordinates and insert size fields
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools fixmate -m \
${OUTPATH}/${f}.mapped.bam ${OUTPATH}/${f}.fixmate.bam --threads 16;
 
echo "Sort and Index BAM $date"
 
# Sort BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools sort \
-o ${OUTPATH}/${f}.fixmate.sorted.bam ${OUTPATH}/${f}.fixmate.bam;
 
# Index BAM file 
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.fixmate.sorted.bam;
 
rm ${OUTPATH}/${f}.mapped.sam
rm ${OUTPATH}/${f}.mapped.bam
 
echo "Mark Duplicates $date"

# mark duplicate alignments in a coordinate sorted file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools markdup -rs \
${OUTPATH}/${f}.fixmate.sorted.bam ${OUTPATH}/${f}.markdup.bam --threads 16;
 
rm ${OUTPATH}/${f}.fixmate.sorted.bam*
 
echo "Index and calculate deapth $date"

# index BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.markdup.bam;
 
# calculate depth
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools depth -a -d 0 \
${OUTPATH}/${f}.markdup.bam > ${OUTPATH}/${f}_samtools_depth.txt;
 
echo "Replace readgroup info $date"

# Replace readgroupinfo. Add file ID to SM field
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools addreplacerg \
-m overwrite_all \
-r ID:1 -r SM:${f} -r PL:bgi \
-o ${OUTPATH}/${f}.markdup.fixedRG.bam \
${OUTPATH}/${f}.markdup.bam;
 
rm ${OUTPATH}/${f}.markdup.bam*
 
# index BAM file
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools index \
${OUTPATH}/${f}.markdup.fixedRG.bam;
 
echo "Call haplotypes (SNP/Indels) $date"

#Calling haplotypes (SNP/Indels)
export _JAVA_OPTIONS=-Djava.io.tmpdir=/home_beegfs/edgars01/Ineta/WGS/tmp
/home_beegfs/edgars01/tools/gatk-4.2.6.1/gatk HaplotypeCaller \
--reference $HG19FASTAPATH \
--input ${OUTPATH}/${f}.markdup.fixedRG.bam \
-ERC GVCF \
-O ${OUTPATH_VCF}/${f}.genome.raw.snps.indels.g.vcf;

date 
echo "Finished mapping: $f"

done