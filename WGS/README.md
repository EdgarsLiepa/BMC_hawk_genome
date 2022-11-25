# Antiobiotic resistance detection

## Description

Assembly of Goshawk genome against reference sequence bAccGen1.1 from short read sequences.

**Reference used** - [bAccGen1.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_929443795.1/) 

## Files

combineFiles.sh - combine raw reads files from sequencer to correspond to sample and add corect name, because sequencing reads from same sample where distributed in a multiple places on a sequecing chip. 

collect_insert_size.sh - CollectInsertSizeMetrics and get read counts from BAM file

coverage.sh - calculate genome coverage with bedtools

deapthAndMapped.sh - Calculate Read depth, coverage, proportion of mapped reads

get_fastqc_report.sh - Create FastQC report from BAM files with marked duplicates

gatk.databse_vanagi - create GenomicsDB from g.vcf files

extractVCF.sh - create VCF file from GenomicsDB

makeDB.sh - create Genomcs DB from VCF files

mapping.sh - map reads to reference, m ate coordinates, sort, index, mark duplicates, calculate deapth, index and call haplotypes. 

index_Reference.sh - index Reference sequence

intervals.list - file for RGI #todo check content

raw_reads - folder with raw read data

raw_reads/starpfaili - folder with intermidiate result files (BAM, SAM, VCF)

samT_stats.sh - gather statistics with *samtools satats* program

WGS_raw_pipeline.sh

## Dependencies

- bwa
- gatk
- samtools
- picard
- bedtools
- fastqc

## TODO

- [ ] Results folder where analysis main results are compiled would be nice.
- [ ] Describe analysis pipeline