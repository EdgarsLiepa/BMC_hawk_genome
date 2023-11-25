# Variant Calling using short read sequences against reference

Author: Edgars Liepa  
email: edgars.liepa@biomed.lu.lv  
Latvian Biomedical Research and Study Center  

# !!! Important

Scripts have defind static file and dependency paths. If you are planing to use any of scripts at this development stage 
beware that you need to change them to corespond to your needs. 

## Description

Assembly of Goshawk genome against reference sequence bAccGen1.1 from short read sequences.

Analysis was done on RTU HPC that uses Cent OS and Torque resource maneger. 

**Reference used** - [bAccGen1.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_929443795.1/) 

### Pipeline

#### [PreProcessing](docs/preprocess.md)

![PreProcessing](docs/preProcessing.png)


![Mapping](docs/mapping.png)

![Variant Calling]()

![Results]()


## Files

combineFiles.sh - combine raw reads files from sequencer to correspond to sample and add corect name, because sequencing reads from same sample where distributed in a multiple places on a sequecing chip. 

collect_insert_size.sh - CollectInsertSizeMetrics and get read counts from BAM file

coverage.sh - calculate genome coverage with bedtools

deapthAndMapped.sh - Calculate Read depth, coverage, proportion of mapped reads

get_fastqc_report.sh - Create FastQC report from BAM files with marked duplicates

gatk.databse_vanagi - create GenomicsDB from g.vcf files

extractVCF.sh - create VCF file from GenomicsDB

extract_VCF_GPU.sh - create VCF file from combined GVCF file using NVIDIA Parabricks GPU acceleration. 

makeDB.sh - create Genomcs DB from VCF files

mapping.sh - map reads to reference, m ate coordinates, sort, index, mark duplicates, calculate deapth, index and call haplotypes. 

prepare_reference.sh - index fasta Reference sequence and create sequence dictionary.

Provide path to refence fasta file.
Acceptable extension names: .fasta, .fna, .ffn, .faa, .frn, .fa

~~~

REFERENCE - path to reference file in fasta format.

Acceptable file extension names: .fasta, .fna, .ffn, .faa, .frn, .fa

# submit job using bash

sh prepare_reference.sh -i /path_to_reference_files/ref.fasta

sh prepare_reference.sh -i /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05_test/GCA_929443795.1_bAccGen1.1_genomic.fasta

# submit job using qsub Torque Resource/Queue Manager

qsub -F "-i ~/path_to_reference_files/ref.fasta" <Path_to_script>/index_reference.sh"

qsub -F "-i /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05_test/GCA_929443795.1_bAccGen1.1_genomic.fasta"  /home_beegfs/edgars01/Ineta/WGS/index_reference.sh

~~~

intervals.list - file for RGI #todo check content

parseFASTQC.scv - Gather results from FASTQC result HTML files and create csv results file.

raw_reads - folder with raw read data

raw_reads/starpfaili - folder with intermidiate result files (BAM, SAM, VCF)

samT_stats.sh - gather statistics with *samtools satats* program

vcf_filtering.sh - filter chromosomes, low deapth and quality. Create statistics. 

qualimap.sh - create BAM file statistics using *qualimap bamqc*

WGS_raw_pipeline.sh

## Dependencies

- bwa
- gatk
- samtools
- picard
- fastqc
- bedtools
    - bedtools genomecov
- NVIDIA parabrics
- vcfstats
- vcftools
- qualimap bamqc

## TODO



- [ ] Results folder where analysis main results are compiled would be nice.
    - [X] Create results and temporary file folders
    - [ ] Make sure that all in scripts paths are leading to results folder
- [ ] Create conda env and add all tools as external libs.
- [ ] Describe analysis pipeline
    - [ ] Create a pipiline describing diagrams
    - [ ] Describe how to use scripts
    - [ ] Describe result files generated
    - [ ] Describe directory structure
- [ ] Add all used programs
- [ ] Should I configure file paths as comand line passable arguments?? 
    - [ ] Move from envoirmental variables in qsub to cmd atguments
- [ ] Kā norādīt skriptiem prog. atrašanās vietu?
    - [ ] Picard
    - [ ] added with load module
        - [ ] VCFtolls.
        - [ ] BCF tools

### For every script file check:

1. Description Header
2. Display job information
3. Input handling 
4. Error handling
5. Programm paths
6. File paths
7. Script runs
8. Output of resources used.
9. Output messages.
10. Added to documentation
11. Added to git

### List of verified scripts:

prepare_reference.sh