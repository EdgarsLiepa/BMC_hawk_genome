#!/bin/bash
#PBS -N PL_BIOR_COVID_samtools
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:10:00
#PBS -j oe
 
/home_beegfs/edgars01/tools/samtools-1.15.1/bin/samtools faidx /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna --fai-idx /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fai 
 