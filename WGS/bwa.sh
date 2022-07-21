#!/bin/bash
#PBS -N PL_BIOR_COVID_GoshHawk
#PBS -l nodes=1:ppn=4
#PBS -l walltime=5:00:00
#PBS -j oe
 
/home/groups/bmc/projects/Innas_Eksomi/nikita/bwa-0.7.17/bwa index /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna
