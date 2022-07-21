#!/bin/bash
#PBS -N PL_BIOR_COVID_picard
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -j oe

# šito vajag lejuplādēt  
java -jar /home_beegfs/edgars01/tools/picard.jar CreateSequenceDictionary \
-R /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.fna \
-O /home_beegfs/edgars01/Ineta/WGS/GosHawkReference-ncbi-genomes-2022-07-05/GCA_929443795.1_bAccGen1.1_genomic.dict