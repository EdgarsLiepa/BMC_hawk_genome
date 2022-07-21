#!/bin/bash

#PBS -N rez_assembly
#PBS -l walltime=336:00:00 
#PBS -l nodes=1:ppn=32,pmem=10g
#PBS -q batch
#PBS -j oe
#PBS -q long


cd $PBS_O_WORKDIR

../../MaSuRCA-4.0.6/bin/masurca masurca_config.txt

./assemble.sh
