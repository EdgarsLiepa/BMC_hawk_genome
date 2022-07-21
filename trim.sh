#!/bin/bash
#PBS -N trimming
#PBS -l walltime=50:00:00
#PBS -l nodes=1:ppn=4
#PBS -j oe
 
export LC_ALL=lv_LV.utf8
export LANG=lv_LV.utf8
 
for f in $(ls /home_beegfs/edgars01/Ineta/combined/*.fq.gz | sed -e 's/_1.fq.gz//' -e 's/_2.fq.gz//' | sort -u)
do
java -jar /home_beegfs/edgars01/trimmomatic/trimmomatic-0.39.jar PE -phred33 ${f}_1.fq.gz ${f}_2.fq.gz ${f}_F_paired.fq.gz ${f}_F_unpaired.fq.gz ${f}_R_paired.fq.gz ${f}_R_unpaired.fq.gz LEADING:30 TRAILING:30 MINLEN:36
done
