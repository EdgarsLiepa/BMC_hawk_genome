#!/bin/bash
#PBS -l nodes=1:ppn=4,mem=10G
#PBS -l walltime=48:59:59
#PBS -j oe
#PBS -N bamqc 
#PBS -A bmc_pl_bior_covid
#PBS -q batch

# -- Displaying job information
echo Running on host `hostname`
echo Time is `date`
echo Current working directory is `pwd`
echo "Node file: $PBS_NODEFILE :"
echo "---------------------"
cat $PBS_NODEFILE
echo "---------------------"
echo Using ${NPROCS} processors across ${NNODES} nodes

FILEPATH="/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM"    
OUTPATH="/home_beegfs/edgars01/Ineta/WGS/Results/qualimap_results" 
qualimapPath="/home_beegfs/edgars01/tools/qualimap_v2.2.1"


cd $PBS_O_WORKDIR

module load java
module load singularity
module load conda

~/tools/nextflow run /home_beegfs/edgars01/Ineta/WGS/qualimap.nf \
	-c /home_beegfs/edgars01/Antibiotic_Resistance/Metagenome/NextFlowPipeline/nextflow.config \
    	-with-report -resume 


# # remove filepath from ls output
# for f in $(ls ${FILEPATH}/*.markdup.fixedRG.bam | sed -e 's/\/home_beegfs\/edgars01\/Ineta\/WGS\/starpfaili\/mappedBAMandSAM\///'  | sed -e 's/.markdup.fixedRG.bam//' | sort -u); do
#     echo "Start: $f - $date"
#     echo "analyze: ${FILEPATH}/$f.markdup.fixedRG.bam "
# 
#     # find path to qualimap 
#     ${qualimapPath}/qualimap bamqc --java-mem-size=30G -bam ${FILEPATH}/${f}.markdup.fixedRG.bam -outdir ${OUTPATH}/$f
#     
#     echo "Done: $f"
# done

echo Time is `date`

### Output job resurces used
/usr/local/bin/qstat -f $PBS_JOBID | /bin/grep -e Job -e resources
