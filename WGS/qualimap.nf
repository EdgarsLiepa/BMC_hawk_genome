#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input_bam_dir = "/home_beegfs/edgars01/Ineta/WGS/starpfaili/mappedBAMandSAM"
params.bam_pattern = "*.markdup.fixedRG.bam"
params.outdir = "/home_beegfs/edgars01/Ineta/WGS/Results/qualimap_results"
params.java_mem = "30G"

process QUALIMAP_BAMQC {
    
    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample_id}/**"
    
    cpus 6
    memory '32 GB'
    time '1d'   
 
    input:
    tuple val(sample_id), path(bam_file)
    
    output:
    path "${sample_id}/**", emit: qc_reports
    
    script:
    """
    echo "Start: ${sample_id} - \$(date)"
    echo "analyze: ${bam_file}"
    
    /home_beegfs/edgars01/tools/qualimap_v2.2.1/qualimap bamqc \\
        --java-mem-size=${params.java_mem} \\
        -bam ${bam_file} \\
        -outdir ${sample_id}
    
    echo "Done: ${sample_id}"
    """
}

workflow {
    // Create channel from BAM files
    bam_ch = Channel
        .fromPath("${params.input_bam_dir}/${params.bam_pattern}")
        .map { bam_file ->
            // Extract sample ID by removing the directory path and file extension
            def sample_id = bam_file.baseName.replaceAll(/\.markdup\.fixedRG$/, '')
            [sample_id, bam_file]
        }
    
    // Run Qualimap bamqc
    QUALIMAP_BAMQC(bam_ch)
}
