#!/usr/bin/env nextflow

process alignment {
    tag "${sample_id}"
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/aligned/", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.align_cores
    memory "${params.align_mem}"

    input:
    tuple val(sample_id), path(fastqs)
    path(adapter_fa)
    
    output:
    val(sample_id), emit: sample_id
    path('aligned/*.bam'), emit: aligned
    
    script:

    """
    alignment_wrapper.py \
    """
}
