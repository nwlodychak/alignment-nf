#!/usr/bin/env nextflow

process gatk {
    tag "${sample_id}"
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/vcf/", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.trimming_cores
    memory "${params.trimming_mem}"