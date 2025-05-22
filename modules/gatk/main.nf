#!/usr/bin/env nextflow

if (params.genome == "human") {


process gatk {
    tag "${sample_id}"

    input:
    set val(name), file(bam) from bams
     
    output:
    set val(name), "${name}.ro.dedup.bam" into preprocessChannel
    file("${name}.ro.dedup.bam.bai") into preprocessBaiChannel

    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/vcf/", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.trimming_cores
    memory "${params.trimming_mem}"


workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}