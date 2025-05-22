#!/usr/bin/env nextflow

// Params
params.genome_reference = null
params.gtf = null
params.genome_dir = "genomes"
params.publish_dir_mode = "copy"
params.alignment_image_uri = "your_container"
params.alignment_image_version = "latest"
params.genome_build_cpus = 8
params.genome_build_memory = "32 GB"

// Logs
log.info """\
         RNA-SEQ ALIGNMENT PIPELINE
         =========================
         genome reference : ${params.genome_reference}
         gtf file         : ${params.gtf}
         output directory : ${params.genome_dir}
         """
         .stripIndent()

// Channels
genome_ch = Channel.fromPath(params.genome_reference, checkIfExists: true)
gtf_ch = Channel.fromPath(params.gtf, checkIfExists: true)

process genome_build {
    tag "Building STAR index for ${genome}"

    input:
    file(genome_reference), file(gtf)
    
    output:
    path genome_reference
    path gtf

    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.genome_dir}", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.genome_build_cpus
    memory params.genome_build_memory
    
    script: 
    """
    export LC_CTYPE=en_US.UTF-8
    export LC_ALL=en_US.UTF-8
    ################################################################################################
    ### GENOME BUILD ###
    ################################################################################################

    mkdir -p ${params.genome_dir}
    STAR --runThreadN ${task.cpus} \
            --runMode genomeGenerate \
            --genomeDir ${params.genome_dir} \
            --genomeFastaFiles ${genome_reference} \
            --sjdbGTFfile ${gtf} \
            --sjdbOverhang 84
    """

workflow {
    genome_build(genome_ch, gtf_ch)
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}