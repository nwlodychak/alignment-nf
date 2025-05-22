#!/usr/bin/env nextflow

// Params
params.platform = null
params.target_dir = "trimmed"
params.sample_id = ""
params.min_len = 100
params.max_len = 300
params.min_avg_qual = 20
params.trim_adapters = true
params.no_porechop = true
params.adapter_fa = "genomes/adapters.fa"
params.threads = 1


// Logs
log.info """\
         RNA-SEQ ALIGNMENT PIPELINE
         =========================
         sample : ${params.sample_id}
         output directory : ${params.trimming_dir}
         """
         .stripIndent()

// Channels
sample_id_ch = Channel.fromPath(params.genome_reference, checkIfExists: true)


// Process
process cutadapt {
    tag "Trimming for ${sample_id}..."

    input:
    tuple val(sample_id), path(fastqs)
    path(adapter_fa)
    
    output:
    val(sample_id), emit: sample_id
    path('trimmed/*.filtered.fastq.gz'), emit: filter_reads

    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/trimmed/", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.trimming_cores
    memory params.genome_build_memory
    
    script:
    trim_str = params.trim_illumina_adapters ? "--trim_illumina_adapters" : ""
    nochop_str = params.nochop ? "--skip_porechop" : ""
    """
    filter_wrapper.py \
        --platform ${params.platform}
        --target_dir ${parmas.target_dir}
        --sample_id ${sample_id}
        --min_len ${parms.min_len}
        --max_len ${params.max_len}
        --min_avg_qual ${params.min_avg_qual}
        --trim_adapters ${trim_str}
        --no_porechop ${nochop_str}
        --adapter_fa ${params.adatper_fa}
        --threads ${cpus}
    """
}

workflow {
    cutadapt(sample_id_ch)
}

workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}