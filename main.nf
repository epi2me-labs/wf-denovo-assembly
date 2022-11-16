#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'

process summariseReads {
    // concatenate fastq and fastq.gz in a dir

    label "wfdenovoassembly"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val("${meta.sample_id}"), path("${meta.sample_id}.stats"), path("${meta.sample_id}.fastq.gz")
    shell:
    """
    fastcat -s "${meta.sample_id}" -r "${meta.sample_id}.stats" -x "${directory}" > "${meta.sample_id}.fastq"
    bgzip "${meta.sample_id}.fastq"
    """
}


process getVersions {
    label "wfdenovoassembly"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    """
}


process getParams {
    label "wfdenovoassembly"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label "wfdenovoassembly"
    input:
        val metadata
        path "seqs.txt"
        path "versions/*"
        path "params.json"
    output:
        path "wf-template-*.html"
    script:
        report_name = "wf-template-" + params.report_name + '.html'
        def metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    report.py $report_name \
        --versions versions \
        seqs.txt \
        --params params.json \
        --metadata metadata.json
    """
}


process deNovo {
    label "wfdenovoassembly"
    cpus params.threads
    input:
        tuple val(sample_id), path("reads.fastq.gz")
    output:
        tuple val(sample_id), path("${sample_id}.draft_assembly.fasta.gz"), path("${sample_id}_flye_stats.tsv")
        
    """
    flye --nano-hq reads.fastq.gz --read-error ${params.read_error} --genome-size ${params.genome_size} --out-dir output --threads "${task.cpus}"
    mv output/assembly.fasta "./${sample_id}.draft_assembly.fasta"
    mv output/assembly_info.txt "./${sample_id}_flye_stats.tsv"
    bgzip "${sample_id}.draft_assembly.fasta"
    """
}


process assemblyStats {
    label "wfdenovoassembly"
    input:
         tuple val(sample_id), path("draft_assembly.fasta.gz"), path("flye_stats.tsv")

    output:
        tuple val(sample_id), path("${sample_id}_quast_report.pdf")
    """
    quast -o quast_output -t $task.cpus "draft_assembly.fasta.gz" || true
    mv "quast_output/report.pdf" "${sample_id}_quast_report.pdf"
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfdenovoassembly"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// workflow module
workflow pipeline {
    take:
        reads
    main:
        summary = summariseReads(reads)
        software_versions = getVersions()
        workflow_params = getParams()
        metadata = reads.map { it -> return it[1] }.toList()
        denovo = deNovo(summary.map{ it -> tuple(it[0], it[2])})
        assembly_quality = assemblyStats(denovo)
        report = makeReport(metadata, summary.map{ it -> it[1]}, software_versions.collect(), workflow_params)
    emit:
        results = summary.map({ it -> tuple(it[1], it[2])})
        .concat(report, workflow_params,
                denovo.map{it -> it[1]},
                assembly_quality.map{it -> it[1]} )
        
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "unclassified":params.analyse_unclassified])

    pipeline(samples)
    output(pipeline.out.results)
} 

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
