#! /usr/bin/env nextflow

println "\n1. Preparing files for analysis.\n"

process combineFastQ {

    script:
    """
    cat ${params.workdir}/*.fastq.gz > ${params.workdir}/${params.sampleID}.fastq.gz
    """
}

println "\n2. Generating sequencing statistics.\n"

process runNanoPlot {

    script:
    """
    NanoPlot --fastq ${params.workdir}/${params.sampleID}.fastq.gz -o ${params.workdir}/NanoPlot -t $params.threads
    """
}

println "\nDone"