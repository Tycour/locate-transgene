#! /usr/bin/env nextflow

println "\nI will NanoPlot $params.reads and store results in $params.outdir using $params.threads CPU threads."

process runNanoPlot {

    script:
    """
    NanoPlot --fastq $params.reads -o $params.outdir -t $params.threads
    """
}
