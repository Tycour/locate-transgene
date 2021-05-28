#! /usr/bin/env nextflow

params.outdir = "NanoPlot"
params.reads = "reads.fasta"

println "I will NanoPlot $params.reads and store results in $params.outdir"
