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

println "\n3. Appending transgene sequence to reference genome file.\n"

process appendSequence2Genome {

    script:
    """
    cat ${params.workdir}/${params.transgene} >> ${params.workdir}/${params.genome}
    """
}

println "\n4. Aligning reads to reference genome.\n"

process mapReads2Genome {

    script:
    """
    minimap2 -a -o ${params.workdir}/${params.sampleID}.sam ${params.workdir}/${params.genome} ${params.workdir}/${params.sampleID}.fastq.gz
    """
}

println "\nTODO - Guide user through IGV process."

println "\nDone"