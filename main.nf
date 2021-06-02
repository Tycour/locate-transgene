#! /usr/bin/env nextflow

sampleID = params.sampleID
readsDir = params.readsDir
transgene = file(params.transgene)
genome = file(params.genome)
threads = params.threads

log.info """\n
L O C A T E   T R A N S G E N E  -  N F    v 0.1
================================
sampleID : $sampleID
readsDir : $readsDir
genome   : $genome
transgene: $transgene
threads  : $threads
\n
"""

process prepareReads {

    output:
    file "allReads.fq.gz" into reads_ch

    """
    cat $readsDir/*.fastq.gz > allReads.fq.gz
    NanoPlot --fastq allReads.fq.gz -o $readsDir/NanoPlot -t $threads
    """
}

process alignSequences {

    input:
    file allReads from reads_ch

    """
    mkdir $readsDir/Alignment
    cat $transgene $genome >> $readsDir/Alignment/AgamP4_${sampleID}.fa
    minimap2 -t $threads -a -o alignment.sam $readsDir/Alignment/AgamP4_${sampleID}.fa $allReads
    samtools sort alignment.sam > $readsDir/Alignment/${sampleID}.sorted.bam
    samtools index $readsDir/Alignment/${sampleID}.sorted.bam $readsDir/Alignment/${sampleID}.sorted.bai -@ $threads
    """

}