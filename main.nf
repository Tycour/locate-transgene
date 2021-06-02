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

process combineReads {

    output:
    file "allReads.fastq.gz" into allReads_ch

    """
    cat $readsDir/*.fastq.gz > allReads.fastq.gz
    NanoPlot --fastq allReads.fastq.gz -o $readsDir/nanoplot -t $threads
    """
}

process filterReads {

    input:
    file allReads from allReads_ch

    output:
    stdout into filtReads_ch

    """
    filtlong --min_length 750 --min_mean_q 70 $allReads
    """
}

process alignSequences {

    input:
    stdin from filtReads_ch

    """
    mkdir $readsDir/Alignment
    cat $transgene $genome >> $readsDir/Alignment/AgamP4_${sampleID}.fa
    minimap2 -t $threads -a -o alignment.sam $readsDir/Alignment/AgamP4_${sampleID}.fa -
    samtools sort alignment.sam > $readsDir/Alignment/${sampleID}.sorted.bam
    samtools index $readsDir/Alignment/${sampleID}.sorted.bam $readsDir/Alignment/${sampleID}.sorted.bai -@ $threads
    """
}