#! /usr/bin/env nextflow

/*
 * 'Locate Transgene-NF' - A Nextflow pipeline for locating randomly inserted transgenes with Nanopore data
 */

/*
 * Define the default parameters
 */

sampleID = params.sampleID
readsDir = params.readsDir
transgene = file(params.transgene)
genome = file(params.genome)
threads = params.threads

log.info """\n
L O C A T E   T R A N S G E N E  -  N F    v 1.0
================================
sampleID : $sampleID
readsDir : $readsDir
genome   : $genome
transgene: $transgene
threads  : $threads
\n
"""

/*
 * Step 1. Combine fastq files and create sequencing statistics with NanoPlot
 */

process combineReads {
    conda 'NanoPlot'

    output:
    file "allReads.fastq.gz" into allReads_ch

    """
    cat $readsDir/*.fastq.gz > allReads.fastq.gz
    NanoPlot --fastq allReads.fastq.gz -o $readsDir/nanoplot -t $threads
    """
}

/*
 * Step 2. Filter reads according to length and quality score
 */

process filterReads {
    conda 'filtlong'

    input:
    file allReads from allReads_ch

    output:
    stdout into filtReads_ch

    """
    filtlong --min_length 750 --min_mean_q 70 $allReads
    """
}

/*
 * Step 3. Align Nanopore reads to modified reference genome (with added transgene contig)
 * Step 4. Sort and index for viewing in preferred Genome Browser (e.g. Integrated Genome Viewer)
 */

process alignSequences {
    conda 'minimap2 samtools'

    input:
    stdin from filtReads_ch

    """
    mkdir -p $readsDir/Alignment
    cat $transgene $genome >> $readsDir/Alignment/AgamP4_${sampleID}.fa
    minimap2 -t $threads -a -o alignment.sam $readsDir/Alignment/AgamP4_${sampleID}.fa -
    samtools sort alignment.sam > $readsDir/Alignment/${sampleID}.sorted.bam
    samtools index $readsDir/Alignment/${sampleID}.sorted.bam $readsDir/Alignment/${sampleID}.sorted.bai -@ $threads
    """
}