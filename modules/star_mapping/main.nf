#!/usr/bin/env nextflow

process STAR_MAPPING {
    label 'process_high'
    container 'ghcr.io/bf528/star:latest'
    publishDir params.outdir

    input:
    tuple val(name), path(reads)
    path(index)

    output:
    tuple val(name), path('*.Aligned.sortedByCoord.out.bam'), emit: bam
    tuple val(name), path('*.Log.final.out'), emit: log

    shell:
    """
    STAR --runThreadN $task.cpus --genomeDir $index --readFilesIn ${reads.join(' ')} --readFilesCommand zcat --outFileNamePrefix $name. --outSAMtype BAM SortedByCoordinate
    """
}