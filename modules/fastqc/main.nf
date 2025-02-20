#!/usr/bin/env nextflow

process FASTQC {
    container 'ghcr.io/bf528/fastqc:latest'
    label 'process_single'
    publishDir params.outdir

    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path('*.zip'), emit: zip
    tuple val(name), path('*.html'), emit: html

    shell:
    """
    fastqc $fastq -t $task.cpus
    """
}