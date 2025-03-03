#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'
    container 'ghcr.io/bf528/multiqc:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*.html'), emit: html

    shell:
    """
    multiqc -f .
    """
}