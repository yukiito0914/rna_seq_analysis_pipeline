#!/usr/bin/env nextflow

process VERSE {
    container 'ghcr.io/bf528/verse:latest'
    label 'process_single'
    publishDir params.outdir, mode: "copy"

    input:
    tuple val(name), path(bam)
    path(gtf)

    output:
    tuple val(name), path("*exon.txt"), emit: counts

    shell:
    """
    verse -a $gtf -o $name $bam -S
    """
}