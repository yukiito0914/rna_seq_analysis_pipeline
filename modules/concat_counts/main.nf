#!/usr/bin/env nextflow

process CONCAT_COUNTS {
    label 'process_single'
    container 'ghcr.io/bf528/pandas:latest'
    publishDir params.outdir, mode: "copy"

    input:
    path(verse_counts)

    output:
    path("combined_counts.txt"), emit: combined

    shell:
    """
    python ${workflow.projectDir}/bin/concat_counts.py -i $verse_counts -o combined_counts.txt
    """
}