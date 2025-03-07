#!/usr/bin/env nextflow

process FILTER_COUNTS {
    label 'process_single'
    container 'ghcr.io/bf528/pandas:latest'
    publishDir params.outdir, mode: "copy"

    input:
    path(counts_matrix)

    output:
    path("filtered_counts.csv"), emit: flltered
    path("filtering_summary.txt"), emit: summary

    shell:
    """
    python ${workflow.projectDir}/bin/filter_counts.py -i $counts_matrix -t ${params.filtering_threshold} -o filtered_counts.csv -s filtering_summary.txt
    """
}