#!/usr/bin/env nextflow

process PARSE_GTF {
    label 'process_single'
    container 'ghcr.io/bf528/biopython:latest'
    publishDir params.outdir

    input:
    path(gtf)

    output:
    path("*.txt"), emit: txt

    shell:
    """
    python ${workflow.projectDir}/bin/parse_gtf.py -i $gtf -o id_genename.txt
    """
}