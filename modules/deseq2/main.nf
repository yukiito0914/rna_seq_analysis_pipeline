#!/usr/bin/env nextflow

process DESEQ2 {
    label 'process_high'
    conda 'envs/deseq2_env.yml'
    publishDir params.outdir, mode: "copy"

    input:
    path(counts_matrix)
    path(id_genename)

    output:
    path("DESeq2_results.csv"), emit: result
    path("DESeq2_top10.csv"), emit: top
    path("DESeq2_summary.txt"), emit: summary
    path("PCA_plot.png")
    path("Heatmap.png")
    path("Top10_Pathways.png")

    shell:
    """
    Rscript ${workflow.projectDir}/bin/deseq2.R $counts_matrix ${params.metadata} DESeq2_results.csv DESeq2_top10.csv DESeq2_summary.txt ${params.deseq2_padj_threshold} $id_genename
    """
}