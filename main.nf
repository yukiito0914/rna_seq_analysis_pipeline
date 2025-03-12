#!/usr/bin/env nextflow

include {FASTQC} from './modules/fastqc'
include {STAR_INDEX} from './modules/star_index'
include {PARSE_GTF} from './modules/parse_gtf'
include {STAR_MAPPING} from './modules/star_mapping'
include {MULTIQC} from './modules/multiqc'
include {VERSE} from './modules/verse'
include {CONCAT_COUNTS} from './modules/concat_counts'
include {FILTER_COUNTS} from './modules/filter_counts'
include {DESEQ2} from './modules/deseq2'

workflow {

    Channel.fromFilePairs(params.reads)
    | set {align_ch}

    Channel.fromFilePairs(params.reads).transpose()
    | set {fastqc_ch}

    // Perform Quality Control
    FASTQC(fastqc_ch)
    // Generate a genome index
    STAR_INDEX(params.genome, params.gtf)
    // Generate a file containing the gene IDs and their corresponding human gene symbols
    PARSE_GTF(params.gtf)
    // Align reads to the genome
    STAR_MAPPING(align_ch, STAR_INDEX.out.index)

    FASTQC.out.zip.map{it[1]}.collect()
    | set {fastqc_out}

    STAR_MAPPING.out.log.map{it[1]}.collect()
    | set {star_log}

    fastqc_out.mix(star_log).flatten().collect()
    | set {multiqc_ch}

    // Perform post-alignment QC
    MULTIQC(multiqc_ch)
    // Quantify alignments to the genome
    VERSE(STAR_MAPPING.out.bam, params.gtf)

    VERSE.out.counts.map{it[1]}.collect()
    | set {concat_ch}

    // Concatenate count outputs into a single matrix
    CONCAT_COUNTS(concat_ch)
    // Filtering the counts matrix
    FILTER_COUNTS(CONCAT_COUNTS.out.txt)
    // Perform differential expression analysis using the filtered counts
    DESEQ2(FILTER_COUNTS.out.flltered, PARSE_GTF.out.txt)

}
