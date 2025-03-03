#!/usr/bin/env nextflow

include {FASTQC} from './modules/fastqc'
include {STAR_INDEX} from './modules/star_index'
include {PARSE_GTF} from './modules/parse_gtf'
include {STAR_MAPPING} from './modules/star_mapping'
include {MULTIQC} from './modules/multiqc'
include {VERSE} from './modules/verse'
include {CONCAT_COUNTS} from './modules/concat_counts'

workflow { 

    Channel.fromFilePairs(params.reads)
    | set {align_ch}

    Channel.fromFilePairs(params.reads).transpose()
    | set {fastqc_ch}

    FASTQC(fastqc_ch)
    STAR_INDEX(params.genome, params.gtf)
    PARSE_GTF(params.gtf)
    STAR_MAPPING(align_ch, STAR_INDEX.out.index)

    FASTQC.out.zip.map{it[1]}.collect()
    | set {fastqc_out}

    STAR_MAPPING.out.log.map{it[1]}.collect()
    | set {star_log}

    fastqc_out.mix(star_log).flatten().collect()
    | set {multiqc_ch}

    MULTIQC(multiqc_ch)
    VERSE(STAR_MAPPING.out.bam, params.gtf)

    VERSE.out.counts.map{it[1]}.collect()
    | set {concat_ch}

    CONCAT_COUNTS(concat_ch)

}
