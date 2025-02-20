#!/usr/bin/env nextflow

include {FASTQC} from './modules/fastqc'
include {STAR_INDEX} from './modules/star_index'
include {PARSE_GTF} from './modules/parse_gtf'

workflow {

    Channel.fromFilePairs(params.reads)
    | set {align_ch}

    Channel.fromFilePairs(params.reads).transpose()
    | set {fastqc_ch}

    FASTQC(fastqc_ch)
    STAR_INDEX(params.genome, params.gtf)
    PARSE_GTF(params.gtf)

}
