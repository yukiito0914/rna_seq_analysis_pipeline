# RNA-seq Differential Expression Analysis Pipeline
This Nextflow pipeline processes RNA-seq data from raw FASTQ files through quality control, alignment, quantification, filtering, and differential expression analysis using DESeq2, followed by gene set enrichment analysis (GSEA) using fgsea and MSigDB canonical pathways.

## Pipeline Overview
FASTQ → QC (FASTQC) → Alignment (STAR) → Quantification (VERSE) → Count Matrix → Filtering → DESeq2 Differential Expression → GSEA → Plots and Reports

## Pipeline Modules
1. FASTQC – Quality control of raw reads  
2. STAR_INDEX – Genome indexing using FASTA and GTF  
3. PARSE_GTF – Gene ID to gene symbol mapping  
4. STAR_MAPPING – Alignment to reference genome  
5. MULTIQC – Aggregated quality control summary  
6. VERSE – Read quantification at the gene level  
7. CONCAT_COUNTS – Merge individual count files  
8. FILTER_COUNTS – Filter low-expression genes  
9. DESEQ2 – Differential expression and GSEA analysis  

## Input Requirements
The following input files and parameters must be provided:

| Parameter                  | Description                                            |
|---------------------------|--------------------------------------------------------|
| `params.reads`            | Path to paired-end FASTQ files                         |
| `params.genome`           | Reference genome in FASTA format                       |
| `params.gtf`              | Gene annotation file in GTF format                     |
| `params.metadata`         | Sample metadata file with `sample` and `condition`     |
| `params.deseq2_padj_threshold` | Adjusted p-value threshold for significance        |

## Running the Pipeline
To execute the pipeline with default parameters:
```bash
nextflow run main.nf
```

## Output Files
The pipeline produces the following outputs under ```results/```:

- FASTQC reports for input FASTQ files
- STAR alignment output and BAM files
- Raw and filtered count matrices
- Combined QC report
- DE results, plots, and enrichment analysis

## Input Format Examples
metadata.csv:
```csv
sample,condition
sample1,control
sample2,treated
...
```