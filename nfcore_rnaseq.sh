#!/bin/bash
#===============================================================================
# NF-CORE RNA-SEQ PIPELINE
# Aligner: STAR + Salmon
# Reference: GRCh37 (hg19) with GENCODE v47lift37
#===============================================================================

# Usage: bash nfcore_rnaseq.sh

# Set paths (modify these for your system)
INPUT_SAMPLESHEET="~/tso-500-@uka-RNA/nfcore_rnaseq_samplesheet.csv"
OUTPUT_DIR="~/tso-500-@uka-RNA/rnaseq_results"
FASTA="~/rnaseq-ref-genome/ref-genome/GRCh37.primary_assembly.genome.fa"
GTF="~/rnaseq-ref-genome/gtf-files/gencode.v47lift37.basic.annotation.gtf"
CONFIG="~/tso-500-@uka-RNA/custom.config"

# Run nf-core/rnaseq pipeline
nextflow run nf-core/rnaseq \
  -r 3.22.2 \
  -profile docker \
  --input ${INPUT_SAMPLESHEET} \
  --outdir ${OUTPUT_DIR} \
  --fasta ${FASTA} \
  --gtf ${GTF} \
  --aligner star_salmon \
  --gencode \
  --max_cpus 16 \
  --max_memory '64 GB' \
  --featurecounts_group_type gene_type \
  --remove_ribo_rna \
  --save_trimmed \
  -c ${CONFIG}

#===============================================================================
# SAMPLESHEET FORMAT (nfcore_rnaseq_samplesheet.csv)
#===============================================================================
# sample,fastq_1,fastq_2,strandedness
# SAMPLE1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,auto
# SAMPLE2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,auto

#===============================================================================
# OUTPUT STRUCTURE
#===============================================================================
# rnaseq_results/
# ├── fastqc/                    # Quality control reports
# ├── trimgalore/                # Trimmed reads
# ├── star_salmon/               # Aligned BAMs + Salmon quant
# │   ├── *.bam                  # Aligned reads
# │   └── salmon.merged.gene_counts.tsv  # Gene counts matrix
# ├── multiqc/                   # Combined QC report
# └── pipeline_info/             # Execution reports
