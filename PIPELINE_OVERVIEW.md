# RNA-seq Meta-analysis Pipeline

## Overview

This pipeline performs comprehensive bulk RNA-seq meta-analysis across multiple studies. It includes quality control, read processing, alignment/quantification, batch effect correction, differential expression analysis, and meta-analysis statistics.

## Pipeline Components

### 1. Quality Control & Preprocessing
- **FastQC**: Quality assessment of raw sequencing reads
- **FASTP**: Read trimming, filtering, and quality control

### 2. Alignment & Quantification
Two alignment strategies are supported:

#### STAR + featureCounts (default)
- **STAR_GENOMEGENERATE**: Generate genome index
- **STAR_ALIGN**: Align reads to reference genome
- **SUBREAD_FEATURECOUNTS**: Count reads mapping to genomic features

#### Salmon (alternative)
- **SALMON_INDEX**: Generate transcript index
- **SALMON_QUANT**: Quantify transcript abundance

### 3. Meta-analysis Modules
- **BATCH_CORRECTION**: Remove batch effects using ComBat-seq or SVA
- **DIFFERENTIAL_EXPRESSION**: Perform differential expression analysis per study using limma-voom/edgeR
- **META_ANALYSIS**: Combine results across studies using random-effects meta-analysis

### 4. Reporting
- **MultiQC**: Comprehensive quality control report

## Key Parameters

### Input
- `--input`: Samplesheet with sample information
- `--sample_info`: Sample metadata including batch and condition information

### Reference Files
- `--fasta`: Reference genome FASTA file
- `--gtf`: Gene annotation GTF file
- `--transcript_fasta`: Transcript sequences (for Salmon)

### Analysis Options
- `--aligner`: Choose between 'star' (default) or 'salmon'
- `--perform_batch_correction`: Enable batch effect correction (default: true)
- `--perform_meta_analysis`: Enable meta-analysis across studies (default: true)

## Output Structure

```
results/
├── fastqc/                    # FastQC reports
├── fastp/                     # Read trimming results
├── star/                      # STAR alignment results (if using STAR)
├── salmon/                    # Salmon quantification (if using Salmon)
├── featurecounts/            # Gene count matrices (if using STAR)
├── batch_correction/         # Batch-corrected count matrices
├── differential_expression/  # DE analysis results per study
├── meta_analysis/           # Meta-analysis results
├── multiqc/                 # Comprehensive QC report
└── pipeline_info/           # Pipeline execution information
```

## Sample Sheet Format

The input samplesheet should contain the following columns:
- `sample`: Sample identifier
- `fastq_1`: Path to R1 FASTQ file
- `fastq_2`: Path to R2 FASTQ file (for paired-end)
- `study`: Study identifier for meta-analysis
- `condition`: Experimental condition
- `batch`: Batch identifier (optional, for batch correction)

## Usage Example

```bash
nextflow run nf-core/rnaseqmeta \\
    --input samplesheet.csv \\
    --sample_info sample_metadata.txt \\
    --fasta genome.fa \\
    --gtf annotations.gtf \\
    --outdir results \\
    --aligner star \\
    --perform_batch_correction true \\
    --perform_meta_analysis true \\
    -profile docker
```

## Meta-analysis Features

1. **Batch Effect Correction**: Uses ComBat-seq or SVA to remove technical batch effects
2. **Cross-study Normalization**: Ensures comparable expression levels across studies
3. **Random Effects Meta-analysis**: Combines effect sizes using metafor package
4. **Heterogeneity Assessment**: Evaluates consistency of effects across studies
5. **Comprehensive Visualization**: Forest plots, volcano plots, PCA plots

This pipeline is designed to handle multiple RNA-seq studies simultaneously and provide robust meta-analysis results for identifying consistently differentially expressed genes across studies.