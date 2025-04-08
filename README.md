# ChIP-seq Analysis Pipeline

This pipeline automates the analysis of ChIP-seq datasets from the ENCODE database, providing a continuous workflow from data acquisition to visualization.

## Overview

This pipeline handles the following steps:
- Automated download of ChIP-seq datasets from ENCODE with multiple input options from the user
- Quality control of sequencing reads
- Alignment to reference genomes
- Peak calling
- Visualization of results

## Requirements

- Python 3.6+
- BWA
- SAMtools
- FastQC
- MultiQC
- MACS2
- deepTools
- PyYAML


## Usage

Run the pipeline using the main script:

```bash
python main.py [ENCODE_URL_OR_FILE] --dir [OUTPUT_DIRECTORY] [OPTIONS]
```

### Command Line Arguments

The pipeline accepts various arguments to customize data download and analysis:

- `url_or_file`: List of ENCODE URLs, files, or accession IDs
- `--dir`: Root directory for downloaded data and analysis results
- `--file-types`: List of file types to download (default: fastq)
- `--assemblies`: List of allowed assemblies (default: all)
- `--encode-access-key-id`: ENCODE access key ID for accessing unpublished data
- `--encode-secret-key`: ENCODE secret key
- `--pooled-rep-only`: Download data from pooled replicates only
- `--dry-run`: Test run without downloading files
- `--max-download`: Maximum number of concurrent downloads

## Pipeline Components

### ENCODE Downloader

The `encode_downloader.py` script handles data retrieval from ENCODE:

- Downloads files based on specified criteria
- Supports authentication for accessing unpublished data
- Creates organized directory structure
- Generates metadata files

### Quality Control

The `quality_control.py` module:
- Runs FastQC on downloaded FASTQ files
- Generates MultiQC reports for easy quality assessment

### Read Alignment

The `read_processor.py` module:
- Aligns reads to reference genome using BWA-mem
- Sorts and indexes BAM files using SAMtools

### Peak Calling

The `peak_caller.py` module:
- Calls peaks using MACS2
- Compares treatment samples against controls

### Visualization

The `visualizer.py` module:
- Creates coverage tracks using bamCoverage
- Generates BigWig files normalized with RPKM

## Example Workflow

1. Download ChIP-seq data for a specific experiment:
   ```bash
   python main.py ENCSR000DKU --dir ./chipseq_analysis
   ```

2. The pipeline will automatically:
   - Download the appropriate FASTQ files
   - Run quality control
   - Align reads to the reference genome
   - Call peaks
   - Generate visualization files

3. Results will be organized in the specified directory structure

## Output Structure

```
output_directory/
├── ENCSRXXXXXX/
    ├──released/reads/fastq
    ├──metadata.json
├── qc_results/           # FastQC results
├── multiqc_report/       # MultiQC summary
├── aligned_*.bam         # Aligned reads
├── peaks/                # Peak calling results
└── visualization/        # Coverage tracks
```