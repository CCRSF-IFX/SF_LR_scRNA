# Single-Cell Analysis Pipeline for ONT Long Reads

This repository houses a comprehensive Snakemake workflow designed for the analysis of single-cell sequencing data from Oxford Nanopore Technologies (ONT). It seamlessly integrates data processing using ONT's wf-single-cell Nextflow pipeline, downstream analysis with standard R packages for single-cell analysis, and utilizes SQANTI3 for detailed isoform classification.

## Pipeline Overview

The pipeline processes single-cell sequencing data through the following major stages:
1. **ONT wf-single-cell**: Utilizes the ONT Nextflow pipeline to process raw reads into matrix files.
2. **Single-Cell Data Analysis**: Implements R packages like Seurat and SingleR for cell clustering and typing.
3. **Isoform Classification and Quantification**: Uses SQANTI3 to classify cell-specific isoforms, followed by merging the quantification data back onto the reads.
4. **Visualization and Reporting**: Generates plots and detailed reports for both cluster analysis and isoform classification.

## Prerequisites

### Software Requirements

- **Snakemake**: Orchestrates and manages the workflow execution.
- **Nextflow**: Required for running the ONT-specific workflows.
- **Singularity**: Ensures reproducibility by containerizing the software environments.
- **R**: A statistical computing environment where several analysis scripts are executed. Ensure the following libraries are installed:
  - Seurat
  - SingleR
  - clusterProfiler
- **Python**: Needed particularly for SQANTI3 and data merging scripts.
- **Conda**: Recommended for managing Python and R environments seamlessly.

### System Requirements

- Access to a computational environment with a substantial amount of memory and CPU resources.
- SLURM or another job scheduler for managing job submissions.

## Installation

Clone this repository to your local machine or computational cluster:

```bash
git clone https://github.com/[YourUsername]/single-cell-ont-analysis.git
cd single-cell-ont-analysis
```

Ensure all dependencies are installed by setting up environments through Conda:

```bash
conda env create -f environment.yml
```

## Configuration

Modify the `config.yaml` file to reflect your project-specific settings including paths, reference genomes, and tool parameters. Sample configurations and environment setup details are also provided for clarity.

## Running the Pipeline

To execute the pipeline, use the following command:

```bash
sbatch submit.sh
```

## Outputs
The pipeline will generate the following outputs:

- Gene expression matrices in .mtx.gz format.
- R objects for downstream analysis (.rds).
- Classified isoforms and comprehensive reports in HTML format from SQANTI3.
- Final merged classification data, linking quantifications to read data.

## Detailed Scripts Description
Scripts located in `scripts` are critical for various stages of the pipeline, including data merging and additional analyses not handled directly by Snakemake.

## For any questions, contact

CCRSF_IFX@nih.gov

## Acknowledgments
- Thanks to EPI2ME Labs for providing the wf-single-cell Nextflow pipeline.
- Gratitude to the developers of Seurat, SingleR, and SQANTI3 for their open-source software which this pipeline depends on.