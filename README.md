[![Snakemake](https://img.shields.io/badge/snakemake-≥8.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/tdayris/fair_fastqc_multiqc/workflows/Tests/badge.svg)](https://github.com/tdayris/fair_fastqc_multiqc/actions?query=branch%3Amain+workflow%3ATests)

Snakemake workflow used to control read qualities over fastq files with FastQC.

This pipeline is meant to be used in teaching sessions or imported in other snakemake-workflows.

## Usage

The usage of this workflow is described in the [Snakemake workflow catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_fastqc_multiqc) it is also available [locally](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/report/usage.rst) on a single page.


## Results

A complete description of the results can be found here in [workflow reports](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/report/results.rst).

## Material and Methods

The tools used in this pipeline are described [here](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/report/material_methods.rst) textually.

```

       ┌────────────────────────┐
       │  Copy / Link / Concat  │
       └─┬───────────────────┬──┘
         │                   │
┌────────▼─┐             ┌───▼────────────┐
│          │             │                │
│  FastQC  │             │  FastQ-Screen  │
│          │             │   (Optional)   │
└─────┬────┘             └─────┬──────────┘
      │                        │
      │                        │
      │     ┌───────────┐      │
      │     │           │      │
      └─────►  MultiQC  ◄──────┘
            │           │
            └───────────┘

```

### QC

| Step       | Wrapper                                                                                        |
| ---------- | ---------------------------------------------------------------------------------------------- |
| FastQC     | [fastqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/fastqc.html)     |
| FastScreen | [fastq-screen](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastq_screen.html) |
| MultiQC    | [multiqc-wrapper](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/wrappers/multiqc.html)   |