This pipeline requires two configuration file:

# `config.yaml`

A standard `Snakemake` configuration, yaml-formatted file containing a list of
all parameters accepted in this workflow:

* `samples`: Path to the file containing link between samples and their fastq file(s)
* `params`: Per-tool list of optional parameters

Example:

```
samples: config/samples.csv

# Optional parameters
params:
    fair_fastqc_multiqc:
        # Optional parameters for multiqc
        multiqc: --module fastqc --zip-data-dir --verbose --no-megaqc-upload --no-ansi --force
```

A complete list of accepted keys is available [in schemas](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/schemas/config.schema.yaml),
with their default value, expected type, and human readable description.


# `samples.csv`

A CSV-formatted text file containing the following mandatory columns:

* `sample_id`: Unique name of the sample
* `upstream_file`: Path to upstream fastq file
* `species`: The species name, according to Ensembl standards.
* `build`: The corresponding genome build, according to Ensembl standards.
* `release`: The corresponding genome release, according to Ensembl standards.
* `downstream_file`: Path to downstream fastq file, leave empty in case of Single ended library.

A complete list of accepted keys is available [in schemas](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/schemas/samples.schema.yaml),
with their default value, expected type, and human readable description.

Example:

```
sample_id,upstream_file,downstream_file,species,build,release
sac_a,data/reads/a.scerevisiae.1.fq,data/reads/a.scerevisiae.2.fq,saccharomyces_cerevisiae,R64-1-1,110
```

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.

# `genomes.csv`

This file is fully optional. When missing, the genome sequences
will be downloaded from Ensembl and indexed.

A CSV-formatted text file containing the following mandatory columns:

* `species`: The species name, according to Ensembl standards
* `build`: The corresponding genome build, according to Ensembl standards
* `release`: The corresponding genome release, according to Ensembl standards

Example:

```
species,build,release
homo_sapiens,GRCh38,105
mus_musculus,GRCm38,99
mus_musculus,GRCm39,110
```

A complete list of accepted keys is available [in schemas](https://github.com/tdayris/fair_fastqc_multiqc/blob/main/workflow/schemas/genomes.schema.yaml),
with their default value, expected type, and human readable description.


Note:

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.
