This pipeline requires two configuration file:

# `config.yaml`

A standard `Snakemake` configuration, yaml-formatted file containing a list of
all parameters accepted in this workflow:

* `samples`: Path to the file containing link between samples and their fastq file(s)

Example:

```
samples: config/samples.csv
```

# `samples.csv`

A CSV-formatted text file containing the following mandatory columns:

* `sample_id`: Unique name of the sample
* `upstream_file`: Path to upstream fastq file
* `downstream_file`: Optional path to downstream fastq file

Example:

```
sample_id,upstream_file,downstream_file,species,build,release
sac_a,data/reads/a.scerevisiae.1.fq,data/reads/a.scerevisiae.2.fq,saccharomyces_cerevisiae,R64-1-1,110
```

While `CSV` format is tested and recommended, this workflow uses python
`csv.Sniffer()` to detect column separator. Tabulation and semicolumn are
also accepted as field separator. Remember that only comma-separator is
tested.
