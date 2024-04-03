# 2.2.6

## Features

* Snakemake wrappers up to 3.7.0

# 2.2.5

## Features:

* Documentation update

# 2.2.4

## Features:

* Documentation update
* Use human readable functions rather than raw lookup
* schemas update

# 2.2.3

## Fix:

* Since copies are done on node and not on login-node, the `on_flamingo()` flag is removed

# 2.2.2

## Fix:

* symlink fully absolute

# 2.2.2

## Fix:

* Missing packages in python environment used for copy
* Fix FileNotFoundError in apptainer when using this pipeline as an external module
* `rsync` syntax error

# 2.2.1

## Fix

* Json schema


# 2.2.0

## Features

* Configuration keys are *all* optional
* Snakemake wrappers update to version 3.5.2

# 2.1.2

## Features:

* Changes in configuration

# 2.1.0

## Features:

* Snakemake wrappers updated to 3.4.1: including an exciting MultiQC update!
* A dedicated MultiQC configuration file is created
* Pipeline containerized

## Fixes:

* Documentation errors
* Report typos

# 2.0.4

## Features:

* Debug logs in `link_or_concat.py`


# 2.0.3

## Fixes:

* AttributeError when empty lookup list is returned

# 2.0.2

## Features:

* Removed slurm partition function
* Job reservation based on input file size
* tempfiles, logs and benchmarks paths reorganized: 
    * `tmp/fair_fastqc_multiqc/{rule_name}/{wildcards}.{extension}`
    * `log/fair_fastqc_multiqc/{rule_name}/{wildcards}.log`
    * `benchmark/fair_fastqc_multiqc/{rule_name}/{wildcards}.tsv`
* Use of `lookup` to find configuration values
* Dag as ascii-art

## Fixes:

* Documentation error


# 2.0.1

## Features:

* Add `slurm_partition` for snakemake plugin executor

## Documentation

* Usage updated

# 2.0.0

Snakemake v8.1+ required

## Features

* New QC: [FastqScreen](https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/fastq_screen.html), available if, and only if a configuration file is provided.
* Snakemake Wrappers update to [3.3.6](https://snakemake-wrappers.readthedocs.io/en/v3.3.6/changelog.html)
* Github actions updated
* Use of lookup and queries instead of hand-made samples.csv parsing
* DAG available on readme

## Known bugs

* Report is unavailable as long as the TBD issue is opened on Snakemake


# 1.0.4

## Features

* More precise configuration file validation
* Labed added to MultiQC to make it more noticable

# 1.0.3

## Features

* Hidden keys in documentation are now described
* MultiQC report name changed in order to identify it more easily
* Example report available
* Disk usage reservation
* Explicit time/memory reservation

## Fixes

* Documentation pointing to fastp and not fastqc
* Typo in snakemake report main page

# 1.0.2

## Fixes

* csv.Sniffer having too much data to define delimiter

# 1.0.1

## Features

* Updated readme
* Updated report content
* Snakemake workflow catalog updates

# 1.0.0

## Features

* Run FastQC on single/pair ended fastq files
* Aggregate reports in MultiQC

## Known bug

* Using Snakemake 8.0.0+ the report is broken. Error is known.
