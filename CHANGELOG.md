# 3.0.0

## Breaking change:

SeqKit is now used once for each file. This is due to very large cohorts taking 
too much space in hot-sotage of computing cluster. A concatenation is being done 
at the end to have a complete table. Pair-ended and Single-ended libraries are
mixed, leading to a change in output file signature, and a breaking change as
a consequence.

## Features:

* SeqKit per sample + Concatenation
* Per-sample statistics are available with genome information, read number, etc.
* FastqScreen is being done against all provided genomes if no FastqScreen
  configuration is provided in configuration file.
* fair_genome_indexer update to 3.10.0

# 2.5.6

## Features:

* Snakemake-wrappers update to 7.0.0
* fair_genome_indexer update to 3.9.8

## Fix:

* FastQC memory has a hard limit set in its memory. We ensure snakemake does not 
  overstep this limit.

# 2.5.5

## Features:

* Snakemake-wrappers update to 5.8.3
* fair_genome_indexer update to 3.9.7
* Remove fastq_screen config file

# 2.5.4

## Features:

* Default configuration update

# 2.5.3

## Fix:

* FastQScreen output key error fixed

# 2.5.2

## Features:

* Only raise warning on untouchable files
* Snakemake wrappers update to 5.6.0
* fair_genome_indexer update to 3.9.5

# 2.5.1

## Fix:

* Allow fastq to be untouchable from current node (iRODS, cold-storage, ...)

# 2.5.0

## Features:

* Seqtk added to check fastq files

# 2.4.3

## Features:

* Snakemake-wrappers update to 5.5.0
* Tasks update

# 2.4.2

## Features:

* Easy snakemake-wrappers update
* Easy conda envs update
* New testing pipeline with additional format checks
* Snakemake wrappers update to 5.3.0
* python, bash and Snakemake environment update

## Docs:

* Citation cff file added

# 2.4.1

## Features:

* More clear error on when a source file does not exists

## Documentation:

* Report update

# 2.4.0

## Features:

* fastqinfo used to validate and gather quality metrics over fastq files and maximum theoretical coverage.
* Seqkit used to gather quality metrics over fastq files
* fastq_utils used to gather quality metrics over fastq files

## Fix:

* Input file concatenation edge case

## Documentation:

* Update for the configuration section

# 2.3.6

## Features:

* Snakemake-wrappers update to 4.5.0

# 2.3.5

## Fix:

* fastq screen parameters

# 2.3.4

## Feature:

* drop compresison dependencies

# 2.3.3

## Feature:

* drop dependency

# 2.3.2

## Fix

* Blacklist fix in fair_genome_indexer

# 2.3.1

## Features

* Documentation udpate

# 2.3.0

## Features

* Snakemake updates now follow the updates from fair_genome_indexer
* Improt functions and globals from fair_genome_indexer

# 2.2.8

## Features

* Update snakemake wrappers to 3.10.2
* Report update
* More readable tmp and logs

# 2.2.7

## Fix

* Fix R2 copy issue

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
