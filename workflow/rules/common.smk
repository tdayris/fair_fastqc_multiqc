import csv
import functools
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from snakemake.common.tbdstring import TBDString
from typing import Any, NamedTuple

snakemake.utils.min_version("8.1.0")

# containerized: "docker://snakemake/snakemake:latest"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.readline())
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.readline())
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


stream_list: list[str] = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    stream=r"|".join(stream_list),


def get_multiqc_report_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, list[str]]:
    """
    Return expected input files for MultiQC report, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Dictionnary of all input files as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, list[str]] = {
        "fastqc_single_ended": collect(
            "results/QC/report_pe/{single_ended_data.sample_id}_fastqc.zip",
            single_ended_data=lookup(
                query="downstream_file != downstream_file", within=samples
            ),
        ),
        "fastqc_pair_ended": collect(
            "results/QC/report_pe/{pair_ended_data.sample_id}.{stream}_fastqc.zip",
            pair_ended_data=lookup(
                query="downstream_file == downstream_file", within=samples
            ),
            stream=stream_list,
        ),
        "fastq_screen_single_ended": collect(
            "tmp/fair_fastqc_multiqc/fastq_screen_single_ended/{single_ended_data.sample_id}.fastq_screen.txt",
            single_ended_data=lookup(
                query="downstream_file != downstream_file", within=samples
            ),
        ),
        "fastq_screen_pair_ended": collect(
            "tmp/fair_fastqc_multiqc/fastq_screen_pair_ended/{pair_ended_data.sample_id}.{stream}.fastq_screen.txt",
            pair_ended_data=lookup(
                query="downstream_file == downstream_file", within=samples
            ),
            stream=stream_list,
        ),
    }

    if not config.get("params", {}).get("fastq_screen", {}).get("fastq_screen_config"):
        del results["fastq_screen_single_ended"]
        del results["fastq_screen_pair_ended"]

    return results
