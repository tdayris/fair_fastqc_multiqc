import csv
import functools
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from snakemake.common.tbdstring import TBDString
from typing import Any

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
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

snakemake_wrappers_version: str = "v3.0.0"


report: "../report/workflows.rst"


stream_list: list[str] = ["1", "2"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    stream=r"|".join(stream_list),


# Memory and time reservation
def get_resources_per_attempt(
    wildcards: snakemake.io.Wildcards, input: snakemake.io.InputFiles, attempt: int = 1, multiplier: int = 1, base: int = 0
) -> int:
    """
    Return the amount of resources needed per GB of input.

    Parameters:
    wildcards  (snakemake.io.Wildcards) : Snakemake signature requires this parameter
    input      (snakemake.io.InputFiles): Snakemake input files
    attempt    (int)                    : The # of times the calling rule has been restarted
    multiplier (int)                    : An arbitrary multiplier
    base       (int)                    : Minimal reservation


    Return:
    (int) The amount of resources needed (mb, minutes, etc)
    """
    return int((multiplier * attempt) + base)


get_2gb_per_attempt = functools.partial(get_resources_per_attempt, multiplier=2048)
get_30min_per_attempt = functools.partial(get_resources_per_attempt, multiplier=30)
get_input_size_per_attempt_plus_1gb = functools.partial(get_resources_per_attempt, base=1024)


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (dict[str, str | None]):
    Sample information
    """
    result: str | None = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_fastqc_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str:
    """
    Return expected input files for FastQC, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their genome

    Return (str):
    Path to a fastq file, as required by FastQC's snakemake-wrapper
    """
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    downstream_file: str | None = sample_data.get("downstream_file")
    if "stream" in wildcards.keys():
        if wildcards.stream == "1":
            return {"fastq": sample_data["upstream_file"]}
        elif wildcards.stream == "2" and downstream_file:
            return {"fastq": sample_data["downstream_file"]}
        raise ValueError("Could not guess which fastq we're talking about")
    return {"fastq": sample_data["upstream_file"]}


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
    results: dict[str, list[str]] = {"fastqc": []}
    datatype: str = "dna"
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        sample_data: dict[str, str | None] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        if sample_data.get("downstream_file"):
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

    return results


def get_fair_fastqc_multiqc_target(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return the expected list of output files at the end of the pipeline

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : Configuration file

    Return (dict[str, List(str)]):
    Dictionnary of expected output files
    """
    return {
        "multiqc": [
            "results/QC/MultiQC_FastQC.html",
            "results/QC/MultiQC_FastQC_data.zip",
        ],
    }
