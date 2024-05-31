import csv
import functools
import os
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from pathlib import Path
from snakemake.common.tbdstring import TBDString
from typing import Any, Callable, NamedTuple

snakemake.utils.min_version("8.2.0")


container: "docker://snakemake/snakemake:v8.5.3"


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


snakemake_wrappers_prefix: str = "v3.10.2"
stream_list: list[str] = ["1", "2"]
tmp: str = f"{os.getcwd()}/tmp"


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    stream=r"|".join(stream_list),


def lookup_config(
    dpath: Callable | str,
    default: str | bool | None = None,
    config: dict[str, Any] = config,
) -> str | bool | None:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value

    Parameters:
    dpath       (str)               : Signature within the configuration file
    default     (str | bool | None) : Default value to return in case of missing value
    config      (dict[str, Any]     : Configuration file

    Return
    (str | bool | None) : The value containted in the configuration file, or provided default value
    """
    value: str | bool | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def get_single_ended_samples(samples: pandas.DataFrame = samples) -> NamedTuple:
    """
    Return the list of single ended samples, as a NameTuple usable
    in expand/collect function

    Parameters:
    samples (pandas.DataFrame): Sample information

    Return: NamedTuple of single ended samples
    """
    return lookup(query="downstream_file != downstream_file", within=samples)


def get_pair_ended_samples(samples: pandas.DataFrame = samples) -> NamedTuple:
    """
    Return the list of pair ended samples, as a NameTuple usable
    in expand/collect function

    Parameters:
    samples (pandas.DataFrame): Sample information

    Return: NamedTuple of pair ended samples
    """
    return lookup(query="downstream_file == downstream_file", within=samples)
