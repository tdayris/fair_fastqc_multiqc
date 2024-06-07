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


snakemake.utils.min_version(fair_genome_indexer.snakemake_min_version)


container: fair_genome_indexer.snakemake_docker_image


# Load and check configuration file
configfile: fair_genome_indexer.default_config_file


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")

samples: pandas.DataFrame = fair_genome_indexer.load_table(sample_table_path)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    genomes: pandas.DataFrame = fair_genome_indexer.load_table(genome_table_path)
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


lookup_config: Callable[str, str] = fair_genome_indexer.lookup_config
snakemake_wrappers_prefix: str = fair_genome_indexer.snakemake_wrappers_prefix
samples_id_tuple: tuple[str] = tuple(samples.sample_id)
stream_tuple: tuple[str] = ("1", "2")
tmp: str = fair_genome_indexer.tmp


wildcard_constraints:
    sample=r"|".join(samples_id_tuple),
    stream=r"|".join(stream_tuple),


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
