import csv
import functools
import os
import pandas
import snakemake
import snakemake.utils

from typing import Any, Callable, NamedTuple


snakemake_min_version: str = "8.13.0"
snakemake.utils.min_version(snakemake_min_version)

snakemake_docker_image: str = "docker://snakemake/snakemake:v8.13.0"


container: snakemake_docker_image


# Load and check configuration file
default_config_file: str = "config/config.yaml"


configfile: default_config_file


snakemake.utils.validate(config, "../schemas/config.schema.yaml")


# Load and check samples properties table
def load_table(path: str) -> pandas.DataFrame:
    """
    Load a table in memory, automatically inferring column separators

    Parameters:
    path (str): Path to the table to be loaded

    Return
    (pandas.DataFrame): The loaded table
    """
    with open(path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
        table_stream.seek(0)

    # Load table
    table: pandas.DataFrame = pandas.read_csv(
        path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )

    # Remove empty lines
    table = table.where(table.notnull(), None)

    return table


def used_genomes(
    genomes: pandas.DataFrame, samples: pandas.DataFrame | None = None
) -> tuple[str]:
    """
    Reduce the number of genomes to download to the strict minimum
    """
    if samples is None:
        return genomes

    return genomes.loc[
        genomes.species.isin(samples.species.tolist())
        & genomes.build.isin(samples.build.tolist())
        & genomes.release.isin(samples.release.tolist())
    ]


def load_genomes(
    path: str | None = None, samples: pandas.DataFrame | None = None
) -> pandas.DataFrame:
    """
    Load genome file, build it if genome file is missing and samples is not None.

    Parameters:
    path    (str)               : Path to genome file
    samples (pandas.DataFrame)  : Loaded samples
    """
    if path is not None:
        genomes: pandas.DataFrame = load_table(path)

        if samples is not None:
            genomes = used_genomes(genomes, samples)
        return genomes

    elif samples is not None:
        return samples[["species", "build", "release"]].drop_duplicates(
            ignore_index=True
        )

    raise ValueError(
        "Provide either a path to a genome file, or a loaded samples table"
    )


# Load and check samples properties tables
try:
    if (samples is None) or samples.empty():
        sample_table_path: str = config.get("samples", "config/samples.csv")
        samples: pandas.DataFrame = load_table(sample_table_path)
except NameError:
    sample_table_path: str = config.get("samples", "config/samples.csv")
    samples: pandas.DataFrame = load_table(sample_table_path)

snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")


# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
try:
    if (genomes is None) or genomes.empty:
        genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)
except NameError:
    genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflows.rst"


snakemake_wrappers_prefix: str = "v3.12.0"
release_tuple: tuple[str] = tuple(set(genomes.release.tolist()))
build_tuple: tuple[str] = tuple(set(genomes.build.tolist()))
species_tuple: tuple[str] = tuple(set(genomes.species.tolist()))
datatype_tuple: tuple[str] = ("dna", "cdna", "all", "transcripts")
gxf_tuple: tuple[str] = ("gtf", "gff3")
id2name_tuple: tuple[str] = ("t2g", "id_to_gene")
tmp: str = f"{os.getcwd()}/tmp"
samples_id_tuple: tuple[str] = tuple(samples.sample_id)
stream_tuple: tuple[str] = ("1", "2")


wildcard_constraints:
    release=r"|".join(release_tuple),
    build=r"|".join(build_tuple),
    species=r"|".join(species_tuple),
    datatype=r"|".join(datatype_tuple),
    gxf=r"|".join(gxf_tuple),
    id2name=r"|".join(id2name_tuple),
    sample=r"|".join(samples_id_tuple),
    stream=r"|".join(stream_tuple),


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query: str = (
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
            wildcards=wildcards
        )
    )

    query_result: str | float = getattr(
        lookup(query=query, within=genomes), key, default
    )
    if query_result != query_result:
        # Then the result of the query is nan
        return default
    return query_result


def get_dna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fasta", default=default, genomes=genomes)


def get_cdna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fasta", default=default, genomes=genomes)


def get_transcripts_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fasta", default=default, genomes=genomes
    )


def select_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fasta(wildcards),
            "cdna": get_cdna_fasta(wildcards),
            "transcripts": get_transcripts_fasta(wildcards),
        },
    )


def get_dna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fai", default=default, genomes=genomes)


def get_cdna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fai", default=default, genomes=genomes)


def get_transcripts_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences index
    """
    default: str = (
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fai", default=default, genomes=genomes
    )


def select_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta index file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fai(wildcards),
            "cdna": get_cdna_fai(wildcards),
            "transcripts": get_transcripts_fai(wildcards),
        },
    )


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GTF formatted)
    """
    default: str = (
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_gff(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation (GFF3 formatted)
    """
    default: str = (
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gff3".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gff3", default=default, genomes=genomes)


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
