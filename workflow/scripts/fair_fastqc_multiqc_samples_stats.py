# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import os.path
import pandas
import tempfile
import typing
from snakemake.shell import shell
import yaml


def read_genome_stats(path: str) -> dict[str, typing.Any]:
    """
    Load genome yaml file in memory
    """
    with open(path, "r") as yaml_stream:
        return yaml.safe_load(yaml_stream)


def read_seqkit_qc(path: str) -> dict[str, typing.Any]:
    """
    Load seqkit in memory
    """
    df = pandas.read_csv(
        path,
        sep=",",
        header=0,
        index_col=0,
    )

    # Normalize sample ids
    df.index = [i.split("/")[-1].replace(".fastq.gz", "") for i in df.index]
    print(df)

    if len(df) == 1:
        # Single ended mode
        return df.to_dict(orient="index")
    else:
        # Pair ended mode
        return {
            str(snakemake.wildcards.sample): {
                "format": str(df.iloc[0].format),
                "type": str(df.iloc[0].type),
                "num_seqs": int(df.iloc[0].num_seqs),
                "sum_len": int(df["sum_len"].sum()),
                "min_len": int(df["min_len"].min()),
                "avg_len": float(df["avg_len"].mean()),
                "max_len": int(df["max_len"].max()),
                "Q1": float(df["Q1"].mean()),
                "Q2": float(df["Q2"].mean()),
                "Q3": float(df["Q3"].mean()),
                "sum_gap": int(df["sum_gap"].sum()),
                "N50": int(df["N50"].max()),
                "N50_num": int(df["N50_num"].max()),
                "Q20(%)": float(df["Q20(%)"].mean()),
                "Q30(%)": float(df["Q30(%)"].mean()),
                "AvgQual": float(df["AvgQual"].mean()),
                "GC(%)": float(df["GC(%)"].mean()),
                "sum_n": int(df["sum_n"].sum()),
            }
        }


stats = read_genome_stats(snakemake.input.yaml)
seqkit = read_seqkit_qc(snakemake.input.stats)
stats.update(seqkit)

with open(snakemake.output.yaml, "w") as yaml_stream:
    yaml_stream.write(yaml.dump(stats))
