# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)
input_files = " ".join(
    [
        snakemake.input.get("r1", ""),
        snakemake.input.get("r2", ""),
        snakemake.input.get("fasta", ""),
    ]
)

shell(
    "bash {snakemake.input.launcher} {input_files} "
    "{extra} > {snakemake.output} {log}"
)
