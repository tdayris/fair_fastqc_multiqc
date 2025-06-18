# coding: utf-8

"""Snakemake wrapper for Librarian"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(
    stdout=True,
    stderr=True,
    append=True,
)

shell(
    "{snakemake.input.launcher} "
    "--local {snakemake.params.extra} "
    "{snakemake.input.fastq} {log}"
)

expected_outputs = {
    "composition_map_svg": "librarian_compositions_map.svg",
    "composition_map_png": "librarian_compositions_map.png",
    "probability_map_svg": "librarian_probability_maps.svg",
    "probability_map_png": "librarian_probability_maps.png",
    "prediction_plot_svg": "librarian_prediction_plot.svg",
    "prediction_plot_png": "librarian_prediction_plot.png",
    "heatmap": "librarian_librarian_heatmap.txt",
    "report": "librarian_Librarian_analysis.html",
}

for snakemake_key, default_path in expected_outputs.items():
    tmp = snakemake.output.get(snakemake_key)
    if tmp:
        shell("mv --verbose {default_path} {tmp} {log}")
