# -*- coding: utf-8 -*-

"""Snakemake wrapper for MultiQC configuration file"""

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2024, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"

import yaml
from typing import Any

mqc_config: dict[str, Any] = {
    "title": "Raw quality control report",
    "subtitle": "Produced on raw fastq recieved from sequencer",
    "intro_text": (
        "This pipeline building this report analyses all samples "
        "according to the same parameters, not taking the "
        "wet-lab experimental design, nor sample organisms "
        "into account."
    ),
    "report_comment": (
        "This report was generated using: "
        "https://github.com/tdayris/fair_fastqc_multiqc"
    ),
    "show_analysis_paths": False,
    "show_analysis_time": False,
    "custom_logo": snakemake.input[0],
    "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
    "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
    "report_header_info": [
        {"Contact E-mail": "bigr@gustaveroussy.fr"},
        {"Applivation type": "Any"},
        {"Project Type": "Quality Control"},
    ],
    "software_versions": {
        "Pipeline": {
            "snakemake": "8.5.3",
            "snakemake-wrappers-utils": "0.6.2",
            "fair_fastqc_multiqc": "2.4.2",
        },
    },
    "disable_version_detection": True,
    "run_modules": [
        "fastqc",
        "fastq_screen",
    ],
    "report_section_order": {
        "fastqc": {"order": 1000},
        "fastq_screen": {"order": "900"},
        "software_versions": {"order": "800"},
    },
}

if snakemake.params["extra"]:
    mqc_config.update(**snakemake.params["extra"])


with open(str(snakemake.output[0]), "w") as out_yaml_stream:
    out_yaml_stream.write(yaml.dump(mqc_config, default_flow_style=False))
