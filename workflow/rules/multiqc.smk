rule fair_fastqc_multiqc_bigr_logo:
    output:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_fastqc_multiqc_bigr_logo/bigr_logo.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_bigr_logo/bigr_logo.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_genome_indexer_wget", default="--verbose"
        ),
        address="https://raw.githubusercontent.com/tdayris/fair_fastqc_multiqc/main/images/bigr_logo.png",
    conda:
        "../envs/bash.yaml"
    shell:
        "wget {params.extra} --output-document {output} {params.address} > {log} 2>&1"


rule fair_fastqc_multiqc_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        temp("tmp/fair_fastqc_multiqc_multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_fastqc_multiqc_multiqc_config.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_multiqc_config.tsv"
    params:
        extra=lambda wildcards, input: lookup_config(
            dpath="params/fair_fastqc_multiqc_multiqc_config",
            default={
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
                "custom_logo": input[0],
                "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
                "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
                "report_header_info": [
                    {"Contact E-mail": "bigr@gustaveroussy.fr"},
                    {"Applivation type": "Any"},
                    {"Project Type": "Quality Control"},
                ],
                "software_versions": {
                    "Quality controls": {
                        "fastqc": "1.12.1",
                        "fastq_screen": "0.15.3",
                        "bowtie2": "2.5.3",
                        "bowtie1": "1.3.1",
                        "multiqc": "1.21.0",
                    },
                    "Pipeline": {
                        "snakemake": "8.5.3",
                        "snakemake-wrappers-utils": "0.6.2",
                        "fair_fastqc_multiqc": "2.2.8",
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
            },
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_fastqc_multiqc_multiqc_config.py"


rule fair_fastqc_multiqc_multiqc_report:
    input:
        fastqc_single_ended=collect(
            "results/QC/report_pe/{single_ended_data.sample_id}_fastqc.zip",
            single_ended_data=get_single_ended_samples(),
        ),
        fastqc_pair_ended=collect(
            "results/QC/report_pe/{pair_ended_data.sample_id}.{stream}_fastqc.zip",
            pair_ended_data=get_pair_ended_samples(),
            stream=stream_tuple,
        ),
        fastq_screen_single_ended=branch(
            condition=use_fqscreen,
            then=collect(
                "tmp/fair_fastqc_multiqc/fastq_screen_single_ended/{single_ended_data.sample_id}.fastq_screen.txt",
                single_ended_data=get_single_ended_samples(),
            ),
            otherwise=[],
        ),
        fastq_screen_pair_ended=branch(
            condition=use_fqscreen,
            then=collect(
                "tmp/fair_fastqc_multiqc_fastq_screen_pair_ended/{pair_ended_data.sample_id}.{stream}.fastq_screen.txt",
                pair_ended_data=get_pair_ended_samples(),
                stream=stream_tuple,
            ),
            otherwise=[],
        ),
        config="tmp/fair_fastqc_multiqc_multiqc_config.yaml",
        logo="tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        report(
            "results/QC/MultiQC_FastQC.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Raw",
            },
        ),
        "results/QC/MultiQC_FastQC_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, input, attempt: attempt
        * 30
        * max(1, int(input.size_mb / 1024)),
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_multiqc/extra",
            default="--verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=lookup_config(
            dpath="params/fair_fastqc_multiqc_multiqc/use_input_file_only",
            default=True,
        ),
    log:
        "logs/fair_fastqc_multiqc_multiqc_report.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_multiqc_report.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"
