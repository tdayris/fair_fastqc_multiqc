rule fair_fastqc_multiqc_bigr_logo:
    output:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 100 + 200,
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
        temp("tmp/fair_fastqc_multiqc_multiqc_config/multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 * attempt * 100,
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
            default=None,
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
        fastq_screen_single_ended=collect(
            "tmp/fair_fastqc_multiqc_fastq_screen_single_ended/{single_ended_data.sample_id}.fastq_screen.txt",
            single_ended_data=get_single_ended_samples(),
        ),
        fastq_screen_pair_ended=collect(
            "tmp/fair_fastqc_multiqc_fastq_screen_pair_ended/{pair_ended_data.sample_id}.{stream}.fastq_screen.txt",
            pair_ended_data=get_pair_ended_samples(),
            stream=stream_tuple,
        ),
        librarian="results/QC/Librarian/librarian_heatmap.txt",
        config="tmp/fair_fastqc_multiqc_multiqc_config/multiqc_config.yaml",
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
        mem_mb=lambda wildcards, attempt: 800 + attempt * 200,
        runtime=lambda wildcards, input, attempt: attempt
        * 30
        * max(1, int(input.size_mb / 1024)),
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_multiqc/extra",
            default="--verbose --no-megaqc-upload --force --no-version-check --no-ai",
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
        "v7.0.0/bio/multiqc"
