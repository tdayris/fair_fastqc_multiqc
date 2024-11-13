rule fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended:
    input:
        expand(
            "tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
            stream=stream_tuple,
            allow_missing=True,
        ),
    output:
        report(
            "results/QC/fastq_info/{sample}.pe.txt",
            caption="../report/fastq_utils_info.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "txt",
                "sample": "{sample}",
                "library": "pair_ended",
            },
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended/{sample}.tsv"
    params:
        command="fastq_info",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended",
            default="",
        ),
    conda:
        "../envs/fastq_utils.yaml"
    script:
        "../scripts/fair_fastqc_multiqc_fastq_utils.py"


use rule fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended as fair_fastqc_multiqc_fastq_utils_fastq_info_single_ended with:
    input:
        "tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
    output:
        report(
            "results/QC/fastq_info/{sample}.se.txt",
            caption="../report/fastq_utils_info.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "txt",
                "sample": "{sample}",
                "library": "single_ended",
            },
        ),
    log:
        "logs/fair_fastqc_multiqc_fastq_utils_fastq_info_single_ended/{sample}.log",
    benchmark:
        "logs/fair_fastqc_multiqc_fastq_utils_fastq_info_single_ended/{sample}.tsv"
    params:
        command="fastq_info",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_utils_fastq_info_pair_ended",
            default="",
        ),
