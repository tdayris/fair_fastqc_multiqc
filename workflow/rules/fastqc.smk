rule fair_fastqc_multiqc_fastqc_pair_ended:
    input:
        sample="tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
    output:
        html=report(
            "results/QC/report_pe/{sample}.{stream}.html",
            caption="../report/fastqc.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "pair_ended",
            },
        ),
        zip="results/QC/report_pe/{sample}.{stream}_fastqc.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, input, attempt: attempt
        * 10
        * max(1, int(input.size_mb / 1024)),
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_fastqc_pair_ended/{sample}.{stream}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastqc_pair_ended/{sample}.{stream}.tsv"
    params:
        extra=lookup_config(dpath="params/fair_fastqc_multiqc_fastqc", default=""),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/fastqc"


use rule fair_fastqc_multiqc_fastqc_pair_ended as fair_fastqc_multiqc_fastqc_single_ended with:
    input:
        sample="tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
    output:
        html=report(
            "results/QC/report_pe/{sample}.html",
            caption="../report/fastqc.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "html",
                "sample": "{sample}",
                "library": "single_ended",
            },
        ),
        zip="results/QC/report_pe/{sample}_fastqc.zip",
    log:
        "logs/fair_fastqc_multiqc_fastqc_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastqc_single_ended/{sample}.tsv"
