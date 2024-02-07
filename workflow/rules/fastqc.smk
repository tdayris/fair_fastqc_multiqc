rule fair_fastqc_multiqc_fastqc_pair_ended:
    input:
        unpack(get_fastqc_fastqscreen_input),
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
        runtime=lambda wildcards, attempt: attempt * 30,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 30),
    log:
        "logs/fastqc/{sample}.{stream}.log",
    benchmark:
        "benchmark/fastqc/{sample}.{stream}.tsv"
    params:
        extra=config.get("params", {}).get("fastqc", ""),
    wrapper:
        "v3.3.6/bio/fastqc"


use rule fair_fastqc_multiqc_fastqc_pair_ended as fair_fastqc_multiqc_fastqc_single_ended with:
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
        "logs/fastqc/{sample}.log",
    benchmark:
        "benchmark/fastqc/{sample}.tsv"
