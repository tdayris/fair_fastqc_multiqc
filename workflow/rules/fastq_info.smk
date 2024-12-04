rule fair_fastqc_multiqc_fastq_info_download:
    output:
        temp("tmp/fair_fastqc_multiqc_fastq_info_download/fastqinfo.sh"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_fastq_info_download.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_info_download.tsv"
    params:
        address="https://raw.githubusercontent.com/raymondkiu/fastq-info/refs/heads/master/bin/fastqinfo-2.0.sh",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_info_download",
            default="--verbose",
        ),
    conda:
        "../envs/bash.yaml"
    shell:
        "wget {params.extra} {params.address} --output-document {output} > {log} 2>&1"


rule fair_fastqc_multiqc_fastq_info_process_pair_ended:
    input:
        launcher="tmp/fair_fastqc_multiqc_fastq_info_download/fastqinfo.sh",
        fasta=lambda wildcards: select_fasta(wildcards),
        r1="tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.1.fastq.gz",
        r2="tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.2.fastq.gz",
    output:
        "results/QC/fastqinfo/{sample}.{species}.{build}.{release}.{datatype}.pe.tsv",
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_fastq_info_process_pair_ended/{species}.{build}.{release}/{sample}.{datatype}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_info_process_pair_ended/{species}.{build}.{release}/{sample}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_info_process_pair_ended",
            default="",
        ),
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/fair_fastqc_multiqc_fastqinfo.py"


use rule fair_fastqc_multiqc_fastq_info_process_pair_ended as fair_fastqc_multiqc_fastq_info_process_single_ended with:
    input:
        launcher="tmp/fair_fastqc_multiqc_fastq_info_download/fastqinfo.sh",
        fasta=lambda wildcards: select_fasta(wildcards),
        r1="tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
    output:
        "results/QC/fastqinfo/{sample}.{species}.{build}.{release}.{datatype}.se.tsv",
    shadow:
        "minimal"
    log:
        "logs/fair_fastqc_multiqc_fastq_info_process_single_ended/{species}.{build}.{release}/{sample}.{datatype}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_info_process_single_ended/{species}.{build}.{release}/{sample}.{datatype}.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_info_process_single_ended",
            default="",
        ),
