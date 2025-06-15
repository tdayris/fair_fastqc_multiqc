rule fair_fastqc_multiqc_seqkit_stats_pair_ended:
    input:
        expand(
            "tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
            stream=stream_tuple,
            allow_missing=True,
        ),
    output:
        stats=temp("tmp/fair_fastqc_multiqc_seqkit_stats_pair_ended/{sample}.tsv"),
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_seqkit_stats_pair_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqkit_stats_pair_ended/{sample}.log"
    params:
        command="stats",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqkit_stats_pair_ended",
            default=" --all --tabular ",
        ),
    wrapper:
        "v7.0.0/bio/seqkit"


use rule fair_fastqc_multiqc_seqkit_stats_pair_ended as fair_fastqc_multiqc_seqkit_stats_single_ended with:
    input:
        ["tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz"],
    output:
        stats=temp("tmp/fair_fastqc_multiqc_seqkit_stats_single_ended/{sample}.tsv"),
    log:
        "logs/fair_fastqc_multiqc_seqkit_stats_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqkit_stats_single_ended/{sample}.log"
    params:
        command="stats",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqkit_stats_single_ended",
            default=" --all --tabular ",
        ),


rule fair_fastqc_multiqc_aggregate_seqkit_reports:
    input:
        pe=expand(
            "tmp/fair_fastqc_multiqc_seqkit_stats_pair_ended/{sample.sample_id}.tsv",
            sample=get_pair_ended_samples(),
        ),
        se=expand(
            "tmp/fair_fastqc_multiqc_seqkit_stats_single_ended/{sample.sample_id}.tsv",
            sample=get_single_ended_samples(),
        ),
    output:
        stats=report(
            "results/QC/SeqKit.Stats.csv",
            caption="../report/seqkit_stats.rst",
            category="Quality Controls",
            subcategory="Raw",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_aggregate_seqkit_reports.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_aggregate_seqkit_reports.tsv"
    params:
        subcommand="cat rows",
        extra="",
    wrapper:
        "v7.0.0/utils/xsv"
