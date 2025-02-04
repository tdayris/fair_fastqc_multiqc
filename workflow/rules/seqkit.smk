rule fair_fastqc_multiqc_seqkit_stats_pair_ended:
    input:
        expand(
            "tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample.sample_id}.{stream}.fastq.gz",
            sample=lookup(query="downstream_file == downstream_file", within=samples),
            stream=stream_tuple,
        ),
    output:
        stats=report(
            "results/QC/SeqKit.Stats.pe.tsv",
            caption="../report/seqkit_stats.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "txt",
                "library": "pair_ended",
            },
        ),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_seqkit_stats_pair_ended.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqkit_stats_pair_ended.log"
    params:
        command="stats",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqkit_stats_pair_ended",
            default="",
        ),
    wrapper:
        "v5.6.0/bio/seqkit"


use rule fair_fastqc_multiqc_seqkit_stats_pair_ended as fair_fastqc_multiqc_seqkit_stats_single_ended with:
    input:
        expand(
            "tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample.sample_id}.fastq.gz",
            sample=lookup(query="downstream_file != downstream_file", within=samples),
        ),
    output:
        stats=report(
            "results/QC/SeqKit.Stats.se.tsv",
            caption="../report/seqkit_stats.rst",
            category="Quality Controls",
            subcategory="Raw",
            labels={
                "report": "txt",
                "library": "single_ended",
            },
        ),
    log:
        "logs/fair_fastqc_multiqc_seqkit_stats_single_ended.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqkit_stats_single_ended.log"
    params:
        command="stats",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqkit_stats_single_ended",
            default="",
        ),
