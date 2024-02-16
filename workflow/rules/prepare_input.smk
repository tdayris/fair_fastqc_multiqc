rule fair_fastqc_multiqc_link_or_concat_single_ended_input:
    output:
        temp(
            "tmp/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.fastq.gz"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024,
        runtime=lambda wildcards, input, attempt: attempt
        * max(1, int(input.size_mb / 1024))
        * 15,
        tmpdir="tmp",
    log:
        "logs/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.tsv"
    params:
        in_files=lambda wildcards: get_fair_fastqc_multiqc_link_or_concat_single_ended_input(
            wildcards
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/link_or_concat.py"


use rule fair_fastqc_multiqc_link_or_concat_single_ended_input as fair_fastqc_multiqc_link_or_concat_pair_ended_input with:
    output:
        temp(
            "tmp/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz"
        ),
    log:
        "logs/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.tsv"
