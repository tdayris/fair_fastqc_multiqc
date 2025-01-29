rule fair_fastqc_multiqc_seqtk_fqchk_pair_ended:
    input:
        "tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
    output:
        "results/QC/Seqtk/{sample}.{stream}.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 500 + attempt * 100,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_seqtk_fqchk_pair_ended/{sample}.{stream}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqtk_fqchk_pair_ended/{sample}.{stream}.tsv"
    params:
        command="fqchk",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqtk_fqchk_pair_ended",
            default="-q0",
        ),
    wrapper:
        "v5.5.0/bio/seqtk"


use rule fair_fastqc_multiqc_seqtk_fqchk_pair_ended as fair_fastqc_multiqc_seqtk_fqchk_single_ended with:
    input:
        "tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
    output:
        "results/QC/Seqtk/{sample}.tsv",
    log:
        "logs/fair_fastqc_multiqc_seqtk_fqchk_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_seqtk_fqchk_single_ended/{sample}.tsv"
    params:
        command="fqchk",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_seqtk_fqchk_single_ended",
            default="-q0",
        ),
