rule fair_fastqc_multiqc_fastq_screen_single_ended:
    input:
        sample="tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
        indexes=select_fastq_screen(files_only=True),
    output:
        txt=temp(
            "tmp/fair_fastqc_multiqc_fastq_screen_single_ended/{sample}.fastq_screen.txt"
        ),
        png="results/QC/fastq_screen/{sample}.fastq_screen.png",
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500,
        runtime=lambda wildcards, input, attempt: attempt
        * 15
        * max(1, int(input.size_mb / 1024)),
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_fastq_screen_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_screen_single_ended/{sample}.tsv"
    params:
        subset=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_screen_subset",
            default=10_000,
        ),
        aligner=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_screen_aligner",
            default="bowtie2",
        ),
        fastq_screen_config=select_fastq_screen(),
    wrapper:
        "v7.0.0/bio/fastq_screen"


use rule fair_fastqc_multiqc_fastq_screen_single_ended as fair_fastqc_multiqc_fastq_screen_pair_ended with:
    input:
        sample="tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
        indexes=select_fastq_screen(files_only=True),
    output:
        txt=temp(
            "tmp/fair_fastqc_multiqc_fastq_screen_pair_ended/{sample}.{stream}.fastq_screen.txt"
        ),
        png="results/QC/fastq_screen/{sample}.{stream}.fastq_screen.png",
    log:
        "logs/fair_fastqc_multiqc_fastq_screen_pair_ended/{sample}.{stream}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_fastq_screen_pair_ended/{sample}.{stream}.tsv"
