rule fair_fastqc_multiqc_fastq_screen_single_ended:
    input:
        sample="tmp/fair_fastqc_multiqc/link_or_concat_single_ended_input/{sample}.fastq.gz",
    output:
        txt=temp(
            "tmp/fair_fastqc_multiqc/fastq_screen_single_ended/{sample}.fastq_screen.txt"
        ),
        tmp="results/QC/fastq_screen/{sample}.fastq_screen.png",
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, input, attempt: attempt
        * 15
        * max(1, int(input.size_mb / 1024)),
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc/fastq_screen_single_ended/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc/fastq_screen_single_ended/{sample}.tsv"
    params:
        subset=dlookup(
            dpath="params/fair_fastqc_multiqc/fastq_screen/subset",
            within=config,
            default=10000,
        ),
        aligner=dlookup(
            dpath="params/fair_fastqc_multiqc/fastq_screen/aligner",
            within=config,
            default="bowtie2",
        ),
        fastq_screen_config=dlookup(
            dpath="params/fair_fastqc_multiqc/fastq_screen/fastq_screen_config",
            within=config,
        ),
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/fastq_screen"


use rule fair_fastqc_multiqc_fastq_screen_single_ended as fair_fastqc_multiqc_fastq_screen_pair_ended with:
    input:
        sample="tmp/fair_fastqc_multiqc/link_or_concat_pair_ended_input/{sample}.{stream}.fastq.gz",
    output:
        txt=temp(
            "tmp/fair_fastqc_multiqc/fastq_screen_pair_ended/{sample}.{stream}.fastq_screen.txt"
        ),
        tmp="results/QC/fastq_screen/{sample}.{stream}.fastq_screen.png",
    log:
        "logs/fair_fastqc_multiqc/fastq_screen_pair_ended/{sample}.{stream}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc/fastq_screen_pair_ended/{sample}.{stream}.tsv"
