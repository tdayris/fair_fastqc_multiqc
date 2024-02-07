rule fair_fastqc_multiqc_fastq_screen_single_ended:
    input:
        unpack(get_fastqc_fastqscreen_input),
    output:
        txt=temp("tmp/fastq_screen/{sample}.fastq_screen.txt"),
        tmp="results/QC/fastq_screen/{sample}.fastq_screen.png",
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 10,
        tmpdir="tmp",
        slurm_partition=lambda wildcards, attempt: get_partition(wildcards, attempt, 10),
    log:
        "logs/fastq_screen/{sample}.log"
    benchmark:
        "benchmark/fastq_screen/{sample}.tsv"
    params:
        aligner=config.get("params", {}).get("fastq_screen", {}).get("aligner", "bowtie2"),
        subset=config.get("params", {}).get("fastq_screen", {}).get("subset", 100_000),
        fastq_screen_config=config.get("params", {}).get("fastq_screen", {}).get("fastq_screen_config"),
    wrapper:
        "v3.3.6/bio/fastq_screen"


use rule fair_fastqc_multiqc_fastq_screen_single_ended as fair_fastqc_multiqc_fastq_screen_pair_ended  with:
    output:
        txt=temp("tmp/fastq_screen/{sample}.{stream}.fastq_screen.txt"),
        tmp="results/QC/fastq_screen/{sample}.{stream}.fastq_screen.png",
    log:
        "logs/fastq_screen/{sample}.{stream}.log"
    benchmark:
        "benchmark/fastq_screen/{sample}.{stream}.tsv"
