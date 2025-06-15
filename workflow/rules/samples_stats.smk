rule fair_fastqc_multiqc_samples_stats:
    input:
        yaml="reference/genomes/{species}.{build}.{release}.{datatype}.statistics.yaml",
        stats="results/QC/SeqKit.Stats.csv",
    output:
        yaml=temp(
            "tmp/fair_fastqc_multiqc_samples_stats/{sample}.{species}.{build}.{release}.{datatype}.statistics.yaml"
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_samples_stats/{sample}.{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_samples_stats/{sample}.{species}.{build}.{release}.{datatype}.tsv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_fastqc_multiqc_samples_stats.py"
