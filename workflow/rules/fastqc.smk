rule fastqc:
    input:
        unpack(get_fastqc_input),
    output:
        html="results/{species}.{build}.{release}/QC/FastQC/{sample}.html",
        zip="results/{species}.{build}.{release}/QC/FastQC/{sample}_fastqc.zip",
    log:
        "log/fastqc/{species}.{build}.{release}/{sample}.log",
    benchmark:
        "benchmark/fastqc/{species}.{build}.{release}/{sample}.tsv",
    params:
        extra=config.get("params", {}).get("fastqc", ""),
    wrapper:
        "v3.3.3/bio/fastqc"