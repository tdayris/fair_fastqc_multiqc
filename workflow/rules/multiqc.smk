rule fair_fastqc_multiqc_multiqc_report:
    input:
        unpack(get_multiqc_report_input),
    output:
        report(
            "results/QC/MultiQC_FastQC.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Raw",
            },
        ),
        "results/QC/MultiQC_FastQC_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * (1024 * 2),
        runtime=lambda wildcards, input, attempt: attempt
        * 30
        * max(1, int(input.size_mb / 1024)),
        tmpdir="tmp",
    params:
        extra=lookup(dpath="params/multiqc", within=config),
        use_input_files_only=True,
    log:
        "logs/fair_fastqc_multiqc/multiqc_report.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc/multiqc_report.tsv"
    wrapper:
        "v3.3.6/bio/multiqc"
