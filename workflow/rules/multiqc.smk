rule fastqc_multiqc_report:
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
        mem_mb=get_2gb_per_attempt,
        runtime=get_30min_per_attempt,
        disk=get_input_size_per_attempt_plus_1gb,
        tmpdir="tmp",
    params:
        extra=config.get("params", {}).get(
            "multiqc",
            "--module fastqc --zip-data-dir --verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=True,
    log:
        "logs/multiqc.log",
    benchmark:
        "benchmark/multiqc.tsv"
    wrapper:
        "v3.3.3/bio/multiqc"
