use rule fair_fastqc_multiqc_fastq_info_download as fair_fastqc_multiqc_download_librarian with:
    output:
        temp("tmp/fair_fastqc_multiqc_download_librarian/librarian.tar.gz"),
    log:
        "logs/fair_fastqc_multiqc_download_librarian.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_download_librarian.tsv"
    params:
        address="https://github.com/DesmondWillowbrook/Librarian/releases/latest/download/librarian.tar.gz",
        extra=lookup_config(
            dpath="params/fair_fastqc_multiqc_fastq_info_download",
            default="--verbose",
        ),


rule fair_fastqc_multiqc_untar_librarian:
    input:
        "tmp/fair_fastqc_multiqc_download_librarian/librarian.tar.gz",
    output:
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/librarian",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian_analysis.Rmd",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/.DS_Store",
        directory(
            "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/compositions_umap/"
        ),
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/librarian_plotting_test_samples_server_220623.R",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/test_library_composition_5.txt",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/prediction_plot.svg",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/compositions_map.svg",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/probability_maps.png",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/Librarian.Rproj",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/.DS_Store",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/.Rhistory",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/librarian_heatmap.txt",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/librarian_offline_analysis.R",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/test_library_predictions.txt",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/probability_maps.svg",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/Librarian_offline_analysis.Rmd",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/compositions_map.png",
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/Librarian/prediction_plot.png",
        directory(
            "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/compositions_umap_results/"
        ),
        "tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/scripts/exec_analysis.sh",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_untar_librarian.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_untar_librarian.tsv"
    params:
        extra=" -xvzf ",
        prefix=subpath(output[0], ancestor=2),
    conda:
        "../envs/librarian.yaml"
    shell:
        "tar --directory {params.prefix} {params.extra} {input} > {log} 2>&1"


rule fair_fastqc_multiqc_librarian_local_mode:
    input:
        launcher="tmp/fair_fastqc_multiqc_download_librarian/librarian_v1.1/librarian",
        fastq=branch(
            lambda wildcards: is_paired(wildcards, samples),
            then="tmp/fair_fastqc_multiqc_link_or_concat_pair_ended_input/{sample}.1.fastq.gz",
            otherwise="tmp/fair_fastqc_multiqc_link_or_concat_single_ended_input/{sample}.fastq.gz",
        ),
    output:
        composition_map_svg="results/QC/Librarian/Composition_map/{sample}_cmap.svg",
        composition_map_png="results/QC/Librarian/Composition_map/{sample}_cmap.png",
        probability_map_svg="results/QC/Librarian/Probability_map/{sample}_pmap.svg",
        probability_map_png="results/QC/Librarian/Probability_map/{sample}_pmap.png",
        prediction_plot_svg="results/QC/Librarian/Prediction_plot/{sample}_pred.svg",
        prediction_plot_png="results/QC/Librarian/Prediction_plot/{sample}_pred.png",
        heatmap=temp(
            "tmp/fair_fastqc_multiqc_librarian_local_mode/{sample}.heatmap.txt"
        ),
        report="results/QC/Librarian/Reports/{sample}.html",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_librarian_local_mode/{sample}.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_librarian_local_mode/{sample}.tsv"
    params:
        extra="",
    conda:
        "../envs/librarian.yaml"
    script:
        "../scripts/fair_fastqc_multiqc_librarian_local_mode.py"


rule fair_fastqc_multiqc_concat_librarian:
    input:
        expand(
            "tmp/fair_fastqc_multiqc_librarian_local_mode/{sample}.heatmap.txt",
            sample=samples.sample_id,
        ),
    output:
        "results/QC/Librarian/librarian_heatmap.csv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1_000,
        runtime=lambda wildcards, attempt: attempt * 15,
        tmpdir=tmp,
    log:
        "logs/fair_fastqc_multiqc_concat_librarian.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_concat_librarian.tsv"
    params:
        subcommand=" cat rows ",
        extra=" --delimiter $'\t'",
    wrapper:
        "v7.0.0/utils/xsv"


use rule fair_fastqc_multiqc_concat_librarian as fair_fastqc_multiqc_format_librarian_for_multiqc with:
    input:
        "results/QC/Librarian/librarian_heatmap.csv",
    output:
        temp("results/QC/Librarian/librarian_heatmap.txt"),
    log:
        "logs/fair_fastqc_multiqc_format_librarian_for_multiqc.log",
    benchmark:
        "benchmark/fair_fastqc_multiqc_format_librarian_for_multiqc.tsv"
    params:
        subcommand="fmt",
        extra=" --out-delimiter $'\t'",
