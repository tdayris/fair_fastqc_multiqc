Matierial and methods
=====================

Quality controls were performed with FastQC_ [#fastqcpaper]_. Possible contaminations
were assessed with FastqScreen_ [#fastqscreenpapaer]_. Quality repord produced during 
both trimming and mapping steps have been aggregated with MultiQC_ [#multiqcpaper]_. 

FastQ files formats were checked using Fastq_utils_ [#fastqutilspaper]_, additional
qualities were gathered using Seqkit_ [#seqkitpaper]_, Fastqinfo_ [#fastqinfopaper]_
and SeqTK_ [#seqtkpaper]_.

The whole pipeline_ [#fair_fastqc_multiqc_quote]_ was powered by Snakemake_ [#snakemakepaper]_,
and relies on fair_genome_indexer_ [#fair_genome_indexer_quote]_ pipeline.


.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#fastqscreenpapaer] Wingett, Steven W., and Simon Andrews. "FastQ Screen: A tool for multi-genome mapping and quality control." F1000Research 7 (2018).
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#fastqutilspaper] Nuno Fonseca, & Jonathan Manning. (2023). nunofonseca/fastq_utils: 0.25.2 (0.25.2). Zenodo. https://doi.org/10.5281/zenodo.7755574
.. [#seqkitpaper] Shen, Wei, et al. "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation." PloS one 11.10 (2016): e0163962.
.. [#fastqinfopaper] Kiu R, fastq-info: compute estimated sequencing depth (coverage) of prokaryotic genomes
.. [#seqtkpaper] Li, Heng. "seqtk Toolkit for processing sequences in FASTA/Q formats." GitHub 767 (2012): 69.
.. [#fair_fastqc_multiqc_quote] Dayris, T. (2024). fair-fastqc-multiqc (Version 2.5.0) [Computer software]. https://github.com/tdayris/fair_fastqc_multiqc
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.
.. [#fair_genome_indexer_quote] Dayris, T. (2024). fair-genome-indexer (Version 3.9.3) [Computer software]. https://github.com/tdayris/fair_genome_indexer


.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_fastqc_multiqc
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_fastqc_multiqc
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/fastqc.html
.. _FastqScreen: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/fastq_screen.html
.. _SeqTK: https://snakemake-wrappers.readthedocs.io/en/v5.6.0/wrappers/seqtk.html
.. _Fastq_utils: https://github.com/nunofonseca/fastq_utils
.. _Seqkit: https://bioinf.shenwei.me/seqkit/
.. _Fastqinfo: https://github.com/raymondkiu/fastq-info
.. _pipeline: https://github.com/tdayris/fair_fastqc_multiqc
.. _fair_genome_indexer: https://github.com/tdayris/fair_genome_indexer


:Authors:
    Thibault Dayris

:Version: 2.5.4 of 2025-03-07
