Matierial and methods
=====================

Quality controls were performed with FastQC_ [#fastqcpaper]_. Possible contaminations
were assessed with FastqScreen_ [#fastqscreenpapaer]_
Quality repord produced during 
both trimming and mapping steps have been aggregated with MultiQC_ [#multiqcpaper]_. The 
whole pipeline was powered by Snakemake_ [#snakemakepaper]_.

.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#fastqscreenpapaer] Wingett, Steven W., and Simon Andrews. "FastQ Screen: A tool for multi-genome mapping and quality control." F1000Research 7 (2018).
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.


.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_fastqc_multiqc
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_fastqc_multiqc
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/fastqc.html
.. _FastqScreen: https://snakemake-wrappers.readthedocs.io/en/v3.4.1/wrappers/fastq_screen.html


:Authors:
    Thibault Dayris

:Version: 2.1.0 of 03/01/2024
