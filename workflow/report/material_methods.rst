Matierial and methods
=====================

Quality controls were performed with FastQC_ [#fastqcpaper]_. Quality repord produced during 
both trimming and mapping steps have been aggregated with MultiQC_ [#multiqcpaper]_. The 
whole pipeline was powered by Snakemake_ [#snakemakepaper]_.

.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.


.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_fastqc_multiqc
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_fastqc_multiqc
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v3.3.3/wrappers/fastqc.html


:Authors:
    Thibault Dayris

:Version: 1.0.0 of 01/11/2024
