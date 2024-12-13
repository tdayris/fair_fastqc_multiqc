Results
=======


Alongside with the report, you may find directories called `reference`,
and `results`.

Reference
---------

You shall find all genome-related files in it. Considering a genome named `XXX`,
the following files are present:

::

    reference/
    ├── XXX.all.vcf
    ├── XXX.cdna.fasta
    ├── XXX.cdna.fasta.fai
    ├── XXX.dna.dict
    ├── XXX.dna.fasta
    ├── XXX.dna.fasta.fai
    └── XXX.gtf


+---------------+-----------------------------+
| Extension     | Content                     |
+===============+=============================+
| `.gtf`        | Genome annotation           |
+---------------+-----------------------------+
| `.fasta`      | Genome sequences            |
+---------------+-----------------------------+
| `.fasta.fai`  | Genome sequences index      |
+---------------+-----------------------------+
| `.dict`       | Genome sequences dictionary |
+---------------+-----------------------------+
| `.vcf`        | Genome known variations     |
+---------------+-----------------------------+

These files are quite volumous and are not embeded in this HTML page. Please
find them directly on file system.


Results
-------

Given an samples called `YYY` and a genome called `XXX`,
the following files are present:

::

    QC
    ├── fastq_info
    │   ├── YYY.pe.txt
    │   └── YYY.se.txt
    ├── fastqinfo
    │   └── YYY.cdna.txt
    ├── MultiQC_FastQC_data.zip
    ├── MultiQC_FastQC.html
    ├── report_pe
    │   └── YYY.html
    ├── report_se
    │   └── YYY.html
    ├── Seqtk
    │   └── YYY.txt
    ├── SeqKit.Stats.pe.tsv
    └── SeqKit.Stats.se.tsv




+---------------+-----------------------+------------------------------------------------------+
| Directory     | File Extension        | Content                                              |
+===============+=======================+======================================================+
| QC            | `MultiQC_data.zip`    | Zipped figures and tables                            |
+               +-----------------------+------------------------------------------------------+
|               | `MultiQC.html`        | Complete quality report, includes all samples        |
+               +-----------------------+------------------------------------------------------+
|               | `SeqKit.Stats.se.tsv` | SeqKit statistics over all single ended samples      |
+               +-----------------------+------------------------------------------------------+
|               | `SeqKit.Stats.pe.tsv` | SeqKit statistics over all pair ended samples        |
+---------------+-----------------------+------------------------------------------------------+
| QC/fastq_info | `YYY.pe.tsv`          | Fastq format checks and QC for a pair ended sample   |
+               +-----------------------+------------------------------------------------------+
|               | `YYY.se.tsv`          | Fastq format checks and QC for a single ended sample |
+---------------+-----------------------+------------------------------------------------------+
| QC/fastqinfo  | `YYY.cdna.txt`        | Fastq format checks and QC relatively to cDNA        |
+               +-----------------------+------------------------------------------------------+
|               | `YYY.dna.txt`         | Fastq format checks and QC relatively to DNA         |
+               +-----------------------+------------------------------------------------------+
|               | `YYY.transcripts.txt` | Fastq format checks and QC relatively to transcripts |
+---------------+-----------------------+------------------------------------------------------+
| QC/report_pe  | `YYY.html`            | Sequence quality report for PE sample `YYY`          |
+---------------+-----------------------+------------------------------------------------------+
| QC/report_se  | `YYY.html`            | Sequence quality report for SE sample `YYY`          |
+---------------+-----------------------+------------------------------------------------------+
+---------------+-----------------------+------------------------------------------------------+
| QC/Seqtk      | `YYY.pe.txt`          | Fastq format checks and QC for a pair ended sample   |
+               +-----------------------+------------------------------------------------------+
|               | `YYY.tsv`             | Fastq format checks and QC for a single ended sample |
+---------------+-----------------------+------------------------------------------------------+
