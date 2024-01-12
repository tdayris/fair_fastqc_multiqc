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
    ├── MultiQC_data.zip
    ├── MultiQC.html
    ├── report_pe
    │   └── YYY.html
    └── report_se
        └── YYY.html



+---------------+---------------------+-----------------------------------------------+
| Directory     | File Extension      | Content                                       |
+===============+=====================+===============================================+
| QC            | `MultiQC_data.zip`  | Zipped figures and tables                     |
+               +---------------------+-----------------------------------------------+
|               | `MultiQC.html`      | Complete quality report, includes all samples |
+---------------+---------------------+-----------------------------------------------+
| QC/report_pe  | `YYY.html`          | Sequence quality report for PE sample `YYY`   |
+---------------+---------------------+-----------------------------------------------+
| QC/report_se  | `YYY.html`          | Sequence quality report for SE sample `YYY`   |
+---------------+---------------------+-----------------------------------------------+