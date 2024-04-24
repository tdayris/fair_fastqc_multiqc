SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR: linter_info.txt format.txt pipeline.txt

MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# List of snakefiles to lint
snakefiles := \
	"../workflow/Snakefile" \
	"../workflow/rules/common.smk" \
	"../workflow/rules/fastqc.smk" \
	"../workflow/rules/fastq_screen.smk" \
	"../workflow/rules/multiqc.smk" \
	"../workflow/rules/prepare_input.smk" \


python_scripts := \
	"../workflow/scripts/link_or_concat.py"

.PHONY: all
all: report.txt


dag.png:
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--rulegraph --forceall | dot -Tpng > dag.png && \
	cp --verbose dag.png ../dag.png


report.txt: pipeline.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile --profile local \
		--report report.zip 2>&1 | tee report.txt


pipeline.txt: linter.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--cores 2 \
		--profile local \
		2>&1 | tee pipeline.txt


linter.txt: format.txt black.txt
	mamba activate test-snakemake-wrappers && \
	snakemake -s ../workflow/Snakefile \
		--lint 2>&1 | tee linter.txt


format.txt:
	mamba activate test-snakemake-wrappers && \
	snakefmt $(snakefiles) 2>&1 | tee format.txt

black.txt:
	mamba activate test-snakemake-wrappers && \
	black $(python_scripts) 2>&1 | tee black.txt


.PHONY: clean
clean:
	rm -rvf linter.txt \
		   format.txt \
		   reference \
		   results \
		   tmp \
		   logs \
		   benchmark \
		   report.html \
		   report.txt \
		   pipeline.txt \
		   genomes.csv \
		   dag.png \
		   report.zip \
		   report \
		   black.txt \
		   .wget-hsts \
		   .conda \
		   .cache \
		   .java \
		   .snakemake
