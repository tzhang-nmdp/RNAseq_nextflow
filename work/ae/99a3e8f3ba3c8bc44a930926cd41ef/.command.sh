#!/bin/bash -euo pipefail
[ ! -f  RAP1_UNINDUCED_REP1.fastq.gz ] && ln -s SRR6357073_1.fastq.gz RAP1_UNINDUCED_REP1.fastq.gz
trim_galore \
    --fastqc_args '-t 2' \
    --cores 1 \
    --gzip \
    RAP1_UNINDUCED_REP1.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
