#!/bin/bash -euo pipefail
[ ! -f  WT_REP2_1.fastq.gz ] && ln -s SRR6357072_1.fastq.gz WT_REP2_1.fastq.gz
[ ! -f  WT_REP2_2.fastq.gz ] && ln -s SRR6357072_2.fastq.gz WT_REP2_2.fastq.gz
trim_galore \
    --fastqc_args '-t 2' \
    --cores 1 \
    --paired \
    --gzip \
    WT_REP2_1.fastq.gz \
    WT_REP2_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
