#!/bin/bash -euo pipefail
samtools faidx genome_gfp.fasta
cut -f 1,2 genome_gfp.fasta.fai > genome_gfp.fasta.sizes

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:CUSTOM_GETCHROMSIZES":
    getchromsizes: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
