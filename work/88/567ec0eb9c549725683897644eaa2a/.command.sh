#!/bin/bash -euo pipefail
fasta2gtf.py \
    -o gfp.gtf \
    -b gene_biotype \
    gfp.fa

cat genome.fasta gfp.fa > genome_gfp.fasta
cat genes.gtf gfp.gtf > genome_gfp.gtf

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:CAT_ADDITIONAL_FASTA":
    python: $(python --version | sed 's/Python //g')
END_VERSIONS
