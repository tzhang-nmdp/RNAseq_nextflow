#!/bin/bash -euo pipefail
gunzip \
    -f \
     \
    gfp.fa.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GUNZIP_ADDITIONAL_FASTA":
    gunzip: $(echo $(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*$//')
END_VERSIONS
