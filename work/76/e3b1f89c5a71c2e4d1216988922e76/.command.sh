#!/bin/bash -euo pipefail
gunzip \
    -f \
     \
    genes.gtf.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:GUNZIP_GTF":
    gunzip: $(echo $(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*$//')
END_VERSIONS
