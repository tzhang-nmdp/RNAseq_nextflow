#!/bin/bash -euo pipefail
bbsplit.sh \
    -Xmx4915M \
    ref_primary=genome_gfp.fasta \
    ref_sarscov2=GCA_009858895.3_ASM985889v3_genomic.200409.fna ref_human=chr22_23800000-23980000.fa \
    path=bbsplit \
    threads=2 \
    build=1

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:BBMAP_BBSPLIT":
    bbmap: $(bbversion.sh | grep -v "Duplicate cpuset")
END_VERSIONS
