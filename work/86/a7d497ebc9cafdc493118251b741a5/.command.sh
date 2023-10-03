#!/bin/bash -euo pipefail
bbsplit.sh \
    -Xmx4915M \
    path=bbsplit \
    threads=2 \
    in=RAP1_UNINDUCED_REP2_trimmed.fq.gz \
    basename=RAP1_UNINDUCED_REP2_%.fastq.gz \
    refstats=RAP1_UNINDUCED_REP2.stats.txt \
    build=1 ambiguous2=all maxindel=150000

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:BBMAP_BBSPLIT":
    bbmap: $(bbversion.sh | grep -v "Duplicate cpuset")
END_VERSIONS
