#!/bin/bash -euo pipefail
bbsplit.sh \
    -Xmx4915M \
    path=bbsplit \
    threads=2 \
    in=WT_REP2_1_val_1.fq.gz in2=WT_REP2_2_val_2.fq.gz \
    basename=WT_REP2_%_#.fastq.gz \
    refstats=WT_REP2.stats.txt \
    build=1 ambiguous2=all maxindel=150000

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:BBMAP_BBSPLIT":
    bbmap: $(bbversion.sh | grep -v "Duplicate cpuset")
END_VERSIONS
