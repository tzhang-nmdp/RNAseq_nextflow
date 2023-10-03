#!/bin/bash -euo pipefail
bbsplit.sh \
    -Xmx4915M \
    path=bbsplit \
    threads=2 \
    in=RAP1_IAA_30M_REP1_1_val_1.fq.gz in2=RAP1_IAA_30M_REP1_2_val_2.fq.gz \
    basename=RAP1_IAA_30M_REP1_%_#.fastq.gz \
    refstats=RAP1_IAA_30M_REP1.stats.txt \
    build=1 ambiguous2=all maxindel=150000

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:BBMAP_BBSPLIT":
    bbmap: $(bbversion.sh | grep -v "Duplicate cpuset")
END_VERSIONS
