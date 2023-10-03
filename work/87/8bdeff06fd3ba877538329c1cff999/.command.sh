#!/bin/bash -euo pipefail
cat input1/SRR6357074_1.fastq.gz input2/SRR6357075_1.fastq.gz > RAP1_UNINDUCED_REP2.merged.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:CAT_FASTQ":
    cat: $(echo $(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*$//')
END_VERSIONS
