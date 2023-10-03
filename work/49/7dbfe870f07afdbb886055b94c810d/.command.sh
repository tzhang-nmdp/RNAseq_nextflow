#!/bin/bash -euo pipefail
printf "%s %s\n" SRR6357076_1.fastq.gz RAP1_IAA_30M_REP1_1.gz SRR6357076_2.fastq.gz RAP1_IAA_30M_REP1_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 2 RAP1_IAA_30M_REP1_1.gz RAP1_IAA_30M_REP1_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
