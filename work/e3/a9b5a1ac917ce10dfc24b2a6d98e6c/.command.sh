#!/bin/bash -euo pipefail
printf "%s %s\n" RAP1_UNINDUCED_REP2.merged.fastq.gz RAP1_UNINDUCED_REP2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 2 RAP1_UNINDUCED_REP2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
