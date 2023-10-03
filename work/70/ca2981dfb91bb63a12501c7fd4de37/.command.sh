#!/bin/bash -euo pipefail
printf "%s %s\n" SRR6357072_1.fastq.gz WT_REP2_1.gz SRR6357072_2.fastq.gz WT_REP2_2.gz | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done
fastqc --quiet --threads 2 WT_REP2_1.gz WT_REP2_2.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_FASTQC_UMITOOLS_TRIMGALORE:FASTQC":
    fastqc: $( fastqc --version | sed -e "s/FastQC v//g" )
END_VERSIONS
