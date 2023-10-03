#!/bin/bash -euo pipefail
salmon quant \
    --geneMap genome_gfp.gtf \
    --threads 2 \
    --libType=A \
    --index salmon \
    -1 WT_REP1.subsampled_R1.fastq.gz -2 WT_REP1.subsampled_R2.fastq.gz \
    --skipQuant \
    -o WT_REP1

if [ -f WT_REP1/aux_info/meta_info.json ]; then
    cp WT_REP1/aux_info/meta_info.json "WT_REP1_meta_info.json"
fi

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:FASTQ_SUBSAMPLE_FQ_SALMON:SALMON_QUANT":
    salmon: $(echo $(salmon --version) | sed -e "s/salmon //g")
END_VERSIONS
