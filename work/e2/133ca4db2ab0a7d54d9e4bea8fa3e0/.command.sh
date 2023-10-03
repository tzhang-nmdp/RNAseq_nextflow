#!/bin/bash -euo pipefail
samtools faidx genome_gfp.fasta
NUM_BASES=`gawk '{sum = sum + $2}END{if ((log(sum)/log(2))/2 - 1 > 14) {printf "%.0f", 14} else {printf "%.0f", (log(sum)/log(2))/2 - 1}}' genome_gfp.fasta.fai`

mkdir star
STAR \
    --runMode genomeGenerate \
    --genomeDir star/ \
    --genomeFastaFiles genome_gfp.fasta \
    --sjdbGTFfile genome_gfp.gtf \
    --runThreadN 2 \
    --genomeSAindexNbases $NUM_BASES \
    --limitGenomeGenerateRAM 6342450944 \


cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:PREPARE_GENOME:STAR_GENOMEGENERATE":
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*$//')
END_VERSIONS
