#!/bin/bash -euo pipefail
STAR \
    --genomeDir star \
    --readFilesIn input1/RAP1_UNINDUCED_REP1_primary.fastq.gz  \
    --runThreadN 2 \
    --outFileNamePrefix RAP1_UNINDUCED_REP1. \
     \
    --sjdbGTFfile genome_gfp.gtf \
    --outSAMattrRGline 'ID:RAP1_UNINDUCED_REP1'  'SM:RAP1_UNINDUCED_REP1'  \
    --quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outSAMstrandField intronMotif



if [ -f RAP1_UNINDUCED_REP1.Unmapped.out.mate1 ]; then
    mv RAP1_UNINDUCED_REP1.Unmapped.out.mate1 RAP1_UNINDUCED_REP1.unmapped_1.fastq
    gzip RAP1_UNINDUCED_REP1.unmapped_1.fastq
fi
if [ -f RAP1_UNINDUCED_REP1.Unmapped.out.mate2 ]; then
    mv RAP1_UNINDUCED_REP1.Unmapped.out.mate2 RAP1_UNINDUCED_REP1.unmapped_2.fastq
    gzip RAP1_UNINDUCED_REP1.unmapped_2.fastq
fi

cat <<-END_VERSIONS > versions.yml
"NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN":
    star: $(STAR --version | sed -e "s/STAR_//g")
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    gawk: $(echo $(gawk --version 2>&1) | sed 's/^.*GNU Awk //; s/, .*$//')
END_VERSIONS
