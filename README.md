# Rnaseq nextflow pipeline 

## Main functions 
1. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html))
2. Sub-sample FastQ files and auto-infer strandedness ([`fq`](https://github.com/stjude-rust-labs/fq), [`Salmon`](https://combine-lab.github.io/salmon/))
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. UMI extraction ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
5. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
6. Removal of genome contaminants ([`BBSplit`](http://seqanswers.com/forums/showthread.php?t=41288))
7. Removal of ribosomal RNA ([`SortMeRNA`](https://github.com/biocore/sortmerna))
8. Choice of multiple alignment and quantification routes:
   1. [`STAR`](https://github.com/alexdobin/STAR) -> [`Salmon`](https://combine-lab.github.io/salmon/)
   2. [`STAR`](https://github.com/alexdobin/STAR) -> [`RSEM`](https://github.com/deweylab/RSEM)
   3. [`HiSAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml) -> **NO QUANTIFICATION**
9. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
10. UMI-based deduplication ([`UMI-tools`](https://github.com/CGATOxford/UMI-tools))
11. Duplicate read marking ([`picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
12. Transcript assembly and quantification ([`StringTie`](https://ccb.jhu.edu/software/stringtie/))
13. Create bigWig coverage files ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
14. Extensive quality control:
    1. [`RSeQC`](http://rseqc.sourceforge.net/)
    2. [`Qualimap`](http://qualimap.bioinfo.cipf.es/)
    3. [`dupRadar`](https://bioconductor.org/packages/release/bioc/html/dupRadar.html)
    4. [`Preseq`](http://smithlabresearch.org/software/preseq/)
    5. [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
15. Pseudoalignment and quantification ([`Salmon`](https://combine-lab.github.io/salmon/) or ['Kallisto'](https://pachterlab.github.io/kallisto/); _optional_)
16. Present QC for raw read, alignment, gene biotype, sample similarity, and strand-specificity checks ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))
First, prepare a samplesheet with your input data that looks as follows:

### prepare the sample input
**samplesheet.csv**:
```csv
sample,fastq_1,fastq_2,strandedness
CONTROL_REP1,{sample_id}_R1_001.fastq.gz,{sample_id}_R2_001.fastq.gz,auto
CONTROL_REP1,{sample_id}_R1_001.fastq.gz,{sample_id}_R2_001.fastq.gz,auto
CONTROL_REP1,{sample_id}_R1_001.fastq.gz,{sample_id}_R2_001.fastq.gz,auto
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end). Rows with the same sample identifier are considered technical replicates and merged automatically. The strandedness refers to the library preparation and will be automatically inferred if set to `auto`.

### run Rnaseq nextflow:
```bash
nextflow run nf-core/rnaseq \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --genome GRCh37 \
    -profile <docker/singularity/.../institute>
```
