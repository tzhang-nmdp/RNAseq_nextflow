/* environment and parameter setup */
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                   } from '../modules'
include { QUALIMAP_RNASEQ             } from '../modules'
include { STRINGTIE_STRINGTIE         } from '../modules'
include { SUBREAD_FEATURECOUNTS       } from '../modules'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { FASTQ_SUBSAMPLE_FQ_SALMON        } from '../subworkflows'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows'
include { FASTQ_FASTQC_UMITOOLS_FASTP      } from '../subworkflows'
include { BAM_MARKDUPLICATES_PICARD        } from '../subworkflows'
include { BAM_RSEQC                        } from '../subworkflows'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../subworkflows'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../subworkflows'

params.input=""
params.indir=""
params.outdir=""

/* work flow */
workflow {

/* input */
ref_genome_index()

seq_alignment()

feature_summary()

multi_qc()

de_analysis()

/* output */

        }
process alignment_index()
    {
    /* input */

    if (aligner=="STAR")
        {
        script '''
         STAR \
            --runThreadN ${n} \
            --runMode genomeGenerate \
            --genomeDir gv_index/ \
            --genomeFastaFiles gencode.v44.transcripts.fa  \
            --sjdbGTFfile  gencode.v44.annotation.gtf \
            --sjdbOverhang 100 \
            --genomeChrBinNbits 10 
        '''
        } else if (aligner=="salmon") {
        script '''
        salmon index \
        -t ${salmon_ref}.fa.gz \
        -i ${salmon_ref}.fa_index
        '''
        }

    }

process alignment_index()
    {
    /* input */

    if (aligner=="STAR")
        {
         script '''
         STAR \
            --runThreadN ${n} \
            --genomeDir ${star_references} \
            --sjdbGTFfile ${eferences_gtf} \ 
            --sjdbOverhang 100 \
            --readFilesIn ${fastq_file} \
            --outFileNamePrefix ${sample_id} \
            --outSAMtype BAM SortedByCoordinate Unsorted
            '''   
        } else if (aligner=="salmon") {
            script '''
            salmon quant \
                -i ${salmon_ref}.fa_index \
                -l A \
                -1 ${fn}/${sample_id}_1.fastq.gz \
                -2 ${fn}/${sample_id}_2.fastq.gz \
                -p 8 --validateMappings -o quants/${sample_id}_quant   
            '''
        } else if (aligner=="salmon_align") {
            script '''
            salmon quant \
                -t ${salmon_ref}.fa \
                -l A \
                -a ${fn}/${sample_id}.bam \
                -o quants/${sample_id}_quant   
            '''
        } 

    /* output */


    }

process feature_summary()
    {
    if (extracter=="RSEM") 
        {
        script '''
        rsem-calculate-expression 
        --bam \
        --no-bam-output -p 10 \
        --paired-end --forward-prob 0 \
        ${sample_id}.bam \
        ${rsem_ref} \
        ${rsem_outdir}
        '''
        } else if (extracter=="tximport") 
        {  
        script '''
        Rscript tximport_summary.R -i  -o 
        '''  
        } else if (extracter=="tximeta") 
        {  
        script '''
        Rscript tximeta_summary.R -i  -o 
        '''  
        }      
    }

process multi_qc()
    {

        script '''
        Rscript multi_qc.R -i  -o 
        '''  
      
    }    

process de_analysis()
    {

        script '''
        Rscript de_analysis.R -i  -o 
        '''          
    }    
