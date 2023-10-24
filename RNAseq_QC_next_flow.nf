#!/usr/local/bin/ nextflow

nextflow.enable.dsl=2
 
// MODULE: Installed directly from nf-core/modules
//
// include { CAT_FASTQ                   } from './modules'
// include { QUALIMAP_RNASEQ             } from './modules'
// include { SALMON_TX2GENE } from './modules'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
// include { FASTQ_FASTQC_UMITOOLS_FASTP } from './subworkflows'

// process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
// docker.enabled = true
// docker.runOptions = '-u $(id -u):$(id -g)'
// process.container = 'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0'
// singularity.enabled = true

/*  Comments are uninterpreted text included with the script.
    They are useful for describing complex parts of the workflow
    or providing useful information such as workflow usage.

    Usage:
       nextflow run wc.nf --input <input_file>

    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash. */
//  Single line comments start with a double slash // and finish on the same line

/*  Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator. */
params.sample_id = "test"
params.refdir = "/home/tzhang/Jupyter_lab_space/git_test/test_pipeline/reference"
params.indir = "/home/tzhang/Jupyter_lab_space/git_test/test_pipeline"
params.outdir = "/home/tzhang/Jupyter_lab_space/git_test/test_pipeline/work"
params.kmer = 31
params.cpus = 8
params.memory = '60.GB'
params.aligner="all"
 
text = """
    This is a multi-line string
    using triple quotes.
    """

//  The default workflow
workflow {

    //  Input data is received through channels
    ref=Channel.value('GRCh38')
    sample_list=['SRR1039509']
    sample_id=Channel.fromList(sample_list)   
    sample_id.view()       
    // read_pair_fq = Channel.fromFilePairs('*_{1,2}.fastq.gz')
    // read_pair_fq.view()     
    refdir = Channel.fromPath(params.refdir)
    refdir.view()
    indir = Channel.fromPath(params.indir)    
    indir.view()    
    outdir = Channel.fromPath(params.outdir)   
    outdir.view()    
    aligner= "all"
    extracter="tximport" 
    printf aligner

    // FASTP(sample_id, outdir)
    // alignment_index(aligner,refdir)

    // alignment_seq(aligner, sample_id, refdir, indir, outdir)

    feature_summary(extracter, sample_id, refdir, indir, outdir)

    // SALMON_TX2GENE()

    // qualimap_qc()

    // deseq2_qc()

/* output */

        }

//process_tuple_io_fastp.nf
process FASTP {
  input:
  val(sample_id)
  val(outdir)
  
  // output:
  // tuple val(sample_id), path("*FP*.fq.gz")
  
  script:
  """
  echo ${outdir}
  fastp \
   -i ${outdir}/${sample_id}_1.fastq.gz \
   -I ${outdir}/${sample_id}_2.fastq.gz \
   -o ${outdir}/${sample_id}_FP_R1.fq.gz \
   -O ${outdir}/${sample_id}_FP_R2.fq.gz
  """
}

process alignment_index()
    {
    input:
    val(aligner)
    path(refdir)

    script: 
    if (aligner=="STAR" )
        {
        """ 
        echo "STAR reference building" && date         
        STAR \\
            --runThreadN 8 \
            --runMode genomeGenerate \
            --genomeDir ${refdir}/star \
            --genomeFastaFiles ${refdir}/GRCh38.primary_assembly.genome.fa  \
            --sjdbGTFfile ${refdir}/gencode.v44.annotation.gtf \
            --sjdbOverhang 100 \
            --genomeChrBinNbits 10 
        """              
        } else if (aligner=="salmon" || aligner=="all" ) {
        """  
        echo "salmon reference building" && date          
        salmon index \
            -t ${refdir}/gencode.v44.transcripts.fa \
            -i ${refdir}/salmon \
            --kmer 31   
        """              
        } else if (aligner=="rsem" || aligner=="all" ) {
        """  
        echo "rsem reference building" && date          
         rsem-prepare-reference \
            --gtf ${refdir}/gencode.v44.annotation.gtf \
            --star \
            --star-path /usr/local/bin \
            -p 8 \
            ${refdir}/GRCh38.primary_assembly.genome.fa \
            ${refdir}/rsem
        """              
        }                         
    }

 
process alignment_seq()
    {
    
    input:
    val(aligner)    
    val(sample_id)
    val(refdir)
    val(indir)
    val(outdir)

    script:   
    if (aligner=="STAR" || aligner=="all" )
        {
        """  
        STAR \
            --runThreadN 4 \
            --quantMode TranscriptomeSAM \
            --genomeDir ${refdir}/star \
            --readFilesIn ${outdir}/${sample_id}_FP_R1.fq.gz \
                        ${outdir}/${sample_id}_FP_R2.fq.gz \
            --outFileNamePrefix ${outdir}/${sample_id} \
            --sjdbGTFfile ${refdir}/gencode.v44.annotation.gtf \
            --sjdbOverhang 100 \
            --outSAMtype BAM SortedByCoordinate Unsorted \
            --readFilesCommand zcat
        """                
        } else if (aligner=="salmon" || aligner=="all" ) {
        """             
            salmon quant \
                -i ${refdir}/salmon \
                -l A \
                -1 ${outdir}/${sample_id}_FP_R1.fq.gz \
                -2 ${outdir}/${sample_id}_FP_R2.fq.gz \
                -p 8 --validateMappings \
                -o ${outdir}/salmon/quants/${sample_id}_quant   
        """                   
        } else if (aligner=="salmon_align" || aligner=="all" ) {
        """             
            salmon quant \
                -t ${refdir}/salmon \
                -l A \
                -a ${outdir}/${sample_id}Aligned.out.bam \
                -o ${outdir}/salmon/quants/${sample_id}_quant 
        """                 
        } 
    }

process feature_summary()
    {
    input:
    val(extracter)    
    val(sample_id)
    val(refdir)
    val(indir)
    val(outdir)

    script: 
    if (extracter=="RSEM" || extracter=="all" ) 
        {
        """            
        rsem-calculate-expression \
        --bam \
        --no-bam-output -p 10 \
        --paired-end --forward-prob 0 \
        ${outdir}/star/${sample_id}
        Aligned.toTranscriptome.out.bam \
        ${refdir}/rsem/rsem \
        ${outdir}/rsem
        """          
        } else if (extracter=="tximport" || extracter=="all" ) {  
        """       
        python bin/salmon_tx2gene.py \
        --gtf ${refdir}/gencode.v44.annotation.gtf \
        --salmon ${outdir}/salmon/quants/SRR1039509_quant \
        --id gene_id \
        --extra gene_name \
        -o ${outdir}/salmon/quants/${sample_id}_salmon_tx2gene.tsv                   
        // Rscript ${indir}/bin/salmon_tximport.r \
        // NULL \
        // ${outdir}/salmon/quants/${sample_id}_quant
        // ${outdir}/salmon/quants/${sample_id}_quant
        """                
        } else if (extracter=="tximeta" || extracter=="all" ) {  
        """             
        Rscript bin/tximeta_summary.R -i  -o 
        """          
        }      
    }

// process qualimap_qc()
//     {
//         script: 
//         """
//             unset DISPLAY
//             mkdir -p tmp
//             export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
//             qualimap \\
//             --java-mem-size=${memory} \\
//             rnaseq \\
//             $args \\
//             -bam ${indir}/${sample_id}.bam \\
//             -gtf $gtf \\
//             -p $strandedness \\
//             $paired_end \\
//             -outdir $prefix
//         """  
      
//     }    

// process deseq2_qc()
//     {
//         script:
//         """
//         Rscript ${indir}/bin/deseq2_qc.r \\
//         --count_file $counts \\
//         --outdir ./ \\
//         --cores ${cpus}
//         """          
//     }    

   
