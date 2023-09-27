#!/usr/local/bin/ nextflow

nextflow.enable.dsl=2

include { INDEX; QUANT; INDEX as SALMON_INDEX  } from './modules/rnaseq-tasks'

process.container = 'quay.io/biocontainers/salmon:1.5.2--h84f40af_0'
docker.enabled = true
docker.runOptions = '-u $(id -u):$(id -g)'
process.container = 'https://depot.galaxyproject.org/singularity/salmon:1.5.2--h84f40af_0'
singularity.enabled = true

/*  Comments are uninterpreted text included with the script.
    They are useful for describing complex parts of the workflow
    or providing useful information such as workflow usage.

    Usage:
       nextflow run wc.nf --input <input_file>

    Multi-line comments start with a slash asterisk /* and finish with an asterisk slash. */
//  Single line comments start with a double slash // and finish on the same line

/*  Workflow parameters are written as params.<parameter>
    and can be initialised using the `=` operator. */
params.input = "test_1.fastq.gz"
params.sleep = 2
text = """
    This is a multi-line string
    using triple quotes.
    """

//  The default workflow
workflow {

    //  Input data is received through channels
    ref=Channel.value('GRCh38')
    chrom_no=Channel.value(['chr1','chr2'])
    chrom_no.view()
    chrom_run=Channel.of('chr1','chr2')
    chrom_run.view()
    test_list=['chr1','chr2']
    chrom_test=Channel.fromList(test_list) 
    chrom_test.view()       
    input_ch = Channel.fromPath(params.input)
    input_ch.view()   
    read_pair_ch = Channel.fromFilePairs('test*_{1,2}.fastq.gz')
    read_pair_ch.view()   
/*    sra_ch =Channel.fromSRA('SRP043510')
    sra_ch.view() */

    /*  The script to execute is called by its process name,
        and input is provided between brackets. */
    NUM_LINES(input_ch)

    /*  Process output is accessed using the `out` channel.
        The channel operator view() is used to print
        process output to the terminal. */
    NUM_LINES.out.view()
}

/*  A Nextflow process block
    Process names are written, by convention, in uppercase.
    This convention is used to enhance workflow readability. */
process NUM_LINES {

    input:
    path read

    output:
    stdout

    script:
    /* Triple quote syntax """, Triple-single-quoted strings may span multiple lines. The content of the string can cross line boundaries without the need to split the string in several pieces and without concatenation or newline escape characters. */
    """
    sleep ${params.sleep}
    printf '${read}'
    gunzip -c ${read} | wc -l
    """
}

//process_python.nf
projectDir="***"

process PYSTUFF {
  script:
  """
  #!/usr/bin/env python
  import gzip

  reads = 0
  bases = 0

  with gzip.open('${projectDir}/data/yeast/reads/ref1_1.fq.gz', 'rb') as read:
      for id in read:
          seq = next(read)
          reads += 1
          bases += len(seq.strip())
          next(read)
          next(read)

  print("reads", reads)
  print("bases", bases)
  """
}

//process_conditional.nf
params.aligner = 'kallisto'
params.transcriptome = "$projectDir/data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
params.kmer = 31

process INDEX {
  script:
  if( params.aligner == 'kallisto' ) {
    """
    echo indexed using kallisto
    kallisto index -i index -k $params.kmer $params.transcriptome
    """
  }  
  else if( params.aligner == 'salmon' ) {
    """
    echo indexed using salmon
    salmon index -t $params.transcriptome -i index --kmer $params.kmer
    """
  }  
  else {
    """
    echo Unknown aligner $params.aligner >> err.log
    """
  }  
}

//process_output_file.nf
methods_ch = channel.of('salmon', 'kallisto')

process METHOD {
  input:
  val x

  output:
  path 'method.txt'

  """
  echo $x > method.txt
  """
}

workflow {
  METHOD(methods_ch)
  // use the view operator to display contents of the channel
  METHOD.out.view({ "Received: $it" })
}

//process_tuple_io_fastp.nf
process FASTP {
  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path("*FP*.fq.gz")
  
  script:
  """
  fastp \
   -i ${reads[0]} \
   -I ${reads[1]} \
   -o ${sample_id}_FP_R1.fq.gz \
   -O ${sample_id}_FP_R2.fq.gz
  """
}

reads_ch = Channel.fromFilePairs('data/yeast/reads/ref1_{1,2}.fq.gz')

workflow {
  FASTP(reads_ch)
  FASTP.out.view()
}

//process_when.nf
nextflow.enable.dsl=2

process CONDITIONAL {
  input:
  val chr

  when:
  chr <= 5

  script:
  """
  echo $chr
  """
}

chr_ch = channel.of(1..22)

workflow {
  CONDITIONAL(chr_ch)
}

//process_directive.nf
nextflow.enable.dsl=2

process PRINTCHR {
  tag "tagging with chr$chr"
  cpus 1
  echo true

  input:
  val chr

  script:
  """
  echo processing chromosome: $chr
  echo number of cpus $task.cpus
  """
}

chr_ch = channel.of(1..22, 'X', 'Y')

workflow {
  PRINTCHR(chr_ch)
}

//process_exercise_pubilishDir.nf
nextflow.enable.dsl=2

process INDEX {
   //add publishDir directive here
   input:
   path transcriptome

   output:
   path "index"

   script:
   """
   salmon index -t $transcriptome -i index
   """
}

params.transcriptome = "data/yeast/transcriptome/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
transcriptome_ch = channel.fromPath(params.transcriptome, checkIfExists: true)

workflow {
  INDEX(transcriptome_ch)
}

//workflow_01.nf
nextflow.enable.dsl=2

process INDEX {
    input:
      path transcriptome
    output:
      path 'index'
    script:
      """
      salmon index -t $transcriptome -i index
      """
}

 process QUANT {
    input:
      each  path(index)
      tuple(val(pair_id), path(reads))
    output:
      path pair_id
    script:
      """
      salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
      """
}

workflow {
    transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz',checkIfExists: true)
    read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz',checkIfExists: true)

    //index process takes 1 input channel as a argument
    index_ch = INDEX(transcriptome_ch)

    //quant channel takes 2 input channels as arguments
    QUANT( index_ch, read_pairs_ch ).view()
}

//workflow_02.nf
nextflow.enable.dsl=2

process INDEX {

  input:
  path transcriptome

  output:
  path 'index', emit: salmon_index

  script:
  """
  salmon index -t $transcriptome -i index
  """
}

process QUANT {
   input:
     each  path(index)
     tuple(val(pair_id), path(reads))
   output:
     path pair_id
   script:
     """
     salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
     """
}

workflow {
  transcriptome_ch = channel.fromPath('data/yeast/transcriptome/*.fa.gz')
  read_pairs_ch = channel.fromFilePairs('data/yeast/reads/*_{1,2}.fq.gz')
  
  //call INDEX process
  INDEX(transcriptome_ch)
  
  //access INDEX object named output
  QUANT(INDEX.out.salmon_index,read_pairs_ch).view()
}

chr_ch = channel
  .of( 1..22, 'X', 'Y' )
  .filter( Number )
  .filter(~/^1.*/)
  .filter { it < 5 }  
  .view()
autosomes_ch =chr_ch.filter( Number )
autosomes_ch.view()

chr = channel
  .of( 'chr1', 'chr2' )
  .map ({ it.replaceAll("chr","") })
chr.view()

ch = channel
    .of( 1, 2, 3, 4 )
    .collect()
    .view()

ch = channel
     .of( ['wt','wt_1.fq'], ['wt','wt_2.fq'], ["mut",'mut_1.fq'], ['mut', 'mut_2.fq'] )
     .groupTuple()
     .view()

ch1 = channel.of( 1,2,3 )
ch2 = channel.of( 'X','Y' )
ch3 = channel.of( 'mt' )
ch4 = ch1.mix(ch2,ch3).view()

reads1_ch = channel
  .of(['wt', 'wt_1.fq'], ['mut','mut_1.fq'])
reads2_ch= channel
  .of(['wt', 'wt_2.fq'], ['mut','mut_2.fq'])
reads_ch = reads1_ch
  .join(reads2_ch)
  .view()

channel
     .of( 'chr1', 'chr2', 'chr3' )
     .into({ ch1; ch2 })

ch1.view({"ch1 emits: $it"})
ch2.view({"ch2 emits: $it"})

ch = channel
    .of(1..22,'X','Y')
    .count()
    .view()

/* Channel.of("val1\tval2\tval3\nval4\tval5\tval6\n")
  .splitCsv(sep: "\t")
  .view()   */
  
   