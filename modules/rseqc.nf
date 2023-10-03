process RSEQC_BAMSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam_stat.txt"), emit: txt
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_stat.py \\
        -i $bam \\
        $args \\
        > ${prefix}.bam_stat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(bam_stat.py --version | sed -e "s/bam_stat.py //g")
    END_VERSIONS
    """
}
process RSEQC_INFEREXPERIMENT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.infer_experiment.txt"), emit: txt
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    infer_experiment.py \\
        -i $bam \\
        -r $bed \\
        $args \\
        > ${prefix}.infer_experiment.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(infer_experiment.py --version | sed -e "s/infer_experiment.py //g")
    END_VERSIONS
    """
}
process RSEQC_INNERDISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*distance.txt"), optional:true, emit: distance
    tuple val(meta), path("*freq.txt")    , optional:true, emit: freq
    tuple val(meta), path("*mean.txt")    , optional:true, emit: mean
    tuple val(meta), path("*.pdf")        , optional:true, emit: pdf
    tuple val(meta), path("*.r")          , optional:true, emit: rscript
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (!meta.single_end) {
        """
        inner_distance.py \\
            -i $bam \\
            -r $bed \\
            -o $prefix \\
            $args \\
            > stdout.txt
        head -n 2 stdout.txt > ${prefix}.inner_distance_mean.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    } else {
        """
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
        END_VERSIONS
        """
    }
}
process RSEQC_JUNCTIONANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.xls")         , emit: xls
    tuple val(meta), path("*.r")           , emit: rscript
    tuple val(meta), path("*.log")         , emit: log
    tuple val(meta), path("*.junction.bed"), optional:true, emit: bed
    tuple val(meta), path("*.Interact.bed"), optional:true, emit: interact_bed
    tuple val(meta), path("*junction.pdf") , optional:true, emit: pdf
    tuple val(meta), path("*events.pdf")   , optional:true, emit: events_pdf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_annotation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args \\
        2> ${prefix}.junction_annotation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_annotation.py --version | sed -e "s/junction_annotation.py //g")
    END_VERSIONS
    """
}
process RSEQC_JUNCTIONSATURATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    tuple val(meta), path("*.r")  , emit: rscript
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    junction_saturation.py \\
        -i $bam \\
        -r $bed \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(junction_saturation.py --version | sed -e "s/junction_saturation.py //g")
    END_VERSIONS
    """
}
process RSEQC_READDISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.read_distribution.txt"), emit: txt
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_distribution.py \\
        -i $bam \\
        -r $bed \\
        > ${prefix}.read_distribution.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_distribution.py --version | sed -e "s/read_distribution.py //g")
    END_VERSIONS
    """
}
process RSEQC_READDUPLICATION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*seq.DupRate.xls"), emit: seq_xls
    tuple val(meta), path("*pos.DupRate.xls"), emit: pos_xls
    tuple val(meta), path("*.pdf")           , emit: pdf
    tuple val(meta), path("*.r")             , emit: rscript
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    read_duplication.py \\
        -i $bam \\
        -o $prefix \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(read_duplication.py --version | sed -e "s/read_duplication.py //g")
    END_VERSIONS
    """
}
process RSEQC_TIN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::rseqc=3.0.1 'conda-forge::r-base>=3.5'"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:3.0.1--py37h516909a_1' :
        'biocontainers/rseqc:3.0.1--py37h516909a_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  bed

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.xls"), emit: xls
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tin.py \\
        -i $bam \\
        -r $bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(tin.py --version | sed -e "s/tin.py //g")
    END_VERSIONS
    """
}

//
// Run RSeQC modules
//
/*
include { RSEQC_BAMSTAT            } from '../../../modules/rseqc/bamstat/main'
include { RSEQC_INNERDISTANCE      } from '../../../modules/rseqc/innerdistance/main'
include { RSEQC_INFEREXPERIMENT    } from '../../../modules/rseqc/inferexperiment/main'
include { RSEQC_JUNCTIONANNOTATION } from '../../../modules/rseqc/junctionannotation/main'
include { RSEQC_JUNCTIONSATURATION } from '../../../modules/rseqc/junctionsaturation/main'
include { RSEQC_READDISTRIBUTION   } from '../../../modules/rseqc/readdistribution/main'
include { RSEQC_READDUPLICATION    } from '../../../modules/rseqc/readduplication/main'
include { RSEQC_TIN                } from '../../../modules/rseqc/tin/main'
*/

workflow BAM_RSEQC {
    take:
    ch_bam_bai    // channel: [ val(meta), [ bam, bai ] ]
    ch_bed        //    file: /path/to/genome.bed
    rseqc_modules //    list: rseqc modules to run

    main:

    ch_versions = Channel.empty()

    ch_bam_bai
        .map { [ it[0], it[1] ] }
        .set { ch_bam }

    //
    // Run RSeQC bam_stat.py
    //
    bamstat_txt = Channel.empty()
    if ('bam_stat' in rseqc_modules) {
        RSEQC_BAMSTAT ( ch_bam )
        bamstat_txt = RSEQC_BAMSTAT.out.txt
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions.first())
    }

    //
    // Run RSeQC inner_distance.py
    //
    innerdistance_distance = Channel.empty()
    innerdistance_freq     = Channel.empty()
    innerdistance_mean     = Channel.empty()
    innerdistance_pdf      = Channel.empty()
    innerdistance_rscript  = Channel.empty()
    if ('inner_distance' in rseqc_modules) {
        RSEQC_INNERDISTANCE ( ch_bam, ch_bed )
        innerdistance_distance = RSEQC_INNERDISTANCE.out.distance
        innerdistance_freq     = RSEQC_INNERDISTANCE.out.freq
        innerdistance_mean     = RSEQC_INNERDISTANCE.out.mean
        innerdistance_pdf      = RSEQC_INNERDISTANCE.out.pdf
        innerdistance_rscript  = RSEQC_INNERDISTANCE.out.rscript
        ch_versions = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions.first())
    }

    //
    // Run RSeQC infer_experiment.py
    //
    inferexperiment_txt = Channel.empty()
    if ('infer_experiment' in rseqc_modules) {
        RSEQC_INFEREXPERIMENT ( ch_bam, ch_bed )
        inferexperiment_txt = RSEQC_INFEREXPERIMENT.out.txt
        ch_versions = ch_versions.mix(RSEQC_INFEREXPERIMENT.out.versions.first())
    }

    //
    // Run RSeQC junction_annotation.py
    //
    junctionannotation_bed          = Channel.empty()
    junctionannotation_interact_bed = Channel.empty()
    junctionannotation_xls          = Channel.empty()
    junctionannotation_pdf          = Channel.empty()
    junctionannotation_events_pdf   = Channel.empty()
    junctionannotation_rscript      = Channel.empty()
    junctionannotation_log          = Channel.empty()
    if ('junction_annotation' in rseqc_modules) {
        RSEQC_JUNCTIONANNOTATION ( ch_bam, ch_bed )
        junctionannotation_bed          = RSEQC_JUNCTIONANNOTATION.out.bed
        junctionannotation_interact_bed = RSEQC_JUNCTIONANNOTATION.out.interact_bed
        junctionannotation_xls          = RSEQC_JUNCTIONANNOTATION.out.xls
        junctionannotation_pdf          = RSEQC_JUNCTIONANNOTATION.out.pdf
        junctionannotation_events_pdf   = RSEQC_JUNCTIONANNOTATION.out.events_pdf
        junctionannotation_rscript      = RSEQC_JUNCTIONANNOTATION.out.rscript
        junctionannotation_log          = RSEQC_JUNCTIONANNOTATION.out.log
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONANNOTATION.out.versions.first())
    }

    //
    // Run RSeQC junction_saturation.py
    //
    junctionsaturation_pdf     = Channel.empty()
    junctionsaturation_rscript = Channel.empty()
    if ('junction_saturation' in rseqc_modules) {
        RSEQC_JUNCTIONSATURATION ( ch_bam, ch_bed )
        junctionsaturation_pdf     = RSEQC_JUNCTIONSATURATION.out.pdf
        junctionsaturation_rscript = RSEQC_JUNCTIONSATURATION.out.rscript
        ch_versions = ch_versions.mix(RSEQC_JUNCTIONSATURATION.out.versions.first())
    }

    //
    // Run RSeQC read_distribution.py
    //
    readdistribution_txt = Channel.empty()
    if ('read_distribution' in rseqc_modules) {
        RSEQC_READDISTRIBUTION ( ch_bam, ch_bed )
        readdistribution_txt = RSEQC_READDISTRIBUTION.out.txt
        ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions.first())
    }

    //
    // Run RSeQC read_duplication.py
    //
    readduplication_seq_xls = Channel.empty()
    readduplication_pos_xls = Channel.empty()
    readduplication_pdf     = Channel.empty()
    readduplication_rscript = Channel.empty()
    if ('read_duplication' in rseqc_modules) {
        RSEQC_READDUPLICATION ( ch_bam )
        readduplication_seq_xls = RSEQC_READDUPLICATION.out.seq_xls
        readduplication_pos_xls = RSEQC_READDUPLICATION.out.pos_xls
        readduplication_pdf     = RSEQC_READDUPLICATION.out.pdf
        readduplication_rscript = RSEQC_READDUPLICATION.out.rscript
        ch_versions = ch_versions.mix(RSEQC_READDUPLICATION.out.versions.first())
    }

    //
    // Run RSeQC tin.py
    //
    tin_txt = Channel.empty()
    if ('tin' in rseqc_modules) {
        RSEQC_TIN ( ch_bam_bai, ch_bed )
        tin_txt     = RSEQC_TIN.out.txt
        ch_versions = ch_versions.mix(RSEQC_TIN.out.versions.first())
    }

    emit:
    bamstat_txt                     // channel: [ val(meta), txt ]

    innerdistance_distance          // channel: [ val(meta), txt ]
    innerdistance_freq              // channel: [ val(meta), txt ]
    innerdistance_mean              // channel: [ val(meta), txt ]
    innerdistance_pdf               // channel: [ val(meta), pdf ]
    innerdistance_rscript           // channel: [ val(meta), r   ]

    inferexperiment_txt             // channel: [ val(meta), txt ]

    junctionannotation_bed          // channel: [ val(meta), bed ]
    junctionannotation_interact_bed // channel: [ val(meta), bed ]
    junctionannotation_xls          // channel: [ val(meta), xls ]
    junctionannotation_pdf          // channel: [ val(meta), pdf ]
    junctionannotation_events_pdf   // channel: [ val(meta), pdf ]
    junctionannotation_rscript      // channel: [ val(meta), r   ]
    junctionannotation_log          // channel: [ val(meta), log ]

    junctionsaturation_pdf          // channel: [ val(meta), pdf ]
    junctionsaturation_rscript      // channel: [ val(meta), r   ]

    readdistribution_txt            // channel: [ val(meta), txt ]

    readduplication_seq_xls         // channel: [ val(meta), xls ]
    readduplication_pos_xls         // channel: [ val(meta), xls ]
    readduplication_pdf             // channel: [ val(meta), pdf ]
    readduplication_rscript         // channel: [ val(meta), r   ]

    tin_txt                         // channel: [ val(meta), txt ]

    versions = ch_versions          // channel: [ versions.yml ]
}

