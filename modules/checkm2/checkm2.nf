process checkm2 {

    tag "${sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0a/0af812c983aeffc99c0fca9ed2c910816b2ddb9a9d0dcad7b87dab0c9c08a16f/data':
        'community.wave.seqera.io/library/checkm2:1.1.0--60f287bc25d7a10d' }"

    publishDir "${params.outdir}/${prefix}", pattern: "${prefix}_checkm2.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(assembly)

    output:
    tuple val(sample_id), path("${prefix}_checkm2.csv")                   , emit: checkm2_output
    tuple val(sample_id), path("${prefix}_checkm2_provenance.yml")        , emit: provenance
    
    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${sample_id}"
    """
    printf -- "- process_name: checkm2\\n"                               >> ${prefix}_checkm2_provenance.yml
    printf -- "  tools:\\n"                                              >> ${prefix}_checkm2_provenance.yml
    printf -- "    - tool_name: checkm2\\n"                              >> ${prefix}_checkm2_provenance.yml
    printf -- "      tool_version: \$(checkm2 --version | tr -d 'v')\\n" >> ${prefix}_checkm2_provenance.yml
    printf -- "      subcommand: predict\\n"                             >> ${prefix}_checkm2_provenance.yml
    printf -- "      parameters:\\n"                                     >> ${prefix}_checkm2_provenance.yml
    printf -- "        - parameter: --database_path\\n"                  >> ${prefix}_checkm2_provenance.yml
    printf -- "          value: ${params.checkm2_db}\\n"                 >> ${prefix}_checkm2_provenance.yml

    mkdir -p tmp

    checkm2 predict \
        --threads ${task.cpus} \
        --database_path ${params.checkm2_db} \
        --input ${assembly} \
        ${args} \
        --output-directory ${prefix}_checkm2_output
    
    cat ${prefix}_checkm2_output/quality_report.tsv | tr \$'\\t' ',' > ${prefix}_checkm2.csv
    """
}