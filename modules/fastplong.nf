process fastplong {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.json", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_RL.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastplong.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastplong.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastplong_provenance.yml"), emit: provenance
    

    script:
    """
    printf -- "- process_name: fastplong\\n"                                          >> ${sample_id}_fastplong_provenance.yml
    printf -- "  tools:\\n"                                                           >> ${sample_id}_fastplong_provenance.yml
    printf -- "    - tool_name: fastplong\\n"                                         >> ${sample_id}_fastplong_provenance.yml
    printf -- "      tool_version: \$(fastplong --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastplong_provenance.yml
    printf -- "      parameters:\\n"                                                  >> ${sample_id}_fastplong_provenance.yml
    printf -- "        - parameter: --cut_tail\\n"                                    >> ${sample_id}_fastplong_provenance.yml
    printf -- "          value: null\\n"                                              >> ${sample_id}_fastplong_provenance.yml

    fastplong \
	-t ${task.cpus} \
	-i ${reads} \
	-o ${sample_id}_RL.trim.fastq.gz \
	-h ${sample_id}_fastplong.html \
	-j ${sample_id}_fastplong.json
    """
}

process fastp_json_to_csv {

    tag { sample_id }

    executor 'local'

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(fastp_json)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.csv")

    script:
    """
    fastp_json_to_csv.py -s ${sample_id} ${fastp_json} > ${sample_id}_fastp.csv
    """
}
