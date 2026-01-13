process fastplong {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastplong.{json,html,csv}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_RL.trim.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastplong.json"), emit: json
    tuple val(sample_id), path("${sample_id}_fastplong.html"), emit: html
    tuple val(sample_id), path("${sample_id}_fastplong.csv"), emit: csv
    tuple val(sample_id), path("${sample_id}_fastplong_provenance.yml"), emit: provenance
    

    script:
    """
    printf -- "- process_name: fastplong\\n"                                          >> ${sample_id}_fastplong_provenance.yml
    printf -- "  tools:\\n"                                                           >> ${sample_id}_fastplong_provenance.yml
    printf -- "    - tool_name: fastplong\\n"                                         >> ${sample_id}_fastplong_provenance.yml
    printf -- "      tool_version: \$(fastplong --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastplong_provenance.yml

    fastplong \
	-t ${task.cpus} \
	-i ${reads} \
	-o ${sample_id}_RL.trim.fastq.gz \
	-h ${sample_id}_fastplong.html \
	-j ${sample_id}_fastplong.json

    fastplong_json_to_csv.py ${sample_id}_fastplong.json -s ${sample_id} > ${sample_id}_fastplong.csv
    """
}

