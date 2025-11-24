process autocycler_estimate_genome_size {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_autocycler_genome_size_estimate.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_autocycler_genome_size_estimate.txt"),             emit: genome_size
    tuple val(sample_id), path("${sample_id}_autocycler_estimate_genome_size_provenance.yml"),  emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_estimate_genome_size\\n"                >> ${sample_id}_autocycler_estimate_genome_size_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_estimate_genome_size_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_estimate_genome_size_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_estimate_genome_size_provenance.yml
    printf -- "      subcommand: helper genome_size\\n"                           >> ${sample_id}_autocycler_estimate_genome_size_provenance.yml

    autocycler helper genome_size \
      --threads ${task.cpus} \
      --reads ${reads[0]} \
      > ${sample_id}_autocycler_genome_size_estimate.txt
    """
}

process autocycler_subsample {

    tag { sample_id }

    input:
    tuple val(sample_id), path(reads), path(genome_size_estimate)

    output:
    tuple val(sample_id), path("${sample_id}_subsampled_reads"),                    emit: subsampled_reads
    tuple val(sample_id), path("${sample_id}_autocycler_subsample.log"),            emit: log
    tuple val(sample_id), path("${sample_id}_autocycler_subsample_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_subsample\\n"                           >> ${sample_id}_autocycler_subsample_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_subsample_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_subsample_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_subsample_provenance.yml
    printf -- "      subcommand: subsample\\n"                                    >> ${sample_id}_autocycler_subsample_provenance.yml

    autocycler subsample \
      --reads ${reads[0]} \
      --out_dir ${sample_id}_subsampled_reads \
      --genome_size \$(cat ${genome_size_estimate}) \
      2>> ${sample_id}_autocycler_subsample.log
    """
}


process autocycler_assemble {

    tag { sample_id + ' / ' + assembler }

    input:
    tuple val(sample_id), path(subsampled_reads), path(genome_size_estimate), val(assembler)

    output:
    tuple val(sample_id), path("${sample_id}_assemblies_${assembler}"),                         emit: assembly_dir
    tuple val(sample_id), path("${sample_id}_autocycler_assemble_${assembler}_provenance.yml"), emit: provenance

    script:
    threads_per_job = task.cpus / 4
    """
    printf -- "- process_name: autocycler_assemble\\n"                            >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "      subcommand: helper ${assembler}\\n"                          >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "      parameters:\\n"                                              >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "        - parameter: min_depth_rel                                 >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml
    printf -- "          value: 0.1                                               >> ${sample_id}_autocycler_assemble_${assembler}_provenance.yml

    mkdir -p ${sample_id}_assemblies_${assembler}
    mkdir -p tmp_01 tmp_02 tmp_03 tmp_04

    for reads_subsample in 01 02 03 04; do
        echo "autocycler helper ${assembler} --threads ${threads_per_job} --dir tmp_\${reads_subsample} --genome_size \$(cat ${genome_size_estimate}) --read_type ${params.read_type} --reads ${sample_id}_subsampled_reads/sample_\${reads_subsample}.fastq --out_prefix ${sample_id}_assemblies_${assembler}/${sample_id}_${assembler}_\${reads_subsample} --min_depth_rel 0.1" >> ${sample_id}_assemblies_${assembler}/jobs.txt
    done

    parallel --jobs 4 \
        --joblog ${sample_id}_assemblies_${assembler}/joblog.tsv \
	--results ${sample_id}_assemblies_${assembler}/logs \
	--timeout "8h" \
	< ${sample_id}_assemblies_${assembler}/jobs.txt


    shopt -s nullglob
    # Give circular contigs from Plassembler extra clustering weight
    if [ ${assembler} = "plassembler" ]; then
        for f in ${sample_id}_assemblies_${assembler}/plassembler*.fasta; do
            sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "\$f"
        done
    fi

    # Give circular contigs from canu extra clustering weight
    if [ ${assembler} = "canu" ]; then
        for f in ${sample_id}_assemblies_${assembler}/canu*.fasta; do
            sed -i 's/^>.*\$/& Autocycler_consensus_weight=2/' "\$f"
        done
    fi

    # Give circular contigs from flye extra clustering weight
    if [ ${assembler} = "flye" ]; then
        for f in ${sample_id}_assemblies_${assembler}/flye*.fasta; do
            sed -i 's/^>.*\$/& Autocycler_consensus_weight=2/' "\$f"
        done
    fi
    shopt -u nullglob

    
    """
}


process autocycler_compress {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_input_assemblies",            mode: 'copy'

    input:
    tuple val(sample_id), path(assembly_dirs)   

    output:
    tuple val(sample_id), path("${sample_id}_autocycler_out"), path("${sample_id}_autocycler.stderr"),  emit: autocycler_out
    tuple val(sample_id), path("${sample_id}_input_assemblies"),                                        emit: input_assemblies
    tuple val(sample_id), path("${sample_id}_autocycler_compress_provenance.yml"),                      emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_compress\\n"                            >> ${sample_id}_autocycler_compress_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_compress_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_compress_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_compress_provenance.yml
    printf -- "      subcommand: compress\\n"                                     >> ${sample_id}_autocycler_compress_provenance.yml

    mkdir ${sample_id}_input_assemblies

    cd ${sample_id}_input_assemblies

    ln -s ../${sample_id}_assemblies_*/*.fasta .

    cd ..

    autocycler compress \
        --threads ${task.cpus} \
        --assemblies_dir ${sample_id}_input_assemblies \
	--autocycler_dir ${sample_id}_autocycler_out \
	2>> ${sample_id}_autocycler.stderr
    """
}


process autocycler_cluster {

    tag { sample_id }

    input:
    tuple val(sample_id), path(autocycler_out), path(autocycler_stderr)

    output:
    tuple val(sample_id), path("${sample_id}_autocycler_out"), path("${sample_id}_autocycler.stderr"),  emit: autocycler_out
    tuple val(sample_id), path("${sample_id}_autocycler_cluster_provenance.yml"),                       emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_cluster\\n"                             >> ${sample_id}_autocycler_cluster_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_cluster_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_cluster_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_cluster_provenance.yml
    printf -- "      subcommand: cluster\\n"                                      >> ${sample_id}_autocycler_cluster_provenance.yml

    autocycler cluster \
        --autocycler_dir ${sample_id}_autocycler_out \
	2>> ${sample_id}_autocycler.stderr
    """
}


process autocycler_trim_resolve {

    tag { sample_id }

    input:
    tuple val(sample_id), path(autocycler_out),  path(autocycler_stderr)

    output:
    tuple val(sample_id), path("${sample_id}_autocycler_out"), path("${sample_id}_autocycler.stderr"),  emit: autocycler_out
    tuple val(sample_id), path("${sample_id}_autocycler_trim_resolve_provenance.yml"),                  emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_trim_resolve\\n"                        >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "      subcommand: trim\\n"                                         >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_trim_resolve_provenance.yml
    printf -- "      subcommand: resolve\\n"                                      >> ${sample_id}_autocycler_trim_resolve_provenance.yml


    for c in ${sample_id}_autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c" 2>> ${sample_id}_autocycler.stderr
        autocycler resolve -c "\$c" 2>> ${sample_id}_autocycler.stderr
    done
    """
}


process autocycler_combine {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_autocycler_out",            mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_autocycler_long.{fa,gfa}",  mode: 'copy'
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_autocycler.stderr",         mode: 'copy'

    input:
    tuple val(sample_id), path(autocycler_out), path(autocycler_stderr)

    output:
    tuple val(sample_id), path("${sample_id}_autocycler_out"), path("${sample_id}_autocycler.stderr"),  emit: autocycler_out
    tuple val(sample_id), path("${sample_id}_autocycler_combine_provenance.yml"),                       emit: provenance

    script:
    """
    printf -- "- process_name: autocycler_combine\\n"                             >> ${sample_id}_autocycler_combine_provenance.yml
    printf -- "  tools:\\n"                                                       >> ${sample_id}_autocycler_combine_provenance.yml
    printf -- "    - tool_name: autocycler\\n"                                    >> ${sample_id}_autocycler_combine_provenance.yml
    printf -- "      tool_version: \$(autocycler --version | cut -d ' ' -f 2)\\n" >> ${sample_id}_autocycler_combine_provenance.yml
    printf -- "      subcommand: combine\\n"                                      >> ${sample_id}_autocycler_combine_provenance.yml

    autocycler combine \
        -a ${sample_id}_autocycler_out \
	-i ${sample_id}_autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa \
	2>> ${sample_id}_autocycler.stderr


    cp ${sample_id}_autocycler_out/consensus_assembly.fasta ${sample_id}_autocycler_long.fa
    cp ${sample_id}_autocycler_out/consensus_assembly.gfa ${sample_id}_autocycler_long.gfa
    """
}
