#!/usr/bin/env nextflow

import java.time.LocalDateTime

nextflow.enable.dsl = 2

include { hash_files }                      from './modules/hash_files.nf'
include { fastplong }                       from './modules/fastplong.nf'
// include { fastp_json_to_csv }            from './modules/fastplong.nf'
include { autocycler_estimate_genome_size } from './modules/autocycler.nf'
include { autocycler_subsample }            from './modules/autocycler.nf'
include { autocycler_assemble }             from './modules/autocycler.nf'
include { autocycler_compress }             from './modules/autocycler.nf'
include { autocycler_cluster }              from './modules/autocycler.nf'
include { autocycler_trim_resolve }         from './modules/autocycler.nf'
include { autocycler_combine }              from './modules/autocycler.nf'
include { reorient_contigs }                from './modules/autocycler.nf'
include { combine_metrics }                 from './modules/autocycler.nf'
include { quast }                           from './modules/assembly_qc.nf'
include { bandage }                         from './modules/assembly_qc.nf'
// include { prokka }                       from './modules/prokka.nf'
// include { bakta }                        from './modules/bakta.nf'
// include { bandage }                      from './modules/long_read_qc.nf'
include { pipeline_provenance }             from './modules/provenance.nf'
include { collect_provenance }              from './modules/provenance.nf'


workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    valid_assemblers = [
        "raven",
	"myloasm",
	"miniasm",
	"flye",
	"metamdbg",
	"necat",
	"nextdenovo",
	"plassembler",
	"canu"
    ].toSet()

    assemblers_list = file(params.assemblers_list).readLines().unique()
    num_assemblers = assemblers_list.size()
    ch_autocycler_assemblers = Channel.fromList(assemblers_list)
    
    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    if (params.samplesheet_input != 'NO_FILE') {
        ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], [it['LONG']]] }
    } else {
        ch_fastq = Channel.fromPath( params.long_reads_search_path ).map{ it -> [it.baseName.split("_")[0], [it]] }
    }

    main:
    ch_sample_ids = ch_fastq.map{ it -> it[0] }
    ch_provenance = ch_fastq.map{ it -> it[0] }

    hash_files(ch_fastq.combine(Channel.of("fastq-input")))

    fastplong(ch_fastq)

    autocycler_estimate_genome_size(ch_fastq)
    ch_genome_size = autocycler_estimate_genome_size.out.genome_size

    autocycler_subsample(ch_fastq.join(ch_genome_size))
    ch_subsampled_reads = autocycler_subsample.out.subsampled_reads

    autocycler_assemble(ch_subsampled_reads.join(ch_genome_size).combine(ch_autocycler_assemblers))
    ch_assembly_dirs = autocycler_assemble.out.assembly_dir.groupTuple(size: num_assemblers)

    autocycler_compress(ch_assembly_dirs)
    ch_autocycler_compress_out = autocycler_compress.out.autocycler_out

    autocycler_cluster(ch_autocycler_compress_out)
    ch_autocycler_cluster_out = autocycler_cluster.out.autocycler_out

    autocycler_trim_resolve(ch_autocycler_cluster_out)
    ch_autocycler_trim_resolve_out = autocycler_trim_resolve.out.autocycler_out

    autocycler_combine(ch_autocycler_trim_resolve_out)
    ch_assembly = autocycler_combine.out.consensus_assembly
    ch_assembly_graph = autocycler_combine.out.consensus_assembly_graph

    reorient_contigs(ch_assembly_graph)
    ch_reoriented_assembly = reorient_contigs.out.assembly
    ch_reoriented_assembly_graph = reorient_contigs.out.assembly_graph

    ch_read_metrics = autocycler_subsample.out.read_metrics
    ch_assembly_metrics = autocycler_combine.out.assembly_metrics
    combine_metrics(ch_read_metrics.join(ch_assembly_metrics))

    quast(ch_reoriented_assembly)

    bandage(ch_reoriented_assembly_graph)

    if (params.prokka) {
	// prokka(ch_assembly)
    }

    if (params.bakta) {
	// bakta(ch_assembly)
    }

    


    //
    // Provenance collection processes
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it -> [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_files.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    

    ch_provenance = ch_provenance.join(fastplong.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    if (params.prokka) {
	// ch_provenance = ch_provenance.join(prokka.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    if (params.bakta) {
	// ch_provenance = ch_provenance.join(bakta.out.provenance).map{ it -> [it[0], it[1] << it[2]] }
    }

    collect_provenance(ch_provenance)
}
