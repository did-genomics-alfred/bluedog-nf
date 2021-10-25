#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Processes to run on reads
include { get_read_prefix_and_type} from './src/channel_helpers.nf'
include { fastqc;multiqc } from './src/processes/read_processes.nf'
include { assemble;assembly_stats;combine_stats } from './src/processes/read_processes.nf'

// Processes to run on assemblies
include { mlst;combine_mlst } from './src/processes/assembly_processes.nf'
include {speciator;combine_speciator} from './src/processes/assembly_processes.nf'
include { kleborate;combine_kleborate } from './src/processes/assembly_processes.nf'

// Utility functions
//include { print_splash } from './src/utilities.nf'
include { check_host;check_arguments;check_boolean_option } from './src/utilities.nf'

// Parameter checks
// need to provide either reads OR assemblies, can't have both or neither
check_arguments(params)
check_host(workflow)

// Require some variables to be boolean
// We must check and change values if needed
run_read_qc = check_boolean_option(params.read_qc, 'read_qc')
run_speciator = check_boolean_option(params.run_speciator, 'run_speciator')
run_mlst = check_boolean_option(params.run_mlst, 'run_mlst')
run_kleborate = check_boolean_option(params.run_kleborate, 'run_kleborate')

// Get reads and seperate into pe and se channels based on prefix
if (params.reads) {
  reads = Channel.fromPath(params.reads).ifEmpty {
      exit 1, "ERROR: did not find any read files with '${params.reads}'"
    }.map {
      get_read_prefix_and_type(it)
    }.branch {
      paired: it[0] == 'pe'
  }
  reads_pe = reads.paired.map { it[1..-1] }.groupTuple()
  // Check that we have the expected number of reads for each prefix in channel and flatten tuple
  reads_pe = reads_pe.map {
    if (it[1].size() != 2) {
      exit 1, "ERROR: didn't get exactly two readsets prefixed with ${it[0]}:\n${it[1]}"
    }
    [it[0], *it[1]]
  }
}

// Get assembly files and create tuple with filename + full file path
if (params.assemblies){
	assemblies = Channel.fromPath(params.assemblies).ifEmpty {
    exit 1, "ERROR: did not find any assembly files with '${params.assemblies}'"
  }.map {
    file -> tuple(file.simpleName, file) }
}

// This workflow is for the assembly steps (includes optional FASTQC on read files)
workflow ASSEMBLE_FROM_READS {

  take:
  read_pairs_ch

  main:

  if( run_read_qc ) {
    fastqc_ch = fastqc(reads_pe)
    multiqc(fastqc_ch.collect())
  }

  //run the assembly process
  assemble(reads_pe)
  stats_ch = assembly_stats(assemble.out.assembly_fasta)
  combine_stats(stats_ch.collect())

  emit:
  assemble.out.assembly_fasta

}

// This is the workflow for analysing assembly files, all steps optional
workflow ANALYSE_ASSEMBLIES {

  take:
  assemblies_ch

  main:
  if ( run_speciator ) {
    species_ch = speciator(assemblies_ch)
    combine_speciator(species_ch.collect())
  }
  if ( run_kleborate ) {
    kleborate_ch = kleborate(assemblies_ch)
    combine_kleborate(kleborate_ch.collect())
  }

  if ( run_mlst ) {
    mlst_ch = mlst(assemblies_ch)
    combine_mlst(mlst_ch.collect())
  }

}

// This is the implicit workflow that calls the others, depending on input
workflow {
  if (params.reads) {
    ASSEMBLE_FROM_READS(reads_pe)
    ANALYSE_ASSEMBLIES(ASSEMBLE_FROM_READS.out)
  }
  if (params.assemblies) {
    ANALYSE_ASSEMBLIES(assemblies)
    //reads_pe.view()
    //assemblies.view()
  }
}
