#!/usr/bin/env nextflow

//ADD ABRITAMR, PROKKA ALTERNATIVE, POSSIBLY SPECIES DETECTOR
//WORK OUT A WAY TO PROVIDE ASSEMBLIES AS INPUT ALTERNATIVE

nextflow.enable.dsl=2

include { get_read_prefix_and_type;get_file_id } from './src/channel_helpers.nf'
include { fastqc;multiqc } from './src/processes/read_processes.nf'
include { assemble;assembly_stats;combine_stats } from './src/processes/read_processes.nf'

include { mlst;combine_mlst } from './src/processes/assembly_processes.nf'
include {speciator} from './src/processes/assembly_processes.nf'

// Get reads and seperate into pe and se channels based on prefix
reads = Channel.fromPath(params.reads).map {
    get_read_prefix_and_type(it)
  }.branch {
    paired: it[0] == 'pe'
}
reads_pe = reads.paired.map { it[1..-1] }.groupTuple()
// Check that we have the expected number of reads for each prefix in pe and se channels and flatten tuple
reads_pe = reads_pe.map {
  if (it[1].size() != 2) {
    exit 1, "ERROR: didn't get exactly two readsets prefixed with ${it[0]}:\n${it[1]}"
  }
  [it[0], *it[1]]
}

assemblies = Channel.fromPath(params.assemblies).map {
  get_file_id(it)}

workflow ASSEMBLE_FROM_READS {
  //get input data
  take:
  read_pairs_ch
  main:
  //fastqc
  if( params.read_qc ) {
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

workflow ANALYSE_ASSEMBLIES {
  take:
  assemblies_ch
  main:
  if (params.run_speciator) {
    species_ch = speciator(assemblies_ch)
    combine_species(species_ch.collect())
  }
  if ( params.run_kleborate ) {
    kleborate_ch = kleborate(assemblies_ch)
    combine_kleborate(kleborate_ch.collect())
  }

  if ( params.run_mlst ) {
    mlst_ch = mlst(assemblies_ch)
    combine_mlst(mlst_ch)
  }

}

//This is the implicit workflow that calls all the others
workflow{
  if (params.reads) {
    ASSEMBLE_FROM_READS(reads_pe)
    ANALYSE_ASSEMBLIES(ASSEMBLE_FROM_READS.out)
  }
  else {
    ANALYSE_ASSEMBLIES(assemblies)
    //reads_pe.view()
    //assemblies.view()
  }
}
