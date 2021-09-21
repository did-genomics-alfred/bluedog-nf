#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { get_read_prefix_and_type } from './src/channel_helpers.nf'

workflow {
  //get input data
  read_pairs = Channel.fromFilePairs(params.reads)

  // Get reads and seperate into pe and se channels based on prefix
  reads = Channel.fromPath(params.reads).ifEmpty {
      exit 1, "ERROR: did not find any read files with '${params.reads}'"
    }.map {
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

  //run the assembly process
  ASSEMBLE(reads_pe)

  //view the output
  //ASSEMBLE.out.view()
}

process ASSEMBLE {
  cache 'lenient'
  publishDir path: {"${params.output_dir}/fasta"},  mode: 'copy', saveAs: {filename -> "${isolate_id}.fasta"}, pattern: '*.fasta'
  publishDir path: {"${params.output_dir}/graph"}, mode: 'copy', saveAs: {filename -> "${isolate_id}.gfa"}, pattern: '*.gfa'
  publishDir path: {"${params.output_dir}/logs"}, mode: 'copy', saveAs: {filename -> "${isolate_id}_unicycler.log"}, pattern: '*.log'

  input:
  tuple val(isolate_id), path(reads_fwd), path(reads_rev)

  output:
  tuple val(isolate_id), path("assembly*"), path("unicycler.log")

  script:
  """
  unicycler -1 ${reads_fwd} -2 ${reads_rev} -o .
  """

}
