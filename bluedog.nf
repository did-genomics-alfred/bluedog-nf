#!/usr/bin/env nextflow

//ADD ABRITAMR, PROKKA ALTERNATIVE, POSSIBLY SPECIES DETECTOR
//WORK OUT A WAY TO PROVIDE ASSEMBLIES AS INPUT ALTERNATIVE

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

  //fastqc
  if( params.read_qc ) {
    fastqc_ch = FASTQC(reads_pe)
    MULTIQC(fastqc_ch.collect())
  }

  //run the assembly process
  assembly_ch = ASSEMBLE(reads_pe)
  stats_ch = STATS(assembly_ch)
  COMBINE(stats_ch.collect())

  if (params.run_kleborate) {
    kleborate_ch = KLEBORATE(assembly_ch)
    COMBINE_KLEB(kleborate_ch.collect())
  }

}

process FASTQC {
  label "short_job"
  cache 'lenient'

  input:
  tuple val(isolate_id), path(reads_fwd), path(reads_rev)

  output:
  path("*.zip")

  script:
  """
  fastqc --noextract $reads_fwd $reads_rev
  """
}

process MULTIQC {
  label "short_job"
  cache 'lenient'
  publishDir path: {"${params.output_dir}"}, mode: 'copy', saveAs: {filename -> "multiQC_report.txt"}

  input:
  path("*")

  output:
  path("multiqc_data/multiqc_fastqc.txt")

  script:
  """
  multiqc --force .
  """
}

process ASSEMBLE {
  cache 'lenient'
  publishDir path: {"${params.output_dir}/fasta"},  mode: 'copy', saveAs: {filename -> "${isolate_id}.fasta"}, pattern: '*.fasta'
  publishDir path: {"${params.output_dir}/graph"}, mode: 'copy', saveAs: {filename -> "${isolate_id}.gfa"}, pattern: 'assembly.gfa'
  publishDir path: {"${params.output_dir}/logs"}, mode: 'copy', saveAs: {filename -> "${isolate_id}_unicycler.log"}, pattern: '*.log'

  input:
  tuple val(isolate_id), path(reads_fwd), path(reads_rev)

  output:
  tuple val(isolate_id), path("assembly.fasta"), path("assembly.gfa"), path("unicycler.log")

  script:
  """
  unicycler -1 ${reads_fwd} -2 ${reads_rev} -o .
  """

}

process STATS {
  label "short_job"
  cache 'lenient'
  publishDir path:{"${params.output_dir}/stats"}, mode: 'copy', saveAs: {filename -> "${isolate_id}_stats.txt"}, pattern: '*.txt'

  input:
  tuple val(isolate_id), path(fasta_file), path(graph_file), path(log_file)

  output:
  path("${isolate_id}_stats.txt")

  script:
  """
  assembly_stats.py -a $fasta_file --id $isolate_id > ${isolate_id}_stats.txt
  """
}

process COMBINE {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}"), mode: 'copy'

  input:
  file("*_stats.txt")

  output:
  file('assembly_stats.txt')

  script:
  """
  awk 'FNR==1 && NR!=1 { while (/^assembly/) getline; } 1 {print}' *_stats.txt > assembly_stats.txt
  """
}

/* Need to comment out the following lines for Kleborate as Kleborate usings the fasta file name
to determine the identifier that goes in the output file.
As we are using the raw assembly.fasta outputs from unicycler, all entries when combined are called
"assembly", which is obviously extremely unhelpful! So not using this for now until we work this out

Additionally, it is not possible to install Kleborate via conda, due to the reliance on Kaptive.
Therefore the bluedog environment on M3 will need to have Kleborate installed manually if we want to include this
*/

/*
process KLEBORATE {
  label "short_job"
  cache 'lenient'

  input:
  tuple val(isolate_id), file(fasta_file), file(graph_file), file(unicycler_log)

  output:
  path("*_results.txt")

  script:
  """
  /Users/jane/Documents/assembly_pipe_dsl2/Kleborate/kleborate-runner.py --resistance -o ${isolate_id}_results.txt -a $fasta_file
  """
}

process COMBINE_KLEB {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}"), mode: 'copy'

  input:
  path("*_results.txt")

  output:
  file("kleborate_results.txt")

  script:
  """
  awk 'FNR==1 && NR!=1 { while (/^strain/) getline; } 1 {print}' *_results.txt > kleborate_results.txt
  """
}
*/
