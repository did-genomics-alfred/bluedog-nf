process mlst {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}/mlst"), mode: 'copy', saveAs: {filename -> "${isolate_id}_mlst.txt"}, pattern: '*_mlst.txt'

  input:
  //tuple val(isolate_id), path(fasta_file), path(graph_file), path(log_file)
  tuple val(isolate_id), path(fasta_file)

  output:
  path("{$isolate_id}_mlst.txt")

  script:
  """
  mlst --label $isolate_id $fasta_file > ${isolate_id}_mlst.txt
  """
}

process combine_mlst {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}"), mode: 'copy'

  input:
  path("*_mlst.txt")

  output:
  file("mlst_results.txt")

  script:
  """
  cat *_mlst.txt > mlst_results.txt
  """
}

process speciator {
  label "short_job"
  cache 'lenien'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  file("species_result.json")

  script:
  """
  singularity run --containall --pwd /bactinspector --bind $PWD/$fasta_file:/bactinspector/input.fasta speciator-3.0.1.sif 2>/dev/null
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
