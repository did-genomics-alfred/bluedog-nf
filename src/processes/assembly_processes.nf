process mlst {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}/mlst"), mode: 'copy', saveAs: {filename -> "${isolate_id}_mlst.txt"}, pattern: '*_mlst.txt'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  path("${isolate_id}_mlst.txt")

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
  cache 'lenient'
  publishDir path:("${params.output_dir}/speciator"), mode: 'copy', pattern: '*_species.txt'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  file("${isolate_id}_species.txt")

  script:
  """
  singularity run --containall --pwd /bactinspector --bind $fasta_file:/bactinspector/input.fasta /projects/js66/software/speciator/speciator-3.0.1.sif > ${isolate_id}_species.json
  parse_speciator.py --json ${isolate_id}_species.json --output ${isolate_id}_species.txt
  """
}

process combine_speciator {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}"), mode: 'copy'

  input:
  path("*_species.txt")

  output:
  file("speciator_results.txt")

  script:
  """
  awk 'FNR==1 && NR!=1 { while (/^isolate/) getline; } 1 {print}' *_species.txt > speciator_results.txt
  """
}

process kleborate {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}/kleborate"), mode: 'copy', saveAs: {filename -> "${isolate_id}_kleborate.txt"}, pattern: '*_kleborate.txt'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  path("*_kleborate.txt")

  script:
  """
  kleborate --resistance -o ${isolate_id}_kleborate.txt -a $fasta_file
  """
}

process combine_kleborate {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}"), mode: 'copy'

  input:
  path("*_kleborate.txt")

  output:
  file("kleborate_results.txt")

  script:
  """
  awk 'FNR==1 && NR!=1 { while (/^strain/) getline; } 1 {print}' *_kleborate.txt > kleborate_results.txt
  """
}

process amrfinder {
  label "short_job"
  cache 'lenient'
  publishDir path:("${params.output_dir}/amrfinder"), mode: 'copy'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  path("*_amr_results.txt")

  script:
  if (params.amr_organism){
    """
    amrfinder -n $fasta_file --plus -o ${params.organism} > ${isolate_id}_amr_results.txt
    """
  }
  else {
    """
    amrfinder -n $fasta_file --plus > ${isolate_id}_amr_results.txt
    """
  }
}
