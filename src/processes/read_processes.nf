process fastqc {
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

process multiqc {
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

process assemble {
  cache 'lenient'
  publishDir path: {"${params.output_dir}/fasta"},  mode: 'copy', saveAs: {filename -> "${isolate_id}.fasta"}, pattern: '*.fasta'
  publishDir path: {"${params.output_dir}/graph"}, mode: 'copy', saveAs: {filename -> "${isolate_id}.gfa"}, pattern: 'assembly.gfa'
  publishDir path: {"${params.output_dir}/logs"}, mode: 'copy', saveAs: {filename -> "${isolate_id}_unicycler.log"}, pattern: '*.log'

  input:
  tuple val(isolate_id), path(reads_fwd), path(reads_rev)

  output:
  tuple val(isolate_id), path("assembly.fasta"), emit: assembly_fasta
  path("assembly.gfa")
  path("unicycler.log")


  script:
  """
  unicycler -1 ${reads_fwd} -2 ${reads_rev} -o .
  """

}

process assembly_stats {
  label "short_job"
  cache 'lenient'
  publishDir path:{"${params.output_dir}/stats"}, mode: 'copy', saveAs: {filename -> "${isolate_id}_stats.txt"}, pattern: '*.txt'

  input:
  tuple val(isolate_id), path(fasta_file)

  output:
  path("${isolate_id}_stats.txt")

  script:
  """
  assembly_stats.py -a $fasta_file --id $isolate_id > ${isolate_id}_stats.txt
  """
}

process combine_stats {
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
