/* Note that these functions are drawn from Stephen Watt's functions in RedDog https://github.com/scwatts/reddog-nf/blob/master/src/utilities.nf */

def print_splash() {
  log.info('--------------------------------------------')
  log.info("""
  ‌‌

  _       _                      _
 | |__   | |  _   _    ___    __| |   ___     __ _
 | '_ \  | | | | | |  / _ \  / _` |  / _ \   / _` |
 | |_) | | | | |_| | |  __/ | (_| | | (_) | | (_| |
 |_.__/  |_|  \__,_|  \___|  \__,_|  \___/   \__, |
                                             |___/
  ‌‌
  """.stripIndent())
  log.info('-------------------------------------------')
}

def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile massive"
  }
}

def check_output_dir(params) {
  // Do not run if output exists and contains files other than the run info directory (which is created by this point)
  output_dir_files = []
  output_dir = file(params.output_dir)
  output_dir.eachFile { output_dir_files.add(it.name) }
  //run_info_dirname = file(params.run_info_dir).simpleName
  //output_dir_files.remove(run_info_dirname)
  if (output_dir_files.size() > 0 && ! params.force) {
    exit 1, "ERROR: output directory '${output_dir}' already exists and contains other files, remove or use --force to overwrite"
  }
}

def check_arguments(params) {
  // Check required input and outputs
  if (! params.reads && ! params.assemblies) {
    exit 1, "ERROR: Either 'reads' or 'assemblies' must be set in nextflow.config"
  }
  if (params.reads && params.assemblies){
    exit 1, "ERROR: Cannot set both 'reads' and 'assemblies' in nextflow.config, must choose one"
  }
  if (! params.output_dir) {
    exit 1, "ERROR: option 'output_dir' must be set in nextflow.configrequired"
  }
}

// For an optional stage param variable, check that it is either a Boolean or String
// If it is a string and either 'true' or 'false', return the boolean equivalent
def check_boolean_option(option, name) {
  if (option.getClass() == java.lang.Boolean) {
    return option
  } else if (option.getClass() == java.lang.String) {
    if (option.toLowerCase() == 'true') {
      return true
    } else if (option.toLowerCase() == 'false') {
      return false
    }
  }
  exit 1, "ERROR: ${name} option must be true or false"
}

def write_param_data_to_run_config() {
  File run_info_fh = new File("${params.run_info_dir}/run_config.tsv")
  if (params.reads){
    run_info_fh.append("reads\t${params.reads}\n")
  }
  if (params.assemblies){
    run_info_fh.append("assemblies\t${params.reads}\n")
  }
  run_info_fh.append("read_qc\t${params.read_qc}\n")
  run_info_fh.append("speciator\t${params.run_speciator}\n")
  run_info_fh.append("mlst\t${params.run_mlst}\n")
  run_info_fh.append("kleborate\t${params.run_kleborate}\n")
}
