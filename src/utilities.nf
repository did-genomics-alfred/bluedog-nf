/* Note that these functions are drawn from Stephen Watt's functions in RedDog https://github.com/scwatts/reddog-nf/blob/master/src/utilities.nf */


def check_host(workflow) {
  // Do not run on MASSIVE unless user specifies profile to use to avoid inadvertently using a local executor
  massive_hostnames = ['m3-login1', 'm3-login2']
  on_massive = massive_hostnames.contains(InetAddress.getLocalHost().getHostName())
  profile_explicit = workflow.commandLine.tokenize(' ').contains('-profile')
  if (on_massive && ! profile_explicit) {
    exit 1, "ERROR: to run on MASSIVE you must explicitly set -profile massive"
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

