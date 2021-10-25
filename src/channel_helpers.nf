// Get the prefix for an input readset and whether it is pe or se
def get_read_prefix_and_type(filepath) {
  // NOTE: nf requires escaping '$'
  regex_pe = "^(.+?)_R?[12](?:_001)?.fastq(?:.gz)?\$"

  java.util.regex.Matcher matcher;
  String read_type;

  if ((matcher = (filepath.getName() =~ /$regex_pe/))) {
    read_type = 'pe'
  } else {
    exit 1, "ERROR: did not find any readsets with the provided glob: ${filepath}"
  }
  return [read_type, matcher.group(1), filepath]
}
