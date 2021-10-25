#!/usr/bin/env python3

import json
from argparse import ArgumentParser

def get_arguments():
    parser = ArgumentParser(description='Parse speciator json file')

    # job submission options
    parser.add_argument('--json', type=str, required=True, help='json file from speciator')
    parser.add_argument('--output', type=str, required=True, help='name of output file')

    return parser.parse_args()

def main():

    args = get_arguments()

    with open(args.json) as f:
        result_dict = json.load(f)

    iso_name = args.json.split('_species.json')[0]

    header = ['isolate', 'speciesName', 'confidence', 'referenceId', 'mashDistance', 'matchingHashes', 'source']

    with open(args.output, 'w') as out:
        out.write('\t'.join(header) + '\n')
        out.write(iso_name + '\t')
        for val in header[1:]:
            out.write(str(result_dict[val]) + '\t')

if __name__ == '__main__':
    main()

