# /usr/bin/env python3

"""
Read BED lines from stdin, filter for a set of known chromosomes, write passing lines and headers to stdout.
This script is used to filter unknown chromosomes from annotation files. For example, if a no-ALT version of the
reference is used, this script will remove all records on ALT chromosomes. Anything not in the FAI file supplied
to the script is dropped.
"""

import argparse
import sys

if __name__ == '__main__':

    # Get command line arguments
    parser = argparse.ArgumentParser('Filter unknown (not in FAI) records from BED (in stdin) and write to stdout.')

    parser.add_argument('-g', '--fai', help='Reference FAI file')

    args = parser.parse_args()

    # Read contigs
    with open(args.fai, 'r') as in_file:
        chrom_set = {
            line.split('\t', 1)[0] for line in [
                line2.strip() for line2 in in_file
            ] if line and not line.startswith('#')
        }

    # Process input
    for line in sys.stdin:
        if line.startswith('#') or line.split('\t', 1)[0] in chrom_set:
            sys.stdout.write(line)
