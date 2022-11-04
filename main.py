#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.1"

# IMPORTS
import sys
import argparse


# FUNCTIONS
def parse_command_line_args():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(
        description="""Calculate the severity score of a given SNP in a given genomic coding sequence 
                       based on conservation in its protein family.\n 
                       (Defaults to sequence and family located in `/testdata` dir)""")

    # Add arguments
    parser.add_argument("location",
                        type=int, metavar="snp_location",
                        help="the location of where the SNP should be introduced")
    parser.add_argument("nucleotide",
                        type=str, metavar="snp_nucleotide",
                        help="the nucleotide the SNP should mutate to")
    parser.add_argument("-s", "--sequence",
                        type=str, metavar="sequence",
                        default="testdata/ZCCHC17_C-Lupus_DNA.fasta",
                        help="a fasta file with a single genomic coding sequence of a protein")
    parser.add_argument("-f", "--family",
                        type=str, metavar="protein_family_file",
                        default="testdata/ZCCHC17_Protein_Family.fasta",
                        help="a multi-fasta file or MSA file of the protein family "
                             "the given genomic coding sequence is from")

    # Parse the arguments and return
    return parser.parse_args()


def main():
    """Main function"""
    return 0


if __name__ == "__main__":
    sys.exit(main())
