#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.1.1"

# IMPORTS
import sys
import argparse
from classes import *


# FUNCTIONS
def parse_command_line_args():
    # Create ArgumentParser object
    parser = argparse.ArgumentParser(
        description="""Calculate the severity score of a given SNP in a given genomic coding sequence 
                       based on amino acid conservation in its protein family.\n 
                       (Defaults to sequence and family located in `/testdata` dir)""")

    # Add arguments
    parser.add_argument("position",
                        type=int, metavar="snp_position",
                        help="the position of where the SNP should be introduced")
    parser.add_argument("nucleotide",
                        type=str, metavar="snp_nucleotide",
                        choices=["A", "a", "T", "t", "C", "c", "G", "g"],
                        help="the nucleotide the SNP should mutate to")
    parser.add_argument("-s", "--sequence",
                        type=str, metavar="sequence",
                        default="testdata/ZCCHC17_C-Lupus_DNA.fasta",
                        help="a fasta file with a single genomic coding sequence of a protein")
    parser.add_argument("-f", "--family",
                        type=str, metavar="protein_family_file",
                        default="testdata/ZCCHC17_Protein_Family.fasta",
                        help="a multi-fasta file with protein sequences like normal or as multi sequence alignment, "
                             "of the protein family the given genomic coding sequence belongs to")
    parser.add_argument("--show_all_frequencies", action='store_true',
                        help="if this argument is used, the frequencies of amino acids for each position will be shown")

    # Parse the arguments and return
    return parser.parse_args()


def main():
    """Main function"""
    # Parse command line arguments
    args = parse_command_line_args()

    # Create MutatedGene object
    mutated_gene = MutatedGene(args.position, args.nucleotide, args.sequence)
    # Create ProteinFamily object
    protein_family = ProteinFamily(args.family)

    # Check severity of MutatedGene in ProteinFamily
    protein_family.check_snp_severity(mutated_gene)

    # If user desires, print frequencies for all positions
    if args.show_all_frequencies:
        protein_family.print_all_frequencies()
    return 0


if __name__ == "__main__":
    sys.exit(main())
