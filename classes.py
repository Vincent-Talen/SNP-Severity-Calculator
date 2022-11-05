#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.3"

# IMPORTS
import sys
from math import ceil
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, MutableSeq
from Bio.Align.Applications import ClustalOmegaCommandline


# CLASSES
class ProteinFamily:
    def __init__(self, protein_family_file):
        self.new_msa_file = "output/msa.fasta"
        self.alignment = self.__create_msa(protein_family_file)

    def __create_msa(self, in_file):
        # Create alignment
        clustalomega_cline = ClustalOmegaCommandline(
            "./clustalo", infile=in_file, outfile=self.new_msa_file,
            verbose=True, auto=True, force=True
        )
        clustalomega_cline()

        # Read alignment from file and return as alignment object
        return AlignIO.read(self.new_msa_file, "fasta")


class MutatedGene:
    def __init__(self, snp_position, snp_nucleotide, sequence_file):
        self.snp_position = snp_position
        self.snp_nucleotide = snp_nucleotide.upper()

        # Load original sequence and create mutated sequence with SNP
        self.original_sequence = self.__read_from_fasta(sequence_file)
        self.mutated_sequence = self.__create_mutated_sequence()

        # Translate both sequences to protein
        self.original_protein = self.original_sequence.translate()
        self.mutated_protein = self.mutated_sequence.translate()

    def __read_from_fasta(self, sequence_file):
        sequence = SeqIO.read(sequence_file, "fasta").seq
        if not 1 <= self.snp_position <= len(sequence):
            print(f"ERROR!\tGiven position ({self.snp_position}) is not a valid position for the given sequence!\n"
                  f"\tPlease give a position between 1 and the length of the given sequence: {len(sequence)}")
            sys.exit(1)
        return sequence

    def __create_mutated_sequence(self):
        old_nucleotide = self.original_sequence[self.snp_position-1].upper()
        
        # Script can be stopped if nucleotide isn't changed
        if old_nucleotide == self.snp_nucleotide:
            print(f"Nucleotide substitution at position {self.snp_position} unsuccessful, "
                  f"nucleotide was already '{old_nucleotide}'!")
            sys.exit(0)

        # Create MutableSeq object and substitute nucleotide at the given position
        mutable_seq = MutableSeq(self.original_sequence)
        mutable_seq[self.snp_position-1] = self.snp_nucleotide
        print(f"Nucleotide substitution at position {self.snp_position} successful! "
              f"Mutation: '{old_nucleotide}' was substituted with '{self.snp_nucleotide}'")
        return Seq(mutable_seq)


class SeverityCalculator:
    def __init__(self):
        pass
