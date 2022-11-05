#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.2"

# IMPORTS
import sys
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq


# CLASSES
class ProteinFamily:
    def __init__(self, protein_family_file):
        self.proteins = self.__read_multi_fasta(protein_family_file)

    @staticmethod
    def __read_multi_fasta(multi_fasta_file):
        return list(SeqIO.parse(multi_fasta_file, "fasta"))


class MutatedGene:
    def __init__(self, snp_location, snp_nucleotide, sequence_file):
        self.snp_location = snp_location
        self.snp_nucleotide = snp_nucleotide
        self.original_sequence = self.__read_from_fasta(sequence_file)

        # Perform stuffs
        self.mutated_sequence = self.__create_mutated_sequence()

    def __read_from_fasta(self, sequence_file):
        sequence = SeqIO.read(sequence_file, "fasta").seq
        if not 1 <= self.snp_location <= len(sequence):
            msg = f"Given location ({self.snp_location}) is not a valid location on given sequence! " \
                  f"Please give a location between 1 and length of given sequence: {len(sequence)}."
            raise IndexError(msg)
        return sequence

    def __create_mutated_sequence(self):
        old_nucleotide = self.original_sequence[self.snp_location-1].upper()
        new_nucleotide = self.snp_nucleotide.upper()
        # Raise ValueError if snp wouldn't create a nucleotide substitution
        if old_nucleotide == self.snp_nucleotide:
            msg = f"Nucleotide substitution at location {self.snp_location} unsuccessful, " \
                  f"nucleotide was already '{old_nucleotide}'!"
            raise ValueError(msg)

        # Create MutableSeq object and substitute nucleotide at position given
        mutable_seq = MutableSeq(self.original_sequence)
        mutable_seq[self.snp_location-1] = new_nucleotide
        print(f"Nucleotide substitution at location {self.snp_location} successful! "
              f"Mutation: '{old_nucleotide}' was substituted with '{new_nucleotide}'.")
        return Seq(mutable_seq)


class SeverityCalculator:
    pass
    def __init__(self):
        pass
