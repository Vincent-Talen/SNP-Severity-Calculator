#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.4"

# IMPORTS
import sys
from math import ceil
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, MutableSeq
from Bio.Align.Applications import ClustalOmegaCommandline


# CLASSES
class MutatedGene:
    def __init__(self, snp_position: int, snp_nucleotide: str, sequence_file: str):
        self.__snp_genomic_pos = snp_position
        self.__snp_nucleotide = snp_nucleotide.upper()
        self.snp_protein_pos = ceil(self.__snp_genomic_pos / 3)

        # Load original sequence and create mutated sequence with SNP
        self.gene_info = dict()
        self.original_sequence = self.__read_from_fasta(sequence_file)
        self.mutated_sequence = self.__create_mutated_sequence()

        # Translate both sequences to protein
        self.original_protein = self.original_sequence.translate()
        self.mutated_protein = self.mutated_sequence.translate()

        # Script can be stopped if synonymous
        self.__synonymous_check()

    def __read_from_fasta(self, sequence_file: str) -> Seq:
        # Load from file and save id, name and description to self.gene_info
        seq_record = SeqIO.read(sequence_file, "fasta")
        self.gene_info["id"] = seq_record.id
        self.gene_info["name"] = seq_record.name
        self.gene_info["description"] = seq_record.description

        # Get sequence and check if given snp position is valid
        sequence = seq_record.seq
        if not 1 <= self.__snp_genomic_pos <= len(sequence):
            print(f"ERROR!\tGiven position ({self.__snp_genomic_pos}) is not a valid position for the given sequence!\n"
                  f"\tPlease give a position between 1 and the length of the given sequence: {len(sequence)}")
            sys.exit(1)
        return sequence

    def __create_mutated_sequence(self) -> Seq:
        old_nucleotide = self.original_sequence[self.__snp_genomic_pos-1].upper()
        # Stop script if nucleotide hasn't changed
        if old_nucleotide == self.__snp_nucleotide:
            print(f"Nucleotide substitution at position {self.__snp_genomic_pos} unsuccessful, "
                  f"nucleotide was already '{old_nucleotide}'!")
            sys.exit(0)

        # Create MutableSeq object and substitute nucleotide at the given position
        mutable_seq = MutableSeq(self.original_sequence)
        mutable_seq[self.__snp_genomic_pos-1] = self.__snp_nucleotide
        print(f"Nucleotide substitution at position {self.__snp_genomic_pos} successful! "
              f"Mutation: '{old_nucleotide}' was substituted with '{self.__snp_nucleotide}'")

        # Return Seq object of mutated sequence
        return Seq(mutable_seq)

    def __synonymous_check(self):
        # Get original and mutated amino acid
        original_aa = self.original_protein[self.snp_protein_pos - 1]
        mutated_aa = self.mutated_protein[self.snp_protein_pos - 1]

class SeverityCalculator:
    def __init__(self):
        pass
        # If the 'mutated' amino acid is the same as the original stop the script
        if original_aa == mutated_aa:
            print(f"Synonymous! SNP mutation did not cause a change "
                  f"of amino acid '{original_aa}' on position {self.snp_protein_pos}")
            sys.exit(0)
        else:
            print(f"Non-Synonymous! SNP mutation caused change of amino acid "
                  f"on position {self.snp_protein_pos}: '{original_aa}' to '{mutated_aa}'!")


class ProteinFamily:
    def __init__(self, protein_family_file: str):
        new_msa_file = "output/msa.fasta"
        self.alignment = self.__create_msa(protein_family_file, new_msa_file)

    @staticmethod
    def __create_msa(in_file: str, out_file: str) -> AlignIO.MultipleSeqAlignment:
        # Create alignment command line task
        clustalomega_cline = ClustalOmegaCommandline(
            "./clustalo", infile=in_file, outfile=out_file,
            verbose=True, auto=True, force=True
        )
        # Perform alignment command line task
        clustalomega_cline()

        # Read alignment from file and return as alignment object
        return AlignIO.read(out_file, "fasta")
