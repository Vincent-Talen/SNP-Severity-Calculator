#!/usr/bin/env python3

#  Copyright (c) 2022. Vincent Talen.
#  Licensed under GPLv3. See LICENSE file.

"""
SNP Severity Calculator
"""

# METADATA
__author__ = "Vincent Talen"
__version__ = "0.5"

# IMPORTS
import sys
from math import ceil
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
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
        self.original_protein = self.original_sequence.translate().rstrip("*")
        self.mutated_protein = self.mutated_sequence.translate().rstrip("*")

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

        # If the 'mutated' amino acid is the same as the original stop the script
        if original_aa == mutated_aa:
            print(f"Synonymous! SNP mutation did not cause a change "
                  f"in amino acid '{original_aa}' on position {self.snp_protein_pos}")
            sys.exit(0)
        elif mutated_aa == "*":
            print(f"Nonsense non-synonymous! SNP mutation caused amino acid "
                  f"'{original_aa}' on position {self.snp_protein_pos} to change to 'STOP'!")
            sys.exit(0)
        else:
            print(f"Missense non-synonymous! SNP mutation caused change of amino acid "
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

    @staticmethod
    def __get_alignment_adjusted_position(aligned_sequence, desired_position) -> int:
        position_without_gaps = 0
        # For each position in the aligned sequence
        for cur_pos, char in enumerate(aligned_sequence):
            # If the character is not a gap it is an amino acid
            if char != '-':
                # If the position without gaps is the desired position; return
                if position_without_gaps == desired_position:
                    return cur_pos
                # If not the desired position, count up one position
                position_without_gaps += 1

    @staticmethod
    def __normalize_freq_dict(freq_dict) -> dict:
        normalized_freq_dict = dict()
        # Calculate minimum value and value range for normalization
        freq_min = list(freq_dict.values())[-1]
        freq_range = list(freq_dict.values())[0] - freq_min
        for aa, freq in freq_dict.items():
            # Perform normalization for each amino acid frequency
            normalized_freq = (freq - freq_min) / freq_range
            normalized_freq_dict[aa] = normalized_freq
        return normalized_freq_dict

    @staticmethod
    def __get_column_freq_dict(current_column) -> dict:
        # Get dictionary with counts per amino acid
        counts_dict = {"A": 0, "R": 0, "N": 0, "D": 0, "C": 0, "E": 0, "Q": 0,
                       "G": 0, "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0,
                       "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0, "-": 0}
        for char in current_column:
            counts_dict[char] += 1

        # Whilst sorting convert counts to frequencies
        sorted_freq_dict = dict()
        for aa, count in sorted(counts_dict.items(), key=lambda item: item[1], reverse=True):
            freq = count / len(current_column)
            sorted_freq_dict[aa] = freq
        return sorted_freq_dict

    def __align_mutated_protein(self, mutated_gene: MutatedGene) -> Seq:
        # Create SeqRecord object for mutated gene protein sequence
        mutated_protein_record = SeqRecord(
            mutated_gene.mutated_protein,
            id=mutated_gene.gene_info["id"] + ".WITH_SNP_MUTATION ",
            name=mutated_gene.gene_info["name"],
            description=mutated_gene.gene_info["description"]
        )

        # Filenames
        combined_fasta_file = "output/combined_proteins.fasta"
        combined_msa_file = "output/combined_msa.fasta"

        # Create list with all records combined and write to file
        all_records = list(self.alignment)
        all_records.insert(0, mutated_protein_record)
        SeqIO.write(all_records, combined_fasta_file, "fasta")

        # Create new MSA and return aligned mutated protein sequence
        combined_msa = self.__create_msa(combined_fasta_file, combined_msa_file)
        return combined_msa[0].seq

    @staticmethod
    def __dict_print_string(given_dict) -> str:
        print_list = list()
        for aa, freq in given_dict.items():
            if freq != 0:
                print_list.append("'{}': {:.2f}".format(aa, freq))
        return ", ".join(print_list)

    def check_snp_severity(self, mutated_gene: MutatedGene):
        # Align the mutated protein
        aligned_protein = self.__align_mutated_protein(mutated_gene)

        # Get adjusted position of snp in aligned sequence
        unaligned_snp_pos = mutated_gene.snp_protein_pos - 1
        aa_position = self.__get_alignment_adjusted_position(aligned_protein, unaligned_snp_pos)

        # Get mutated amino acid and protein family frequency dict for affected column
        mutated_aa = aligned_protein[aa_position]
        affected_column = self.alignment[:, aa_position]
        freq_dict = self.__get_column_freq_dict(affected_column)

        # Print score/probability is the normalized frequency
        normalized_freq_dict = self.__normalize_freq_dict(freq_dict)
        print(f"SNP mutation has position {aa_position+1} in the MSA with "
              f"the protein family consisting of {len(self.alignment)} genes")
        print(f"\tAmino acid frequencies for position: " + self.__dict_print_string(freq_dict))
        print(f"\tSeverity score: {normalized_freq_dict[mutated_aa]} (<0.05 is deleterious)")

    def print_all_frequencies(self):
        print("\nShowing amino acid frequencies for all positions:")
        # For each column/position in the alignment
        for i in range(self.alignment.get_alignment_length()):
            # Get the frequency dict for that column
            current_column = self.alignment[:, i]
            freq_dict = self.__get_column_freq_dict(current_column)

            # Print stuff
            print(f"\tMSA Position {i + 1} ->\t" + self.__dict_print_string(freq_dict))
