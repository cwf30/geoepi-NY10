#!/usr/bin/env python
# This script filters a FASTA file to include only sequences with IDs present in a separate list of sequence IDs.

import sys
import os
from Bio import SeqIO

# Parse the input FASTA file and load all sequences into a list
with open(sys.argv[1], "r") as file:
    sequences = list(SeqIO.parse(file, "fasta"))  # Read sequences in FASTA format

# Read the sequence IDs from the second input file
with open(sys.argv[2], "r") as sequenceIDs_file:
    sequenceIDs = [line.strip() for line in sequenceIDs_file]  # Strip whitespace from each line

# Generate the output file name by prefixing "subset_" to the input FASTA file name
fasta_file = os.path.basename(sys.argv[1])  # Extract the base name of the input FASTA file
output_file_name = "subset_" + fasta_file  # Prefix "subset_" to the file name

# Write a new FASTA file containing only the filtered sequences
with open(output_file_name, "w") as output_file:
    for sequence in sequences:
        if sequence.id in sequenceIDs:  # Check if the sequence ID is in the provided list
            SeqIO.write(sequence, output_file, "fasta")  # Write the sequence to the output file
