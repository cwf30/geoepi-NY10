#!/usr/bin/env python
# This script splits a nucleotide or amino acid alignment into separate files for each gene partition.

import sys
from Bio import SeqIO

# Define the nucleotide partitions (1-based coordinates)
partitions_nt = {
    'C': (1, 369), 'prM': (370, 645), 'M': (646, 870), 
    'E': (871, 2373), 'NS1': (2374, 3429), 'NS2A': (3430, 4122),
    'NS2B': (4123, 4515), 'NS3': (4516, 6372), 'NS4A': (6373, 6819),
    'NS4B': (6820, 7584), 'NS5': (7585, 10314)
}

# Convert nucleotide partitions to amino acid partitions (1-based coordinates)
partitions_aa = {
    k: ((start-1)//3 + 1, end//3) for k, (start, end) in partitions_nt.items()
}

# Determine the partition type and file extension based on the input file
if 'fna' in sys.argv[1]:  # If input is a nucleotide file
    partitions = partitions_nt
    ext = 'fna'
else:  # If input is an amino acid file
    partitions = partitions_aa
    ext = 'faa'

# Parse the input sequences from the FASTA file
with open(sys.argv[1], "r") as file:
    sequences = list(SeqIO.parse(file, "fasta"))  # Read the sequences into a list

# Create separate files for each gene partition
for gene, (start, end) in partitions.items():
    with open(f"{gene}_partition.{ext}", "w") as file:
        # Extract the gene-specific portion of each sequence and write it to the gene's file
        for sequence in sequences:
            gene_sequence = sequence[start-1:end]  # Extract the relevant range (0-based indexing)
            SeqIO.write(gene_sequence, file, "fasta")  # Write the gene sequence to the file
