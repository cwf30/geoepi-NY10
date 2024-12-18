#!/usr/bin/env python
# This script processes input FASTA files, simplifying and unifying the format of sequence IDs and saving the results to a new file.

import sys

# Define the name of the output file for sanitized sequences
sanitized_file = 'sanitized_sequences.fasta'

# Open the output file for writing
with open(sanitized_file, 'w') as sanitized:
    # Iterate through each FASTA file provided as a command-line argument
    for fasta_file in sys.argv[1:]:
        # Read and parse the input FASTA file
        with open(fasta_file, 'r') as fasta:
            fasta_data = fasta.read()
            sequences = fasta_data.split('>')[1:]  # Split into individual sequences

        # Process each sequence
        for sequence in sequences:
            # Extract and simplify the sequence ID
            sequence_lines = sequence.strip().split('\n')
            sequence_id = sequence_lines[0].split(' ')[0].split('.')[0].split('_')[0]
            sequence_data = '\n'.join(sequence_lines[1:])  # Combine sequence lines

            # Write sanitized sequence to the output file
            sanitized.write(f'>{sequence_id}\n{sequence_data}\n')
