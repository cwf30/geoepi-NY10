#!/usr/bin/env python
# This script filters FASTA sequences based on genome quality and logs the results.
# Sequences passing the quality threshold are saved to an output file.

import sys

# Define the name of the output file for sequences that pass the quality filter
output_file = "quality_sequences.fasta"

# Function to check the quality of a genome sequence
def quality_check(description, genome, quality_threshold, output_file):
    """
    Evaluates a genome sequence against a quality threshold based on its 'N' content.
    - Sequences shorter than 10,200 bases are automatically excluded.
    - Sequences with a percentage of 'N's below the threshold are written to the output file.
    - Logs the result (PASS/FAIL and N percentage) to a log file.

    Parameters:
        description (str): The FASTA header line for the sequence.
        genome (str): The genome sequence.
        quality_threshold (float): Maximum allowable percentage of 'N's for a sequence to pass.
        output_file (str): Path to the output file for sequences that pass.
    """
    # Calculate the total length of the genome
    total_length = len(genome)
    # Exclude sequences shorter than 10,200 bases
    if total_length < 10_200:
        return
    # Calculate the percentage of 'N' characters in the genome
    n_count = genome.count('N')
    n_percentage = (n_count / total_length) * 100
    # Determine if the sequence passes the quality filter
    status = 'PASS' if n_percentage < quality_threshold else 'FAIL'

    # If the sequence passes, write it to the output file
    if n_percentage < quality_threshold:
        with open(output_file, 'a') as f:
            f.write(description + genome)

    # Log the quality check result to a log file
    with open('quality_filter_log.txt', 'a') as i:
        i.write(f'{description.strip()} - {status} - {n_percentage:.2f}% N\'s\n')

# Main function to process the input FASTA file
def main(input_file, threshold, output_file):
    """
    Reads a FASTA file, evaluates each sequence for quality, and applies the quality filter.

    Parameters:
        input_file (str): Path to the input FASTA file.
        threshold (float): Maximum allowable percentage of 'N's for a sequence to pass.
        output_file (str): Path to the output file for sequences that pass.
    """
    with open(input_file, 'r') as f:
        genome = ''  # Initialize variable to store the genome sequence
        description = ''  # Initialize variable to store the FASTA header
        for line in f:
            if line.startswith('>'):  # Detect a new FASTA sequence header
                # Perform quality check on the previous sequence if it exists
                if genome:
                    quality_check(description, genome, threshold, output_file)
                # Update the description with the current header
                description = line
                genome = ''  # Reset the genome sequence for the new entry
            else:
                # Append the current line to the genome sequence
                genome += line
        # Perform quality check on the last sequence in the file
        if genome:
            quality_check(description, genome, threshold, output_file)

# Entry point of the script
if __name__ == '__main__':
    # Set the quality threshold (default to 5% if not provided)
    arg2 = 5 if len(sys.argv) < 3 else float(sys.argv[2])
    # Call the main function with the input file, threshold, and output file
    main(sys.argv[1], arg2, output_file)
