#!/usr/bin/env python
# This script filters a CSV file to include only rows corresponding to a list of sequence IDs provided in a separate file.

import sys
import os
import csv

# Parse the input CSV file containing sequence data
with open(sys.argv[1], "r") as file:
    reader = csv.reader(file)
    header = next(reader)  # Extract the header row
    sequences = list(reader)  # Store the remaining rows as a list

# Read the sequence IDs from the second input file
with open(sys.argv[2], "r") as sequenceIDs_file:
    sequenceIDs = [line.strip() for line in sequenceIDs_file]  # Strip whitespace from each line

# Generate the output file name based on the input CSV file name
csv_file = os.path.basename(sys.argv[1])  # Extract the base name of the input CSV file
output_file_name = "subset_" + csv_file  # Prefix the base name with "subset_"

# Create a new CSV file containing only the filtered sequences
with open(output_file_name, "w") as output_file:
    writer = csv.writer(output_file)
    writer.writerow(header)  # Write the header row to the output file
    for sequence in sequences:
        if sequence[0] in sequenceIDs:  # Check if the sequence ID is in the provided list
            writer.writerow(sequence)  # Write the row to the output file
