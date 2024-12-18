#!/usr/bin/env python
# This script processes a multiple sequence alignment of amino acid sequences
# and compares them to a reference sequence to identify polymorphisms. The output includes:
# - A sparse CSV file indicating presence/absence of polymorphisms.
# - A dense CSV file detailing specific amino acid substitutions.
# - An optional subset CSV file focusing on user-specified polymorphisms.

import sys
import csv
from Bio import SeqIO 

# Input arguments
AA_Alignment = sys.argv[1]  # Path to amino acid alignment file
path_to_reference = sys.argv[2]  # Path to reference sequence file

# Check if a subset of mutations is provided
subset_mutations = False
if sys.argv[3] != "NO_FILE":
    AA_polymorphisms = sys.argv[3]
    subset_mutations = True

# Output file names
output_CSV_sparse = "polymorphisms_sparse.csv"
output_CSV_dense = "polymorphisms_dense.csv"
output_polymorphism_subset = "key_polymorphisms_sparse.csv"

# Define NY99 protein boundaries (1-based indices)
NY_99_proteins = {
    "C": (1, 123), "prM": (124, 215), "M": (216, 290), "E": (291, 791), 
    "NS1": (792, 1143), "NS2A": (1144, 1374), "NS2B": (1375, 1505),
    "NS3": (1506, 2124), "NS4A": (2125, 2273), "NS4B": (2274, 2528),
    "NS5": (2529, 3438)
}

# Read amino acid sequences from alignment
AA_sequences = {}
with open(AA_Alignment, "r") as file:
    sequences = SeqIO.parse(file, "fasta")
    for sequence in sequences:
        AA_sequences[sequence.id] = sequence.seq  # Store sequences by accession ID

# Read the reference sequence
reference_sequence = ""
with open(path_to_reference, "r") as file:
    sequences = SeqIO.parse(file, "fasta")
    for sequence in sequences:
        reference_sequence = str(sequence.seq)  # Use the first sequence as the reference
        break

# Read subset mutations if provided
if subset_mutations:
    with open(AA_polymorphisms, "r") as file:
        subset_of_mutations = [line.strip() for line in file]

# Initialize data structures
mutations_sparse = set()
mutations_dense = set()
mutation_data = {
    seq_id: {'mutations_dense': {}, 'mutations_sparse': []} for seq_id in AA_sequences.keys()
}

# First pass: Identify mutations and populate sparse data
for sequence_id, sequence in AA_sequences.items():
    for i, (seq, ref) in enumerate(zip(sequence, reference_sequence)):
        if seq == ref or seq == "X" or ref == "X" or seq == "-":
            continue  # Skip matches, ambiguous amino acids, or gaps

        # Determine the protein and relative position of the mutation
        for protein, positions in NY_99_proteins.items():
            if positions[0] <= i + 1 <= positions[1]:
                relative_position = i - positions[0] + 2
                mut_sparse = f'{protein}-{ref}{relative_position}{seq}'  # E.g., "E-V159A"
                mutation_data[sequence_id]['mutations_sparse'].append(mut_sparse)
                mutations_sparse.add(mut_sparse)

                mut_dense = f'{protein}-{ref}{relative_position}'  # E.g., "E-V159"
                mutations_dense.add(mut_dense)

# Second pass: Populate dense data
for sequence_id, sequence in AA_sequences.items():
    for i, (seq, ref) in enumerate(zip(sequence, reference_sequence)):
        for protein, positions in NY_99_proteins.items():
            if positions[0] <= i + 1 <= positions[1]:
                relative_position = i - positions[0] + 2
                mut_dense = f'{protein}-{ref}{relative_position}'
                if mut_dense not in mutations_dense:
                    continue
                mutation_data[sequence_id]['mutations_dense'][mut_dense] = seq

# Write dense output file
sorted_mutations = sorted(list(mutations_dense), key=lambda x: (x.split("-")[0], int(x.split("-")[1][1:])))
with open(output_CSV_dense, "w") as file:
    writer = csv.writer(file)
    writer.writerow(['accession'] + sorted_mutations)
    for sequence_id, mutations in mutation_data.items():
        row = [sequence_id] + [mutations['mutations_dense'].get(mut, '-') for mut in sorted_mutations]
        writer.writerow(row)

# Write sparse output file
sorted_mutations = sorted(list(mutations_sparse), key=lambda x: (x.split("-")[0], int(x.split("-")[1][1:-1])))
with open(output_CSV_sparse, "w") as file:
    writer = csv.writer(file)
    writer.writerow(['accession'] + sorted_mutations)
    for sequence_id, mutations in mutation_data.items():
        row = [sequence_id] + ['1' if mut in mutations['mutations_sparse'] else '0' for mut in sorted_mutations]
        writer.writerow(row)

# Exit if no subset file is provided
if not subset_mutations:
    sys.exit(0)

# Write subset output file
filtered_mutation_data = {
    sequence_id: data for sequence_id, data in mutation_data.items()
    if any(mut in data['mutations_sparse'] for mut in subset_of_mutations)
}
filtered_sorted_mutations = [mut for mut in sorted_mutations if mut in subset_of_mutations]
with open(output_polymorphism_subset, "w") as file:
    writer = csv.writer(file)
    writer.writerow(['accession'] + filtered_sorted_mutations)
    for sequence_id, data in filtered_mutation_data.items():
        row = [sequence_id] + ['1' if mut in data['mutations_sparse'] else '0' for mut in filtered_sorted_mutations]
        writer.writerow(row)
