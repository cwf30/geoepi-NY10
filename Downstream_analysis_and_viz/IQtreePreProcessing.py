import pandas as pd
from datetime import datetime
from Bio import SeqIO

# -------------------- Configuration -------------------- #
# File paths
METADATA_PATH = '/path/to/metadata.csv'
ALIGNMENT_INPUT_PATH = "/path/to/masked_alignment_polyproteins.fna"
ALIGNMENT_OUTPUT_PATH = "/path/to/masked_alignment_polyproteins_no_dups.fna"
DATE_FILE_PATH = "/path/to/date_file.txt"

# Sequence trimming length
TRIM_LENGTH = 9504
# -------------------------------------------------------- #

# -------------------- Step 1: Trim Sequences -------------------- #
trimmed_records = []
sequence_dict = {}

for record in SeqIO.parse(ALIGNMENT_INPUT_PATH, "fasta"):
    trimmed_seq = record.seq[:TRIM_LENGTH]
    seq_str = str(trimmed_seq)
    if seq_str not in sequence_dict:
        sequence_dict[seq_str] = record
        trimmed_records.append(record)

# Write the non-redundant trimmed sequences to the output file
SeqIO.write(trimmed_records, ALIGNMENT_OUTPUT_PATH, "fasta")

# Identify removed duplicate IDs
all_ids = set(record.id for record in SeqIO.parse(ALIGNMENT_INPUT_PATH, "fasta"))
unique_ids = set(record.id for record in trimmed_records)

# -------------------- Step 2: Build Date File -------------------- #
# Read metadata
metadata_df = pd.read_csv(METADATA_PATH)

# Filter metadata for non-redundant sequences
filtered_metadata = metadata_df[
    (metadata_df['accession'].isin(unique_ids))
].copy()

# Function to format the collection date
def format_date(date_str):
    if len(date_str) == 4:
        return date_str  # Year only
    try:
        return datetime.strptime(date_str, '%d-%b-%Y').strftime('%Y-%m-%d')
    except ValueError:
        return date_str  # Return original if format is unexpected

# Apply date formatting
filtered_metadata['formatted_date'] = filtered_metadata['collection_date'].apply(format_date)

# Write the accession and formatted date to the date file
filtered_metadata[['accession', 'formatted_date']].to_csv(
    DATE_FILE_PATH,
    sep='\t',
    header=False,
    index=False
)

