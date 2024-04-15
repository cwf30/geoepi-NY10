# Genomic Data Cleaning and Processing

This notebook uses a metadata file containing information about WNV genomes (for all genomes with county level location data available publically as of Feb 2024) to aggregate all the relevant genomes into various formats needed for downstream analysis. In this notebook I will be:

1) Renaming sequences to a common format (accession only) for easier cross-referencing with metadata file later
2) Filtering sequences by 95% 'N' to exclude low-quality sequences
3) Extracting the polyprotein sequences and translating them
4) Saving the nucleotide and amino acid sequence data to multiFASTA files


```python
import os
import csv

#All genomes are in fasta_folder_path, but within multiple folders
fasta_folder_path = '/Users/cwcf/Documents/ORISE/data/Genomic'
#metadata path
csv_file_path = '/Users/cwcf/Documents/ORISE/data/Data_products/data_summary/metadata.csv'

#list of output files and their paths produced by this notebook:

#contains all genomes that are described in the metadata file
genomes_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/all_WNV_feb_2024.fna'
#same genomes above, but only the nucleotide sequence for the polyprotein
poly_protein_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/polyprotein_feb_2024.fna'
#as above, but translated using standard genetic code
poly_protein_translated_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/polyprotein_feb_2024.faa'

#contains all genomes in metadata file, and with <5% ambiguous base calls
high_quality_genomes_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/high_quality_WNV_feb_2024.fna'
#same genomes above, but only the nucleotide sequence for the polyprotein
high_quality_poly_protein_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/high_quality_polyprotein_feb_2024.fna'
#as above, but translated using standard genetic code
high_quality_poly_protein_translated_path = '/Users/cwcf/Documents/ORISE/data/Data_products/multiFASTA/high_quality_polyprotein_feb_2024.faa'

```


```python
# Recursive function to find all fasta files in the folder
def find_fasta_files(folder_path):
    fasta_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.fasta') or file.endswith('.fa'):
                fasta_files.append(os.path.join(root, file))
    return fasta_files
```


```python
all_genome_paths = find_fasta_files(fasta_folder_path)
print(f'Number of fasta files found: {len(all_genome_paths)}')
```

    Number of fasta files found: 1264



```python
# open the metadata file and extract the accession numbers
accession_numbers = []
with open(csv_file_path, 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    for row in csv_reader:
        accession_numbers.append(row['accession'].strip())
print(f'Number of accession numbers: {len(accession_numbers)}')
```

    Number of accession numbers: 2174


## Cross-Referencing Genomes With Metadata

Within each fasta file are multiple genomes that we need to go through and check against those that we have good metadata for.

We can just remove genomes if they are not in the metadata file, but to do that we need to open each FASTA file:


```python
genomes_with_metadata = []
accession_with_metadata = []
for fasta_file in all_genome_paths:
    with open(fasta_file, 'r') as fasta:
        fasta_data = fasta.read()
        sequences = fasta_data.split('>')[1:]
        for sequence in sequences:            
            #strip everything from the sequence description except the accesion number. 
            sequence_id = sequence.strip().split('\n')[0].split(' ')[0].split('.')[0].split('_')[0]
            # Look for the sequence's accession number in the metadata file
            if any(accession == sequence_id for accession in accession_numbers):
                #change the sequence id in sequence to sequence_id
                sequence = sequence.replace(sequence.strip().split('\n')[0], sequence_id)
                genomes_with_metadata.append(sequence)
                accession_with_metadata.append(sequence_id)
accession_with_metadata = list(set(accession_with_metadata))

print(f'Number of genomes with metadata: {len(accession_with_metadata)}\n')          
print(f'Sample data from genomes_with_metadata: \n{genomes_with_metadata[0][:60]}...')

```

    Number of genomes with metadata: 2174
    
    Sample data from genomes_with_metadata: 
    KX547621
    AGTGTTTGTGAGGATTAACAACAATTAACACAGTGCGAGCTGTTTCTTAGC...


### We have genomes for every row in our metadata file, and can now save these to all_WNV_feb_2024.fna


```python
with open(genomes_path, 'w') as output_file:
    all_genomes_fasta_data = ''.join(['>' + sequence for sequence in genomes_with_metadata])
    output_file.write(all_genomes_fasta_data)
```

## Quality Filtering

If a genome is longer than 10.4kb and has fewer than than 5% ambiguous base calls ('N'), we want to exclude from almost all analyses. let's remove anythign that doesn't meet this threshold. We will also add in a new column for the metadata to denote that the genome is low quality.


```python
high_quality_genomes = []
high_quality_accessions = []
for genome in genomes_with_metadata:
    total_length = len(genome)
    if total_length < 10_400:
        continue
    n_count = genome.count('N')
    n_percentage = (n_count / total_length) * 100

    if n_percentage < 5:
        high_quality_genomes.append(genome)
        high_quality_accessions.append(genome.split('\n')[0].strip())

print(f'number of high quality genomes that are in our metadata file: {len(high_quality_accessions)}')
```

    number of high quality genomes that are in our metadata file: 1838


Now to update the metadata file with sequence quality:


```python
# Open the CSV file
with open(csv_file_path, 'r') as file:
    reader = csv.reader(file)
    rows = list(reader)

# Strip leading and trailing spaces from each value
for row in rows:
    for i in range(len(row)):
        row[i] = row[i].strip()

# Add a new column called '<5% N' if it doesn't already exist
header = rows[0]
if '<5% N' not in header:
    header.append('<5% N')
    # Update the rows with the new column values
    for row in rows[1:]:
        row.append('')  # Add an empty value for the new column
else:
     # Replace the values in the existing column with 1
    column_index = header.index('<5% N')
    accession_index = header.index('accession')
    for row in rows[1:]:
        if row[accession_index] in high_quality_accessions:
            row[column_index] = 1
        else:
            row[column_index] = 0

# Write the updated rows back to the CSV file
with open(csv_file_path, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(rows)
```

### Saving high quality genomes to high_quality_WNV_feb_2024.fna


```python
with open(high_quality_genomes_path, 'w') as output_file:
    high_quality_genomes_fasta_data = ''.join(['>' + sequence for sequence in high_quality_genomes])
    output_file.write(high_quality_genomes_fasta_data)
```

## Polyprotein Annotation and Translation

WNV has a single polyprotein that is cleaved into individual proteins. The polyprotein is ca. 3400 AA long, and we can use this to find the coding sequence and extract it from the genome


```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```


```python
def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((seq_len-frame) // 3) #Multiple of three
            trans = str(nuc[frame:frame+length].translate(trans_table))
            proteins = trans.split("*")
            prev_stop = 0
            for pro in proteins:
                if len(pro) >= min_protein_length:
                    location = pro.find("M")
                    if location != -1:
                        start = frame + (prev_stop + location) * 3
                        pro = pro[location:]
                        end = start + len(pro) * 3
                        answer.append((start, end, len(pro)))
                prev_stop += len(pro) + 1  # +1 to account for the stop codon
    answer.sort(key = lambda x: -x[2]) # Sorted by protein length
    return answer

def find_orfs_in_all_frames(sequence, min_protein_length=3000):
    table = 1 # Standard Genetic Code
    orfs = find_orfs_with_trans(sequence, table, min_protein_length)
    return orfs
```


```python
def WNV_polyprotein_search(input_multiFASTA, min_protein_length=100):
    nt_polyproteins = []
    aa_polyproteins = []
    longest_seq = None

    for record in SeqIO.parse(input_multiFASTA, "fasta"):
        orfs = find_orfs_in_all_frames(record.seq, 3000)
        if orfs:
            longest_orf = orfs[0]  # The first ORF in the list is the longest
            start, end, _ = longest_orf
            nt_seq = record.seq[start:end]
            aa_seq = Seq(nt_seq).translate(to_stop=True)
            nt_polyproteins.append(SeqRecord(nt_seq, id=record.id, description=""))
            aa_polyproteins.append(SeqRecord(aa_seq, id=record.id, description=""))
            if longest_seq is None or len(nt_seq) > len(longest_seq):
                longest_seq = nt_seq
                longest_id = record.id
    return(nt_polyproteins, aa_polyproteins)
```


```python
nt_polyproteins, aa_polyproteins = WNV_polyprotein_search(genomes_path)
high_quality_nt_polyproteins, high_quality_aa_polyproteins = WNV_polyprotein_search(high_quality_genomes_path)

print(f'Of {len(genomes_with_metadata)} unfiltered genomes, {len(nt_polyproteins)} had intact polyprotein sequences (longer than 3k AA)')
print(f'Of {len(high_quality_genomes)} high quality genomes, {len(high_quality_nt_polyproteins)} had intact polyprotein sequences (longer than 3k AA)')
```

    Of 2175 unfiltered genomes, 2143 had intact polyprotein sequences (longer than 3k AA)
    Of 1838 high quality genomes, 1830 had intact polyprotein sequences (longer than 3k AA)


### Saving polyprotein sequences to FASTA files

polyprotein_feb_2024.fna,  
polyprotein_feb_2024.faa,  
high_quality_polyprotein_feb_2024.fna,  
high_quality_polyprotein_feb_2024.faa


```python
SeqIO.write(nt_polyproteins, poly_protein_path, "fasta")
SeqIO.write(aa_polyproteins, poly_protein_translated_path, "fasta")
SeqIO.write(high_quality_nt_polyproteins, high_quality_poly_protein_path, "fasta")
SeqIO.write(high_quality_aa_polyproteins, high_quality_poly_protein_translated_path, "fasta")
```




    1830



fin.
