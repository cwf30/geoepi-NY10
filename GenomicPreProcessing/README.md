# West Nile Virus Genome Processing Pipeline Documentation

## Overview

This pipeline processed **West Nile Virus (WNV) genomes** for downstream analysis in Fautt et al. 2024. It sanitizes sequence IDs, filters sequences based on quality, extracts polyproteins, aligns sequences, partitions genomes into gene regions, calls non-synonymous polymorphisms, and optionally subsets output files by sequence IDs or specific mutations of interest.

------------------------------------------------------------------------

## File Structure

The pipeline expects the following directory structure:

```         
project/
├── Process_WNV_Genomes.nf      # The Nextflow pipeline file
|
├── bin/                        # Directory containing all scripts necessary for the pipeline
|
├── data/
│   ├── query_genomes/          # Directory containing input genome FASTA files
│   │   ├── genome1.fna
│   │   ├── genome2.fasta
│   │   └── ...
│   ├── reference_genome/       # Directory containing the reference genome
│   │   └── NY99.faa
│   ├── seqIDs.txt              # (Optional) File containing sequence IDs for subsetting
│   │   
│   └── mutations.txt           # (Optional) File containing mutations of interest
│       
└── output/                     # Directory where pipeline outputs are stored
```

### Input Files

1.  **Query Genomes** (Required):
    -   Location: `data/query_genomes/**.{fa,fasta,fna}`
    -   Format: Nucleotide FASTA files.
    -   Description: The genomes to be processed by the pipeline.
2.  **Reference Genome** (Required):
    -   Location: `data/reference_genome/NY99.faa`
    -   Format: Amino acid FASTA file.
    -   Description: Used for polymorphism calling.
3.  **Mutation File** (Optional):
    -   Location: `data/mutations.txt`
    -   Format: Plain text file, one non-synonymous mutation per line (e.g. NS2A-R188K).
    -   Description: Used for outputting polymorphism and alignment files for only subsetted sequences.
4.  **Sequence ID File** (Optional):
    -   Location: `data/seqIDs.txt`
    -   Format: Plain text file, one sequence ID per line. Expects ID to be simplified sequence name (i.e. 'MH170276' but not 'MH170276.1')
    -   Description: Used for outputting polymorphism and alignment files for only subsetted sequences.

------------------------------------------------------------------------

## Inputs and Options

1.  **Inputs**:
    -   `--fasta_path`: Path to query genome files. Default: `data/query_genomes/**.{fa,fasta,fna}`
    -   `--reference_genome`: Path to the reference AA polyprotein sequence for polymorphism calling. Default: `data/reference_genome/NY99.faa`
    -   `--mut`: Path to the mutations file (optional). Default: `NO_FILE`
    -   `--seqIDs`: Path to the sequence IDs file for subsetting (optional). Default: `NO_FILE`
2.  **Default Setup**: If the file structure matches the expected layout and filenames (e.g., `NY99.faa` as the reference genome), the pipeline runs without specifying any parameters.

------------------------------------------------------------------------

## Outputs

The pipeline generates the following outputs:

1.  **Sanitized Sequences**:
    -   File: `output/sanitized_sequences.fasta`
    -   Description: Query genomes in a single FASTA file with sanitized sequence IDs (sequence names stripped of everything after the first period or underscore.)
2.  **Quality-Filtered Sequences**:
    -   File: `output/quality_sequences.fasta`
    -   Description: Sequences that are longer than 10,200bp and contain less than 5% 'N'
3.  **Extracted Polyproteins**:
    -   Files:
        -   `output/polyproteins.fna`: Nucleotide sequences.
        -   `output/polyproteins.faa`: Amino acid sequences.
    -   Description: Extracted polyprotein regions from the genomes.
4.  **Alignments**:
    -   File: `output/alignment_<polyprotein_file>`
    -   Description: Alignments of polyprotein sequences using MAFFT.
5.  **Masked Alignments**:
    -   File: `output/masked_<alignment_file>`
    -   Description: Gaps in alignments masked using ClipKit.
6.  **Gene Partitions**:
    -   Directory: `output/partitions/`
    -   Files: Gene-specific partitions, e.g., `C_partition.fna`, `E_partition.fna`.
7.  **Polymorphism Data**:
    -   Files:
        -   `output/polymorphisms_sparse.csv`: Sparse representation of polymorphisms.
        -   `output/polymorphisms_dense.csv`: Dense representation of polymorphisms.
8.  **Subset Data** (Optional):
    -   Files:
        -   `output/subset/subset_<fasta_file>`: Subset of FASTA sequences.
        -   `output/subset/subset_<CSV_file>`: Subset of polymorphism data.

------------------------------------------------------------------------

## Running the Pipeline

### Default Command

If the input file structure matches the default layout, you can run the pipeline without additional parameters:

``` bash
nextflow run Process_WNV_Genomes.nf
```

### Custom Parameters

Use the following parameters to customize inputs:

``` bash
nextflow run Process_WNV_Genomes.nf \
  --fasta_path 'custom_path/query_genomes/**.fasta' \
  --reference_genome 'custom_path/reference_genome/custom_ref.faa' \
  --mut 'custom_path/mutations/mutations.txt' \
  --seqIDs 'custom_path/sequence_IDs/custom_seqIDs.txt'
```
