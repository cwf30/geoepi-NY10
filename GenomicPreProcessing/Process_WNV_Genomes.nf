#!/usr/bin/env nextflow
// Nextflow pipeline for processing genomic sequences: sanitizing IDs, filtering for quality, extracting polyproteins, 
// aligning, masking, partitioning, and calling polymorphisms, with optional sequence subsetting.

params.fasta_path = 'data/query_genomes/**.{fa,fasta,fna}'         // Path to query genomes
params.reference_genome = 'data/reference_genome/NY99.faa'         // Reference genome
params.mut = 'NO_FILE'                                             // Optional mutations file, can be used to subset genomes that contain any of the mutations
params.seqIDs = 'NO_FILE'                                          // Optional sequence IDs file for subsetting by sequence ID

// Process: Sanitize Sequence IDs
process sanitize_ID {
    publishDir 'output/', mode: 'copy' // Output sanitized sequences to 'output/' directory

    input:
    path fasta_files // Input FASTA files

    output:
    path 'sanitized_sequences.fasta' // Output sanitized sequences in FASTA format

    script:
    """
    sanitize_seqID.py $fasta_files
    """
}

// Process: Filter sequences based on quality
process quality_filter {
    publishDir 'output/', mode: 'copy' // Output filtered sequences to 'output/' directory

    input:
    path fasta_files // Input sanitized FASTA files

    output:
    path 'quality_sequences.fasta' // Output high-quality sequences in FASTA format

    script:
    """
    quality_filter.py $fasta_files
    """
}

// Process: Extract polyproteins from sequences
process extract_polyprotein {
    publishDir 'output/', mode: 'copy' // Output extracted polyproteins to 'output/' directory

    input:
    path high_quality_reads // Input high-quality sequences

    output:
    path 'polyproteins.*' // Output nucleotide and amino acid polyproteins

    script:
    """
    extract_polyproteins.py $high_quality_reads
    """
}

// Process: Align sequences using MAFFT
process align {
    input:
    path polyprotein_multifasta // Input polyprotein multi-FASTA file

    output:
    path "alignment_${polyprotein_multifasta}" // Output alignment file

    script:
    """
    mafft --auto --maxiterate 1000 ${polyprotein_multifasta} > alignment_${polyprotein_multifasta}
    """
}

// Process: Mask gaps in alignments using clipKit
process mask_gaps {
    publishDir 'output/', mode: 'copy' // Output masked alignments to 'output/' directory

    input:
    path alignment // Input alignment file

    output:
    path "masked_${alignment}" // Output masked alignment file

    script:
    """
    clipKit ${alignment} -m gappy -o masked_${alignment}
    """
}

// Process: Partition genomes into gene regions
process find_partitions {
    publishDir 'output/partitions', mode: 'copy' // Output partitioned genomes to 'output/partitions/'

    input:
    path alignment // Input masked alignment

    output:
    path "*partition*" // Output gene partitions in separate files

    script:
    """
    partition_genomes.py ${alignment}
    """
}

// Process: Call polymorphisms from the alignment
process polymorphism_caller {
    publishDir 'output/', mode: 'copy' // Output polymorphism data to 'output/' directory

    input:
    path polyprotein_alignment // Input polyprotein alignment
    path reference_genome      // Input reference genome
    path mut                   // Optional mutations file

    output:
    path "*polymorphisms*.csv" // Output polymorphism data as CSV

    script:
    """
    call_polymorphisms.py ${polyprotein_alignment} ${reference_genome} ${mut}
    """
}

// Process: Subset sequences from a FASTA file
process subset_sequences_from_fasta {
    publishDir 'output/subset', mode: 'copy' // Output subsets to 'output/subset/'

    when:
    params.seqIDs != 'NO_FILE' // Only execute if sequence IDs are provided

    input:
    tuple path(fasta_file), path(seqIDs) // Input FASTA file and sequence IDs

    output:
    path "subset_${fasta_file}" // Output subset FASTA file

    script:
    """
    subset_sequences_from_fasta.py $fasta_file $seqIDs
    """
}

// Process: Subset sequences from a CSV file
process subset_sequences_from_CSV {
    publishDir 'output/subset', mode: 'copy' // Output subsets to 'output/subset/'

    when:
    params.seqIDs != 'NO_FILE' // Only execute if sequence IDs are provided

    input:
    tuple path(CSV_file), path(seqIDs) // Input CSV file and sequence IDs

    output:
    path "subset_${CSV_file}" // Output subset CSV file

    script:
    """
    subset_sequences_from_CSV.py $CSV_file $seqIDs
    """
}

// Workflow definition
workflow {
    // Initialize input channels
    Channel
        .fromPath(params.fasta_path)
        .collectFile(name: 'fasta_files', newLine: true)
        .set { fasta_files }

    Channel
        .fromPath(params.reference_genome)
        .set { reference_genome }

    Channel
        .fromPath(params.seqIDs)
        .set { seqIDs }

    Channel
        .fromPath(params.mut)
        .set { mut }

    // Sanitize sequence IDs
    sanitize_ID(fasta_files).collectFile(name: 'sanitized_fasta_files', newLine: true)
        .set { sanitized_fasta_files }

    // Filter sequences for quality
    quality_filter(sanitized_fasta_files)

    // Extract polyproteins
    extract_polyprotein(quality_filter.out)

    extract_polyprotein.out.flatten().set { polyprotein_multifasta }

    // Align sequences
    align(polyprotein_multifasta)

    // Mask gaps in alignments
    mask_gaps(align.out)

    // Partition genomes into gene regions
    find_partitions(mask_gaps.out)

    // Branch alignments by file type (nucleotide or amino acid)
    mask_gaps.out.flatten().branch {
        aa: it.name.endsWith('.faa')
        nt: it.name.endsWith('.fna')
    }.set { branched_alignments }

    // Call polymorphisms
    polymorphism_caller(branched_alignments.aa, reference_genome, mut)
        .flatten()
        .combine(seqIDs)
        .set { polymorphism_caller_out }

    // Collect outputs for subsetting
    mask_gaps.out.flatten().set { mask_gaps_out_collected }
    find_partitions.out.flatten().set { find_partitions_out_collected }
    polyprotein_multifasta.concat(mask_gaps_out_collected, find_partitions_out_collected)
        .combine(seqIDs)
        .set { subset_fasta_inputs }

    // Subset sequences if desired
    subset_sequences_from_fasta(subset_fasta_inputs)
    subset_sequences_from_CSV(polymorphism_caller_out)
}
