#!/usr/bin/env python
# Script to identify and extract nucleotide and amino acid polyproteins from a multi-FASTA input
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Function to find open reading frames (ORFs) and translate them
def find_orfs(sequence, trans_table=1, min_protein_length=3000):
    """
    Identifies ORFs in a nucleotide sequence, translates them, and filters by minimum protein length.

    Parameters:
        sequence (Seq): Input nucleotide sequence.
        trans_table (int): Translation table (default: 1 for the standard genetic code).
        min_protein_length (int): Minimum protein length to consider (default: 3000).

    Returns:
        list of tuples: Each tuple contains (start, end, protein_length), sorted by protein length.
    """
    orfs = []
    seq_len = len(sequence)
    # Check both forward (+1) and reverse complement (-1) strands
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):  # Check all three reading frames
            length = 3 * ((seq_len - frame) // 3)  # Ensure length is a multiple of 3
            trans = str(nuc[frame:frame + length].translate(trans_table))  # Translate to protein sequence
            proteins = trans.split("*")  # Split at stop codons
            prev_stop = 0
            for pro in proteins:
                if len(pro) >= min_protein_length:  # Filter by protein length
                    location = pro.find("M")  # Look for the first start codon (Methionine)
                    if location != -1:
                        start = frame + (prev_stop + location) * 3  # Calculate nucleotide start
                        pro = pro[location:]  # Trim protein to start from the first Methionine
                        end = start + len(pro) * 3  # Calculate nucleotide end
                        orfs.append((start, end, len(pro)))
                prev_stop += len(pro) + 1  # Increment to account for the stop codon
    # Sort ORFs by protein length in descending order
    orfs.sort(key=lambda x: -x[2])
    return orfs


# Function to search for WNV polyproteins in a multi-FASTA input
def WNV_polyprotein_search(input_multiFASTA):
    """
    Extracts the longest nucleotide and amino acid polyproteins from each sequence in a multi-FASTA file.

    Parameters:
        input_multiFASTA (str): Path to the input multi-FASTA file.

    Returns:
        tuple: Two lists containing SeqRecords for nucleotide and amino acid polyproteins.
    """
    nt_polyproteins = []  # List to store nucleotide polyproteins
    aa_polyproteins = []  # List to store amino acid polyproteins
    longest_seq = None  # Tracks the longest nucleotide sequence

    for record in SeqIO.parse(input_multiFASTA, "fasta"):
        # Identify ORFs in the sequence
        orfs = find_orfs(record.seq, min_protein_length=3000)
        if orfs:
            longest_orf = orfs[0]  # Select the longest ORF
            start, end, _ = longest_orf
            nt_seq = record.seq[start:end]  # Extract the nucleotide sequence of the ORF
            aa_seq = Seq(nt_seq).translate(to_stop=True)  # Translate to amino acid sequence
            # Store the sequences as SeqRecord objects
            nt_polyproteins.append(SeqRecord(nt_seq, id=record.id, description=""))
            aa_polyproteins.append(SeqRecord(aa_seq, id=record.id, description=""))
            # Update the longest sequence if applicable
            if longest_seq is None or len(nt_seq) > len(longest_seq):
                longest_seq = nt_seq
                longest_id = record.id
    return nt_polyproteins, aa_polyproteins

# Process the input multi-FASTA file and save the polyproteins to output files
nt_polyproteins, aa_polyproteins = WNV_polyprotein_search(sys.argv[1])
SeqIO.write(nt_polyproteins, 'polyproteins.fna', "fasta")  # Write nucleotide polyproteins to a FASTA file
SeqIO.write(aa_polyproteins, 'polyproteins.faa', "fasta")  # Write amino acid polyproteins to a FASTA file
