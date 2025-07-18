from Bio import SeqIO


def check_optimised_sequences(reference_file, optimised_files):
    """
    Check optimised sequences against reference.
    Compares one or more optimised sequence files against a reference file.
    """
    for optimised_file in optimised_files:
        # Read the reference sequence
        reference_seq = SeqIO.read(reference_file, "fasta")
        
        # Read the optimised sequence
        optimised_seq = SeqIO.read(optimised_file, "fasta")
        
        # Compare translation of sequences
        if reference_seq.seq.translate() != optimised_seq.seq.translate():
            raise ValueError(f"Optimised sequence {optimised_file} does not match reference sequence {reference_file} after translation.")

    print(f"All optimised sequences match the reference sequence {reference_file} after translation.")