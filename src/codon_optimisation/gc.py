
from Bio import SeqIO


def check_gc_content(seq: str) -> float:
    """Calculate GC content of a sequence.
    GC content is the percentage of nucleotides in the sequence that are either G or C.
    """
    seq_record = SeqIO.read(seq, "fasta")
    gc_content = (seq_record.seq.count('G') + seq_record.seq.count('C')) / len(seq_record.seq) * 100
    return gc_content