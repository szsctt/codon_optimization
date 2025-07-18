# shared functions

import re
import sys

from Bio import SeqIO

import typer
from typing_extensions import Annotated

def codon_optimise(
        seq: Annotated[str, typer.Argument(help="Fasta file containing sequence to codon optimse")],
        type: Annotated[str, typer.Option(help="Type of sequence: 'dna' or 'protein'", default='dna')],
        before: Annotated[str, typer.Option(help="Sequence to add before the codon optimised sequence", default='')],
        after: Annotated[str, typer.Option(help="Sequence to add after the codon optimised sequence", default='')],
        avoid: Annotated[list[str], typer.Option(
            help="List of sequences to avoid in the codon optimised sequence. "
                 "Use 'N' for any nucleotide, e.g. 'GGCC' or 'GG.NN'. "
                 "If multiple sequences are provided, they will be avoided in the order "
                 "they are provided. If a sequence is still present after all codons "
                 "have been changed, an error will be raised. ", default=None)],
        output_file: Annotated[str, typer.Option(help="Output file to write the codon optimised sequence to", default='codon_optimised.fasta')]
        

):
    """
    Reverse translate a protein sequence to nucleotides, using the most 
    frequent codon for each amino acid.
    Avoid using any sequences in the avoid list.  If they are present,
    change the codons involved in the avoid sequence, starting with
    the first codon, then the second, etc.  If we reach the end and the avoid
    sequence is still present, go back to the first codon and change it to the 
    next most frequent codon, then the second, etc.  If we reach the end and
    the avoid sequence is still present, raise an error
    """
    
    # read fasta file with biopython into dict - assume fasta format
    seq_dict = SeqIO.to_dict(SeqIO.parse(seq, 'fasta'))
    if len(seq_dict) != 1:
        raise ValueError("Fasta file must contain exactly one sequence")
    name = next(iter(seq_dict.keys()))
    seq = next(iter(seq_dict.values())).seq

    if type == 'dna':
        # for DNA, we first translate to protein, then reverse translate
        seq = seq.upper()
        assert all(i in 'ACGTN' for i in seq)
        # translate to protein
        seq = str(seq.translate())
        opt = reverse_translate(seq, before=before, after=after, avoid=process_avoid_seqs(avoid))
    
    elif type == 'protein':
        seq = seq.upper()
        assert all(i in 'ACDEFGHIKLMNPQRSTVWY*' for i in seq) # this will fail for seuqnces that don't contain all of these amino acids
        opt = reverse_translate(str(seq.seq), before=before, after=after, avoid=process_avoid_seqs(avoid))
    else:
        raise ValueError(f"Unknown type: {type}")
    
    # write to file
    if output_file is None or output_file == '-':
        out = sys.stdout
    else:
        out = open(output_file, 'w')
    with out as f:
        f.write(f'>{name}_codon_optimised\n{opt}\n')
    # close the file if it was opened
    if out is not sys.stdout:
        out.close()



# https://www.genscript.com/tools/codon-frequency-table
codon_table = {'F': ['TTC', 'TTT'],
               'L': ['CTG', 'CTC', 'TTG', 'CTT', 'TTA', 'CTA'],
               'Y': ['TAC', 'TAT'],
               '*': ['TGA', 'TAA', 'TAG'],
               'H': ['CAC', 'CAT'],
               'Q': ['CAG', 'CAA'],
               'I': ['ATC', 'ATT', 'ATA'],
               'M': ['ATG'],
               'N': ['AAC', 'AAT'],
               'K': ['AAG', 'AAA'],
               'V': ['GTG', 'GTC', 'GTT', 'GTA'],
               'D': ['GAC', 'GAT'],
               'E': ['GAG', 'GAA'],
               'S': ['AGC', 'TCC', 'TCT', 'TCA', 'AGT', 'TCA'],
               'C': ['TGC', 'TGT'],
               'W': ['TGG'],
               'P': ['CCC', 'CCT', 'CCA', 'CCG'],
               'R': ['CGG', 'AGA', 'AGG', 'CGC', 'CGA', 'CGT'],
               'T': ['ACC', 'ACA', 'ACT', 'ACG'],
               'A': ['GCC', 'GCT', 'GCA', 'GCG'],
               'G': ['GGC', 'GGA', 'GGG', 'GGT']
               }

def process_avoid_seqs(avoid):

    if avoid is None:
        return []
    avoid = [i.replace("N", '.') for i in avoid]
    avoid = [re.compile(i) for i in avoid]
    return avoid

def reverse_translate(prot_seq, before='', after='', avoid=None):
    '''
    Reverse translate a protein sequence to nucleotides, using the most 
    frequent codon for each amino acid.
    Avoid using any sequences in the avoid list.  If they are present,
    change the codons involved in the avoid sequence, starting with
    the first codon, then the second, etc.  If we reach the end and the avoid
    sequence is still present, go back to the first codon and change it to the 
    next most frequent codon, then the second, etc.  If we reach the end and
    the avoid sequence is still present, raise an error
    '''
    assert all(len(i) == 1 for i in prot_seq)

    # initially try with most frequent codon for each amino acid
    codon_usage = [(aa, 0) for aa in prot_seq]
    dna_seq = aa_to_nt(codon_usage, before, after)
    
    offset = len(before)
    
    # check that we don't have any sequences from the avoid list
    while (avoid_pos := check_avoid(dna_seq, avoid)) != -1:
        # change one codon to a lower frequency one
        codon_usage = change_codon(codon_usage, avoid_pos, offset = offset)
        # check for avoid sequence again
        dna_seq = aa_to_nt(codon_usage, before, after)

    return dna_seq
    
def change_codon(codon_usage, avoid_pos, offset = 0):
    """
    Change one codon to a lower frequency codon
    Change the codon with the highest frequency that can be changed
    A codon can't be changed if it is already the lowest frequency codon
    """
    # if nothing to avoid, don't change anything
    if avoid_pos == -1:
        return codon_usage

    # get start and end of avoid sequence
    avoid_start, avoid_end = avoid_pos

    # find the codons that intersect with the avoid sequence
    codon_starts = [offset + i*3 for i in range(len(codon_usage))]
    avoid_codons = [i for i in range(len(codon_starts)) 
                    if check_overlap(codon_starts[i], codon_starts[i]+3, avoid_start, avoid_end)]

    # if we can't change any codons, raise an error
    if check_all_last_codon(codon_usage, avoid_codons):
        raise ValueError('Cannot avoid sequence')

    # get codon to change
    idx_c = get_high_freq_codon(codon_usage, avoid_codons)

    # change codon
    codon_usage[idx_c] = (codon_usage[idx_c][0], codon_usage[idx_c][1] + 1)

    return codon_usage

def get_high_freq_codon(codon_usage, avoid_codons):
    """
    Return the index of the highest frequency codon that can be changed
    From the indices in avoid_codons
    """

    # remove index of any codons that can't be changed
    avoid_codons = [i for i in avoid_codons if not check_last_codon(codon_usage[i])]

    if len(avoid_codons) == 0:
        raise ValueError('Cannot avoid sequence')

    # get the codon with the highest frequency (lowest index)
    codon_to_change = min(avoid_codons, key = lambda x: codon_usage[x][1])

    return codon_to_change
    
def check_last_codon(codon_usage):

    assert len(codon_usage) == 2
    assert codon_usage[0] in codon_table

    # index can't be higher than the number of codons for that amino acid
    assert codon_usage[1] < len(codon_table[codon_usage[0]])
    assert codon_usage[1] >= 0

    # check if we are using the last available codon for one amino acid
    if codon_usage[1] == len(codon_table[codon_usage[0]]) - 1:
        return True
    return False
        
def check_all_last_codon(codon_usage, avoid_codons):
    """
    Check if codon_usage specificeies the last available codon for each
    amino acid in avoid_codons
    e.g. for codon_usage = [(M, 0)] and avoid_codons = [0], we are using the 
    last available codon for M (since this is only one), so return True
    """
    assert len(avoid_codons) > 0
    assert all(i in range(len(codon_usage)) for i in avoid_codons)
    return all(check_last_codon(codon_usage[i]) for i in avoid_codons)

def check_overlap(a_start, a_stop, b_start, b_stop):
    """
    Check for overlap between two intervals
    Intervals are half-open: [start, stop)
    so check_overlap(0, 10, 10, 15) is False
    """
    if a_start >= b_stop:
        return False
    if b_start >= a_stop:
        return False
    return True
    
def aa_to_nt(codon_usage, before='', after=''):

    nt = ''
    for aa, codon in codon_usage:
        nt += codon_table[aa][codon]

    return before + nt + after

def check_avoid(dna_seq, avoid):
    if avoid is None:
        return -1
    for a in avoid:
        if (match := a.search(dna_seq)) is not None:
            return match.start(), match.end()
    return -1

def check_consecutive(t, n):

    if n == 0:
        return False

    assert n > 0
    assert len(t) >= n
    # check that there are no more than n consecutive amino acids from the same group
    for i in range(len(t)-n+1):
        slice_set = set(t[i:i+n])
        if len(slice_set) == 1:
            return True
    return False


if __name__ == '__main__':
    main()