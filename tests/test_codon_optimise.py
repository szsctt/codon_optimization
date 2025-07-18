import re
import itertools
import pytest
from Bio.Seq import Seq

from codon_optimisation.codon_optimise import (check_consecutive, reverse_translate, check_overlap, 
                       check_all_last_codon, aa_to_nt, check_avoid, change_codon,
                       get_high_freq_codon, check_last_codon)

def test_check_consecutive_1():

    assert check_consecutive([1,2,3,4,5], 0) is False

def test_check_consecutive_2():

    assert check_consecutive([1,2,3,4,5], 1) is True

def test_check_consecutive_3():

    assert check_consecutive([1,2,3,4,5], 2) is False

def test_check_consecutive_4():

    assert check_consecutive([1,1,3,4,5], 2) is True

def test_check_consecutive_5():

    with pytest.raises(AssertionError):
        check_consecutive([1,2,3,4,5], -1)

def test_check_consecutive_6():

    with pytest.raises(AssertionError):
        check_consecutive([1,2,3,4,5], 6)


def test_reverse_translate_1():

    assert reverse_translate('M') == 'ATG'    

def test_reverse_translate_2():

    assert reverse_translate('L') == 'CTG'

def test_reverse_translate_3():

    assert reverse_translate('V') == 'GTG'

def test_reverse_translate_4():

    assert reverse_translate('IG') == 'ATCGGC'

def test_reverse_translate_5():

    peptide = 'MIGVDEQTAGC'
    rev = reverse_translate(peptide)
    assert peptide == str(Seq(rev).translate())

def test_reverse_translate_6():
    
    # can't avoid ATG, only one codon for M
    peptide = 'MAG'
    avoid = [re.compile('ATG')]
    with pytest.raises(ValueError):
        reverse_translate(peptide, avoid=avoid)

def test_reverse_translate_7():

    # change the first codon for I
    peptide = 'MAG'
    avoid = [re.compile('GCC')]
    nucl = reverse_translate(peptide, avoid=avoid)
    assert nucl == 'ATGGCTGGC'
    assert str(Seq(nucl).translate()) == peptide

def test_reverse_translate_8():

    # multiple occurences of avoid codon - can't avoid
    peptide = 'MAG'
    avoid = [re.compile('GGC')]
    with pytest.raises(ValueError):
        reverse_translate(peptide, avoid=avoid)

def test_reverse_translate_9():

    # have to change second and thrid codons
    peptide = 'MAGG'
    avoid = [re.compile('GGCC'), re.compile('CTGGC')]
    nucl = reverse_translate(peptide, avoid=avoid)
    assert nucl == 'ATGGCTGGAGGC'

def test_reverse_translate_10():

    # try a large number of peptides
    peptides = itertools.product('MILVAGDE', repeat=5)
    avoid = [re.compile('GGCC'),re.compile('CTGGC')]
    res = [reverse_translate(pep, avoid=avoid) for pep in peptides]

    # check that all peptides are translated correctly
    assert all(str(Seq(nucl).translate()) == pep for nucl, pep in zip(res, peptides))
    
    # check that all peptides don't contain avoid sequences
    for nucl in res:
        assert all(a not in nucl for a in ('GGCC', 'CTGGC'))

def test_reverse_translate_11():

    # check with avoid sequence with ambiguity
    peptide = 'MAG'
    avoid = [re.compile('GG.CGG')]
    nucl = reverse_translate(peptide, avoid=avoid)
    assert nucl == 'ATGGCTGGC'


def test_check_overlap_1():

    assert check_overlap(1, 2, 3, 4) is False

def test_check_overlap_2():

    assert check_overlap(1, 3, 2, 4) is True

def test_check_overlap_3():

    assert check_overlap(1, 4, 2, 3) is True

def test_check_overlap_4():

    assert check_overlap(3, 4, 1, 2) is False

def test_check_overlap_5():

    assert check_overlap(1, 3, 3, 4) is False

def test_check_overlap_6():

    assert check_overlap(3, 4, 1, 3) is False

def test_check_overlap_7():

    assert check_overlap(1, 3, 1, 3) is True

def test_check_overlap_8():

    assert check_overlap(1, 3, 1, 4) is True

def test_check_overlap_9():

    assert check_overlap(1, 4, 2, 4) is True

def test_check_all_last_codon_1():

    # only one codon for M
    assert check_all_last_codon([('M', 0)], avoid_codons=[0]) is True

def test_check_all_last_codon_2():

    # doesn't make sense if there are no codons to avoid
    with pytest.raises(AssertionError):
        check_all_last_codon([('M', 0)], avoid_codons=[])

def test_check_all_last_codon_3():

    # if avoid codons specifies something outstide of codons range, should raise error
    with pytest.raises(AssertionError):
        check_all_last_codon([('M', 0)], avoid_codons=[1])

def test_check_all_last_codon_4():

    # there are more than two codons for I
    assert check_all_last_codon([('M', 0), ('I', 1)], avoid_codons=[1]) is False

def test_check_all_last_codon_5():

    # there are more than two codons for I
    assert check_all_last_codon([('M', 0), ('I', 1)], avoid_codons=[0, 1]) is False

def test_check_all_last_codon_6():

    # M has only one codon, I has three
    assert check_all_last_codon([('M', 0), ('I', 2)], avoid_codons=[0, 1]) is True

def test_change_codon_1():

    # if avoid_pos is -1, should return same codons
    codons = [('M', 0), ('I', 0)]
    changed = change_codon(codons, -1)
    assert changed == codons

def test_change_codon_2():

    # if avoid_pos isn't a length 2 tuple, should raise error
    codons = [('M', 0), ('I', 0)]
    with pytest.raises(TypeError):
        change_codon(codons, 2)
    
def test_change_codon_3():

    # if we can't change any codons, should raise error
    codons = [('M', 0), ('I', 0)]
    with pytest.raises(ValueError):
        change_codon(codons, (0,3))

def test_change_codon_4():

    # change the first isoleucine because it is the first one that
    # intersects with avoid coords and can be changed
    codons = [('M', 0), ('I', 0)]
    changed = change_codon(codons, (3, 6))
    assert changed == [('M', 0), ('I', 1)]

def test_change_codon_5():

    # change the first isoleucine because it is the first one that 
    # intersects with avoid coords and can be changed
    codons = [('M', 0), ('I', 0), ('I', 0)]
    changed = change_codon(codons, (3, 9))
    assert changed == [('M', 0), ('I', 1), ('I', 0)]

def test_change_codon_6():

    # change the second isoleucine because it has the lowest frequency
    # of the two isoleucines that intersect with avoid coords
    codons = [('M', 0), ('I', 1), ('I', 0)]
    changed = change_codon(codons, (3, 9))
    assert changed == [('M', 0), ('I', 1), ('I', 1)]

def test_change_codon_7():

    # testing with an offset
    codons = [('M', 0), ('I', 1), ('I', 0)]
    changed = change_codon(codons, (4, 10), offset=1)
    assert changed == [('M', 0), ('I', 1), ('I', 1)]

def test_change_codon_8():

    # coords cross boundries of two codons
    codons = [('M', 0), ('I', 1), ('I', 0)]
    changed = change_codon(codons, (0, 5), offset=1)
    assert changed == [('M', 0), ('I', 2), ('I', 0)]

def test_change_codon_9():

    # avoid coords don't intersect with any codons
    codons = [('M', 0), ('I', 1), ('I', 0)]
    with pytest.raises(AssertionError):
        change_codon(codons, (10, 20))



def test_get_high_freq_codon_1():

    # methione can't be changed so has to be isoleucine
    codons = [('M', 0), ('I', 0)]
    assert get_high_freq_codon(codons, [0, 1]) == 1

def test_get_high_freq_codon_2():

    # methione can't be changed, first isoleucine has lower frequency so change second one
    codons = [('M', 0), ('I', 1), ('I', 0)]
    assert get_high_freq_codon(codons, [0, 1, 2]) == 2

def test_get_high_freq_codon_3():
    
    # methione can't be changed, has to be first isoleucine
    codons = [('M', 0), ('I', 1), ('I', 0)]
    assert get_high_freq_codon(codons, [0, 1]) == 1

def test_get_high_freq_codon_4():

    codons = [('M', 0), ('I', 0), ('I', 1)]
    assert get_high_freq_codon(codons, [0, 1, 2]) == 1

def test_get_high_freq_codon_5():

    # change change methionine so raise error
    codons = [('M', 0)]
    with pytest.raises(ValueError):
        get_high_freq_codon(codons, [0])

def test_get_high_freq_codon_6():

    # more codons than possible
    codons = [('I', 5)]
    with pytest.raises(AssertionError):
        get_high_freq_codon(codons, [0])

def test_check_last_codon_1():

    assert check_last_codon(('M', 0)) is True

def test_check_last_codon_2():

    assert check_last_codon(('I', 0)) is False

def test_check_last_codon_3():

    assert check_last_codon(('I', 2)) is True

def test_check_last_codon_4():

    with pytest.raises(AssertionError):
        check_last_codon(('I', 3))

def test_check_last_codon_5():

    with pytest.raises(AssertionError):
        check_last_codon(('I', -1))

def test_check_last_codon_6():

    with pytest.raises(AssertionError):
        check_last_codon(('B', 3))

def test_check_last_codon_7():

    with pytest.raises(AssertionError):
        check_last_codon([('M', 0), ('I', 0)])

def test_aa_to_nt_1():

    assert aa_to_nt([('M', 0)]) == 'ATG'

def test_aa_to_nt_2():

    with pytest.raises(IndexError):
        aa_to_nt([('M', 1)])

def test_aa_to_nt_3():

    assert aa_to_nt([('M', 0,), ('L', 0)]) == 'ATGCTG'

def test_aa_to_nt_4():

    assert aa_to_nt([('M', 0), ('L', 0)], before = 'A', after = 'T') == 'AATGCTGT'

def test_aa_to_nt_5():

    assert aa_to_nt([], before = 'A', after = 'T') == 'AT'

def test_aa_to_nt_6():

    assert aa_to_nt([('M', 0), ('L', 2)], before = 'A', after = 'T') == 'AATGTTGT'

def test_aa_to_nt_7():

    nt = aa_to_nt([('T', 0), ('L', 2), ('A', 0), ('M', 0)])
    assert str(Seq(nt).translate()) == 'TLAM'

def test_aa_to_nt_8():

    nt = aa_to_nt([('T', 0), ('L', 2), ('*', 0), ('M', 0)])
    assert str(Seq(nt).translate()) == 'TL*M'

def test_check_avoid_1():

    avoid = [re.compile('ATG')]
    assert check_avoid('ATGCTG', avoid) == (0, 3)

def test_check_avoid_2():

    avoid = [re.compile('ATG'), re.compile('CTG')]
    assert check_avoid('ATGCTG', avoid) == (0, 3)

def test_check_avoid_3():

    avoid = [re.compile('AA')]
    assert check_avoid('ATGCTG', avoid) == -1

def test_check_avoid_4():

    avoid = [re.compile('AA'), re.compile('ATG')]
    assert check_avoid('ATGCTG', avoid) == (0, 3)

def test_check_avoid_5():

    avoid = [re.compile('AAAAAAAAAAA')]
    assert check_avoid('ATGCTG', avoid) == -1