from functools import partial

from bioutils.normalize import NormalizationMode, normalize

import pytest


sequence = "CCCCCCCCACACACACACTAGCAGCAGCA"

    
normalize_seq = partial(normalize, sequence=sequence)
normalize_trim = partial(normalize_seq, mode=NormalizationMode.TRIMONLY)
normalize_left = partial(normalize_seq, mode=NormalizationMode.LEFTSHUFFLE)
normalize_right = partial(normalize_seq, mode=NormalizationMode.RIGHTSHUFFLE)
normalize_expand = partial(normalize_seq, mode=NormalizationMode.EXPAND)
normalize_vcf = partial(normalize_seq, mode=NormalizationMode.VCF)


def test_trim():
    assert ((25, 25), ('', 'AC')) == normalize_trim(interval=(22,25), alleles=(None, "AGCAC"))
    assert ((24, 25), ('C', '', 'CAC')) == normalize_trim(interval=(22,25), alleles=(None, "AG", "AGCAC"))
    assert ((23, 24), ('G', '', 'GCA')) == normalize_trim(interval=(22,25), alleles=(None, "AC", "AGCAC"))
    assert ((22, 24), ('AG', 'G', 'AGCA')) == normalize_trim(interval=(22,25), alleles=(None, "GC", "AGCAC"))


def test_anchor():
    assert ((23, 25), ('GC', '')) == normalize_trim(interval=(22,25), alleles=(None, "A"), anchor_length=0)
    assert ((22, 26), ('AGCA', 'AA')) == normalize_trim(interval=(22,25), alleles=(None, "A"), anchor_length=1)
    assert ((21, 27), ('CAGCAG', 'CAAG')) == normalize_trim(interval=(22,25), alleles=(None, "A"), anchor_length=2)

    # off the left
    assert ((1, 1), ('', 'C')) == normalize_trim(interval=(1,1), alleles=(None, "C"), anchor_length=0)
    assert ((0, 2), ('CC', 'CCC')) == normalize_trim(interval=(1,1), alleles=(None, "C"), anchor_length=1)
    assert ((0, 3), ('CCC', 'CCCC')) == normalize_trim(interval=(1,1), alleles=(None, "C"), anchor_length=2)

    # off the right
    assert ((28, 28), ('', 'C')) == normalize_trim(interval=(28,28), alleles=(None, "C"), anchor_length=0)
    assert ((27, 29), ('CA', 'CCA')) == normalize_trim(interval=(28,28), alleles=(None, "C"), anchor_length=1)
    assert ((26, 29), ('GCA', 'GCCA')) == normalize_trim(interval=(28,28), alleles=(None, "C"), anchor_length=2)


def test_error_distinct():
    """Must have at least two distinct allele sequences (incl. ref) to normalize"""
    with pytest.raises(ValueError):
        normalize_trim(interval=(22,25), alleles=(None, "AGC"))

def test_error_ref_allele():
    "First allele is ref allele and must be None"
    with pytest.raises(ValueError):
        normalize_trim(interval=(22,25), alleles=("foo", "AGC"))
    

#def test_partial():
#def test_multiallele():
#def test_bounds_honored():
#def test_input_alleles_not_modified
