from functools import partial

import pytest

from bioutils.normalize import NormalizationMode, normalize

sequence = "CCCCCCCCACACACACACTAGCAGCAGCA"

normalize_seq = partial(normalize, sequence=sequence)
normalize_trim = partial(normalize_seq, mode=NormalizationMode.TRIMONLY)
normalize_trim_no_shuffle = partial(normalize_seq, mode=None, trim=True)
normalize_no_trim_no_shuffle = partial(normalize_seq, mode=None, trim=False)
normalize_left = partial(normalize_seq, mode=NormalizationMode.LEFTSHUFFLE)
normalize_right = partial(normalize_seq, mode=NormalizationMode.RIGHTSHUFFLE)
normalize_expand = partial(normalize_seq, mode=NormalizationMode.EXPAND)
normalize_vcf = partial(normalize_seq, mode=NormalizationMode.VCF)
normalize_left_no_trim = partial(normalize_seq, mode=NormalizationMode.LEFTSHUFFLE, trim=False)
normalize_right_no_trim = partial(normalize_seq, mode=NormalizationMode.RIGHTSHUFFLE, trim=False)
normalize_expand_no_trim = partial(normalize_seq, mode=NormalizationMode.EXPAND, trim=False)
normalize_vcf_no_trim = partial(normalize_seq, mode=NormalizationMode.VCF, trim=False)


@pytest.mark.parametrize("normalize_fn", [normalize_trim, normalize_trim_no_shuffle])
def test_trim(normalize_fn):
    """Should trim common prefix and suffix when trim=True."""
    assert ((25, 25), ("", "AC")) == normalize_fn(interval=(22, 25), alleles=(None, "AGCAC"))
    assert ((24, 25), ("C", "", "CAC")) == normalize_fn(interval=(22, 25), alleles=(None, "AG", "AGCAC"))
    assert ((23, 24), ("G", "", "GCA")) == normalize_fn(interval=(22, 25), alleles=(None, "AC", "AGCAC"))
    assert ((22, 24), ("AG", "G", "AGCA")) == normalize_fn(interval=(22, 25), alleles=(None, "GC", "AGCAC"))


@pytest.mark.parametrize("normalize_fn", [normalize_trim, normalize_trim_no_shuffle])
def test_anchor(normalize_fn):
    assert ((23, 25), ("GC", "")) == normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=0)
    assert ((22, 26), ("AGCA", "AA")) == normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=1)
    assert ((21, 27), ("CAGCAG", "CAAG")) == normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=2)

    # off the left
    assert ((1, 1), ("", "C")) == normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=0)
    assert ((0, 2), ("CC", "CCC")) == normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=1)
    assert ((0, 3), ("CCC", "CCCC")) == normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=2)

    # off the right
    assert ((28, 28), ("", "C")) == normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=0)
    assert ((27, 29), ("CA", "CCA")) == normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=1)
    assert ((26, 29), ("GCA", "GCCA")) == normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=2)


def test_trinuc():
    """LEFTSHUFFLE, RIGHTSHUFFLE, EXPAND normalization for trinucleotide"""
    # 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
    #  C C C C C C C C A C A C A C A C A C T A G C A G C A G C A T
    #                                             ^ [22,22): ['', 'AGC']   | Starting alleles
    #      LEFTSHUFFLE                      ^ [19,19): ['', 'AGC']
    #      RIGHTSHUFFLE                                         ^ [29,29): ['', 'GCA']
    #      EXPAND                           ^-------------------^ [19,29): ['AGCAGCAGCA', 'AGCAGCAGCAGCA']
    assert ((19, 19), ("", "AGC")) == normalize_left(interval=(22, 22), alleles=(None, "AGC"))
    assert ((29, 29), ("", "GCA")) == normalize_right(interval=(22, 22), alleles=(None, "AGC"))
    assert ((19, 29), ("AGCAGCAGCA", "AGCAGCAGCAGCA")) == normalize_expand(interval=(22, 22), alleles=(None, "AGC"))


def test_bounds():
    """ensure that bounds are honored"""
    assert ((20, 24), ("GCAG", "GCAGCAG")) == normalize_expand(
        interval=(22, 22), alleles=(None, "AGC"), bounds=(20, 24)
    )


def test_no_trim_no_shuffle():
    """Should not trim or shuffle when mode=None, trim=False."""
    assert ((22, 25), ("AGC", "AGC")) == normalize_no_trim_no_shuffle(interval=(22, 25), alleles=(None, "AGC"))
    assert ((22, 25), ("AGC", "AGCT")) == normalize_no_trim_no_shuffle(interval=(22, 25), alleles=(None, "AGCT"))


def test_shuffle_no_trim():
    """Should shuffle but not trim when mode!=None and trim=False."""
    assert ((19, 22), ("AGC", "AGC")) == normalize_left_no_trim(interval=(22, 25), alleles=(None, "AGC"))
    assert ((26, 29), ("GCA", "GCA")) == normalize_right_no_trim(interval=(22, 25), alleles=(None, "AGC"))
    assert ((19, 29), ("AGCAGCAGCA", "AGCAGCAGCA")) == normalize_expand_no_trim(
        interval=(22, 25), alleles=(None, "AGC")
    )


# TODO: def test_multiallele():


def test_mode_string():
    "test that mode as string is accepted"
    _normalize = partial(normalize_seq, interval=(28, 28), alleles=(None, "C"))
    vcf_out = ((26, 27), ("G", "GC"))
    assert vcf_out != _normalize(), "not VCF output by default"
    assert vcf_out == _normalize(mode="VCF"), "mode as string recognized"


def test_input_alleles_not_modified():
    """ensure that alleles list is not modified"""
    alleles = (None, "AGCAC")
    normalize_trim(interval=(22, 25), alleles=alleles)
    assert (None, "AGCAC") == alleles


@pytest.mark.parametrize("normalize_fn", [normalize_trim, normalize_trim_no_shuffle])
def test_error_distinct(normalize_fn):
    """Must have at least two distinct allele sequences (incl. ref) to normalize"""
    with pytest.raises(ValueError):
        normalize_fn(interval=(22, 25), alleles=(None, "AGC"))


def test_error_ref_allele():
    "First allele is ref allele and must be None"
    with pytest.raises(ValueError):
        normalize_trim(interval=(22, 25), alleles=("foo", "AGC"))


def test_error_vcf_mode_no_trim():
    """Should raise error when mode=VCF, trim=False."""
    with pytest.raises(ValueError) as exc_info:
        normalize_vcf_no_trim(interval=(22, 25), alleles=(None, "AGC"))
    assert str(exc_info.value) == "May not disable trimming with VCF normalization mode"
