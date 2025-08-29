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
    assert normalize_fn(interval=(22, 25), alleles=(None, "AGCAC")) == ((25, 25), ("", "AC"))
    assert normalize_fn(interval=(22, 25), alleles=(None, "AG", "AGCAC")) == (
        (24, 25),
        ("C", "", "CAC"),
    )
    assert normalize_fn(interval=(22, 25), alleles=(None, "AC", "AGCAC")) == (
        (23, 24),
        ("G", "", "GCA"),
    )
    assert normalize_fn(interval=(22, 25), alleles=(None, "GC", "AGCAC")) == (
        (22, 24),
        ("AG", "G", "AGCA"),
    )


@pytest.mark.parametrize("normalize_fn", [normalize_trim, normalize_trim_no_shuffle])
def test_anchor(normalize_fn):
    assert normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=0) == (
        (23, 25),
        ("GC", ""),
    )
    assert normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=1) == (
        (22, 26),
        ("AGCA", "AA"),
    )
    assert normalize_fn(interval=(22, 25), alleles=(None, "A"), anchor_length=2) == (
        (21, 27),
        ("CAGCAG", "CAAG"),
    )

    # off the left
    assert normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=0) == (
        (1, 1),
        ("", "C"),
    )
    assert normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=1) == (
        (0, 2),
        ("CC", "CCC"),
    )
    assert normalize_fn(interval=(1, 1), alleles=(None, "C"), anchor_length=2) == (
        (0, 3),
        ("CCC", "CCCC"),
    )

    # off the right
    assert normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=0) == (
        (28, 28),
        ("", "C"),
    )
    assert normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=1) == (
        (27, 29),
        ("CA", "CCA"),
    )
    assert normalize_fn(interval=(28, 28), alleles=(None, "C"), anchor_length=2) == (
        (26, 29),
        ("GCA", "GCCA"),
    )


def test_trinuc():
    """LEFTSHUFFLE, RIGHTSHUFFLE, EXPAND normalization for trinucleotide"""
    # 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
    #  C C C C C C C C A C A C A C A C A C T A G C A G C A G C A T
    #                                             ^ [22,22): ['', 'AGC']   | Starting alleles
    #      LEFTSHUFFLE                      ^ [19,19): ['', 'AGC']
    #      RIGHTSHUFFLE                                         ^ [29,29): ['', 'GCA']
    #      EXPAND                           ^-------------------^ [19,29): ['AGCAGCAGCA', 'AGCAGCAGCAGCA']
    assert normalize_left(interval=(22, 22), alleles=(None, "AGC")) == ((19, 19), ("", "AGC"))
    assert normalize_right(interval=(22, 22), alleles=(None, "AGC")) == ((29, 29), ("", "GCA"))
    assert normalize_expand(interval=(22, 22), alleles=(None, "AGC")) == (
        (19, 29),
        ("AGCAGCAGCA", "AGCAGCAGCAGCA"),
    )


def test_bounds():
    """Ensure that bounds are honored"""
    assert normalize_expand(interval=(22, 22), alleles=(None, "AGC"), bounds=(20, 24)) == (
        (20, 24),
        ("GCAG", "GCAGCAG"),
    )


def test_no_trim_no_shuffle():
    """Should not trim or shuffle when mode=None, trim=False."""
    assert normalize_no_trim_no_shuffle(interval=(22, 25), alleles=(None, "AGC")) == (
        (22, 25),
        ("AGC", "AGC"),
    )
    assert normalize_no_trim_no_shuffle(interval=(22, 25), alleles=(None, "AGCT")) == (
        (22, 25),
        ("AGC", "AGCT"),
    )


def test_shuffle_no_trim():
    """Should shuffle but not trim when mode!=None and trim=False."""
    assert normalize_left_no_trim(interval=(22, 25), alleles=(None, "AGC")) == (
        (19, 22),
        ("AGC", "AGC"),
    )
    assert normalize_right_no_trim(interval=(22, 25), alleles=(None, "AGC")) == (
        (26, 29),
        ("GCA", "GCA"),
    )
    assert normalize_expand_no_trim(interval=(22, 25), alleles=(None, "AGC")) == (
        (19, 29),
        ("AGCAGCAGCA", "AGCAGCAGCA"),
    )


# TODO: def test_multiallele():


def test_mode_string():
    "Test that mode as string is accepted"
    _normalize = partial(normalize_seq, interval=(28, 28), alleles=(None, "C"))
    vcf_out = ((26, 27), ("G", "GC"))
    assert vcf_out != _normalize(), "not VCF output by default"
    assert vcf_out == _normalize(mode="VCF"), "mode as string recognized"


def test_input_alleles_not_modified():
    """Ensure that alleles list is not modified"""
    alleles = (None, "AGCAC")
    normalize_trim(interval=(22, 25), alleles=alleles)
    assert alleles == (None, "AGCAC")


@pytest.mark.parametrize("normalize_fn", [normalize_trim, normalize_trim_no_shuffle])
def test_error_distinct(normalize_fn):
    """Must have at least two distinct allele sequences (incl. ref) to normalize"""
    with pytest.raises(ValueError):  # noqa: PT011
        normalize_fn(interval=(22, 25), alleles=(None, "AGC"))


def test_error_ref_allele():
    "First allele is ref allele and must be None"
    with pytest.raises(ValueError):  # noqa: PT011
        normalize_trim(interval=(22, 25), alleles=("foo", "AGC"))


def test_error_vcf_mode_no_trim():
    """Should raise error when mode=VCF, trim=False."""
    with pytest.raises(ValueError) as exc_info:  # noqa: PT011
        normalize_vcf_no_trim(interval=(22, 25), alleles=(None, "AGC"))
    assert str(exc_info.value) == "May not disable trimming with VCF normalization mode"
