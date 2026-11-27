import itertools

import pytest

from bioutils.sequences import TranslationTable, translate_cds


def test_translate_examples():
    """test for standard translation table"""

    assert translate_cds("ATGCGA") == "MR"
    assert translate_cds("AUGCGA") == "MR"
    assert translate_cds(None) is None
    assert translate_cds("") == ""
    with pytest.raises(ValueError):
        translate_cds("AUGCG")

    assert translate_cds("AUGCG", full_codons=False) == "M*"
    assert translate_cds("ATGTAN") == "MX"
    assert translate_cds("CCN") == "P"
    assert translate_cds("TRA") == "*"
    assert translate_cds("TTNTA", full_codons=False) == "X*"
    assert translate_cds("CTB") == "L"
    assert translate_cds("AGM") == "X"
    assert translate_cds("GAS") == "X"
    assert translate_cds("CUN") == "L"
    with pytest.raises(ValueError):
        translate_cds("AUGCGQ")


def test_translate_selenoproteins():
    """unit test for sec codon"""
    assert translate_cds("AUGTGATAA") == "M**"
    assert translate_cds("AUGTGATAA", translation_table=TranslationTable.standard) == "M**"
    assert translate_cds("AUGTGATAA", translation_table=TranslationTable.selenocysteine) == "MU*"
    assert (
        translate_cds(
            "AUGTGATA",
            translation_table=TranslationTable.selenocysteine,
            full_codons=False,
        )
        == "MU*"
    )

    with pytest.raises(ValueError):
        translate_cds("AUGTGATA", translation_table=TranslationTable.selenocysteine)


def test_translate_vertebrate_mitochondrial():
    """unit test for vertebrate mitochondrial codons"""
    assert translate_cds("AUGTGATAA") == "M**"
    assert translate_cds("ATATGAAGGAGA", translation_table=TranslationTable.vertebrate_mitochondrial) == "MW**"
    assert (
        translate_cds(
            "ATAAG",
            translation_table=TranslationTable.vertebrate_mitochondrial,
            full_codons=False,
        )
        == "M*"
    )

    with pytest.raises(ValueError):
        translate_cds("ATAAG", translation_table=TranslationTable.vertebrate_mitochondrial)


def test_translate_degenerate_codons_selenocysteine():
    """TRA expands to TAA (*) and TGA (U), so it is ambiguous in the sec table."""
    assert translate_cds("TRA", translation_table=TranslationTable.selenocysteine) == "X"


def test_translate_degenerate_codons_vertebrate_mitochondrial():
    """Degenerate codons reflect NCBI genetic code table 2 (AGA/AGG=*, ATA=M, TGA=W)."""
    vmito = TranslationTable.vertebrate_mitochondrial
    assert translate_cds("AGR", translation_table=vmito) == "*"
    assert translate_cds("ATR", translation_table=vmito) == "M"
    assert translate_cds("TGR", translation_table=vmito) == "W"
    assert translate_cds("ATY", translation_table=vmito) == "I"
    for ambiguous in ("ATH", "ATM", "ATW", "MGA", "MGG", "MGR", "TRA"):
        assert translate_cds(ambiguous, translation_table=vmito) == "X"


@pytest.mark.parametrize(
    "table",
    [
        TranslationTable.standard,
        TranslationTable.selenocysteine,
        TranslationTable.vertebrate_mitochondrial,
    ],
)
def test_degenerate_codons_match_iupac_expansion(table):
    """Every IUPAC triplet translates to the shared AA iff all expansions agree."""
    iupac = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "AG",
        "Y": "CT",
        "S": "GC",
        "W": "AT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
        "N": "ACGT",
    }
    base = {
        "".join(t): translate_cds("".join(t), translation_table=table)
        for t in itertools.product("ACGT", repeat=3)
    }
    for c1, c2, c3 in itertools.product(iupac, repeat=3):
        codon = c1 + c2 + c3
        aas = {base[b1 + b2 + b3] for b1 in iupac[c1] for b2 in iupac[c2] for b3 in iupac[c3]}
        expected = next(iter(aas)) if len(aas) == 1 else "X"
        assert translate_cds(codon, translation_table=table) == expected
