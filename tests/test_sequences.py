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
    assert (
        translate_cds("AUGTGATAA", translation_table=TranslationTable.standard) == "M**"
    )
    assert (
        translate_cds("AUGTGATAA", translation_table=TranslationTable.selenocysteine)
        == "MU*"
    )
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
