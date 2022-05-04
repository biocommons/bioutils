# -*- coding: utf-8 -*-
"""Simple functions and lookup tables for nucleic acid and amino acid sequences.
"""

import logging
import re
from enum import Enum
from string import ascii_lowercase

_logger = logging.getLogger(__name__)

aa3_to_aa1_lut = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Xaa": "X",
    "Ter": "*",
    "Sec": "U",
}

aa1_to_aa3_lut = {v: k for k, v in aa3_to_aa1_lut.items()}

dna_to_aa1_lut = {  # NCBI standard translation table
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAT": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGT": "S",
    "ATA": "I",
    "ATC": "I",
    "ATG": "M",
    "ATT": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAT": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAT": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TAA": "*",
    "TAC": "Y",
    "TAG": "*",
    "TAT": "Y",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "TGA": "*",
    "TGC": "C",
    "TGG": "W",
    "TGT": "C",
    "TTA": "L",
    "TTC": "F",
    "TTG": "L",
    "TTT": "F",
    # degenerate codons
    "AAR": "K",
    "AAY": "N",
    "ACB": "T",
    "ACD": "T",
    "ACH": "T",
    "ACK": "T",
    "ACM": "T",
    "ACN": "T",
    "ACR": "T",
    "ACS": "T",
    "ACV": "T",
    "ACW": "T",
    "ACY": "T",
    "AGR": "R",
    "AGY": "S",
    "ATH": "I",
    "ATM": "I",
    "ATW": "I",
    "ATY": "I",
    "CAR": "Q",
    "CAY": "H",
    "CCB": "P",
    "CCD": "P",
    "CCH": "P",
    "CCK": "P",
    "CCM": "P",
    "CCN": "P",
    "CCR": "P",
    "CCS": "P",
    "CCV": "P",
    "CCW": "P",
    "CCY": "P",
    "CGB": "R",
    "CGD": "R",
    "CGH": "R",
    "CGK": "R",
    "CGM": "R",
    "CGN": "R",
    "CGR": "R",
    "CGS": "R",
    "CGV": "R",
    "CGW": "R",
    "CGY": "R",
    "CTB": "L",
    "CTD": "L",
    "CTH": "L",
    "CTK": "L",
    "CTM": "L",
    "CTN": "L",
    "CTR": "L",
    "CTS": "L",
    "CTV": "L",
    "CTW": "L",
    "CTY": "L",
    "GAR": "E",
    "GAY": "D",
    "GCB": "A",
    "GCD": "A",
    "GCH": "A",
    "GCK": "A",
    "GCM": "A",
    "GCN": "A",
    "GCR": "A",
    "GCS": "A",
    "GCV": "A",
    "GCW": "A",
    "GCY": "A",
    "GGB": "G",
    "GGD": "G",
    "GGH": "G",
    "GGK": "G",
    "GGM": "G",
    "GGN": "G",
    "GGR": "G",
    "GGS": "G",
    "GGV": "G",
    "GGW": "G",
    "GGY": "G",
    "GTB": "V",
    "GTD": "V",
    "GTH": "V",
    "GTK": "V",
    "GTM": "V",
    "GTN": "V",
    "GTR": "V",
    "GTS": "V",
    "GTV": "V",
    "GTW": "V",
    "GTY": "V",
    "MGA": "R",
    "MGG": "R",
    "MGR": "R",
    "TAR": "*",
    "TAY": "Y",
    "TCB": "S",
    "TCD": "S",
    "TCH": "S",
    "TCK": "S",
    "TCM": "S",
    "TCN": "S",
    "TCR": "S",
    "TCS": "S",
    "TCV": "S",
    "TCW": "S",
    "TCY": "S",
    "TGY": "C",
    "TRA": "*",
    "TTR": "L",
    "TTY": "F",
    "YTA": "L",
    "YTG": "L",
    "YTR": "L",
}

# translation table for selenocysteine
dna_to_aa1_sec = dna_to_aa1_lut.copy()
dna_to_aa1_sec["TGA"] = "U"

complement_transtable = bytes.maketrans(b"ACGT", b"TGCA")


def aa_to_aa1(seq):
    """Coerces string of 1- or 3-letter amino acids to 1-letter representation.

    Args:
        seq (str): An amino acid sequence.

    Returns:
        str: The sequence as one of 1-letter amino acids.

    Examples:
        >>> aa_to_aa1("CATSARELAME")
        'CATSARELAME'

        >>> aa_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
        'CATSARELAME'

        >>> aa_to_aa1(None)
    """

    if seq is None:
        return None
    return aa3_to_aa1(seq) if looks_like_aa3_p(seq) else seq


def aa_to_aa3(seq):
    """Coerces string of 1- or 3-letter amino acids to 3-letter representation.

    Args:
        seq (str): An amino acid sequence.

    Returns:
        str: The sequence as one of 3-letter amino acids.

    Examples:
        >>> aa_to_aa3("CATSARELAME")
        'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

        >>> aa_to_aa3("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
        'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

        >>> aa_to_aa3(None)
    """

    if seq is None:
        return None
    return aa1_to_aa3(seq) if not looks_like_aa3_p(seq) else seq


def aa1_to_aa3(seq):
    """Converts string of 1-letter amino acids to 3-letter amino acids.

    Should only be used if the format of the sequence is known; otherwise use ``aa_to_aa3()``.

    Args:
        seq (str): An amino acid sequence as 1-letter amino acids.

    Returns:
        str: The sequence as 3-letter amino acids.

    Raises:
        KeyError: If the sequence is not of 1-letter amino acids.

    Examples:
        >>> aa1_to_aa3("CATSARELAME")
        'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

        >>> aa1_to_aa3(None)
    """

    if seq is None:
        return None
    return "".join(aa1_to_aa3_lut[aa1] for aa1 in seq)


def aa3_to_aa1(seq):
    """Converts string of 3-letter amino acids to 1-letter amino acids.

    Should only be used if the format of the sequence is known; otherwise use ``aa_to_aa1()``.

    Args:
        seq (str): An amino acid sequence as 3-letter amino acids.

    Returns:
        str: The sequence as 1-letter amino acids.

    Raises:
        KeyError: If the sequence is not of 3-letter amino acids.

    Examples:
        >>> aa3_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
        'CATSARELAME'

        >>> aa3_to_aa1(None)
    """

    if seq is None:
        return None
    return "".join(
        aa3_to_aa1_lut[aa3] for aa3 in [seq[i : i + 3] for i in range(0, len(seq), 3)]
    )


def complement(seq):
    """Retrieves the complement of a sequence.

    Args:
        seq (str): A nucleotide sequence.

    Returns:
        str: The complement of the sequence.

    Examples:
        >>> complement("ATCG")
        'TAGC'

        >>> complement(None)
    """

    if seq is None:
        return None
    return seq.translate(complement_transtable)


def elide_sequence(s, flank=5, elision="..."):
    """Trims the middle of the sequence, leaving the right and left flanks.

    Args:
        s (str): A sequence.
        flank (int, optional): The length of each flank. Defaults to five.
        elision (str, optional): The symbol used to represent the part trimmed. Defaults to '...'.

        Returns:
            str: The sequence with the middle replaced by ``elision``.

    Examples:
        >>> elide_sequence("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        'ABCDE...VWXYZ'

        >>> elide_sequence("ABCDEFGHIJKLMNOPQRSTUVWXYZ", flank=3)
        'ABC...XYZ'

        >>> elide_sequence("ABCDEFGHIJKLMNOPQRSTUVWXYZ", elision="..")
        'ABCDE..VWXYZ'

        >>> elide_sequence("ABCDEFGHIJKLMNOPQRSTUVWXYZ", flank=12)
        'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

        >>> elide_sequence("ABCDEFGHIJKLMNOPQRSTUVWXYZ", flank=12, elision=".")
        'ABCDEFGHIJKL.OPQRSTUVWXYZ'
    """

    elided_sequence_len = flank + flank + len(elision)
    if len(s) <= elided_sequence_len:
        return s
    return s[:flank] + elision + s[-flank:]


def looks_like_aa3_p(seq):
    """Indicates whether a string looks like a 3-letter AA string.

    Args:
        seq (str): A sequence.

    Returns:
        bool: Whether the string is of the format of a 3-letter AA string.
    """
    return (
        seq is not None
        and (len(seq) % 3 == 0)
        and (len(seq) == 0 or seq[1] in ascii_lowercase)
    )


def normalize_sequence(seq):
    """Converts sequence to normalized representation for hashing.

    Essentially, removes whitespace and asterisks, and uppercases the string.

    Args:
        seq (str): The sequence to be normalized.

    Returns:
        str: The sequence as a string of uppercase letters.

    Raises:
        RuntimeError: If the sequence contains non-alphabetic characters (besides '*').

    Examples:
        >>> normalize_sequence("ACGT")
        'ACGT'

        >>> normalize_sequence("  A C G T * ")
        'ACGT'

        >>> normalize_sequence("ACGT1")
        Traceback (most recent call last):
        ...
        RuntimeError: Normalized sequence contains non-alphabetic characters
    """

    nseq = re.sub(r"[\s\*]", "", seq).upper()
    m = re.search("[^A-Z]", nseq)
    if m:
        _logger.debug("Original sequence: " + seq)
        _logger.debug("Normalized sequence: " + nseq)
        _logger.debug("First non-[A-Z] at {}".format(m.start()))
        raise RuntimeError("Normalized sequence contains non-alphabetic characters")
    return nseq


def reverse_complement(seq):
    """Converts a sequence to its reverse complement.

    Args:
        seq (str): A nucleotide sequence.

    Returns:
        str: The reverse complement of the sequence.

    Examples:
        >>> reverse_complement("ATCG")
        'CGAT'

        >>> reverse_complement(None)
    """

    if seq is None:
        return None
    return "".join(reversed(complement(seq)))


def replace_t_to_u(seq):
    """Replaces the T's in a sequence with U's.

    Args:
        seq (str): A nucleotide sequence.

    Returns:
        str: The sequence with the T's replaced by U's.

    Examples:
        >>> replace_t_to_u("ACGT")
        'ACGU'

        >>> replace_t_to_u(None)
    """

    if seq is None:
        return None
    return seq.replace("T", "U").replace("t", "u")


def replace_u_to_t(seq):
    """Replaces the U's in a sequence with T's.

    Args:
        seq (str): A nucleotide sequence.

    Returns:
        str: The sequence with the U's replaced by T's.

    Examples:
        >>> replace_u_to_t("ACGU")
        'ACGT'

        >>> replace_u_to_t(None)
    """

    if seq is None:
        return None
    return seq.replace("U", "T").replace("u", "t")


class StrEnum(str, Enum):
    """utility class"""

    def __str__(self) -> str:
        return str.__str__(self)

    def __repr__(self) -> str:
        return str.__repr__(self)


class TranslationTable(StrEnum):
    """An enum that controls switching between standard and selenocysteine translation tables."""

    standard = "standard"
    selenocysteine = "sec"


def translate_cds(
    seq, full_codons=True, ter_symbol="*", translation_table=TranslationTable.standard
):
    """Translates a DNA or RNA sequence into a single-letter amino acid sequence.

    Args:
        seq (str): A nucleotide sequence.
        full_codons (bool, optional): If ``True``, forces sequence to have length
            that is a multiple of 3 and raises an error otherwise.
            If False, ``ter_symbol`` will be added as the last amino acid.
            This corresponds to biopython's behavior of padding the last codon with ``N``s.
            Defaults to ``True``.
        ter_symbol (str, optional): Placeholder for the last amino acid if
            sequence length is not divisible by three and ``full_codons`` is False.
            Defaults to ``'*'``
        translation_table (TranslationTable, optional): One of the options from the TranslationTable. It indicates
            which codon to amino acid translation table to use.
            By default we will use the standard translation table for humans. To enable translation for selenoproteins,
            the TranslationTable.selenocysteine table can get used

    Returns:
        str: The corresponding single letter amino acid sequence.

    Raises:
        ValueError: If ``full_codons`` and the sequence is not a multiple of three.
        ValueError: If a codon is undefined in the table.

    Examples:
        >>> translate_cds("ATGCGA")
        'MR'

        >>> translate_cds("AUGCGA")
        'MR'

        >>> translate_cds(None)

        >>> translate_cds("")
        ''

        >>> translate_cds("AUGCG")
        Traceback (most recent call last):
        ...
        ValueError: Sequence length must be a multiple of three

        >>> translate_cds("AUGCG", full_codons=False)
        'M*'

        >>> translate_cds("ATGTAN")
        'MX'

        >>> translate_cds("CCN")
        'P'

        >>> translate_cds("TRA")
        '*'

        >>> translate_cds("TTNTA", full_codons=False)
        'X*'

        >>> translate_cds("CTB")
        'L'

        >>> translate_cds("AGM")
        'X'

        >>> translate_cds("GAS")
        'X'

        >>> translate_cds("CUN")
        'L'

        >>> translate_cds("AUGCGQ")
        Traceback (most recent call last):
        ...
        ValueError: Codon CGQ at position 4..6 is undefined in codon table
    """

    if seq is None:
        return None

    if len(seq) == 0:
        return ""

    if full_codons and len(seq) % 3 != 0:
        raise ValueError("Sequence length must be a multiple of three")

    if translation_table == TranslationTable.standard:
        trans_table = dna_to_aa1_lut
    elif translation_table == TranslationTable.selenocysteine:
        trans_table = dna_to_aa1_sec
    else:
        raise ValueError("Unsupported translation table {}".format(translation_table))
    seq = replace_u_to_t(seq)
    seq = seq.upper()

    protein_seq = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i : i + 3]
        try:
            aa = trans_table[codon]
        except KeyError:
            # if this contains an ambiguous code, set aa to X, otherwise, throw error
            iupac_ambiguity_codes = "BDHVNUWSMKRYZ"
            if any(
                iupac_ambiguity_code in codon
                for iupac_ambiguity_code in iupac_ambiguity_codes
            ):
                aa = "X"
            else:
                raise ValueError(
                    "Codon {} at position {}..{} is undefined in codon table".format(
                        codon, i + 1, i + 3
                    )
                )
        protein_seq.append(aa)

    # check for trailing bases and add the ter symbol if required
    if not full_codons and len(seq) % 3 != 0:
        protein_seq.append(ter_symbol)

    return "".join(protein_seq)


## <LICENSE>
## Copyright 2014 Bioutils Contributors
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>
