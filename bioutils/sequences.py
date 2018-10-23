# -*- coding: utf-8 -*-
"""simple functions and lookup tables for nucleic acid and amino acid sequences"""

from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import re

import six

from six.moves import range

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
    "Sec": "",
}

aa1_to_aa3_lut = {v: k for k, v in six.iteritems(aa3_to_aa1_lut)}

dna_to_aa1_lut = {      # NCBI standard translation table
    'AAA': 'K',
    'AAC': 'N',
    'AAG': 'K',
    'AAT': 'N',
    'ACA': 'T',
    'ACC': 'T',
    'ACG': 'T',
    'ACT': 'T',
    'AGA': 'R',
    'AGC': 'S',
    'AGG': 'R',
    'AGT': 'S',
    'ATA': 'I',
    'ATC': 'I',
    'ATG': 'M',
    'ATT': 'I',
    'CAA': 'Q',
    'CAC': 'H',
    'CAG': 'Q',
    'CAT': 'H',
    'CCA': 'P',
    'CCC': 'P',
    'CCG': 'P',
    'CCT': 'P',
    'CGA': 'R',
    'CGC': 'R',
    'CGG': 'R',
    'CGT': 'R',
    'CTA': 'L',
    'CTC': 'L',
    'CTG': 'L',
    'CTT': 'L',
    'GAA': 'E',
    'GAC': 'D',
    'GAG': 'E',
    'GAT': 'D',
    'GCA': 'A',
    'GCC': 'A',
    'GCG': 'A',
    'GCT': 'A',
    'GGA': 'G',
    'GGC': 'G',
    'GGG': 'G',
    'GGT': 'G',
    'GTA': 'V',
    'GTC': 'V',
    'GTG': 'V',
    'GTT': 'V',
    'TAA': '*',
    'TAC': 'Y',
    'TAG': '*',
    'TAT': 'Y',
    'TCA': 'S',
    'TCC': 'S',
    'TCG': 'S',
    'TCT': 'S',
    'TGA': '*',
    'TGC': 'C',
    'TGG': 'W',
    'TGT': 'C',
    'TTA': 'L',
    'TTC': 'F',
    'TTG': 'L',
    'TTT': 'F',
}


if six.PY2:                     # pragma: no cover
    # flake8: noqa

    import string
    from string import ascii_lowercase
    complement_transtable = {ord(f): ord(t) for f, t in zip("ACGT", "TGCA")}

elif six.PY3:                   # pragma: no cover
    # flake8: noqa

    from string import ascii_lowercase
    complement_transtable = bytes.maketrans(b"ACGT", b"TGCA")


def aa_to_aa1(seq):
    """coerce string of 1- or 3-letter amino acids to 1-letter

    >>> aa_to_aa1("CATSARELAME")
    'CATSARELAME'

    >>> aa_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
    'CATSARELAME'

    >>> aa_to_aa1(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return aa3_to_aa1(seq) if looks_like_aa3_p(seq) else seq


def aa_to_aa3(seq):
    """coerce string of 1- or 3-letter amino acids to 3-letter

    >>> aa_to_aa3("CATSARELAME")
    'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    >>> aa_to_aa3("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
    'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    >>> aa_to_aa3(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return aa1_to_aa3(seq) if not looks_like_aa3_p(seq) else seq


def aa1_to_aa3(seq):
    """convert string of 1-letter amino acids to 3-letter amino acids

    >>> aa1_to_aa3("CATSARELAME")
    'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    >>> aa1_to_aa3(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return "".join(aa1_to_aa3_lut[aa1] for aa1 in seq)


def aa3_to_aa1(seq):
    """convert string of 3-letter amino acids to 1-letter amino acids

    >>> aa3_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
    'CATSARELAME'

    >>> aa3_to_aa1(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return "".join(aa3_to_aa1_lut[aa3]
                   for aa3 in [seq[i:i + 3] for i in range(0, len(seq), 3)])


def complement(seq):
    """return the complement of a sequence

    >>> complement("ATCG")
    'TAGC'

    >>> complement(None)

    """

    if seq is None:
        return None
    seq = to_unicode(seq)
    return seq.translate(complement_transtable)


def elide_sequence(s, flank=5, elision="..."):
    """trim a sequence to include the left and right flanking sequences of
    size `flank`, with the intervening sequence elided by `elision`.

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
    """string looks like a 3-letter AA string"""
    return (seq is not None and (len(seq) % 3 == 0) and
            (len(seq) == 0 or seq[1] in ascii_lowercase))


def normalize_sequence(seq):
    """return normalized representation of sequence for hashing

    This really means ensuring that the sequence is represented as a
    binary blob and removing whitespace and asterisks and uppercasing.

    >>> normalize_sequence("ACGT")
    'ACGT'

    >>> normalize_sequence("  A C G T * ")
    'ACGT'

    >>> normalize_sequence("ACGT1")
    Traceback (most recent call last):
    ...
    RuntimeError: Normalized sequence contains non-alphabetic characters

    """

    assert isinstance(seq, six.text_type)
    nseq = re.sub("[\s\*]", "", seq).upper()
    m = re.search("[^A-Z]", nseq)
    if m:
        _logger.debug("Original sequence: " + seq)
        _logger.debug("Normalized sequence: " + nseq)
        _logger.debug("First non-[A-Z] at {}".format(m.start()))
        raise RuntimeError("Normalized sequence contains non-alphabetic characters")
    return nseq


def reverse_complement(seq):
    """return the reverse complement of a sequence

    >>> reverse_complement("ATCG")
    'CGAT'

    >>> reverse_complement(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return "".join(reversed(complement(seq)))


def replace_t_to_u(seq):
    """
    >>> replace_t_to_u("ACGT")
    'ACGU'

    >>> replace_t_to_u(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return seq.replace("T", "U").replace("t", "u")


def replace_u_to_t(seq):
    """
    >>> replace_u_to_t("ACGU")
    'ACGT'

    >>> replace_u_to_t(None)

    """
    if seq is None:
        return None
    seq = to_unicode(seq)
    return seq.replace("U", "T").replace("u", "t")


def to_unicode(s):
    return s if isinstance(s, six.text_type) else s.decode("ASCII")


def to_ascii(s):
    return s if isinstance(s, six.binary_type) else s.encode("ASCII")


def translate_cds(seq, full_codons=True, ter_symbol="*"):
    """translate a DNA or RNA sequence into a single-letter amino acid sequence
    using the standard translation table

    If full_codons is True, a sequence whose length isn't a multiple of three
    generates a ValueError; else an 'X' will be added as the last amino acid.
    This matches biopython's behaviour when padding the last codon with 'N's.

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

    seq = replace_u_to_t(seq)
    seq = seq.upper()

    protein_seq = list()
    for i in range(0, len(seq) - len(seq) % 3, 3):
        try:
            aa = dna_to_aa1_lut[seq[i:i + 3]]
        except KeyError:
            raise ValueError("Codon {} at position {}..{} is undefined in codon table".format(
                seq[i:i + 3], i+1, i+3))
        protein_seq.append(aa)

    # check for trailing bases and add the ter symbol if required
    if not full_codons and len(seq) % 3 != 0:
        protein_seq.append(ter_symbol)

    return ''.join(protein_seq)


# legacy equivalents
_looks_like_aa3_p = looks_like_aa3_p
_to_unicode = to_unicode
_to_binary = to_ascii
to_binary = to_ascii



## <LICENSE>
## Copyright 2014 Bioutils Contributors (https://bitbucket.org/biocommons/bioutils)
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
