# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals
"""simple functions and lookup tables for nucleic acid and amino acid sequences"""

import six

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


def _to_unicode(s):
    return s if isinstance(s, six.text_type) else s.decode()


if six.PY2:  # pragma: no cover

    import string
    from string import ascii_lowercase
    complement_transtable = {ord(f): ord(t) for f, t in zip("ACGT", "TGCA")}

elif six.PY3:  # pragma: no cover

    from string import ascii_lowercase
    complement_transtable = bytes.maketrans(b"ACGT", b"TGCA")


def complement(seq):
    """return the complement of a sequence

    >>> complement("ATCG")
    'TAGC'

    >>> complement(None)

    """

    if seq is None:
        return None
    seq = _to_unicode(seq)
    return seq.translate(complement_transtable)


def reverse_complement(seq):
    """return the reverse complement of a sequence

    >>> reverse_complement("ATCG")
    'CGAT'

    >>> reverse_complement(None)

    """
    if seq is None:
        return None
    seq = _to_unicode(seq)
    return "".join(reversed(complement(seq)))


def replace_t_to_u(seq):
    """
    >>> replace_t_to_u("ACGT")
    'ACGU'

    >>> replace_t_to_u(None)

    """
    if seq is None:
        return None
    seq = _to_unicode(seq)
    return seq.replace("T", "U").replace("t", "u")


def replace_u_to_t(seq):
    """
    >>> replace_u_to_t("ACGU")
    'ACGT'

    >>> replace_u_to_t(None)

    """
    if seq is None:
        return None
    seq = _to_unicode(seq)
    return seq.replace("U", "T").replace("u", "t")


def aa3_to_aa1(seq):
    """convert string of 3-letter amino acids to 1-letter amino acids

    >>> aa3_to_aa1("CysAlaThrSerAlaArgGluLeuAlaMetGlu")
    'CATSARELAME'

    >>> aa3_to_aa1(None)

    """
    if seq is None:
        return None
    seq = _to_unicode(seq)
    return "".join(aa3_to_aa1_lut[aa3]
                   for aa3 in [seq[i:i + 3] for i in range(0, len(seq), 3)])


def aa1_to_aa3(seq):
    """convert string of 1-letter amino acids to 3-letter amino acids

    >>> aa1_to_aa3("CATSARELAME")
    'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    >>> aa1_to_aa3(None)

    """
    if seq is None:
        return None
    seq = _to_unicode(seq)
    return "".join(aa1_to_aa3_lut[aa1] for aa1 in seq)


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
    seq = _to_unicode(seq)
    return aa3_to_aa1(seq) if __looks_like_aa3_p(seq) else seq


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
    seq = _to_unicode(seq)
    return aa1_to_aa3(seq) if not __looks_like_aa3_p(seq) else seq


def __looks_like_aa3_p(seq):
    """string looks like a 3-letter AA string"""
    return (seq is not None and (len(seq) % 3 == 0) and
            (len(seq) == 0 or seq[1] in ascii_lowercase))


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
