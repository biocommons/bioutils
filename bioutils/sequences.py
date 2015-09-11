# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

"""simple functions and lookup tables for nucleic acid and amino acid sequences"""

import string


aa3_to_aa1_lut = {
    'Ala': 'A',    'Arg': 'R',    'Asn': 'N',    'Asp': 'D',
    'Cys': 'C',    'Gln': 'Q',    'Glu': 'E',    'Gly': 'G',
    'His': 'H',    'Ile': 'I',    'Leu': 'L',    'Lys': 'K',
    'Met': 'M',    'Phe': 'F',    'Pro': 'P',    'Ser': 'S',
    'Thr': 'T',    'Trp': 'W',    'Tyr': 'Y',    'Val': 'V',
    'Xaa': 'X',    'Ter': '*',    'Sec': 'U',
}
aa1_to_aa3_lut = {v: k for k, v in aa3_to_aa1_lut.iteritems()}

complement_transtable_str = string.maketrans('ACGT', 'TGCA')
complement_transtable_uni = dict(zip(map(ord, u'ACGT'), u'TGCA'))



def complement(s):
    """return the complement of a sequence

    >>> complement(b'ATCG')
    'TAGC'

    >>> complement(u'ATCG')
    u'TAGC'

    """
    if s is None:
        return None
    if isinstance(s, unicode):
        return s.translate(complement_transtable_uni)
    return s.translate(complement_transtable_str)


def reverse_complement(s):
    """return the reverse complement of a sequence

    >>> reverse_complement(b'ATCG')
    'CGAT'

    >>> reverse_complement(u'ATCG')
    u'CGAT'

    """
    return None if s is None else b''.join(reversed(complement(s)))


def replace_t_to_u(s):
    return s.replace('T', 'U').replace('t', 'u') if s else s


def replace_u_to_t(s):
    return s.replace('U', 'T').replace('u', 't') if s else s


def aa3_to_aa1(s):
    """convert string of 3-letter amino acids to 1-letter amino acids

    >>> aa3_to_aa1('CysAlaThrSerAlaArgGluLeuAlaMetGlu')
    u'CATSARELAME'

    """
    return None if s is None else ''.join([aa3_to_aa1_lut[aa3] for aa3 in
                                          [s[i:i + 3] for i in range(0, len(s), 3)]])


def aa1_to_aa3(s):
    """convert string of 1-letter amino acids to 3-letter amino acids

    >>> aa1_to_aa3('CATSARELAME')
    u'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    """
    return None if s is None else ''.join([aa1_to_aa3_lut[aa1] for aa1 in s])


def aa_to_aa1(s):
    """coerce string of 1- or 3-letter amino acids to 1-letter

    >>> aa_to_aa1('CATSARELAME')
    u'CATSARELAME'

    >>> aa_to_aa1('CysAlaThrSerAlaArgGluLeuAlaMetGlu')
    u'CATSARELAME'

    """
    return aa3_to_aa1(s) if __looks_like_aa3_p(s) else s


def aa_to_aa3(s):
    """coerce string of 1- or 3-letter amino acids to 3-letter

    >>> aa_to_aa3('CATSARELAME')
    u'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    >>> aa_to_aa3('CysAlaThrSerAlaArgGluLeuAlaMetGlu')
    u'CysAlaThrSerAlaArgGluLeuAlaMetGlu'

    """
    return aa1_to_aa3(s) if not __looks_like_aa3_p(s) else s


def __looks_like_aa3_p(s):
    "string looks like a 3-letter AA string"
    return (
        s is not None
        and (len(s) % 3 == 0)
        and (len(s) == 0 or s[1] in string.lowercase)
    )




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
