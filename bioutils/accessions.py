# -*- coding: utf-8 -*-
"""simple routines to deal with accessions, identifiers, etc.

"""

from __future__ import absolute_import, division, print_function, unicode_literals

from six import iteritems


def prepend_chr(chr):
    """prefix chr with 'chr' if not present

    Users are strongly encouraged to NOT use this function. Added a
    'chr' prefix means that you're using a name that is not consistent
    with authoritative assembly records.

    >>> prepend_chr('22')
    'chr22'

    >>> prepend_chr('chr22')
    'chr22'

    """
    return chr if chr[0:3] == 'chr' else 'chr' + chr


def strip_chr(chr):
    """remove 'chr' prefix if it exists

    >>> strip_chr('22')
    '22'

    >>> strip_chr('chr22')
    '22'

    """
    return chr[3:] if chr[0:3] == 'chr' else chr


def chr22XY(c):
    """force to name from 1..22, 23, 24, X, Y, M 
    to in chr1..chr22, chrX, chrY, chrM
    str or ints accepted

    >>> chr22XY('1')
    'chr1'
    >>> chr22XY(1)
    'chr1'
    >>> chr22XY('chr1')
    'chr1'
    >>> chr22XY(23)
    'chrX'
    >>> chr22XY(24)
    'chrY'
    >>> chr22XY("X")
    'chrX'
    >>> chr22XY("23")
    'chrX'
    >>> chr22XY("M")
    'chrM'

    """
    c = str(c)
    if c[0:3] == 'chr':
        c = c[3:]
    if c == '23': c = 'X'
    if c == '24': c = 'Y'
    return 'chr' + c


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
