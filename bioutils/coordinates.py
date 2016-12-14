# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals
"""provides utilities for interconverting between coordinate systems
especially as used by the hgvs code.  The three systems are:
 
                  : A : C : G : T : A : C :
  human/hgvs  h   :-3 :-2 :-1 : 1 : 2 : 3 :
  continuous  c   :-2 :-1 : 0 : 1 : 2 : 3 :
  interbase   i  -3  -2  -1   0   1   2   3

Human/hgvs coordinates are the native coordinates used by the HGVS
recommendations. The coordinates are 1-based, inclusive, and refer to
the nucleotides; there is no 0.

Continuous coordinates are similar to hgvs coordinates, but adds 1 to
all negative values so that there is no discontinuity between -1 and 1
(as there is with HGVS).

Interbase coordinates refer to the zero-width junctions between
nucleotides.  The main advantage of interbase coordinates is that
there are no corner cases in the specification of intervals used for
insertions and deletions as there is with numbering systems that refer
to nucleotides themselves.  Numerically, interbase intervals are
0-based, left-closed, and right-open.  Beacuse referring to a single
interbase coordinate is not particularly meaningful, interbase
coordinates are always passed as start,end pairs.

Because it's easy to confuse these coordinates in code, _h, _c, and _i
suffixes are often used to clarify variables.

For code clarity, this module provides functions that interconvert
*intervals* specified in each of the coordinate systems.

"""

PLUS_STRAND = 1
MINUS_STRAND = -1


def strand_pm_to_int(s):
    """
    >>> strand_pm_to_int('+')
    1
    >>> strand_pm_to_int('-')
    -1
    >>> strand_pm_to_int('arglefargle')

    """
    return PLUS_STRAND if s == '+' else MINUS_STRAND if s == '-' else None


def strand_int_to_pm(i):
    """
    >>> strand_int_to_pm(1)
    '+'
    >>> strand_int_to_pm(-1)
    '-'
    >>> strand_int_to_pm(42)

    """
    return '+' if i == PLUS_STRAND else '-' if i == MINUS_STRAND else None


strand_pm = strand_int_to_pm

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
