# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

from .assemblies import get_assembly

import warnings
warnings.warn("GRCh37-specific components of bioutils.accessions are deprecated; see bioutils.assemblies")


def chr22XY(c):
    """force to name in chr1..chr22, chrX, chrY, chrM"""
    if c[0:3] == 'chr':
        c = c[3:]
    if c == '23': c = 'X'
    if c == '24': c = 'Y'
    return 'chr'+c


def prepend_chr(chr):
    """prefix chr with 'chr' if not present

    >>> prepend_chr('22')
    u'chr22'

    >>> prepend_chr('chr22')
    u'chr22'

    """
    return chr if chr[0:3] == 'chr' else 'chr' + chr


def strip_chr(chr):
    """remove 'chr' prefix if it exists

    >>> strip_chr('22')
    u'22'

    >>> strip_chr('chr22')
    u'22'

    """
    return chr[3:] if chr[0:3] == 'chr' else chr        



############################################################################
# Deprecated stuff
# This is all GRCh37 specific

_grch37 = get_assembly("GRCh37.p13")
primary_assembly_accessions = {
    'GRCh37': {s['refseq_ac'] for s in _grch37['sequences'] if s['refseq_ac'].startswith('NC_')},
    }
NC_to_chr_dict = {
    s['refseq_ac']: str(s['name'])
    for s in _grch37['sequences']
    if s['refseq_ac'].startswith('NC_')
    }
NC_to_chr = NC_to_chr_dict
chr_to_NC_dict = {v: k for k, v in NC_to_chr_dict.iteritems()}
chr_to_NC = chr_to_NC_dict
def chr_to_nc(s):
    return chr_to_NC_dict[s]
chr_size = {
    prepend_chr(s['name']): s['length']
    for s in _grch37['sequences']
    if s['refseq_ac'].startswith('NC_')
    }



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
