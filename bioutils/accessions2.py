# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals

import re

import recordtype


class ContigInfo(recordtype.recordtype('ContigInfo',
                                       field_names = ['ac', 'len', 'md5', 'name', 'descr'],
                                       default = None)):
    pass


contig_info = {
    # This data was pulled from several sources. Ideally, we'd build/identify an automated process for this
    'GRCh37': [
        ContigInfo('NC_000001.10'  , 249250621 ,  '65fd477c90391caf96ac6acd936325e6' , 'chr1'  , None),
        ContigInfo('NC_000002.11'  , 243199373 ,  'c998a4be509c7d9db4bea67bab6f3ff4' , 'chr2'  , None),
        ContigInfo('NC_000003.11'  , 198022430 ,  'fb432c7e55d43c4b24610de1f71ceb66' , 'chr3'  , None),
        ContigInfo('NC_000004.11'  , 191154276 ,  'eeab3bc1789bf45593d463bfb1f87434' , 'chr4'  , None),
        ContigInfo('NC_000005.9'   , 180915260 ,  '0f105ea7d3e52a358e6a95240f0e7163' , 'chr5'  , None),
        ContigInfo('NC_000006.11'  , 171115067 ,  '094b3e53ac3332068c326631fa8d56b1' , 'chr6'  , None),
        ContigInfo('NC_000007.13'  , 159138663 ,  '8cd03bbd73290cfdcdf20313d08d2328' , 'chr7'  , None),  
        ContigInfo('NC_000008.10'  , 146364022 ,  '5bc6ccb4932a7fbf89e1a2b590590244' , 'chr8'  , None),
        ContigInfo('NC_000009.11'  , 141213431 ,  '6083f0d2206564c2bdd502f62a5ec453' , 'chr9'  , None),
        ContigInfo('NC_000010.10'  , 135534747 ,  'de25534493509079ebc13d454046fe81' , 'chr10' , None),
        ContigInfo('NC_000011.9'   , 135006516 ,  '3e8a88b2ef4af16238d960650d0d1d5a' , 'chr11' , None),
        ContigInfo('NC_000012.11'  , 133851895 ,  '8aa780adec82a3dfdb2051a7c217393c' , 'chr12' , None),
        ContigInfo('NC_000013.10'  , 115169878 ,  '44f0badb334fecaa379c048cdf12fe00' , 'chr13' , None),
        ContigInfo('NC_000014.8'   , 107349540 ,  '32099c663397885e63e3267e55160aa1' , 'chr14' , None),
        ContigInfo('NC_000015.9'   , 102531392 ,  '3bfa56ac1d4aa240f03e0a260f821871' , 'chr15' , None),
        ContigInfo('NC_000016.9'   ,  90354753 ,  '86dbc200b8cade10a9ae575abc27b09d' , 'chr16' , None),
        ContigInfo('NC_000017.10'  ,  81195210 ,  '9fe3cf8faf16177d4a399f4b3087865d' , 'chr17' , None),
        ContigInfo('NC_000018.9'   ,  78077248 ,  '8b24694abeeb3e1b4bece7aba928be92' , 'chr18' , None),
        ContigInfo('NC_000019.9'   ,  63025520 ,  '49d2652a2c77be977389d1888db438e8' , 'chr20' , None),
        ContigInfo('NC_000020.10'  ,  59128983 ,  'fc5729e25602436f16af7a6538f9399c' , 'chr19' , None),
        ContigInfo('NC_000021.8'   ,  51304566 ,  'e9522b06dc2e59ad4ddbf9f64d36eba2' , 'chr22' , None),
        ContigInfo('NC_000022.10'  ,  48129895 ,  '0547fd82e81ed1d821a3675be36bb396' , 'chr21' , None),
        ContigInfo('NC_000023.10'  , 155270560 ,  '0e89b6ab3065727482ba8b7437926fbf' , 'chrX'  , None),
        ContigInfo('NC_000024.9'   ,  59373566 ,  '92fce73ff4bc1b7146790293ce1d6554' , 'chrY'  , None),
        ContigInfo(None            ,    16571  ,                                None , 'chrM'  , None),
        ContigInfo(None            ,    40103  ,                                None , 'chr11_gl000202_random', None),
        ContigInfo(None            ,  1680828  ,                                None , 'chr17_ctg5_hap1', None),
        ContigInfo(None            ,    37498  ,                                None , 'chr17_gl000203_random', None),
        ContigInfo(None            ,    81310  ,                                None , 'chr17_gl000204_random', None),
        ContigInfo(None            ,   174588  ,                                None , 'chr17_gl000205_random', None),
        ContigInfo(None            ,    41001  ,                                None , 'chr17_gl000206_random', None),
        ContigInfo(None            ,     4262  ,                                None , 'chr18_gl000207_random', None),
        ContigInfo(None            ,    92689  ,                                None , 'chr19_gl000208_random', None),
        ContigInfo(None            ,   159169  ,                                None , 'chr19_gl000209_random', None),
        ContigInfo(None            ,   106433  ,                                None , 'chr1_gl000191_random', None),
        ContigInfo(None            ,   547496  ,                                None , 'chr1_gl000192_random', None),
        ContigInfo(None            ,    27682  ,                                None , 'chr21_gl000210_random', None),
        ContigInfo(None            ,   590426  ,                                None , 'chr4_ctg9_hap1', None),
        ContigInfo(None            ,   189789  ,                                None , 'chr4_gl000193_random', None),
        ContigInfo(None            ,   191469  ,                                None , 'chr4_gl000194_random', None),
        ContigInfo(None            ,  4622290  ,                                None , 'chr6_apd_hap1', None),
        ContigInfo(None            ,  4795371  ,                                None , 'chr6_cox_hap2', None),
        ContigInfo(None            ,  4610396  ,                                None , 'chr6_dbb_hap3', None),
        ContigInfo(None            ,  4683263  ,                                None , 'chr6_mann_hap4', None),
        ContigInfo(None            ,  4833398  ,                                None , 'chr6_mcf_hap5', None),
        ContigInfo(None            ,  4611984  ,                                None , 'chr6_qbl_hap6', None),
        ContigInfo(None            ,  4928567  ,                                None , 'chr6_ssto_hap7', None),
        ContigInfo(None            ,   182896  ,                                None , 'chr7_gl000195_random', None),
        ContigInfo(None            ,    38914  ,                                None , 'chr8_gl000196_random', None),
        ContigInfo(None            ,    37175  ,                                None , 'chr8_gl000197_random', None),
        ContigInfo(None            ,    90085  ,                                None , 'chr9_gl000198_random', None),
        ContigInfo(None            ,   169874  ,                                None , 'chr9_gl000199_random', None),
        ContigInfo(None            ,   187035  ,                                None , 'chr9_gl000200_random', None),
        ContigInfo(None            ,    36148  ,                                None , 'chr9_gl000201_random', None),
        ContigInfo(None            ,   166566  ,                                None , 'chrUn_gl000211', None),
        ContigInfo(None            ,   186858  ,                                None , 'chrUn_gl000212', None),
        ContigInfo(None            ,   164239  ,                                None , 'chrUn_gl000213', None),
        ContigInfo(None            ,   137718  ,                                None , 'chrUn_gl000214', None),
        ContigInfo(None            ,   172545  ,                                None , 'chrUn_gl000215', None),
        ContigInfo(None            ,   172294  ,                                None , 'chrUn_gl000216', None),
        ContigInfo(None            ,   172149  ,                                None , 'chrUn_gl000217', None),
        ContigInfo(None            ,   161147  ,                                None , 'chrUn_gl000218', None),
        ContigInfo(None            ,   179198  ,                                None , 'chrUn_gl000219', None),
        ContigInfo(None            ,   161802  ,                                None , 'chrUn_gl000220', None),
        ContigInfo(None            ,   155397  ,                                None , 'chrUn_gl000221', None),
        ContigInfo(None            ,   186861  ,                                None , 'chrUn_gl000222', None),
        ContigInfo(None            ,   180455  ,                                None , 'chrUn_gl000223', None),
        ContigInfo(None            ,   179693  ,                                None , 'chrUn_gl000224', None),
        ContigInfo(None            ,   211173  ,                                None , 'chrUn_gl000225', None),
        ContigInfo(None            ,    15008  ,                                None , 'chrUn_gl000226', None),
        ContigInfo(None            ,   128374  ,                                None , 'chrUn_gl000227', None),
        ContigInfo(None            ,   129120  ,                                None , 'chrUn_gl000228', None),
        ContigInfo(None            ,    19913  ,                                None , 'chrUn_gl000229', None),
        ContigInfo(None            ,    43691  ,                                None , 'chrUn_gl000230', None),
        ContigInfo(None            ,    27386  ,                                None , 'chrUn_gl000231', None),
        ContigInfo(None            ,    40652  ,                                None , 'chrUn_gl000232', None),
        ContigInfo(None            ,    45941  ,                                None , 'chrUn_gl000233', None),
        ContigInfo(None            ,    40531  ,                                None , 'chrUn_gl000234', None),
        ContigInfo(None            ,    34474  ,                                None , 'chrUn_gl000235', None),
        ContigInfo(None            ,    41934  ,                                None , 'chrUn_gl000236', None),
        ContigInfo(None            ,    45867  ,                                None , 'chrUn_gl000237', None),
        ContigInfo(None            ,    39939  ,                                None , 'chrUn_gl000238', None),
        ContigInfo(None            ,    33824  ,                                None , 'chrUn_gl000239', None),
        ContigInfo(None            ,    41933  ,                                None , 'chrUn_gl000240', None),
        ContigInfo(None            ,    42152  ,                                None , 'chrUn_gl000241', None),
        ContigInfo(None            ,    43523  ,                                None , 'chrUn_gl000242', None),
        ContigInfo(None            ,    43341  ,                                None , 'chrUn_gl000243', None),
        ContigInfo(None            ,    39929  ,                                None , 'chrUn_gl000244', None),
        ContigInfo(None            ,    36651  ,                                None , 'chrUn_gl000245', None),
        ContigInfo(None            ,    38154  ,                                None , 'chrUn_gl000246', None),
        ContigInfo(None            ,    36422  ,                                None , 'chrUn_gl000247', None),
        ContigInfo(None            ,    39786  ,                                None , 'chrUn_gl000248', None),
        ContigInfo(None            ,    38502  ,                                None , 'chrUn_gl000249', None),
    ],
}   
        

def chr_to_nc(s):
    return chr_to_NC_dict[s]

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
## LEGACY

primary_assembly_accessions = {
    'GRCh37': { ci.ac:None for ci in contig_info['GRCh37'] if ci.ac },
    }   
        
NC_to_chr_dict = {
    ci.ac:strip_chr(ci.name) for ci in contig_info['GRCh37'] if re.match('chr(\d+|X|Y)',ci.name)
   }
NC_to_chr = NC_to_chr_dict
        
chr_to_NC_dict = {v: k for k, v in NC_to_chr_dict.iteritems()}
chr_to_NC = chr_to_NC_dict
        
chr_size_by_assembly = {
    'GRCh37': {ci.name:ci.len for ci in contig_info['GRCh37']},
    }

chr_size = chr_size_by_assembly['GRCh37']






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
