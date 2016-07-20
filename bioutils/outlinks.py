# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals

from .accessions import chr_to_NC, strip_chr


def url_for_slice(c, s, e):
    url_fmt = (
        'http://www.ncbi.nlm.nih.gov/projects/sviewer/?id={nc}'
        '&noslider=1&tracks=[key:sequence_track,name:Sequence,'
        'display_name:Sequence,category:Sequence,annots:Sequence'
        ',ShowLabel:false][key:gene_model_track,name:Genes---Unnamed'
        ',display_name:Genes,category:Genes,annots:Unnamed,Options'
        ':MergeAll,SNPs:false,CDSProductFeats:false,'
        'ShowLabelsForAllFeatures:false,HighlightMode:2][key:'
        'alignment_track,name:Alignments---NG%20Alignments,'
        'display_name:NG%20Alignments,category:Alignments,'
        'annots:NG%20Alignments,Layout:Adaptive1000,StatDisplay'
        ':-1,LinkMatePairAligns:true,Color:true,AlignedSeqFeats'
        ':false,Label:true][key:alignment_track,name:Alignments'
        '---Refseq%20Alignments,display_name:Refseq%20Alignments'
        ',category:Alignments,annots:Refseq%20Alignments,Layout'
        ':Adaptive1000,StatDisplay:-1,LinkMatePairAligns:true'
        ',Color:true,AlignedSeqFeats:false,Label:true]&mk=&color'
        '=0&label=0&decor=0&spacing=0&v={start}:{end}&c='
        '33cccc&gflip=false&select='
        'gi|224589812-0390e571-039b5700-0103-dee8c202-ffea8d58'
    )
    return url_fmt.format(nc=chr_to_NC[strip_chr(c)], start=str(s), end=str(e))


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
