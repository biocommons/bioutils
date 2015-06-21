# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals

"""Provides dictionaries of genome assembly data as provided by
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/*.assembly.txt

Assemblies are stored in json files with the package.  Those files are
built with sbin/assembly-to-json, also in this package.

"""


import pkg_resources

import json


_assy_dir = '_data/assemblies'
_assy_path_fmt = _assy_dir + '/' + '{name}.json'


def _get_assembly_names():
    """[INTERNAL ONLY] return list of available assemblies
    
    >>> assy_names = _get_assembly_names()

    >>> 'GRCh37.p13' in assy_names
    True

    """
    
    return [n.replace(".json", "")
            for n in pkg_resources.resource_listdir(__name__, _assy_dir)
            if n.endswith(".json")
            ]

def _get_assembly(name):
    """[INTERNAL ONLY] read a single assembly by name, returning a dictionary of assembly data

    """

    rsrc_strm = pkg_resources.resource_stream(__name__, _assy_path_fmt.format(name=name))
    return json.load(rsrc_strm)


def get_assemblies(names=[]):
    """read specified assemblies, or all if none specified, returning a
    dictionary of name: assembly-info
    
    >>> assemblies = get_assemblies(names=['GRCh37.p13'])
    >>> assemblies.keys()
    [u'GRCh37.p13']

    >>> assemblies = get_assemblies()
    >>> u'GRCh38.p2' in assemblies.keys()
    True

    >>> assy = assemblies['GRCh37.p13']

    >>> assy['description']
    u'Genome Reference Consortium Human Build 37 patch release 13 (GRCh37.p13)'

    >>> assy['submitter']
    u'Genome Reference Consortium'

    >>> assy['date']
    u'2013-06-28'

    >>> len(assy['sequences'])
    297

    >>> rec = [r for r in assy['sequences'] if r['genbank_ac'] == 'CM000663.1'][0]
    >>> rec['refseq_ac'], rec['length'], rec['name'], rec['aliases']
    (u'NC_000001.10', 249904550, u'1', [u'chr1'])

    >>> rec = [r for r in assy['sequences'] if r['genbank_ac'] == 'GL000198.1'][0]
    >>> rec['refseq_ac'], rec['length'], rec['name'], rec['aliases']
    (u'NT_113914.1', 90085, u'HSCHR9_RANDOM_CTG1', [u'chr9_gl000198_random'])

    """

    if names == []:
        names = _get_assembly_names()
    return {a['name']: a for a in (_get_assembly(n) for n in names)}


if __name__ == "__main__":
    import doctest
    doctest.testmod()
