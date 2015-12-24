# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals

"""Provides dictionaries of genome assembly data as provided by
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/*.assembly.txt

Assemblies are stored in json files with the package in
_data/assemblies/.  (Those files are built with sbin/assembly-to-json,
also in this package.)

"""


import json
import pkg_resources

_assy_dir = '_data/assemblies'
_assy_path_fmt = _assy_dir + '/' + '{name}.json'


def get_assembly_names():
    """return list of available assemblies
    
    >>> assy_names = get_assembly_names()

    >>> 'GRCh37.p13' in assy_names
    True

    """
    
    return [n.replace(".json", "")
            for n in pkg_resources.resource_listdir(__name__, _assy_dir)
            if n.endswith(".json")
            ]


def get_assembly(name):
    """read a single assembly by name, returning a dictionary of assembly data

    >>> assy = get_assembly('GRCh37.p13')

    >>> assy['name']
    u'GRCh37.p13'

    >>> assy['description']
    u'Genome Reference Consortium Human Build 37 patch release 13 (GRCh37.p13)'

    >>> assy['refseq_ac']
    u'GCF_000001405.25'

    >>> assy['genbank_ac']
    u'GCA_000001405.14'

    >>> len(assy['sequences'])
    297

    >>> import pprint
    >>> pprint.pprint(assy['sequences'][0])
    {u'aliases': [u'chr1'],
     u'assembly_unit': u'Primary Assembly',
     u'genbank_ac': u'CM000663.1',
     u'length': 249250621,
     u'name': u'1',
     u'refseq_ac': u'NC_000001.10'}
    
    """

    rsrc_strm = pkg_resources.resource_stream(__name__, _assy_path_fmt.format(name=name))
    return json.load(rsrc_strm)


def get_assemblies(names=[]):
    """read specified assemblies, or all if none specified, returning a
    dictionary of assembly-name: assembly.  See get_assembly()
    function for the structure of assembly data.

    >>> assemblies = get_assemblies(names=['GRCh37.p13'])
    >>> assemblies.keys()
    [u'GRCh37.p13']
    >>> assy = assemblies['GRCh37.p13']

    >>> assemblies = get_assemblies()
    >>> u'GRCh38.p2' in assemblies.keys()
    True

    """

    if names == []:
        names = get_assembly_names()
    return {a['name']: a for a in (get_assembly(n) for n in names)}


if __name__ == "__main__":
    import doctest
    doctest.testmod()
