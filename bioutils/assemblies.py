# -*- coding: utf-8 -*-, flake8: noqa
from __future__ import absolute_import, division, print_function, unicode_literals

"""Provides dictionaries of genome assembly data as provided by
ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/*.assembly.txt

Assemblies are stored in json files with the package in
_data/assemblies/.  (Those files are built with sbin/assembly-to-json,
also in this package.)

Definitions:

* accession (ac): symbol used to refer to a sequence (e.g., NC_000001.10)
* name: human-label (e.g., '1', 'MT', 'HSCHR6_MHC_APD_CTG1') that
  refers to a sequence, unique within some domain (e.g., GRCh37.p10)
* chromosome (chr): subset of names that refer to chromosomes 1..22, X, Y, MT
* aliases: list of other names; uniqueness unknown

.. note:: Some users prefer using a 'chr' prefix for chromosomes and
some don't.  Some prefer upper case and others prefer lower.  This
rift is unfortunate and creates unnecessary friction in sharing data.
You say TO-my-to and I say TO-mah-to doesn't apply here.  This code
favors using the authoritative names exactly as defined in the
assembly records.  Users are encouraged to use sequence names
verbatim, without prefixes or case changes.

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


def make_name_ac_map(assy_name, primary_only=False):
    """make map from sequence name to accession for given assembly name

    >>> grch38p5_name_ac_map = make_name_ac_map('GRCh38.p5')
    >>> grch38p5_name_ac_map['1']
    u'NC_000001.11'

    """
    return {s['name']: s['refseq_ac']
            for s in get_assembly(assy_name)['sequences']
            if (not primary_only or _is_primary(s))
            }


def make_ac_name_map(assy_name, primary_only=False):
    """make map from accession (str) to sequence name (str) for given assembly name

    >>> grch38p5_ac_name_map = make_ac_name_map('GRCh38.p5')
    >>> grch38p5_ac_name_map['NC_000001.11']
    u'1'

    """

    return {s['refseq_ac']: s['name']
            for s in get_assembly(assy_name)['sequences']
            if (not primary_only or _is_primary(s))
            }


def prepend_chr(chr):
    """prefix chr with 'chr' if not present

    Users are strongly encouraged to NOT use this function. Added a
    'chr' prefix means that you're using a name that is not consistent
    with authoritative assembly records.

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
# Internal functions

def _is_primary(s):
    return s['assembly_unit'] == 'Primary Assembly'


if __name__ == "__main__":
    import doctest
    doctest.testmod()
