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

    return [
        n.replace(".json", "")
        for n in pkg_resources.resource_listdir(__name__, _assy_dir)
        if n.endswith(".json")
    ]


def get_assembly(name):
    """read a single assembly by name, returning a dictionary of assembly data

    >>> assy = get_assembly('GRCh37.p13')

    >>> assy['name']
    'GRCh37.p13'

    >>> assy['description']
    'Genome Reference Consortium Human Build 37 patch release 13 (GRCh37.p13)'

    >>> assy['refseq_ac']
    'GCF_000001405.25'

    >>> assy['genbank_ac']
    'GCA_000001405.14'

    >>> len(assy['sequences'])
    297

    >>> import pprint
    >>> pprint.pprint(assy['sequences'][0])
    {'aliases': ['chr1'],
     'assembly_unit': 'Primary Assembly',
     'genbank_ac': 'CM000663.1',
     'length': 249250621,
     'name': '1',
     'refseq_ac': 'NC_000001.10',
     'relationship': '=',
     'sequence_role': 'assembled-molecule'}
    """

    js = pkg_resources.resource_string(
        __name__, _assy_path_fmt.format(name=name))
    return json.loads(js.decode("utf-8"))


def get_assemblies(names=[]):
    """read specified assemblies, or all if none specified, returning a
    dictionary of assembly-name: assembly.  See get_assembly()
    function for the structure of assembly data.

    >>> assemblies = get_assemblies(names=['GRCh37.p13'])
    >>> assy = assemblies['GRCh37.p13']

    >>> assemblies = get_assemblies()
    >>> 'GRCh38.p2' in assemblies
    True

    """

    if names == []:
        names = get_assembly_names()
    return {a['name']: a for a in (get_assembly(n) for n in names)}


def make_name_ac_map(assy_name, primary_only=False):
    """make map from sequence name to accession for given assembly name

    >>> grch38p5_name_ac_map = make_name_ac_map('GRCh38.p5')
    >>> grch38p5_name_ac_map['1']
    'NC_000001.11'

    """
    return {
        s['name']: s['refseq_ac']
        for s in get_assembly(assy_name)['sequences']
        if (not primary_only or _is_primary(s))
    }


def make_ac_name_map(assy_name, primary_only=False):
    """make map from accession (str) to sequence name (str) for given assembly name

    >>> grch38p5_ac_name_map = make_ac_name_map('GRCh38.p5')
    >>> grch38p5_ac_name_map['NC_000001.11']
    '1'

    """

    return {
        s['refseq_ac']: s['name']
        for s in get_assembly(assy_name)['sequences']
        if (not primary_only or _is_primary(s))
    }


############################################################################
# Internal functions


def _is_primary(s):
    """return True if this sequence record is part of the primary assembly

    >>> _is_primary({'assembly_unit': 'Primary Assembly'})
    True
    
    >>> _is_primary({'assembly_unit': 'Something else entirely'})
    False

    """
    return s['assembly_unit'] == 'Primary Assembly'
